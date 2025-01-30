import logging
import numpy as np
import pandas as pd
import pyomo.environ as pm
import pypsa
from functions import CRF

class Network:
    """
    A class representing a Network.

    Attributes
    ----------
    type : string
        type of fuel demand. - maybe plant type instead of fuel type
    generators : dictionary
        contains generator types with their potential and maximum capacity.
    n :
        network. Default is None
    Methods
    -------
    set_network(demand_profile, times, country_series):
        sets up the network.
    set_generators_in_network(country_series):
        sets provided generator in the network.
    _create_override_components():
        set up new component attributes as required.
    """
    def __init__(self, type, generators):
        """
        
        """
        self.type = type
        self.generators = generators
        self.n = None

    def set_network(self, demand_profile, times, country_series):
        '''
        Sets up the network.

        Parameters
        ----------
        demand_profile : pandas DataFrame
            hourly dataframe of hydrogen demand in kg.
        times : xarray DataArray
            1D dataarray with timestamps for wind and solar potential.
        country_series: pandas Series
            interest rate and lifetime information.
        '''
        # Set the time values for the network
        if self.n == None:
            self.n = pypsa.Network(override_component_attrs=self._create_override_components())

        if self.type == "hydrogen":
            self.n.set_snapshots(times)

            # Import the design of the H2 plant into the network
            self.n.import_from_csv_folder("parameters/basic_h2_plant")

            # Import demand profile
            # Note: All flows are in MW or MWh, conversions for hydrogen done using HHVs. Hydrogen HHV = 39.4 MWh/t
            self.n.add('Load',
                'Hydrogen demand',
                bus = 'Hydrogen',
                p_set = demand_profile['Demand']/1000*39.4,  # Should clean up these parameters or at least comment re units
                       # we probably want to make HHV an input in the end via snakemake or a file so this can be generalised
                )

            for item in [self.n.links, self.n.stores, self.n.storage_units]:
                item.capital_cost = item.capital_cost * CRF(country_series['Plant interest rate'], 
                                                            country_series['Plant lifetime (years)'])
        elif self.type == "ammonia":
            # Set the time values for the network
            self.n.set_snapshots(demand_profile.index)
            demand_profile['weights'] = 8760 / len(self.n.snapshots)
            self.n.snapshot_weightings = demand_profile['weights']

            # Import the design of the NH3 plant into the network
            self.n.import_from_csv_folder("parameters/basic_nh3_plant")

            # Import demand profile
            # Note: All flows are in MW or MWh, conversions for hydrogen done using HHVs. Hydrogen HHV = 39.4 MWh/t
            # Note: All flows are in MW or MWh, conversions for ammonia done using HHVs. Ammonia HHV = 6.25 MWh/t
            # hydrogen_demand = pd.read_excel(demand_path,index_col = 0) # Excel file in kg hydrogen, convert to MWh
            self.n.add('Load',
                'Ammonia demand',
                bus='Ammonia',
                p_set=demand_profile['Demand'].to_numpy() / 1000 * 6.25,
                )
            
            for item in [self.n.links, self.n.stores]:
                item.capital_cost = item.capital_cost * CRF(country_series['Plant interest rate'],
                                                            country_series['Plant lifetime (years)'])
            self.n.links.loc['HydrogenCompression', 'marginal_cost'] = 0.0001  # Just stops pointless cycling through storage

    def set_generators_in_network(self, country_series):
        '''
        Sets provided generator in the network.

        Parameters
        ----------
        country_series: pandas Series
            interest rate and lifetime information.
        '''
        # Send the weather data to the model
        for gen, gen_list in self.generators.items():
            self.n.generators_t.p_max_pu[gen] = gen_list[0]

            # specify maximum capacity based on land use
            self.n.generators.loc[gen,'p_nom_max'] = gen_list[1] # same as above

            # specify technology-specific and country-specific WACC and lifetime here
            self.n.generators.loc[gen,'capital_cost'] = self.n.generators.loc[gen,'capital_cost']\
                * CRF(country_series[f'{gen} interest rate'], country_series[f'{gen} lifetime (years)'])

    def _create_override_components(self):
        # I assume this is just so that we can have hydrogen and power both as buses? Hmm
        """Set up new component attributes as required"""
        # Modify the capacity of a link so that it can attach to 2 buses.
        override_component_attrs = pypsa.descriptors.Dict(
            {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
        )

        override_component_attrs["Link"].loc["bus2"] = [
            "string",
            np.nan,
            np.nan,
            "2nd bus",
            "Input (optional)",
        ]
        override_component_attrs["Link"].loc["efficiency2"] = [
            "static or series",
            "per unit",
            1.0,
            "2nd bus efficiency",
            "Input (optional)",
        ]
        override_component_attrs["Link"].loc["p2"] = [
            "series",
            "MW",
            0.0,
            "2nd bus output",
            "Output",
        ]
        
        return override_component_attrs
    
def nh3_pyomo_constraints(n, snapshots):
    """Includes a series of additional constraints which make the ammonia plant work as needed:
    i) Battery sizing
    ii) Ramp hard constraints down (Cannot be violated)
    iii) Ramp hard constraints up (Cannot be violated)
    iv) Ramp soft constraints down
    v) Ramp soft constraints up
    (iv) and (v) just softly suppress ramping so that the model doesn't 'zig-zag', which looks a bit odd on operation.
    Makes very little difference on LCOA. """
    timestep = 3 # -- is this the same as freq?
    # The battery constraint is built here - it doesn't need a special function because it doesn't depend on time
    n.model.battery_interface = pm.Constraint(
        rule=lambda model: n.model.link_p_nom['BatteryInterfaceIn'] ==
                        n.model.link_p_nom['BatteryInterfaceOut'] /
                        n.links.efficiency["BatteryInterfaceOut"])

    # Constrain the maximum discharge of the H2 storage relative to its size
    time_step_cycle = 4/8760*timestep*0.5  # Factor 0.5 for 3 hour time step, 0.5 for oversized storage
    n.model.cycling_limit = pm.Constraint(
        rule=lambda model: n.model.link_p_nom['BatteryInterfaceOut'] ==
                        n.model.store_e_nom['CompressedH2Store'] * time_step_cycle)

    # The HB Ramp constraints are functions of time, so we need to create some pyomo sets/parameters to represent them.
    n.model.t = pm.Set(initialize=n.snapshots)
    n.model.HB_max_ramp_down = pm.Param(initialize=n.links.loc['HB'].ramp_limit_down)
    n.model.HB_max_ramp_up = pm.Param(initialize=n.links.loc['HB'].ramp_limit_up)

    # Using those sets/parameters, we can now implement the constraints...
    logging.warning('Pypsa has been overridden - Ramp rates on NH3 plant are included')
    n.model.NH3_pyomo_overwrite_ramp_down = pm.Constraint(n.model.t, rule=_nh3_ramp_down)
    n.model.NH3_pyomo_overwrite_ramp_up = pm.Constraint(n.model.t, rule=_nh3_ramp_up)
    # n.model.NH3_pyomo_penalise_ramp_down = pm.Constraint(n.model.t, rule=_penalise_ramp_down)
    # n.model.NH3_pyomo_penalise_ramp_up = pm.Constraint(n.model.t, rule=_penalise_ramp_up)

def _nh3_ramp_down(model, t):
    """Places a cap on how quickly the ammonia plant can ramp down"""
    # if t == 0:
    timestep = 3 # -- is this the same as freq?
    if t == model.t.at(1):

        old_rate = model.link_p['HB', model.t.at(-1)]
    else:
        # old_rate = model.link_p['HB', t - 1]
        old_rate = model.link_p['HB', t - pd.Timedelta(timestep, unit = 'H')]

    return old_rate - model.link_p['HB', t] <= \
        model.link_p_nom['HB'] * model.HB_max_ramp_down
    # Note 20 is the UB of the size of the ammonia plant; essentially if x = 0 then the constraint is not active


def _nh3_ramp_up(model, t):
    """Places a cap on how quickly the ammonia plant can ramp down"""
    # if t == 0:
    timestep = 3 # -- is this the same as freq?
    if t == model.t.at(1):
        old_rate = model.link_p['HB', model.t.at(-1)]
    else:
        # old_rate = model.link_p['HB', t - 1]
        old_rate = model.link_p['HB', t - pd.Timedelta(timestep, unit = 'H')]


    return model.link_p['HB', t] - old_rate <= \
        model.link_p_nom['HB'] * model.HB_max_ramp_up()