import numpy as np
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

        if self.type == "Hydrogen":
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
        elif self.type == "Ammonia":
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
            self.n.generators_t.p_max_pu[f'{gen}'] = gen_list[0]

            # specify maximum capacity based on land use
            self.n.generators.loc[f'{gen}','p_nom_max'] = gen_list[1] # same as above

            # specify technology-specific and country-specific WACC and lifetime here
            self.n.generators.loc[f'{gen}','capital_cost'] = self.n.generators.loc[f'{gen}','capital_cost']\
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