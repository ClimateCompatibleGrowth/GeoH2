import pandas as pd
import numpy as np

def CRF(interest,lifetime):
    '''
    Calculates the capital recovery factor of a capital investment.

    Parameters
    ----------
    interest : float
        interest rate.
    lifetime : float or integer
        lifetime of asset.

    Returns
    -------
    CRF : float
        present value factor.
    '''
    interest = float(interest)
    lifetime = float(lifetime)

    CRF = (((1 + interest)**lifetime) * interest)/(((1 + interest)**lifetime) - 1)
    
    return CRF

def trucking_costs(transport_state, distance, quantity, interest, transport_params_filepath):
    '''
    Calculates the annual cost of transporting hydrogen by truck.

    Parameters
    ----------
    transport_state : string
        state hydrogen is transported in, one of '500 bar', 'LH2', 'LOHC', or 'NH3'.
    distance : float
        distance between hydrogen production site and demand site.
    quantity : float
        annual amount of hydrogen to transport.
    interest : float
        interest rate on capital investments.
    transport_params_filepath : string
        path to transport_parameters.xlsx file
        
    Returns
    -------
    annual_costs : float
        annual cost of hydrogen transport with specified method.
    '''
    daily_quantity = quantity/365

    transport_params = pd.read_excel(transport_params_filepath,
                                         sheet_name = transport_state,
                                         index_col = 'Parameter'
                                         ).squeeze('columns')

    average_truck_speed = transport_params['Average truck speed (km/h)']
    working_hours = transport_params['Working hours (h/day)']
    diesel_price = transport_params['Diesel price (euros/L)']
    costs_for_driver = transport_params['Costs for driver (euros/h)']
    working_days = transport_params['Working days (per year)']
    max_driving_dist = transport_params['Max driving distance (km/a)']

    spec_capex_truck = transport_params['Spec capex truck (euros)']
    spec_opex_truck = transport_params['Spec opex truck (% of capex/a)']
    diesel_consumption = transport_params['Diesel consumption (L/100 km)']
    truck_lifetime = transport_params['Truck lifetime (a)']

    spec_capex_trailor = transport_params['Spec capex trailer (euros)']
    spec_opex_trailor =transport_params['Spec opex trailer (% of capex/a)']
    net_capacity = transport_params['Net capacity (kg H2)']
    trailor_lifetime = transport_params['Trailer lifetime (a)']
    loading_unloading_time = transport_params['Loading unloading time (h)']

    # Calculate deliveries needed per day
    amount_deliveries_needed = daily_quantity/net_capacity
    # Calculate how many deliveries each truck can do per day
    deliveries_per_truck = working_hours/(loading_unloading_time +
                                          (2 * distance/average_truck_speed))
    # Deliveries per day / Deliveries per truck = Trucks per day
    # In the lines below, the 0.5 is put in place in order to round up so full demand is met
    trailors_needed = round((amount_deliveries_needed/
                             deliveries_per_truck) + 0.5)
    total_drives_day = round(amount_deliveries_needed + 0.5) # not in ammonia calculation
    if transport_state == 'NH3':
        trucks_needed = trailors_needed
    else:
        trucks_needed = max(round((total_drives_day * 2 * distance *
                                        working_days/max_driving_dist) + 0.5),
                                            trailors_needed)
    # Get the capex of all the trucks and trailors needed
    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    # Get fuel costs and wages
    if amount_deliveries_needed < 1:
        # In the lines below, 365 refers to  days to spread it over the year
        # and 100 is there because diesel_consumption is in liters/100km
        fuel_costs = (amount_deliveries_needed * 2 *
                        distance * 365/100) * diesel_consumption * diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed) *
                    2 + loading_unloading_time) * working_days * costs_for_driver
    else:
        fuel_costs = (round(amount_deliveries_needed + 0.5) *
                        2 * distance * 365/100) * diesel_consumption * diesel_price
        wages = round(amount_deliveries_needed + 0.5) * (
                    (distance/average_truck_speed) * 2 + loading_unloading_time
                    ) * working_days * costs_for_driver
    # Get total annual costs including capex, fuel costs, wages
    annual_costs = (capex_trucks * CRF(interest, truck_lifetime) + capex_trailor *
                        CRF(interest, trailor_lifetime)) +\
                            capex_trucks * spec_opex_truck + capex_trailor *\
                                spec_opex_trailor + fuel_costs + wages
    
    return annual_costs


def h2_conversion_stand(final_state, quantity, electricity_costs, heat_costs, interest,
                        conversion_params_filepath):
    # Leader to go through and add comments to further explain the uses of load and unload and standard condition
    '''
    Calculates the annual cost and electricity and heating demand for converting 
    hydrogen to a given state

    Parameters
    ----------
    final_state : string
        final state to convert hydrogen to, one of 'standard condition', '500 bar',
        'LH2', 'LOHC_load', 'LOHC_unload', 'NH3_load', or 'NH3_unload'.
    quantity : float
        annual quantity of hydrogen to convert in kg.
    electricity_costs : float
        unit price for electricity.
    heat_costs : float
        unit costs for heat.
    interest : float
        interest rate applicable to hydrogen converter investments.
    conversion_params_filepath : string
        path to conversion parameters excel sheet.

    Returns
    -------
    elec_demand : float
        annual electricity demand.
    heat_demand : float
        annual heat demand.
    annual_costs : float
        annual hydrogen conversion costs.

    '''
    
    daily_throughput = quantity/365
    
    if final_state == 'standard condition':
        elec_demand = 0 
        heat_demand = 0
        annual_costs = 0 
    else:
        conversion_params = pd.read_excel(conversion_params_filepath,
                                             sheet_name = final_state,
                                             index_col = 'Parameter'
                                             ).squeeze('columns')
        
        if final_state == '500 bar':
            cp = conversion_params['Heat capacity']
            tein = conversion_params['Input temperature (K)']
            pein = conversion_params['Input pressure (bar)']
            k = conversion_params['Isentropic exponent']
            n_isentrop = conversion_params['Isentropic efficiency']
                        
            compressor_lifetime = conversion_params['Compressor lifetime (a)']
            capex_coef = conversion_params['Compressor capex coefficient (euros per kilograms H2 per day)']
            opex_compressor = conversion_params['Compressor opex (% capex)']
            # I don't know what the 500 means - we should probably just assign that to a variable that's named
            # I can look up what it is but like for readability
            # Can Leander suggest a comment to add here to reference the below formula/source?
            elec_demand_per_kg_h2 = (cp * tein * (((500/pein)**((k - 1)/k)) - 1))/n_isentrop
            elec_demand = elec_demand_per_kg_h2 * quantity
            heat_demand = 0 
            # I don't know what the 0.6038 is, same as above
            # Can Leander suggest a comment to add here to reference the below formula/source?
            capex_compressor = capex_coef * ((daily_throughput)**0.6038)

            # Leander: please add comment here to explain below
            annual_costs = (capex_compressor * CRF(interest, compressor_lifetime)) +\
                                (capex_compressor * opex_compressor) +\
                                    elec_demand * electricity_costs +\
                                        heat_demand * heat_costs
            
        elif final_state == 'LH2':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            capex_quadratic_coef = conversion_params['Capex quadratic coefficient (euros (kg H2)-2)']
            capex_linear_coef = conversion_params['Capex linear coefficient (euros per kg H2)']
            capex_constant = conversion_params['Capex constant (euros)']
            opex_liquid_plant = conversion_params['Opex (% of capex)']
            liquid_plant_lifetime = conversion_params['Plant lifetime (a)']
            
            heat_demand = 0
            elec_demand = electricity_unit_demand * quantity
            # Maybe we should put some comments in referring to equations documented somewhere ... hmm
            capex_liquid_plant = capex_quadratic_coef * (daily_throughput**2) +\
                                    capex_linear_coef * daily_throughput +\
                                        capex_constant

            # Leander: please add comment here to explain below
            annual_costs = (capex_liquid_plant * CRF(interest, liquid_plant_lifetime)) +\
                                (capex_liquid_plant * opex_liquid_plant) +\
                                    elec_demand * electricity_costs +\
                                        heat_demand * heat_costs
            
        elif final_state == 'LOHC_load':
            # In this conversion you "load" the hydrogen molecule to a carrier liquid.
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coef = conversion_params['Capex coefficient (euros per kilograms H2 per year)']
            opex_hydrogenation = conversion_params['Opex (% of capex)']
            hydrogenation_lifetime = conversion_params['Hydrogenation lifetime (a)']
            costs_carrier = conversion_params['Carrier costs (euros per kg carrier)']
            ratio_carrier = conversion_params['Carrier ratio (kg carrier: kg hydrogen)']
            
            elec_demand = electricity_unit_demand * quantity 
            heat_demand = heat_unit_demand * quantity              
            capex_hydrogenation = capex_coef * quantity

            # why are daily carrier costs included in net present value calculation?
            # Leander: please add comment here to explain below
            annual_costs = (capex_hydrogenation + costs_carrier *
                                ratio_carrier * daily_throughput) *\
                                    CRF(interest, hydrogenation_lifetime) +\
                                        capex_hydrogenation * opex_hydrogenation +\
                                            elec_demand * electricity_costs +\
                                                heat_demand * heat_costs
            
        elif final_state == 'LOHC_unload':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coef = conversion_params['Capex coefficient (euros per kilograms H2 per year)']
            opex_dehydrogenation = conversion_params['Opex (% of capex)']
            dehydrogenation_lifetime = conversion_params['Hydrogenation lifetime (a)']
            
            elec_demand = electricity_unit_demand * quantity 
            heat_demand = heat_unit_demand * quantity
            capex_dehydrogenation = capex_coef * quantity
            
            # Leander: please add comment here to explain below
            annual_costs = (capex_dehydrogenation *
                                CRF(interest, dehydrogenation_lifetime)) +\
                                    (capex_dehydrogenation * opex_dehydrogenation) +\
                                        elec_demand * electricity_costs +\
                                            heat_demand * heat_costs
            
        elif final_state == 'NH3_load':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coefficient = conversion_params['Capex coefficient (euros per annual g H2)']
            opex_NH3_plant = conversion_params['Opex (% of capex)']
            NH3_plant_lifetime = conversion_params['Plant lifetime (a)']
            
            
            elec_demand = electricity_unit_demand * quantity
            heat_demand = heat_unit_demand * quantity
            capex_NH3_plant = capex_coefficient * quantity

            # Leander: please add comment here to explain below
            annual_costs = capex_NH3_plant * CRF(interest, NH3_plant_lifetime) +\
                                capex_NH3_plant * opex_NH3_plant +\
                                    elec_demand * electricity_costs +\
                                        heat_demand * heat_costs
            
        elif final_state == 'NH3_unload':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coefficient = conversion_params['Capex coefficient (euros per hourly g H2)']
            opex_NH3_plant = conversion_params['Opex (% of capex)']
            NH3_plant_lifetime = conversion_params['Plant lifetime (a)']
            
            elec_demand = electricity_unit_demand * quantity
            heat_demand = heat_unit_demand * quantity

            # Again - I have no idea what those factors are... 365 days/year, 24 hours/day, but why 1000 and 0.7451?
            # Can Leander suggest a comment to add here to reference the below formula/source?
            capex_NH3_plant = capex_coefficient * ((quantity/1000/365/24) ** 0.7451)    

            # Leander: please add comment here to explain below
            annual_costs = capex_NH3_plant *\
                                CRF(interest, NH3_plant_lifetime) +\
                                    capex_NH3_plant * opex_NH3_plant +\
                                        elec_demand * electricity_costs +\
                                            heat_demand * heat_costs
        else:
            raise NotImplementedError(f'Conversion costs for {final_state} not currently supported.')
        
    return elec_demand, heat_demand, annual_costs

def cheapest_trucking_strategy(final_state, quantity, distance, 
                                elec_costs, heat_costs, interest,
                                conversion_params_filepath, transport_params_filepath):
    # Leader to go through and add comments
    '''
    Calculates the lowest-cost state to transport hydrogen by truck

    Parameters
    ----------
    final_state : string
        final state for hydrogen demand.
    quantity : float
        annual demand for hydrogen in kg.
    distance : float
        distance to transport hydrogen.
    elec_costs : float
        cost per kWh of electricity for that country.
    heat_costs : float
        cost per kWh of heat.
    interest : float
        interest on conversion and trucking capital investments (not including roads).
    conversion_params_filepath : string
        path to conversion parameters excel sheet.
    transport_params_filepath : string
        path to transport parameters excel sheet. 
    
    Returns
    -------
    costs_per_unit : float
        storage, conversion, and transport costs for the cheapest trucking option.
    cheapest_option : string
        the lowest-cost state in which to transport hydrogen by truck.

    '''

    # I am confused by this entire function to be honest. Struggling to follow the logic
    # Different between _load and _unload is fuzzy
    if final_state == '500 bar':
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('500 bar', distance, quantity, interest, transport_params_filepath)
    elif final_state == 'NH3':
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('500 bar',distance,quantity,interest,transport_params_filepath) +\
                    h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]
    else:  
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('500 bar',distance,quantity,interest,transport_params_filepath) +\
                    h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]
    
    if final_state == 'LH2':
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('LH2',distance, quantity,interest,transport_params_filepath)
    elif final_state == 'NH3':
        # Should these ones be LH2 in first two lines???
        dist_costs_lh2 = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('500 bar',distance,quantity,interest,transport_params_filepath) +\
                    h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]
    else:
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('LH2',distance, quantity,interest,transport_params_filepath) +\
                    h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]
    
    if final_state == 'NH3':
        dist_costs_nh3 = h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('NH3',distance, quantity, interest,transport_params_filepath)
        dist_costs_lohc = h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('LOHC',distance, quantity, interest,transport_params_filepath) +\
                    h2_conversion_stand('LOHC_unload', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                        h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]
    else:
        dist_costs_nh3 = h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('NH3',distance, quantity,interest,transport_params_filepath) +\
                    h2_conversion_stand('NH3_unload', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                        h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]
        dist_costs_lohc = h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                trucking_costs('LOHC',distance, quantity,interest,transport_params_filepath) +\
                    h2_conversion_stand('LOHC_unload', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2] +\
                        h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]

    lowest_cost = np.nanmin([dist_costs_500bar, dist_costs_lh2, dist_costs_lohc, dist_costs_nh3])
    
    # Allocating cheapest option
    if dist_costs_500bar == lowest_cost:
        cheapest_option = '500 bar'
    elif dist_costs_lh2 == lowest_cost:
        cheapest_option = 'LH2'
    elif dist_costs_lohc == lowest_cost:
        cheapest_option = 'LOHC'
    elif dist_costs_nh3 == lowest_cost: 
         cheapest_option = 'NH3'
    
    costs_per_unit = lowest_cost/quantity
    
    return costs_per_unit, cheapest_option

    
    
def cheapest_pipeline_strategy(final_state, quantity, distance, 
                                elec_costs, heat_costs,interest, 
                                conversion_params_filepath,
                                pipeline_params_filepath,
                                elec_cost_grid = 0.):
    # Leader to go through and add comments
    '''
    Calculates the lowest-cost way to transport hydrogen via pipeline

    Parameters
    ----------
    final_state : string
        final state for hydrogen demand.
    quantity : float
        annual demand for hydrogen in kg.
    distance : float
        distance to transport hydrogen.
    elec_costs : float
        cost per kWh of electricity for that country.
    heat_costs : float
        cost per kWh of heat.
    interest : float
        interest on pipeline capital investments.
    conversion_params_filepath: string
        path to conversion parameters excel sheet.
    pipeline_params_filepath : string
        path to pipeline parameters excel sheet.
    elec_cost_grid : float
        grid electricity costs that pipeline compressors pay. Default 0.

    Returns
    -------
    costs_per_unit : float
        storage, conversion, and transport costs for the cheapest option.
    cheapest_option : string
        the lowest-cost state in which to transport hydrogen by truck.

    '''
    if final_state == 'NH3':
        dist_costs_pipeline = pipeline_costs(distance,quantity, elec_cost_grid, pipeline_params_filepath, interest)[0] +\
                h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]  
    else:
        dist_costs_pipeline = pipeline_costs(distance, quantity, elec_cost_grid, pipeline_params_filepath, interest)[0] +\
                h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath)[2]

    costs_per_unit = dist_costs_pipeline/quantity
    cheapest_option = pipeline_costs(distance, quantity, elec_cost_grid, pipeline_params_filepath, interest)[1] 

    return costs_per_unit, cheapest_option


# Only new pipelines
def pipeline_costs(distance, quantity, elec_cost, pipeline_params_filepath, interest):
    '''
    Calculates the annualized cost of building a pipeline.

    Parameters
    ----------
    distance : float
        distance from production site to demand site in km.
    quantity : float
        annual quantity of hydrogen demanded in kg.
    elec_cost : float
        price of electricity along pipeline in euros.
    pipeline_params_filepath: string
        path to conversion parameters excel sheet.
    interest : float
        interest rate on capital investments.

    Returns
    -------
    annual_costs : float
        annual costs for pipeline.
    string
        size of pipeline to build.

    '''
    all_parameters = pd.read_excel(pipeline_params_filepath,
                                   sheet_name='All',
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    opex = all_parameters['Opex (% of capex)']
    availability = all_parameters['Availability']
    lifetime_pipeline = all_parameters['Pipeline lifetime (a)']
    lifetime_compressors = all_parameters['Compressor lifetime (a)']
    electricity_demand = all_parameters['Electricity demand (kWh/kg*km)']
    large_max_capacity = all_parameters['Large pipeline max capacity (GW)']
    med_max_capacity = all_parameters['Medium pipeline max capacity (GW)']
    sml_max_capacity = all_parameters['Small pipeline max capcity (GW)']

    # 33.333 (kWh/kg) is the energy density of hydrogen, 8760 are hours are in the year.
    large_max_throughput = (((large_max_capacity * (10**6))/33.333)) * 8760 * availability
    med_max_throughput = (((med_max_capacity * (10**6))/33.333)) * 8760 * availability
    sml_max_throughput = (((sml_max_capacity * (10**6))/33.333)) * 8760 * availability

    if quantity <= sml_max_throughput:
        pipeline_type = 'Small'
        
    elif quantity > sml_max_throughput and quantity <= med_max_throughput:
        pipeline_type = 'Medium'
    
    elif quantity > med_max_throughput and quantity <= large_max_throughput:
        pipeline_type = 'Large'

    else:
        return np.nan,'No Pipeline big enough'
    
    pipeline_parameters = pd.read_excel(pipeline_params_filepath,
                                   sheet_name=pipeline_type,
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    capex_pipeline = pipeline_parameters['Pipeline capex (euros)']
    capex_compressor = pipeline_parameters['Compressor capex (euros)']
    
    capex_annual = ((capex_pipeline * distance) *
                        CRF(interest, lifetime_pipeline)) +\
                            ((capex_compressor * distance) *\
                                CRF(interest, lifetime_compressors))
    opex_annual = opex * (capex_pipeline + capex_compressor) * distance
    electricity_costs = electricity_demand * distance * quantity * elec_cost

    annual_costs = capex_annual + opex_annual + electricity_costs

    return annual_costs, f"{pipeline_type} Pipeline"