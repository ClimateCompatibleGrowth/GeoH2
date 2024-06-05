
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

    CRF = (((1+interest)**lifetime)*interest)/(((1+interest)**lifetime)-1)
    return CRF

def trucking_costs(transport_state, distance, quantity, interest, transport_excel_path):
    '''
    calculates the annual cost of transporting hydrogen by truck.

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
    excel_path : string
        path to transport_parameters.xlsx file
        
    Returns
    -------
    annual_costs : float
        annual cost of hydrogen transport with specified method.
    '''
    daily_quantity = quantity/365

    transport_parameters = pd.read_excel(transport_excel_path,
                                         sheet_name = transport_state,
                                         index_col = 'Parameter'
                                         ).squeeze('columns')

    average_truck_speed = transport_parameters['Average truck speed (km/h)']
    working_hours = transport_parameters['Working hours (h/day)']
    diesel_price = transport_parameters['Diesel price (euros/L)']
    costs_for_driver = transport_parameters['Costs for driver (euros/h)']
    working_days = transport_parameters['Working days (per year)']
    max_driving_dist = transport_parameters['Max driving distance (km/a)']

    spec_capex_truck = transport_parameters['Spec capex truck (euros)']
    spec_opex_truck = transport_parameters['Spec opex truck (% of capex/a)']
    diesel_consumption = transport_parameters['Diesel consumption (L/100 km)']
    truck_lifetime = transport_parameters['Truck lifetime (a)']

    spec_capex_trailor = transport_parameters['Spec capex trailer (euros)']
    spec_opex_trailor =transport_parameters['Spec opex trailer (% of capex/a)']
    net_capacity = transport_parameters['Net capacity (kg H2)']
    trailor_lifetime = transport_parameters['Trailer lifetime (a)']
    loading_unloading_time = transport_parameters['Loading unloading time (h)']


    amount_deliveries_needed = daily_quantity/net_capacity
    deliveries_per_truck = working_hours/(loading_unloading_time+(2*distance/average_truck_speed))
    trailors_needed = round((amount_deliveries_needed/deliveries_per_truck)+0.5,0)
    total_drives_day = round(amount_deliveries_needed+0.5,0) # not in ammonia calculation
    if transport_state == 'NH3':
        trucks_needed = trailors_needed
    else:
        trucks_needed = max(round((total_drives_day*2*distance*working_days/max_driving_dist)+0.5,0),trailors_needed)

    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    if amount_deliveries_needed < 1:
        fuel_costs = (amount_deliveries_needed*2*distance*365/100)*diesel_consumption*diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver
    
    else:
        fuel_costs = (round(amount_deliveries_needed+0.5)*2*distance*365/100)*diesel_consumption*diesel_price
        wages = round(amount_deliveries_needed+0.5) * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver

    annual_costs = (capex_trucks*CRF(interest,truck_lifetime)+capex_trailor*CRF(interest,trailor_lifetime))\
        + capex_trucks*spec_opex_truck + capex_trailor*spec_opex_trailor + fuel_costs + wages
    return annual_costs


def h2_conversion_stand(final_state, quantity, electricity_costs, heat_costs, interest,
                        conversion_excel_path):
    '''
    calculates the annual cost and electricity and heating demand for converting 
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
    conversion_excel_path: string
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
    
    if final_state != 'standard condition':
        conversion_parameters = pd.read_excel(conversion_excel_path,
                                             sheet_name = final_state,
                                             index_col = 'Parameter'
                                             ).squeeze('columns')

    if final_state == 'standard condition':
        elec_demand = 0 
        heat_demand = 0
        annual_costs = 0 
        return elec_demand, heat_demand, annual_costs

    elif final_state == '500 bar':
        cp = conversion_parameters['Heat capacity']
        Tein = conversion_parameters['Input temperature (K)']
        pein = conversion_parameters['Input pressure (bar)']
        k = conversion_parameters['Isentropic exponent']
        n_isentrop = conversion_parameters['Isentropic efficiency']
                    
        compressor_lifetime = conversion_parameters['Compressor lifetime (a)']
        capex_coefficient = conversion_parameters['Compressor capex coefficient (euros per kilograms H2 per day)']
        opex_compressor = conversion_parameters['Compressor opex (% capex)']

        elec_demand_per_kg_h2 = (cp*Tein*(((500/pein)**((k-1)/k))-1))/n_isentrop
        elec_demand = elec_demand_per_kg_h2 * quantity
        heat_demand = 0 

        capex_compressor = capex_coefficient * ((daily_throughput)**0.6038)

        annual_costs = (capex_compressor*CRF(interest,compressor_lifetime))\
            + (capex_compressor*opex_compressor)\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LH2':

        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        capex_quadratic_coefficient = conversion_parameters['Capex quadratic coefficient (euros (kg H2)-2)']
        capex_linear_coefficient = conversion_parameters['Capex linear coefficient (euros per kg H2)']
        capex_constant = conversion_parameters['Capex constant (euros)']
        opex_liquid_plant = conversion_parameters['Opex (% of capex)']
        liquid_plant_lifetime = conversion_parameters['Plant lifetime (a)']
        
        heat_demand = 0
        elec_demand = electricity_unit_demand * quantity
        capex_liquid_plant = capex_quadratic_coefficient *(daily_throughput**2)\
            +capex_linear_coefficient*daily_throughput\
                +capex_constant

        annual_costs = (capex_liquid_plant*CRF(interest,liquid_plant_lifetime))\
            + (capex_liquid_plant*opex_liquid_plant)\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs
        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LOHC_load':
 
        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per kilograms H2 per year)']
        opex_hydrogenation = conversion_parameters['Opex (% of capex)']
        hydrogenation_lifetime = conversion_parameters['Hydrogenation lifetime (a)']
        costs_carrier = conversion_parameters['Carrier costs (euros per kg carrier)']
        ratio_carrier = conversion_parameters['Carrier ratio (kg carrier: kg hydrogen)']
        
        elec_demand = electricity_unit_demand * quantity 
        heat_demand = heat_unit_demand * quantity              
        capex_hydrogenation = capex_coefficient * quantity

        # why are daily carrier costs included in net present value calculation?
        annual_costs = (capex_hydrogenation+costs_carrier*ratio_carrier*daily_throughput)*CRF(interest, hydrogenation_lifetime)\
            + capex_hydrogenation*opex_hydrogenation\
                + elec_demand * electricity_costs \
                    + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LOHC_unload':

        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per kilograms H2 per year)']
        opex_dehydrogenation = conversion_parameters['Opex (% of capex)']
        dehydrogenation_lifetime = conversion_parameters['Hydrogenation lifetime (a)']
        
        elec_demand = electricity_unit_demand * quantity 
        heat_demand = heat_unit_demand * quantity
        capex_dehydrogenation = capex_coefficient * quantity
        
        annual_costs = (capex_dehydrogenation*CRF(interest, dehydrogenation_lifetime))\
            + (capex_dehydrogenation*opex_dehydrogenation)\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'NH3_load':

        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per annual g H2)']
        opex_NH3_plant = conversion_parameters['Opex (% of capex)']
        NH3_plant_lifetime =conversion_parameters['Plant lifetime (a)']
        
        
        elec_demand = electricity_unit_demand * quantity
        heat_demand = heat_unit_demand * quantity
        capex_NH3_plant = capex_coefficient * quantity

        annual_costs = capex_NH3_plant*CRF(interest,NH3_plant_lifetime)\
            + capex_NH3_plant*opex_NH3_plant\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs
            
        return elec_demand, heat_demand, annual_costs
    
    elif final_state == 'NH3_unload':

        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per hourly g H2)']
        opex_NH3_plant = conversion_parameters['Opex (% of capex)']
        NH3_plant_lifetime = conversion_parameters['Plant lifetime (a)']
        
        elec_demand = electricity_unit_demand * quantity
        heat_demand = heat_unit_demand * quantity

        capex_NH3_plant = capex_coefficient * ((quantity/1000/365/24) ** 0.7451)    

        annual_costs = capex_NH3_plant*CRF(interest,NH3_plant_lifetime) + capex_NH3_plant*opex_NH3_plant \
            + elec_demand * electricity_costs + heat_demand*heat_costs
            
        return elec_demand, heat_demand, annual_costs

    else:
        raise NotImplementedError(f'Conversion costs for {final_state} not currently supported.')

def cheapest_trucking_strategy(final_state, quantity, distance, 
                                elec_costs, heat_costs, interest,
                                conversion_excel_path, transport_excel_path,
                                elec_costs_demand, elec_cost_grid = 0.):
    '''
    calculates the lowest-cost state to transport hydrogen by truck

    Parameters
    ----------
    final_state : string
        final state for hydrogen demand.
    quantity : float
        annual demand for hydrogen in kg.
    distance : float
        distance to transport hydrogen.
    elec_costs : float
        cost per kWh of electricity at hydrogen production site.
    heat_costs : float
        cost per kWh of heat.
    interest : float
        interest on conversion and trucking capital investments (not including roads).
    conversion_excel_path: string
        path to conversion parameters excel sheet.
    elec_costs_demand : float
        cost per kWh of electricity at hydrogen demand site.
    elec_cost_grid : float
        grid electricity costs that pipeline compressors pay. Default 0.
    
    Returns
    -------
    costs_per_unit : float
        storage, conversion, and transport costs for the cheapest trucking option.
    cheapest_option : string
        the lowest-cost state in which to transport hydrogen by truck.

    '''

    if final_state == '500 bar':
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)
    elif final_state == 'NH3':
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]
    else:  
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]
    if final_state == 'LH2':
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('LH2',distance, quantity,interest,transport_excel_path) 
    elif final_state == 'NH3':
        dist_costs_lh2 = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]
    else:
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('LH2',distance, quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]
    if final_state == 'NH3':
        dist_costs_nh3 = h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('NH3',distance, quantity, interest,transport_excel_path)
        dist_costs_lohc = h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('LOHC',distance, quantity, interest,transport_excel_path)\
                    + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]\
                        + h2_conversion_stand('NH3_load', quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]
    else:
        dist_costs_nh3 = h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('NH3',distance, quantity,interest,transport_excel_path)\
                    + h2_conversion_stand('NH3_unload', quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]\
                        + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]
        dist_costs_lohc = h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest, conversion_excel_path)[2]\
                + trucking_costs('LOHC',distance, quantity,interest,transport_excel_path)\
                    + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]\
                        + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]

    lowest_cost = np.nanmin([dist_costs_500bar, dist_costs_lh2, dist_costs_lohc, dist_costs_nh3])
    
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
                                conversion_excel_path,
                                pipeline_excel_path,
                                elec_costs_demand,
                                elec_cost_grid = 0.):
    '''
    calculates the lowest-cost way to transport hydrogen via pipeline

    Parameters
    ----------
    final_state : string
        final state for hydrogen demand.
    quantity : float
        annual demand for hydrogen in kg.
    distance : float
        distance to transport hydrogen.
    elec_costs : float
        cost per kWh of electricity at hydrogen production site.
    heat_costs : float
        cost per kWh of heat.
    interest : float
        interest on pipeline capital investments.
    conversion_excel_path: string
        path to conversion parameters excel sheet.
    elec_costs_demand : float
        cost per kWh of electricity at hydrogen demand site.
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
        dist_costs_pipeline = pipeline_costs(distance,quantity,elec_cost_grid, pipeline_excel_path, interest)[0]\
                + h2_conversion_stand(final_state+'_load', quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]  
    else:
        dist_costs_pipeline = pipeline_costs(distance,quantity,elec_cost_grid,pipeline_excel_path,interest)[0]\
                + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest, conversion_excel_path)[2]

    costs_per_unit = dist_costs_pipeline/quantity
    cheapest_option = pipeline_costs(distance, quantity, elec_cost_grid, pipeline_excel_path, interest)[1] 

    return costs_per_unit, cheapest_option


#Only new pipelines
def pipeline_costs(distance, quantity, elec_cost, pipeline_excel_path, interest):
    '''
    calculates the annualized cost of building a pipeline.

    Parameters
    ----------
    distance : float
        distance from production site to demand site in km.
    quantity : float
        annual quantity of hydrogen demanded in kg.
    elec_cost : float
        price of electricity along pipeline in euros.
    pipeline_excel_path: string
        path to conversion parameters excel sheet.
    interest : float
        interest rate on capital investments.

    Returns
    -------
    float
        annual costs for pipeline.
    string
        size of pipeline to build

    '''
    all_parameters = pd.read_excel(pipeline_excel_path,
                                   sheet_name='All',
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    opex = all_parameters['Opex (% of capex)']
    availability = all_parameters['Availability']
    lifetime_pipeline = all_parameters['Pipeline lifetime (a)']
    lifetime_compressors = all_parameters['Compressor lifetime (a)']
    electricity_demand = all_parameters['Electricity demand (kWh/kg*km)']
    max_capacity_big = all_parameters['Large pipeline max capacity (GW)']
    max_capacity_med = all_parameters['Medium pipeline max capacity (GW)']
    max_capacity_sml = all_parameters['Small pipeline max capcity (GW)']

    max_throughput_big = (((max_capacity_big*(10**6))/33.333))*8760*availability
    max_throughput_med = (((max_capacity_med*(10**6))/33.333))*8760*availability
    max_throughput_sml = (((max_capacity_sml*(10**6))/33.333))*8760*availability

    if quantity <= max_throughput_sml:
        pipeline_type = 'Small'
        
    elif quantity > max_throughput_sml and quantity <= max_throughput_med:
        pipeline_type = 'Medium'
    
    elif quantity > max_throughput_med and quantity <= max_throughput_big:
        pipeline_type = 'Large'

    else:
        return np.nan,'No Pipeline big enough'
    
    pipeline_parameters = pd.read_excel(pipeline_excel_path,
                                   sheet_name=pipeline_type,
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    capex_pipeline = pipeline_parameters['Pipeline capex (euros)']
    capex_compressor = pipeline_parameters['Compressor capex (euros)']
    
    capex_annual = ((capex_pipeline*distance)*CRF(interest,lifetime_pipeline))\
        + ((capex_compressor*distance)*CRF(interest,lifetime_compressors))
    opex_annual = opex*(capex_pipeline+capex_compressor)*distance
    electricity_costs = electricity_demand * distance * quantity * elec_cost

    annual_costs = capex_annual + opex_annual + electricity_costs

    return annual_costs, f"{pipeline_type} Pipeline"
