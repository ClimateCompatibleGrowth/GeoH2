
# from gurobipy import *
# import gurobipy as gp
# from gurobipy import GRB
import pandas as pd
# import random
import numpy as np
# import xlsxwriter as xw

# !!! put this conversion in an excel sheet or use USD for all calculations
usd_to_euro = 0.95

# RBF function

def RBF(interest,lifetime):
    '''
    Calculates the present value factor of a capital investment.

    Parameters
    ----------
    interest : float
        interest rate.
    lifetime : float or integer
        lifetime of asset.

    Returns
    -------
    rbf : float
        present value factor.

    '''
    interest = float(interest)
    lifetime = float(lifetime)

    rbf = (((1+interest)**lifetime)-1)/(((1+interest)**lifetime)*interest)
    return rbf


# Minimum function
# replace this with numpy nanmin

# def cheapest(option1,option2):
#     if option1<option2:
#         return option1
#     elif option1 == 'Nan':
#         return option2
#     elif option2 == 'Nan':
#         return option1
#     else:
#         return option2

#Transportation annual costs

###!!! experimental transport cost function for all transport states
transport_excel_path = "Data/transport_parameters.xlsx"

def trucking_costs(transport_state, distance, quantity, interest, excel_path):
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

    transport_parameters = pd.read_excel(excel_path,
                                         sheet_name = transport_state,
                                         index_col = 'Parameter'
                                         ).squeeze('columns')

    average_truck_speed = transport_parameters['Average truck speed (km/h)']                #km/h
    working_hours = transport_parameters['Working hours (h/day)']                     #h/day
    diesel_price = transport_parameters['Diesel price (euros/L)']                    #€/l
    costs_for_driver = transport_parameters['Costs for driver (euros/h)']                  #€/h
    working_days = transport_parameters['Working days (per year)']                      #per year
    max_driving_dist = transport_parameters['Max driving distance (km/a)']               #km/a Maximum driving distance per truck per year

    spec_capex_truck = transport_parameters['Spec capex truck (euros)']               #€
    spec_opex_truck = transport_parameters['Spec opex truck (% of capex/a)']                  #% of CAPEX/a
    diesel_consumption = transport_parameters['Diesel consumption (L/100 km)']                 #l/100km
    truck_lifetime = transport_parameters['Truck lifetime (a)']                      #a

    spec_capex_trailor = transport_parameters['Spec capex trailer (euros)']
    spec_opex_trailor =transport_parameters['Spec opex trailer (% of capex/a)']
    net_capacity = transport_parameters['Net capacity (kg H2)']                     #kgh2
    trailor_lifetime = transport_parameters['Trailer lifetime (a)']                   #a
    loading_unloading_time = transport_parameters['Loading unloading time (h)']            #hours 


    # max_day_dist = max_driving_dist/working_days
    amount_deliveries_needed = daily_quantity/net_capacity
    deliveries_per_truck = working_hours/(loading_unloading_time+(2*distance/average_truck_speed))
    trailors_needed = round((amount_deliveries_needed/deliveries_per_truck)+0.5,0)
    total_drives_day = round(amount_deliveries_needed+0.5,0) # not in ammonia calculation
    if transport_state == 'NH3': #!!! double checking if this is needed with Leander
        trucks_needed = trailors_needed
    else:
        trucks_needed = max(round((total_drives_day*2*distance*working_days/max_driving_dist)+0.5,0),trailors_needed)

    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    # capex_total = capex_trailor + capex_trucks
    # this if statement seems suspect to me-- how can a fractional number of deliveries be completed?
    if amount_deliveries_needed < 1:
        fuel_costs = (amount_deliveries_needed*2*distance*365/100)*diesel_consumption*diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver
    
    else:
        fuel_costs = (round(amount_deliveries_needed+0.5)*2*distance*365/100)*diesel_consumption*diesel_price
        wages = round(amount_deliveries_needed+0.5) * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver

    annual_costs = (capex_trucks/RBF(interest,truck_lifetime)+capex_trailor/RBF(interest,trailor_lifetime))\
        + capex_trucks*spec_opex_truck + capex_trailor*spec_opex_trailor + fuel_costs + wages
    return annual_costs

#Transformation of hydrogen: Standard condition is pressure = 25 bar and T = 298.15 K
#There are four options: Standard condition, 500 bar, LH2, LOHC and NH3
#quantity per year in kg
#energy demand yearly
#default opex 2% of capex




conversion_excel_path = "Data/conversion_parameters.xlsx"

def h2_conversion_stand(final_state, quantity, electricity_costs, heat_costs, interest):
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
        
        #specific investment costs derived from: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/760479/H2_supply_chain_evidence_-_publication_version.pdf
            
        compressor_lifetime = conversion_parameters['Compressor lifetime (a)']
        capex_coefficient = conversion_parameters['Compressor capex coefficient (euros per kilograms H2 per day)']
        opex_compressor = conversion_parameters['Compressor opex (% capex)']
        
        # compressor_lifetime = 15
        # capex_coefficient = 40035
        # opex_compressor = 0.08                                                  #Annahme 8% von Capex jährlich Assuming 8% of Capex annually
        
        # cp = 14200*((2.777778)*(10**-7))                                        #Wärmekapazität von h2 = 14200 J/kgK Heat capacity of h2 = 14200 J/kgK
        # Tein = 298.15
        # pein = 25
        # k = 1.402                                                               # Isentropenexponent isentropic exponent
        # n_isentrop = 0.8                                                        #Isentroper Wirkungsgrad Isentropic efficiency
        
        elec_demand_per_kg_h2 = (cp*Tein*(((500/pein)**((k-1)/k))-1))/n_isentrop          #kWh/kgh2
        elec_demand = elec_demand_per_kg_h2 * quantity
        heat_demand = 0 

        capex_compressor = capex_coefficient * ((daily_throughput)**0.6038)

        annual_costs = (capex_compressor/RBF(interest,compressor_lifetime))\
            + (capex_compressor*opex_compressor)\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LH2':
        # electricity_unit_demand = 9.93
        # capex_quadratic_coefficient = -0.0002
        # capex_linear_coefficient = 1781.9
        # capex_constant = 3*(10**7)
        # opex_liquid_plant = 0.06
        # liquid_plant_lifetime = 20
        
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

        annual_costs = (capex_liquid_plant/RBF(interest,liquid_plant_lifetime))\
            + (capex_liquid_plant*opex_liquid_plant)\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs
        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LOHC_load':
        # electricity_unit_demand = 0.35
        # heat_unit_demand = 0
        # capex_coefficient = 0.84
        # opex_hydrogenation = 0.04
        # hydrogenation_lifetime = 25
        # costs_carrier = 2                                                       #€/kg_carrier       https://www.hydrogenious.net/index.php/en/2020/07/21/lohc-global-hydrogen-opportunity/#:~:text=carbon%20hydrogen%20transportation.-,Source%3A%20Hydrogenious%20LOHC%20Technologies.,hydrogen%20cost%20of%20%242.54%2Fkg
        # ratio_carrier = 16.1                                                    #kg_carrier:kg_h2   https://www.tvt.tf.fau.eu/files/2019/03/lohc-lkw_bericht_final_teil_1.pdf
        
        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per kilograms H2 per year)']
        opex_hydrogenation = conversion_parameters['Opex (% of capex)']
        hydrogenation_lifetime = conversion_parameters['Hydrogenation lifetime (a)']
        costs_carrier = conversion_parameters['Carrier costs (euros per kg carrier)']                                                      #€/kg_carrier       https://www.hydrogenious.net/index.php/en/2020/07/21/lohc-global-hydrogen-opportunity/#:~:text=carbon%20hydrogen%20transportation.-,Source%3A%20Hydrogenious%20LOHC%20Technologies.,hydrogen%20cost%20of%20%242.54%2Fkg
        ratio_carrier = conversion_parameters['Carrier ratio (kg carrier: kg hydrogen)']
        
        elec_demand = electricity_unit_demand * quantity 
        heat_demand = heat_unit_demand * quantity                               #https://www.hydrogenious.net/index.php/en/hydrogenious-3/lohc-technology/                             
        capex_hydrogenation = capex_coefficient * quantity

        # why are daily carrier costs included in net present value calculation?
        annual_costs = (capex_hydrogenation+costs_carrier*ratio_carrier*daily_throughput)/RBF(interest, hydrogenation_lifetime)\
            + capex_hydrogenation*opex_hydrogenation\
                + elec_demand * electricity_costs \
                    + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LOHC_unload':
        # electricity_unit_demand = 0.35
        # heat_unit_demand = 12
        # capex_coefficient = 2.46
        # opex_dehydrogenation = 0.04
        # dehydrogenation_lifetime = 25
        
        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per kilograms H2 per year)']
        opex_dehydrogenation = conversion_parameters['Opex (% of capex)']
        dehydrogenation_lifetime = conversion_parameters['Hydrogenation lifetime (a)']
        
        elec_demand = electricity_unit_demand * quantity 
        heat_demand = heat_unit_demand * quantity
        capex_dehydrogenation = capex_coefficient * quantity
        
        annual_costs = (capex_dehydrogenation/RBF(interest, dehydrogenation_lifetime))\
            + (capex_dehydrogenation*opex_dehydrogenation)\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'NH3_load':
        # electricity_unit_demand = 3.109
        # heat_unit_demand = 0
        # capex_coefficient = 4.25379
        # opex_NH3_plant = 0.015
        # NH3_plant_lifetime = 25
        
        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per annual g H2)']
        opex_NH3_plant = conversion_parameters['Opex (% of capex)']
        NH3_plant_lifetime =conversion_parameters['Plant lifetime (a)']
        
        
        elec_demand = electricity_unit_demand * quantity                                         #see excel
        heat_demand = heat_unit_demand * quantity
        capex_NH3_plant = capex_coefficient * quantity

        annual_costs = capex_NH3_plant/RBF(interest,NH3_plant_lifetime)\
            + capex_NH3_plant*opex_NH3_plant\
                + elec_demand * electricity_costs\
                    + heat_demand*heat_costs
            
        return elec_demand, heat_demand, annual_costs
    
    elif final_state == 'NH3_unload':
        # electricity_unit_demand = 4.2
        # heat_unit_demand = 0
        # capex_coefficient = 17262450
        # opex_NH3_plant = 0.02
        # NH3_plant_lifetime = 25
        
        electricity_unit_demand = conversion_parameters['Electricity demand (kWh per kg H2)']
        heat_unit_demand = conversion_parameters['Heat demand (kWh per kg H2)']
        capex_coefficient = conversion_parameters['Capex coefficient (euros per hourly g H2)']
        opex_NH3_plant = conversion_parameters['Opex (% of capex)']
        NH3_plant_lifetime = conversion_parameters['Plant lifetime (a)']
        
        elec_demand = electricity_unit_demand * quantity                                         #see excel
        heat_demand = heat_unit_demand * quantity

        capex_NH3_plant = capex_coefficient * ((quantity/1000/365/24) ** 0.7451)    

        annual_costs = capex_NH3_plant/RBF(interest,NH3_plant_lifetime) + capex_NH3_plant*opex_NH3_plant \
            + elec_demand * electricity_costs + heat_demand*heat_costs
            
        return elec_demand, heat_demand, annual_costs

    else:
        print('Conversion costs for {} not currently supported.'.format(final_state))

def cheapest_transport_strategy(final_state, quantity, distance, 
                                elec_costs, heat_costs,interest, 
                                elec_costs_demand, days_storage,
                                elec_cost_grid = 0., pipeline = True):
    '''
    calculates the lowest-cost way to transport hydrogen, either in different states
    by truck or via pipeline if allowed

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
        interest on capital investments.
    elec_costs_demand : float
        cost per kWh of electricity at hydrogen demand site.
    days_storage : float
        number of days of storage to build at production and demand sites.
    elec_cost_grid : float
        grid electricity costs that pipeline compressors pay. Default 0.
    pipeline : boolean
        If True, building a pipeline is an option. Default True.
        
    Returns
    -------
    costs_per_unit : float
        storage, conversion, and transport costs for the cheapest option.
    cheapest_option : string
        the lowest-cost state in which to transport hydrogen by truck.

    '''

    storage_costs_500bar = storage_costs('500 bar',quantity,days_storage,interest)\
        + storage_costs(final_state,quantity,days_storage,interest)
    storage_costs_lohc = storage_costs('LOHC',quantity,days_storage,interest)\
        + storage_costs(final_state,quantity,days_storage,interest)
    storage_costs_lh2 = storage_costs('LH2',quantity,days_storage,interest)\
        + storage_costs(final_state,quantity,days_storage,interest)
    storage_costs_nh3 = storage_costs('NH3',quantity,days_storage,interest)\
        + storage_costs(final_state,quantity,days_storage,interest)
    storage_costs_pipeline = storage_costs(final_state,quantity,days_storage,interest)


    if final_state == '500 bar':
        dist_costs_500bar = storage_costs_500bar\
            + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)
    elif final_state == 'NH3':
        dist_costs_500bar = storage_costs_500bar\
            + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest)[2]
    else:  
        dist_costs_500bar = storage_costs_500bar\
            + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest)[2]
    if final_state == 'LH2':
        dist_costs_lh2 =  storage_costs_lh2\
            + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('LH2',distance, quantity,interest,transport_excel_path) 
    elif final_state == 'NH3':
        dist_costs_lh2 = storage_costs_500bar\
            + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('500 bar',distance,quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest)[2]
    else:
        dist_costs_lh2 =  storage_costs_lh2\
            + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('LH2',distance, quantity,interest,transport_excel_path)\
                    + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest)[2]
    if final_state == 'NH3':
        dist_costs_nh3 = storage_costs_nh3 \
            + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('NH3',distance, quantity, interest,transport_excel_path) 
        dist_costs_lohc = storage_costs_lohc\
            + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('LOHC',distance, quantity, interest,transport_excel_path)\
                    + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interest)[2]\
                        + h2_conversion_stand('NH3_load', quantity, elec_costs_demand, heat_costs, interest)[2]
    else:
        dist_costs_nh3 = storage_costs_nh3\
            + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('NH3',distance, quantity,interest,transport_excel_path)\
                    + h2_conversion_stand('NH3_unload', quantity, elec_costs_demand, heat_costs, interest)[2]\
                        + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest)[2]
        dist_costs_lohc = storage_costs_lohc\
            + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest)[2]\
                + trucking_costs('LOHC',distance, quantity,interest,transport_excel_path)\
                    + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interest)[2]\
                        + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest)[2]
    if pipeline == True:
        if final_state == 'NH3':
            dist_costs_pipeline = storage_costs_pipeline\
                + pipeline_costs(distance,quantity,elec_cost_grid,interest)[0]\
                    + h2_conversion_stand(final_state+'_load', quantity, elec_costs_demand, heat_costs, interest)[2]  
        else:
            dist_costs_pipeline = storage_costs_pipeline\
                + pipeline_costs(distance,quantity,elec_cost_grid,interest)[0]\
                    + h2_conversion_stand(final_state, quantity, elec_costs_demand, heat_costs, interest)[2]
    else:
        dist_costs_pipeline = np.nan

    lowest_cost = np.nanmin([dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline])
    
    if dist_costs_500bar == lowest_cost:
        cheapest_option = '500 bar'
    elif dist_costs_lh2 == lowest_cost:
        cheapest_option = 'LH2'
    elif dist_costs_lohc == lowest_cost:
        cheapest_option = 'LOHC'
    elif dist_costs_nh3 == lowest_cost: 
         cheapest_option = 'NH3'
    else:
         cheapest_option = pipeline_costs(distance,quantity,elec_cost_grid,interest)[1] 
    
    costs_per_unit = lowest_cost/quantity
    
    return costs_per_unit, cheapest_option

#Only new pipelines
pipeline_excel_path = "Data/pipeline_parameters.xlsx"

def pipeline_costs(distance,quantity,elec_cost,interest):
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
    # opex = 0.0125                                                           #% of capex per year
    # availability = 0.95                                                     #Assumption
    # lifetime_pipeline = 42.5
    # lifetime_compressors = 24
    # electricity_demand = 0.000613819                                        #kWh/kg*km

    # capex_pipeline_big = 2.8 * 1000000                                      #Mio. €/km
    # capex_compression_big = 0.62 * 1000000                                  #Mio. €/km
    # max_capacity_big = 13                                                   #GW

    # capex_pipeline_med = 2.2 * 1000000                                       #Mio. €/km
    # capex_compression_med = 0.31 * 1000000                                   #Mio. €/km
    # max_capacity_med = 4.7                                                   #GW

    # capex_pipeline_sml = 1.5 * 1000000                                       #Mio. €/km
    # capex_compression_sml = 0.09 * 1000000                                   #Mio. €/km
    # max_capacity_sml = 1.2                                                   #GW
    
    max_throughput_big = (((max_capacity_big*(10**6))/33.333))*8760*availability         #kg/year   
    max_throughput_med = (((max_capacity_med*(10**6))/33.333))*8760*availability          #kg/year   
    max_throughput_sml = (((max_capacity_sml*(10**6))/33.333))*8760*availability         #kg/year   

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
    
    capex_annual = ((capex_pipeline*distance)/RBF(interest,lifetime_pipeline))\
        + ((capex_compressor*distance)/RBF(interest,lifetime_compressors))
    opex_annual = opex*(capex_pipeline+capex_compressor)*distance
    electricity_costs = electricity_demand * distance * quantity * elec_cost

    annual_costs = capex_annual + opex_annual + electricity_costs

    return annual_costs, f"{pipeline_type} Pipeline"

def storage_costs(state,quantity,storagedays,interest):
    # !!! put these quantities in an excel sheet

    opex = 0.02
    lifetime_storage = 20 

    if state == '500 bar':
        capex_storage = 2160*((quantity*storagedays/365)**-0.146)
    
    elif state == 'LH2':
        capex_storage = 37.24

    elif state == 'LOHC':
        capex_storage = 3.192
    
    elif state == 'NH3':
        capex_storage = 8.512
    

    capex_storage = capex_storage * (quantity*storagedays/365)
    annual_costs = (capex_storage/RBF(interest,lifetime_storage)) + opex * capex_storage

    return annual_costs

    # !!! put these quantities in an excel sheet-- what are they for???

demand = [100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000, 200000000, 500000000, 1000000000]

distance = list(range(100, 1550, 50))

total_list = []

#for i in demand:
#    row_list = []
#    for j in distance:
#        row_list.append(str(cheapest_dist_option_pipeline('500 bar', i,j,0.04,0.08,0.08,0.04,0.04)[0]))
    
#    total_list.append(row_list)

#print(total_list)

array = np.array(total_list)

storage = ['500 bar', 'LH2', 'LOHC', 'NH3']

#for i in storage:
#    print(storage_costs(i,10000000,3,0.08)/(10000000))


#print(array)

#workbook = xw.Workbook('arrays.xlsx')
#worksheet = workbook.add_worksheet()


#row = 0

#for col, data in enumerate(array):
#    worksheet.write_column(row, col, data)

#workbook.close()


#qu = 1000000

#print((h2_conversion_stand('LOHC_load', qu, 0.04, 0.08, 0.08)[2] + h2_conversion_stand('LOHC_unload', qu, 0.04, 0.08, 0.08)[2] + transport_lohc(700,qu,0.08) + h2_conversion_stand('500 bar', qu, 0.04, 0.03, 0.08)[2])/qu)

#print(( h2_conversion_stand('500 bar', qu, 0.04, 0.08, 0.08)[2] + transport_500bar(700,qu,0.08))/qu)


#transport_500bar(100,10000,0.08)
#print(transport_500bar(200,10000000,0.08)/10000000 )
#print(distance)
