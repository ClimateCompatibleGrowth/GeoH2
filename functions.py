
from gurobipy import *
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import random
import numpy as np
import xlsxwriter as xw

usd_to_euro = 0.95

# RBF function

def RBF(zins,lt):

    zins = float(zins)
    lt = float(lt)

    rbf = (((1+zins)**lt)-1)/(((1+zins)**lt)*zins)
    return rbf


# Minimum function

def cheapest(option1,option2):
    if option1<option2:
        return option1
    elif option1 == 'Nan':
        return option2
    elif option2 == 'Nan':
        return option1
    else:
        return option2

#Transportation annual costs

def transport_500bar(distance,quantity,interest):

    quantity = quantity/365

    average_truck_speed = 70                #km/h
    working_hours = 24                      #h/day
    diesel_price = 1.5                      #€/l
    costs_for_driver = 2.85                   #€/h
    working_days = 365                      #per year
    max_driving_dist = 160000               #km/a Maximum driving distance per truck per year

    spec_capex_truck = 160000               #€
    spec_opex_truck = 0.12                  #% of CAPEX/a
    diesel_consumption = 35                 #l/100km
    truck_lifetime = 8                      #a

    spec_capex_trailor = 660000
    spec_opex_trailor = 0.02
    net_capacity = 1100                     #kgh2
    trailor_lifetime = 12                   #a
    loading_unloading_time = 1.5            #hours 


    max_day_dist = max_driving_dist/working_days
    amount_deliveries_needed = quantity/net_capacity
    deliveries_per_truck = working_hours/(loading_unloading_time+(2*distance/average_truck_speed))
    trailors_needed = round((amount_deliveries_needed/deliveries_per_truck)+0.5,0)
    total_drives_day = round(amount_deliveries_needed+0.5,0)
    trucks_needed = max(round((total_drives_day*2*distance*working_days/max_driving_dist)+0.5,0),trailors_needed)

    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    capex_total = capex_trailor + capex_trucks

    if amount_deliveries_needed < 1:
        fuel_costs = (amount_deliveries_needed*2*distance*365/100)*diesel_consumption*diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver
    
    else:
        fuel_costs = (round(amount_deliveries_needed+0.5)*2*distance*365/100)*diesel_consumption*diesel_price
        wages = round(amount_deliveries_needed+0.5) * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver



    annual_costs = (capex_trucks/RBF(interest,truck_lifetime)+capex_trailor/RBF(interest,trailor_lifetime)) + capex_trucks*spec_opex_truck + capex_trailor*spec_opex_trailor + fuel_costs + wages
    return annual_costs

def transport_lh2(distance,quantity,interest):

    quantity = quantity/365

    average_truck_speed = 70                #km/h
    working_hours = 24                      #h/day
    diesel_price = 1.5                      #€/l
    costs_for_driver = 2.85                 #€/h
    working_days = 365
    max_driving_dist = 160000               #km/a Maximum driving distance per truck per year

    spec_capex_truck = 160000               #€
    spec_opex_truck = 0.12                  #% of CAPEX/a
    diesel_consumption = 35                 #l/100km
    truck_lifetime = 8                      #a

    spec_capex_trailor = 860000
    spec_opex_trailor = 0.02
    net_capacity = 4300                     #kgh2
    trailor_lifetime = 12                   #a
    loading_unloading_time = 3              #hours 


    max_day_dist = max_driving_dist/working_days
    amount_deliveries_needed = quantity/net_capacity
    deliveries_per_truck = working_hours/(loading_unloading_time+(2*distance/average_truck_speed))
    trailors_needed = round((amount_deliveries_needed/deliveries_per_truck)+0.5,0)
    total_drives_day = round(amount_deliveries_needed+0.5,0)
    trucks_needed = max(round((total_drives_day*2*distance*working_days/max_driving_dist)+0.5,0),trailors_needed)

    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    capex_total = capex_trailor + capex_trucks

    if amount_deliveries_needed < 1:
        fuel_costs = (amount_deliveries_needed*2*distance*365/100)*diesel_consumption*diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver
    
    else:
        fuel_costs = (round(amount_deliveries_needed+0.5)*2*distance*365/100)*diesel_consumption*diesel_price
        wages = round(amount_deliveries_needed+0.5) * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver


    annual_costs = (capex_trucks/RBF(interest,truck_lifetime)+capex_trailor/RBF(interest,trailor_lifetime)) + capex_trucks*spec_opex_truck + capex_trailor*spec_opex_trailor + fuel_costs + wages
    return annual_costs

def transport_lohc(distance,quantity,interest):

    quantity = quantity/365

    average_truck_speed = 70                #km/h
    working_hours = 24                      #h/day
    diesel_price = 1.5                      #€/l
    costs_for_driver = 2.85                 #€/h
    working_days = 365
    max_driving_dist = 160000               #km/a Maximum driving distance per truck per year

    spec_capex_truck = 160000               #€
    spec_opex_truck = 0.12                  #% of CAPEX/a
    diesel_consumption = 35                 #l/100km
    truck_lifetime = 8                      #a

    spec_capex_trailor = 150000
    spec_opex_trailor = 0.02
    net_capacity = 1800                     #kgh2
    trailor_lifetime = 12                   #a
    loading_unloading_time = 1.5            #hours 


    max_day_dist = max_driving_dist/working_days
    amount_deliveries_needed = quantity/net_capacity
    deliveries_per_truck = working_hours/(loading_unloading_time+(2*distance/average_truck_speed))
    trailors_needed = round((amount_deliveries_needed/deliveries_per_truck)+0.5,0)
    total_drives_day = round(amount_deliveries_needed+0.5,0)
    trucks_needed = max(round((total_drives_day*2*distance*working_days/max_driving_dist)+0.5,0),trailors_needed)

    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    capex_total = capex_trailor + capex_trucks

    if amount_deliveries_needed < 1:
        fuel_costs = (amount_deliveries_needed*2*distance*365/100)*diesel_consumption*diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver
    
    else:
        fuel_costs = (round(amount_deliveries_needed+0.5)*2*distance*365/100)*diesel_consumption*diesel_price
        wages = round(amount_deliveries_needed+0.5) * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver


    annual_costs = (capex_trucks/RBF(interest,truck_lifetime)+capex_trailor/RBF(interest,trailor_lifetime)) + capex_trucks*spec_opex_truck + capex_trailor*spec_opex_trailor + fuel_costs + wages
    return annual_costs

def transport_NH3(distance,quantity,interest):

    quantity = quantity/365

    average_truck_speed = 70                #km/h
    working_hours = 24                      #h/day
    diesel_price = 1.5                      #€/l
    costs_for_driver = 20                   #€/h
    working_days = 365      

    spec_capex_truck = 160000               #€
    spec_opex_truck = 0.12                  #% of CAPEX/a
    diesel_consumption = 35                 #l/100km
    truck_lifetime = 12                     #a

    spec_capex_trailor = 190000
    spec_opex_trailor = 0.02
    net_capacity = 2600                     #kgh2
    trailor_lifetime = truck_lifetime       #a
    loading_unloading_time = 1.5            #hours 


    amount_deliveries_needed = quantity/net_capacity
    deliveries_per_truck = working_hours/(loading_unloading_time+(2*distance/average_truck_speed))
    trucks_trailors_needed = round((amount_deliveries_needed/deliveries_per_truck)+0.5,0)

    capex_trucks = trucks_trailors_needed * spec_capex_truck
    capex_trailor = trucks_trailors_needed * spec_capex_trailor
    capex_total = capex_trailor + capex_trucks

    if amount_deliveries_needed < 1:
        fuel_costs = (amount_deliveries_needed*2*distance*365/100)*diesel_consumption*diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver
    
    else:
        fuel_costs = (round(amount_deliveries_needed+0.5)*2*distance*365/100)*diesel_consumption*diesel_price
        wages = round(amount_deliveries_needed+0.5) * ((distance/average_truck_speed)*2+loading_unloading_time) * working_days * costs_for_driver

    annual_costs = (capex_trucks/RBF(interest,truck_lifetime)+capex_trailor/RBF(interest,trailor_lifetime)) + capex_trucks*spec_opex_truck + capex_trailor*spec_opex_trailor + fuel_costs + wages
    return annual_costs


#Transformation of hydrogen: Standard condition is pressure = 25 bar and T = 298.15 K
#There are four options: Standard condition, 500 bar, LH2, LOHC and NH3
#quantity per year in kg
#energy demand yearly
#default opex 2% of capex

def h2_conversion_stand(final_state, quantity, electricity_costs, heat_costs, interestrate):
    
    h2_throughput = quantity/365/24

    if final_state == 'standard condition':
        elec_demand = 0 
        heat_demand = 0
        annual_costs = 0 
        return elec_demand, heat_demand, annual_costs

    elif final_state == '500 bar':
            
        cp = 14200*((2.777778)*(10**-7))                                        #Wärmekapazität von h2 = 14200 J/kgK
        Tein = 298.15
        pein = 25
        k = 1.402                                                               # Isentropenexponent 
        n_isentrop = 0.8                                                        #Isentroper Wirkungsgrad
        elec_demand = (cp*Tein*(((500/pein)**((k-1)/k))-1))/n_isentrop          #kWh/kgh2
        elec_demand = elec_demand * quantity

        heat_demand =0 

            #specific investment costs derived from: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/760479/H2_supply_chain_evidence_-_publication_version.pdf
            
        compressor_lifetime = 15
        capex_compressor = 40035 * ((h2_throughput*24)**0.6038)
        opex_compressor = 0.08                                                  #Annahme 8% von Capex jährlich

        annual_costs = (capex_compressor/RBF(interestrate,compressor_lifetime)) + (capex_compressor*opex_compressor) + elec_demand * electricity_costs + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LH2':

        elec_demand = 9.93 * quantity
        heat_demand = 0

        capex_liquid_plant = -0.0002*((h2_throughput*24)**(2))+1781.9*(h2_throughput*24)+3*(10**7)
        opex_liquid_plant = 0.06
        liquid_plant_lifetime = 20

        annual_costs = (capex_liquid_plant/RBF(interestrate,liquid_plant_lifetime)) + (capex_liquid_plant*opex_liquid_plant) + elec_demand * electricity_costs + heat_demand*heat_costs
        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LOHC_load':

        elec_demand = 0.35 * quantity 
        heat_demand = 0 * quantity                                                 #https://www.hydrogenious.net/index.php/en/hydrogenious-3/lohc-technology/                             

        capex_hydrogenation = 0.84 * quantity
        opex_hydrogenation = 0.04
        hydrogenation_lifetime = 25
        costs_carrier = 2                                                       #€/kg_carrier       https://www.hydrogenious.net/index.php/en/2020/07/21/lohc-global-hydrogen-opportunity/#:~:text=carbon%20hydrogen%20transportation.-,Source%3A%20Hydrogenious%20LOHC%20Technologies.,hydrogen%20cost%20of%20%242.54%2Fkg
        ratio_carrier = 16.1                                                    #kg_carrier:kg_h2   https://www.tvt.tf.fau.eu/files/2019/03/lohc-lkw_bericht_final_teil_1.pdf

        annual_costs = ((capex_hydrogenation+(costs_carrier*ratio_carrier*h2_throughput*24))/RBF(interestrate, hydrogenation_lifetime)) + (capex_hydrogenation*opex_hydrogenation) + elec_demand * electricity_costs + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'LOHC_unload':

        elec_demand = 0.35 * quantity 
        heat_demand = 12 * quantity

        capex_dehydrogenation = 2.46 * quantity
        opex_dehydrogenation = 0.04
        dehydrogenation_lifetime = 25

        annual_costs = (capex_dehydrogenation/RBF(interestrate, dehydrogenation_lifetime)) + (capex_dehydrogenation*opex_dehydrogenation) + elec_demand * electricity_costs + heat_demand*heat_costs

        return elec_demand, heat_demand, annual_costs

    elif final_state == 'NH3_load':

        elec_demand = 2.809 * quantity + 0.3 * quantity                                           #see excel
        heat_demand = 0

        capex_NH3_plant = 4253.79 * quantity/1000
        opex_NH3_plant = 0.015
        NH3_plant_lifetime = 25

        annual_costs = capex_NH3_plant/RBF(interestrate,NH3_plant_lifetime) + capex_NH3_plant*opex_NH3_plant + elec_demand * electricity_costs + heat_demand*heat_costs
            
        return elec_demand, heat_demand, annual_costs
    elif final_state == 'NH3_unload':

        elec_demand = 4.2 * quantity                                         #see excel
        heat_demand = 0

        capex_NH3_plant = 18.171 * (((quantity/1000)/365/24) ** 0.7451) * 0.95 * 1000000   
        opex_NH3_plant = 0.02
        NH3_plant_lifetime = 25

        annual_costs = capex_NH3_plant/RBF(interestrate,NH3_plant_lifetime) + capex_NH3_plant*opex_NH3_plant + elec_demand * electricity_costs + heat_demand*heat_costs
            
        return elec_demand, heat_demand, annual_costs

def cheapest_dist_option(final_state, quantity, dist, elec_costs, heat_costs, interestrate, elec_costs_demand, days_storage):
    
    if final_state == '500 bar':
        dist_costs_500bar = storage_costs('500 bar',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lohc = storage_costs('LOHC',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lh2 = storage_costs('LH2',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_nh3 = storage_costs('NH3',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)

        dist_costs_500bar = dist_costs_500bar + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interestrate)[2] + transport_500bar(dist,quantity,interestrate)
        dist_costs_lh2 =  dist_costs_lh2 + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lh2(dist, quantity,interestrate) + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_lohc = dist_costs_lohc + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lohc(dist, quantity,interestrate) + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_nh3 = dist_costs_nh3 + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_NH3(dist, quantity,interestrate) + h2_conversion_stand('NH3_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]

        if dist_costs_500bar == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = '500 bar'
        elif dist_costs_lh2 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = 'LH2'
        elif dist_costs_lohc == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = 'LOHC'
        else:
            cheapest_option = 'NH3'

        return min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3)/quantity, cheapest_option

    elif final_state == 'LH2':
        dist_costs_500bar = storage_costs('500 bar',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lohc = storage_costs('LOHC',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lh2 = storage_costs('LH2',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_nh3 = storage_costs('NH3',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)

        dist_costs_500bar = dist_costs_500bar + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interestrate)[2] + transport_500bar(dist,quantity,interestrate) + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2]
        dist_costs_lh2 =  dist_costs_lh2 + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lh2(dist, quantity,interestrate) 
        dist_costs_lohc = dist_costs_lohc + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lohc(dist, quantity,interestrate) + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('LH2', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_nh3 = dist_costs_nh3 + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_NH3(dist, quantity,interestrate) + h2_conversion_stand('NH3_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('LH2', quantity, elec_costs_demand, heat_costs, interestrate)[2]

        if dist_costs_500bar == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = '500 bar'
        elif dist_costs_lh2 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = 'LH2'
        elif dist_costs_lohc == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = 'LOHC'
        else:
            cheapest_option = 'NH3'

        return min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3)/quantity, cheapest_option
    
    elif final_state == 'NH3':
        dist_costs_500bar = storage_costs('500 bar',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lohc = storage_costs('LOHC',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lh2 = storage_costs('LH2',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_nh3 = storage_costs('NH3',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)

        dist_costs_500bar = dist_costs_500bar + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interestrate)[2] + transport_500bar(dist,quantity,interestrate) + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2]
        dist_costs_lh2 =  dist_costs_lh2 + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lh2(dist, quantity, interestrate) + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2]
        dist_costs_lohc = dist_costs_lohc + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lohc(dist, quantity, interestrate) + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('NH3_load', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_nh3 = dist_costs_nh3 + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_NH3(dist, quantity, interestrate) 

        if dist_costs_500bar == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = '500 bar'
        elif dist_costs_lh2 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = 'LH2'
        elif dist_costs_lohc == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3):
            cheapest_option = 'LOHC'
        else:
            cheapest_option = 'NH3'

        return min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3)/quantity, cheapest_option

#Only new pipelines
def transport_pipeline(distance,quantity,elec_cost,interest):

    opex = 0.0125                                                           #% of capex per year
    availability = 0.95                                                     #Assumption
    lifetime_pipeline = 42.5
    lifetime_compressors = 24
    electricity_demand = 0.000613819                                        #kWh/kg*km

    capex_pipeline_big = 2.8 * 1000000                                      #Mio. €/km
    capex_compression_big = 0.62 * 1000000                                  #Mio. €/km
    max_capacity_big = 13                                                   #GW
    max_throughput_big = (((max_capacity_big*(10**6))/33.333))*8760*availability         #kg/year   

    capex_pipeline_med = 2.2 * 1000000                                       #Mio. €/km
    capex_compression_med = 0.31 * 1000000                                   #Mio. €/km
    max_capacity_med = 4.7                                                   #GW
    max_throughput_med = (((max_capacity_med*(10**6))/33.333))*8760*availability          #kg/year   

    capex_pipeline_sml = 1.5 * 1000000                                       #Mio. €/km
    capex_compression_sml = 0.09 * 1000000                                   #Mio. €/km
    max_capacity_sml = 1.2                                                   #GW
    max_throughput_sml = (((max_capacity_sml*(10**6))/33.333))*8760*availability         #kg/year   

    if quantity < max_throughput_sml:

        capex_annual = ((capex_pipeline_sml*distance)/RBF(interest,lifetime_pipeline)) + ((capex_compression_sml*distance)/RBF(interest,lifetime_compressors))
        opex_annual = opex*(capex_pipeline_sml+capex_compression_sml)*distance
        electricity_costs = electricity_demand * distance * quantity * elec_cost

        annual_costs = capex_annual + opex_annual + electricity_costs

        return annual_costs, "Small Pipeline"
    
    elif quantity > max_throughput_sml and quantity < max_throughput_med:

        capex_annual = ((capex_pipeline_med*distance)/RBF(interest,lifetime_pipeline)) + ((capex_compression_med*distance)/RBF(interest,lifetime_compressors))
        opex_annual = opex*(capex_pipeline_med+capex_compression_med)*distance
        electricity_costs = electricity_demand * distance * quantity

        annual_costs = capex_annual + opex_annual + electricity_costs * elec_cost

        return annual_costs, "Medium Pipeline"

    elif quantity > max_throughput_med and quantity < max_throughput_big:

        capex_annual = ((capex_pipeline_big*distance)/RBF(interest,lifetime_pipeline)) + ((capex_compression_big*distance)/RBF(interest,lifetime_compressors))
        opex_annual = opex*(capex_pipeline_big+capex_compression_big)*distance
        electricity_costs = electricity_demand * distance * quantity * elec_cost

        annual_costs = capex_annual + opex_annual + electricity_costs

        return annual_costs, "Big Pipeline"

    else:
        
        return 'No Pipeline big enough'

def cheapest_dist_option_pipeline(final_state, quantity, dist, elec_costs, heat_costs, interestrate, elec_costs_demand, elec_cost_grid, days_storage):
    
    if final_state == '500 bar':
        dist_costs_500bar = storage_costs('500 bar',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lohc = storage_costs('LOHC',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lh2 = storage_costs('LH2',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_nh3 = storage_costs('NH3',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_pipeline = 0 + storage_costs(final_state,quantity,days_storage,interestrate)

        dist_costs_500bar = dist_costs_500bar + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interestrate)[2] + transport_500bar(dist,quantity,interestrate)
        dist_costs_lh2 =  dist_costs_lh2 + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lh2(dist, quantity,interestrate) + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_lohc = dist_costs_lohc + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lohc(dist, quantity,interestrate) + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_nh3 = dist_costs_nh3 + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_NH3(dist, quantity,interestrate) + h2_conversion_stand('NH3_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_pipeline = dist_costs_pipeline + transport_pipeline(dist,quantity,elec_cost_grid,interestrate)[0] + h2_conversion_stand('500 bar', quantity, elec_costs_demand, heat_costs, interestrate)[2]

        if dist_costs_500bar == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = '500 bar'
        elif dist_costs_lh2 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = 'LH2'
        elif dist_costs_lohc == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = 'LOHC'
        elif dist_costs_nh3 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline): 
            cheapest_option = 'NH3'
        else:
            cheapest_option = transport_pipeline(dist,quantity,elec_cost_grid,interestrate)[1] 

        return min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline)/quantity, cheapest_option

    elif final_state == 'LH2':
        dist_costs_500bar = storage_costs('500 bar',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lohc = storage_costs('LOHC',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lh2 = storage_costs('LH2',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_nh3 = storage_costs('NH3',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_pipeline = 0 + storage_costs(final_state,quantity,days_storage,interestrate)

        dist_costs_500bar = dist_costs_500bar + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interestrate)[2] + transport_500bar(dist,quantity,interestrate) + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2]
        dist_costs_lh2 =  dist_costs_lh2 + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lh2(dist, quantity,interestrate) 
        dist_costs_lohc = dist_costs_lohc + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lohc(dist, quantity,interestrate) + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('LH2', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_nh3 = dist_costs_nh3 + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_NH3(dist, quantity,interestrate) + h2_conversion_stand('NH3_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('LH2', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_pipeline = dist_costs_pipeline + transport_pipeline(dist,quantity,elec_cost_grid,interestrate)[0] + h2_conversion_stand('LH2', quantity, elec_costs_demand, heat_costs, interestrate)[2]


        if dist_costs_500bar == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = '500 bar'
        elif dist_costs_lh2 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = 'LH2'
        elif dist_costs_lohc == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = 'LOHC'
        elif dist_costs_nh3 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline): 
            cheapest_option = 'NH3'
        else:
            cheapest_option = transport_pipeline(dist,quantity,elec_cost_grid,interestrate)[1] 

        return min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline)/quantity, cheapest_option
    

    elif final_state == 'NH3':
        dist_costs_500bar = storage_costs('500 bar',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lohc = storage_costs('LOHC',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_lh2 = storage_costs('LH2',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_nh3 = storage_costs('NH3',quantity,days_storage,interestrate) + storage_costs(final_state,quantity,days_storage,interestrate)
        dist_costs_pipeline = 0 + storage_costs(final_state,quantity,days_storage,interestrate)


        dist_costs_500bar = dist_costs_500bar + h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interestrate)[2] + transport_500bar(dist,quantity,interestrate) + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2]
        dist_costs_lh2 =  dist_costs_lh2 + h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lh2(dist, quantity,interestrate) + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2]
        dist_costs_lohc = dist_costs_lohc + h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_lohc(dist, quantity,interestrate) + h2_conversion_stand('LOHC_unload', quantity, elec_costs_demand, heat_costs, interestrate)[2] + h2_conversion_stand('NH3_load', quantity, elec_costs_demand, heat_costs, interestrate)[2]
        dist_costs_nh3 = dist_costs_nh3 + h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interestrate)[2] + transport_NH3(dist, quantity,interestrate) 
        dist_costs_pipeline = dist_costs_pipeline + transport_pipeline(dist,quantity,elec_cost_grid,interestrate)[0] + h2_conversion_stand('NH3_load', quantity, elec_costs_demand, heat_costs, interestrate)[2]

        if dist_costs_500bar == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = '500 bar'
        elif dist_costs_lh2 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = 'LH2'
        elif dist_costs_lohc == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline):
            cheapest_option = 'LOHC'
        elif dist_costs_nh3 == min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline): 
            cheapest_option = 'NH3'
        else:
            cheapest_option = transport_pipeline(dist,quantity,elec_cost_grid,interestrate)[1] 

        return min(dist_costs_500bar,dist_costs_lh2,dist_costs_lohc,dist_costs_nh3,dist_costs_pipeline)/quantity, cheapest_option

def storage_costs(state,quantity,storagedays,interest):
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
