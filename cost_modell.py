# from calendar import c, prcal
# from turtle import color, distance
import geopandas as gpd
# from matplotlib import markers
# from numpy import minimum
import pandas as pd
import matplotlib.pyplot as plt                 #see: https://geopandas.org/en/stable/docs/user_guide/mapping.html for plotting
from functions import *
from cmath import nan, pi
from shapely.geometry import Point
import shapely.geometry
import shapely.wkt
import geopy.distance
import PySimpleGUI as sg 
import math
from xlsxwriter import Workbook

elec_tech_data = pd.read_excel("Data/technology_parameters.xlsx",sheet_name= 'Electricity')
elec_tech = elec_tech_data.set_index('Technology').agg(list, axis=1).to_dict()

ely_tech_data = pd.read_excel("Data/technology_parameters.xlsx",sheet_name= 'Electrolysis')
ely_para = ely_tech_data.set_index('Parameter').agg(float, axis=1).to_dict()

wind_tech_data = pd.read_excel("Data/technology_parameters.xlsx",sheet_name='Wind turbine')
wind_para = wind_tech_data.set_index('Parameter').agg(float, axis=1).to_dict()

infra_data = pd.read_excel("Data/technology_parameters.xlsx",sheet_name='Infra')
infra_para = infra_data.set_index('Infrastructure').agg(list, axis=1).to_dict()

global_data = pd.read_excel("Data/technology_parameters.xlsx",sheet_name='Global')
global_para = global_data.set_index('Parameter').agg(list, axis=1).to_dict()

sg.theme('Reddit')

#button_menu = ['File',['Hydrogen state at destination', ['500 bar', 'Liquid H2']]]


# How many demand scenarios
layout1 = [[sg.Text('Amount of demand scenarios'), sg.InputText()],
           [sg.Submit(), sg.Exit()]]  

window1 = sg.Window('Demand scenarios', layout1, size=(500,100),font=("Helvetica", 15))

event, values1 = window1.read(close=True)   

# Where are demand locations?
amount_list = []

for i in range(int(values1[0])):
    demand_center = str(i+1)
    amount_list.append([sg.Text('Demand location '+demand_center +":",font=("Helvetica", 12)),sg.Text('Lat',font=("Helvetica", 12)), sg.InputText(default_text='-1.286', size=(10,10),font=("Helvetica", 12)),sg.Text('Lon',font=("Helvetica", 12)), sg.InputText(default_text='36.817', size=(10,10),font=("Helvetica", 12))])
    amount_list.append([sg.Text('Hydrogen demand '+demand_center +" [kg]:",font=("Helvetica", 12)), sg.InputText(default_text='10000000',font=("Helvetica", 12))])
    amount_list.append([sg.Text("Hydrogen state at destination",font=("Helvetica", 12)), sg.Combo(['500 bar','LH2','NH3'], default_value='500 bar',font=("Helvetica", 12))])

    amount_list.append([sg.Text()])


layout2 = [
           [sg.Column(amount_list, size=(500,500), scrollable=True, vertical_scroll_only=True)],
           [sg.Submit(), sg.Exit()]
           ]

window2 = sg.Window('Demand locations').Layout(layout2)

event, values2 = window2.read(close=True)   


#Possible adjustments in
layout3 = [
                 [sg.Text("Allow grid construction"), sg.Radio('Yes', 'group 3', key='Grid construction', default= True), sg.Radio('No', 'group 3')],
                 [sg.Text("Allow road construction"), sg.Radio('Yes', 'group 4', key='Road construction', default= True), sg.Radio('No', 'group 4')],
                 [sg.Text("Allow pipeline construction"), sg.Radio('Yes', 'group 5', default= True, key='Pipeline construction'), sg.Radio('No', 'group 5')],
                 [sg.Text("Show costs of Ammonia at destination (only chose yes, if state of destination is NH3)"), sg.Radio('Yes', 'group 2', key='ammonia_map'), sg.Radio('No', 'group 2',default= True)],
                 
                 [sg.Text('Wanted Graphs')],
                 [sg.Text('Locally prodcution')],
                 [sg.Checkbox('H2 Production costs [€/kg]', key = 'h2_prod_costs'), sg.Checkbox('LCOE [€/MWh]', key = 'cheapest_elec_cost')],
                 [sg.Checkbox('H2 capacity [kt/a]', key = 'h2_potential'), sg.Checkbox('Power generation capacity [TWh/a]', key = 'power_potential')],

                 [sg.Text('Demand coverage')],        
                 [sg.Checkbox('Cheapest production locations', key = 'prod_loc'), sg.Checkbox('Demand coverage', key = 'demand_cov')],
       


                 [sg.Submit(), sg.Exit()]]      


window3 = sg.Window('Coordinates of demand', layout3, size=(750,500),font=("Helvetica", 15))

#while True:
#    event, values = window.read()    
#    if event (sg.WIN_CLOSED, 'EXIT'):
#        break

#window.close

event, values3 = window3.read(close=True)    

demand_center_list = []


for i in range(int(values1[0])):
    
    demand_center_list.append([float(values2[(4*i)]), float(values2[(4*i+1)]), float(values2[(4*i+2)]), values2[(4*i+3)]])

#Fixed parameter declaration

# water_spec_cost = 1.2                       # €/m3
# h2_en_den = 33.33                           #kWh/kgh2
# days_of_storage = 3
# interest = 0.08

water_spec_cost = global_para['Water specific cost (euros/m3)'][0] # €/m3
h2_en_den = global_para['H2 energy density (kWh/ kg H2)'][0]       #kWh/kgh2
days_of_storage = global_para['Storage duration (days)'][0]
interest = global_para['Interest rate'][0]

pv_lifetime = elec_tech['PV'][2]
pv_opex = elec_tech['PV'][1]                #€/kWp*a
pv_capex = elec_tech['PV'][0]               #€/kWp

wind_lifetime = elec_tech['Wind'][2]
wind_opex = elec_tech['Wind'][1]            #€/kWp*a
wind_capex = elec_tech['Wind'][0]           #€/kW

ely_capex = ely_para['CAPEX']               # €/kW
ely_stack_replacement = ely_para['Stack']   # €/kW
ely_opex = ely_para['OPEX']                 # % CAPEX/a
ely_lt = ely_para['Lifetime']               # a
ely_eff = ely_para['Efficiency']                                              
ely_water = ely_para['Water_cons']          #liter/kg
flh_pv = int(ely_para['Fullload_pv'])            #h
flh_wind = int(ely_para['Fullload_wind'])        #h



grid_capex = infra_para['Grid'][0]               
elec_trans_costs = infra_para['Grid'][1]                #€/MWh PLatzhalter see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8661478/pdf/main.pdf
grid_lifetime = infra_para['Grid'][2]                   #years PLATZHALTER

road_capex_long = infra_para['Long road'][0]            #€/km from John Hine, converted to Euro (Assumed earth to paved road)
road_capex_short = infra_para['Short road'][0]          #€/km for raods < 10 km, from John Hine, converted to Euro (Assumed earth to paved road)
road_opex = infra_para['Short road'][1]                 #€/km/year from John Hine, converted to Euro (Assumed earth to paved road)
road_lifetime = infra_para['Short road'][2]             #years, assumption



#Data Input: Hexagon file
hex = gpd.read_file('Data/hex_final.geojson')

print(hex['theo_turbines'][1077],hex['theo_pv'][1077])

#PV electricty costs
pv_elec_cost = []
pv_ely_ratio = []


for i in range(len(hex)):
    pv_hourly_output = []
    input_ely_pv = []

    pv_hourly_output = []
    input_ely_pv = []

    pv_elec_cost_hex = (((pv_capex/RBF(interest,pv_lifetime))+pv_opex)/hex['pv'][i]/365) * 1000
    pv_elec_cost.append(pv_elec_cost_hex)
    #!!! replace this hardcoded number (where does it come from?)
    pv_ely_ratio.append(0.5916057743717814)


hex['pv_ely_ratio'] = pv_ely_ratio
hex['pv_elec_cost'] = pv_elec_cost

#Wind electricty costs
wind_elec_cost = []
wind_yearly_output = []
wind_ely_ratio = []
wind_ely_size = []

cp = wind_para['cp']                            #Coefficient of performance wind turbine 
den_air = wind_para['air density']              #Air density in kg/m3
d_rot = wind_para['rotor diameter']             #Diameter of rotor in kg/m3
field_eff = wind_para['field efficiency'] 
availability = wind_para['availability'] 
p_turb = wind_para['Power turbine'] 
start_speed = wind_para['Start']                #[m/s]
switchoff_speed = wind_para['Switch off']       #[m/s]

for i in range(len(hex)):

    wind_hourly_output = []
    input_ely_wind = []
    
    v_m = hex['wind'][i]                       # [m/s] mean wind speed
    A = v_m *(2/(pi**0.5))
    annual_power_output = 0

    #for p in range((int(switchoff_speed))):

    #    if (p+1) < start_speed:
    #        annual_power_output = annual_power_output + 0
        
    #    else:
    #        windspeed_yearly_output = 0.5 * cp * den_air * ((82**2)*pi/4) * ((p+1) ** 3) * field_eff * availability
    #        windspeed_probability = (1-math.exp(-(((p+1+0.5)/A)**2)))-(1-math.exp(-(((p-1-0.5)/A)**2)))
    #        annual_power_output = annual_power_output + windspeed_probability*windspeed_yearly_output*8760/1000



    for p in range((int(switchoff_speed))):

        if (p+1) < start_speed:
            annual_power_output = annual_power_output + 0
            windspeed_probability = (1-math.exp(-(((p+1+0.5)/A)**2)))-(1-math.exp(-(((p+1-0.5)/A)**2)))

            for h in range(round(8760*windspeed_probability)):
                wind_hourly_output.append(0)

        elif (p+1) == switchoff_speed:
            annual_power_output = annual_power_output + 0
            windspeed_probability = (1-math.exp(-(((p+1+0.5)/A)**2)))-(1-math.exp(-(((p+1-0.5)/A)**2)))

            for h in range(round(8760*windspeed_probability)):
                wind_hourly_output.append(0)
        
        else:
            #!!! where does the 82 come from?
            windspeed_output = 0.5 * cp * den_air * ((82**2)*pi/4) * ((p+1) ** 3) * field_eff * availability
            windspeed_probability = (1-math.exp(-(((p+1+0.5)/A)**2)))-(1-math.exp(-(((p+1-0.5)/A)**2)))

            if (windspeed_output/1000) > p_turb:
                annual_power_output = annual_power_output + windspeed_probability*p_turb*8760
                for h in range(round(8760*windspeed_probability)):
                    wind_hourly_output.append(p_turb)
            else:
                annual_power_output = annual_power_output + windspeed_probability*windspeed_output*8760/1000

                for h in range(round(8760*windspeed_probability)):
                    wind_hourly_output.append(windspeed_output/1000)

    wind_yearly_output.append(annual_power_output)
    turb_out = annual_power_output / p_turb
    wind_elec_cost_hex = ((wind_capex/RBF(interest,wind_lifetime)+wind_opex)/turb_out) * 1000
    #!!! where does this cutoff cost come from?
    if wind_elec_cost_hex > 150:
        wind_elec_cost_hex = nan
    wind_elec_cost.append(wind_elec_cost_hex)

    wind_hourly_output.sort(reverse=True)
    for o in wind_hourly_output:
        if o >= wind_hourly_output[flh_wind]:
            input_ely_wind.append(wind_hourly_output[flh_wind])
        else:
            input_ely_wind.append(o)

    wind_ely_ratio.append(sum(input_ely_wind)/annual_power_output)

    #ratio_wind_flh = wind_hourly_output[flh_wind]/p_turb
    #wind_ely_size.append(ratio_wind_flh)



hex['wind_ely_ratio'] = wind_ely_ratio
hex['wind_output'] = wind_yearly_output
hex['wind_elec_cost'] = wind_elec_cost

#!!! this is where I will need to input optimized cost based on wind-solar mix
#Cheapest electricity costs
cheapest_elec_cost = []
cheapest_elec_tech = []

for i in range(len(hex)):

    cheapest_elec_cost.append(cheapest(hex['wind_elec_cost'][i], hex['pv_elec_cost'][i]))
    if hex['wind_elec_cost'][i] == nan:
        cheapest_elec_tech.append('PV')
    elif hex['wind_elec_cost'][i]<hex['pv_elec_cost'][i]:
        cheapest_elec_tech.append('Wind')
    else:
        cheapest_elec_tech.append('PV')



hex['cheapest_elec_tech'] = cheapest_elec_tech
hex['cheapest_elec_cost'] = cheapest_elec_cost


#calculating power potential

elec_power_pot = []
h2_potential = []



for i in range(len(hex)):
    if hex['cheapest_elec_tech'][i] == 'PV' and hex['theo_pv'][i] >= 1:
        elec_power_pot.append((hex['theo_pv'][i]*hex['pv'][i]*1000)/1000000)
        

    elif hex['cheapest_elec_tech'][i] == 'Wind' and hex['theo_turbines'][i] >= 1:
        elec_power_pot.append((hex['theo_turbines'][i]*hex['wind_output'][i])/1000000)


    else:
        elec_power_pot.append(nan)


hex['power_potential'] = elec_power_pot
# !!! seems like energy density parameter is hardcoded here        
for i in range(len(hex)):
    if hex['cheapest_elec_tech'][i] == 'PV' and hex['theo_pv'][i] > 0:
        h2_potential.append(hex['power_potential'][i]*hex['pv_ely_ratio'][i]*ely_eff/33.33)
    
    elif hex['cheapest_elec_tech'][i] == 'Wind' and hex['theo_turbines'][i] > 0:
        h2_potential.append(hex['power_potential'][i]*hex['wind_ely_ratio'][i]*ely_eff/33.33)

    else:
        h2_potential.append(nan)

hex['h2_potential'] = h2_potential
        


    
#Transmission line: Possibility to trasnmit electricity over grid in other hexagons 


if values3['Grid construction'] == True:

    elec_cost_to_connect = []
    for i in range(len(hex)):
        if hex['grid_dist'][i] != 0:
            #!!! where are these numbers coming from?
            elec_cost_to_connect.append(hex['cheapest_elec_cost'][i]+((hex['grid_dist'][i]*grid_capex/RBF(interest,grid_lifetime))/(2000*8760*0.95*0.9)))
        else:
            elec_cost_to_connect.append(hex['cheapest_elec_cost'][i])

    
    elec_cost_at_grid = []
    for i in range(len(hex)):
        if hex['grid_dist'][i] == 0:
            elec_cost_at_grid.append(hex['cheapest_elec_cost'][i])
    
    
    cheapest_elec_cost_grid = []
    if min(elec_cost_at_grid) == min(elec_cost_to_connect):

        for i in range(len(hex)):
            if hex['grid_dist'][i] == 0:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], (min(elec_cost_at_grid)+elec_trans_costs)))
                if (min(elec_cost_at_grid)+elec_trans_costs) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'
            else:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], (min(elec_cost_at_grid)+elec_trans_costs+((hex['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime)))))
                #cheapest_elec_cost_grid.append(cheapest_elec_cost[i])
                if (min(elec_cost_at_grid)+elec_trans_costs+((hex['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime))) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'
    else:
        for i in range(len(hex)):
            if hex['grid_dist'][i] == 0:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i],min(elec_cost_to_connect)+elec_trans_costs))
                if (min(elec_cost_to_connect)+elec_trans_costs) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'

            elif i == min(range(len(elec_cost_to_connect)), key=elec_cost_to_connect.__getitem__):
                cheapest_elec_cost_grid.append(cheapest_elec_cost[i])
            else:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], (min(elec_cost_to_connect)+elec_trans_costs+((hex['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime)))))
                if (min(elec_cost_to_connect)+elec_trans_costs+((hex['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime))) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'
else:
    cheapest_elec_cost_grid = []
    elec_cost_at_grid = []

    for i in range(len(hex)):
        if hex['grid_dist'][i] == 0:
            elec_cost_at_grid.append(hex['cheapest_elec_cost'][i])

    for i in range(len(hex)):
        if hex['grid_dist'][i] == 0:
            cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], (min(elec_cost_at_grid)+elec_trans_costs)))
            if (min(elec_cost_at_grid)+elec_trans_costs) < cheapest_elec_cost[i]:
                cheapest_elec_tech[i] = 'Grid'
        else:
            cheapest_elec_cost_grid.append(cheapest_elec_cost[i])
            #cheapest_elec_cost_grid.append(cheapest_elec_cost[i])


hex['cheapest_elec_cost'] = cheapest_elec_cost_grid


#for i in range(len(hex)):

#    if cheapest_elec_cost_grid[i] != cheapest_elec_cost[i]:
#        cheapest_elec_tech[i] = 'Grid'


hex['cheapest_elec_tech'] = cheapest_elec_tech


h2o_costs_dom_water_bodies = []
h2o_costs_ocean = []
h2o_costs = []
# !!! add water costs to excel sheet
electricity_demand_h2o_treatment = 0.4                                  #kWh/          https://betterbuildingssolutioncenter.energy.gov/sites/default/files/Primer%20on%20energy%20efficiency%20in%20water%20and%20wastewater%20plants_0.pdf
electricity_demand_h2o_ocean_treatment = 3.7                            #kWh/m3     //https://www.pnas.org/doi/epdf/10.1073/pnas.1902335116, https://essay.utwente.nl/78100/1/Antonyan%2C%20M.%201817078%20_openbaar.pdf, 
water_transport_costs = 0.1                                             #€/100km/3
water_spec_cost = 1.25                                                  #€/m3

for i in range(len(hex)):
    h2o_costs_dom_water_bodies.append(((water_spec_cost + (water_transport_costs/100)*min(hex['waterbody_dist'][i],hex['waterway_dist'][i]) + electricity_demand_h2o_treatment*(hex['cheapest_elec_cost'][i]/1000)*ely_water)/1000))

for i in range(len(hex)):
    h2o_costs_ocean.append(((water_spec_cost + (water_transport_costs/100)*hex['ocean_dist'][i] + electricity_demand_h2o_ocean_treatment*(hex['cheapest_elec_cost'][i]/1000)*ely_water)/1000))

for i in range(len(hex)):
    h2o_costs.append(min(h2o_costs_dom_water_bodies[i],h2o_costs_ocean[i]))

h2_prod_costs = []


for i in range(len(hex)):
    if hex['cheapest_elec_tech'][i] == 'PV' and hex['theo_pv'][i] >= 1:
        ely_costs = ((((ely_capex+ely_stack_replacement)/RBF(interest,ely_lt))/(flh_pv))*(h2_en_den/ely_eff))*(1 + ely_opex)
        h2_prod_costs.append(((hex['cheapest_elec_cost'][i]/1000)* (h2_en_den/ely_eff)) + ely_costs + h2o_costs[i]*ely_water)
    
    elif hex['cheapest_elec_tech'][i] == 'Wind' and hex['theo_turbines'][i] >= 1:
        ely_costs = ((((ely_capex+ely_stack_replacement)/RBF(interest,ely_lt))/(flh_wind))*(h2_en_den/ely_eff))*(1 + ely_opex)
        h2_prod_costs.append(((hex['cheapest_elec_cost'][i]/1000)* (h2_en_den/ely_eff)) + ely_costs + h2o_costs[i]*ely_water)
    
    elif hex['cheapest_elec_tech'][i] == 'Grid':
        ely_costs = ((((ely_capex+ely_stack_replacement)/RBF(interest,ely_lt))/(flh_wind))*(h2_en_den/ely_eff))*(1 + ely_opex)
        h2_prod_costs.append(((hex['cheapest_elec_cost'][i]/1000)* (h2_en_den/ely_eff)) + ely_costs + h2o_costs[i]*ely_water)        

    else: 
        h2_prod_costs.append(nan)


hex['h2_prod_costs'] = h2_prod_costs



#print(hex)
#print(max(h2_prod_costs))
#print(min(h2o_costs_dom_water_bodies))

#distance to demand: Mombassa coordinates:  39.66007221012109, -4.039286378400124 info: flipped from google maps

distance_dict = {}
total_demand = 0

for d in range(len(demand_center_list)):


    total_demand = total_demand + demand_center_list[d][2]

    lat = demand_center_list[d][0]
    lon = demand_center_list[d][1]

    distance_to_demand = []

    for i in range(len(hex)):
        
        #Method takes to long so far
        #supply centre
        #poly = shapely.wkt.loads(str(hex['geometry'][i]))
        #center = poly.centroid
        #coords_1 = (-4.039286378400124, 39.66007221012109)
        #coords_1 = (lat,lon)
        #coords_2 = (center.y, center.x)

        #url = "http://router.project-osrm.org/route/v1/driving/" + str(coords_2[1]) +','+ str(coords_2[0]) +';'+ str(coords_1[1]) +','+ str(coords_1[0]) + '?overview=false'

        #r = requests.get(url)
        # then you load the response using the json libray
        # by default you get only one alternative so you access 0-th element of the `routes`
        #routes = json.loads(r.content)
        #route_1 = json.loads(r.content)["routes"][0]

        #dist = route_1['distance']/1000

        #print(i)


        poly = shapely.wkt.loads(str(hex['geometry'][i]))
        center = poly.centroid
        #coords_1 = (-4.039286378400124, 39.66007221012109)
        coords_1 = (lat,lon)
        coords_2 = (center.y, center.x)
        dist = geopy.distance.geodesic(coords_1, coords_2).km

        #print(dist)

        distance_to_demand.append(dist)


    #determine elec_cost at demand to determine potential energy costs

    demand_location = Point(float(lon),float(lat))
    demand_fid = 0

    for i in range(len(hex)):
        if hex['geometry'][i].contains(demand_location) == True:
            demand_fid = i

    elec_costs_at_demand = float(hex['cheapest_elec_cost'][demand_fid])/1000


    hydrogen_quantity = demand_center_list[d][2]
    h2_costs_incl_conversion = []
    h2_costs_to_demand = []
    road_construction_costs = []
    transport_type = []

    for i in range(len(hex)):
        if hex['road_dist'][i]==0:
            road_construction_costs.append(0)
        elif hex['road_dist'][i]!=0 and hex['road_dist'][i]<10:
            road_construction_costs.append(((hex['road_dist'][i]*road_capex_short)/(RBF(interest,road_lifetime)))+(hex['road_dist'][i]*road_opex))
        else:
            road_construction_costs.append(((hex['road_dist'][i]*road_capex_long)/(RBF(interest,road_lifetime)))+(hex['road_dist'][i]*road_opex))


    
        
    if values3['Pipeline construction'] == True: 
        if values3['Road construction'] == True:

            for i in range(len(hex)):


                if i == demand_fid:
                    if demand_center_list[d][3] == 'NH3':
                    # !!! where are the 0.03 values coming from?
                        h2_costs_incl_conversion.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3]+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        h2_costs_to_demand.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3]+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        transport_type.append('None')

                    else:

                        h2_costs_incl_conversion.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        h2_costs_to_demand.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        transport_type.append('None')


                else:
                
                    if demand_center_list[d][3] == '500 bar':
                        h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)+h2_prod_costs[i]+cheapest_dist_option_pipeline('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[0])

                        transport_state = cheapest_dist_option_pipeline('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,min(cheapest_elec_cost_grid)/1000,days_of_storage)[1]
                        transport_type.append(transport_state)


                        if  transport_state != 'NH3' and transport_state != 'LOHC' and transport_state != "Small Pipeline" and transport_state != "Medium Pipeline" and transport_state != "Large Pipeline":
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        elif transport_state == "Small Pipeline" or transport_state == "Medium Pipeline" or transport_state == "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i])
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                

                    if demand_center_list[d][3] == 'LH2':
                        h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)+h2_prod_costs[i]+cheapest_dist_option_pipeline('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[0])

                        transport_state = cheapest_dist_option_pipeline('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,min(cheapest_elec_cost_grid)/1000,days_of_storage)[1]
                        transport_type.append(transport_state)


                        if  transport_state != 'NH3' and transport_state != 'LOHC' and transport_state != "Small Pipeline" and transport_state != "Medium Pipeline" and transport_state != "Large Pipeline":
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        elif transport_state == "Small Pipeline" or transport_state == "Medium Pipeline" or transport_state == "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i])
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                    if demand_center_list[d][3] == 'NH3':
                        h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)+h2_prod_costs[i]+cheapest_dist_option_pipeline('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[0])

                        transport_state = cheapest_dist_option_pipeline('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,min(cheapest_elec_cost_grid)/1000,days_of_storage)[1]
                        transport_type.append(transport_state)
                            
                        if  transport_state != 'NH3' and transport_state != 'LOHC' and transport_state != "Small Pipeline" and transport_state != "Medium Pipeline" and transport_state != "Large Pipeline":
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        elif transport_state == "Small Pipeline" or transport_state == "Medium Pipeline" or transport_state == "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i])
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                

        else: 
            for i in range(len(hex)):

                if i == demand_fid:

                    h2_costs_incl_conversion.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                    h2_costs_to_demand.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                    transport_type.append('None')

                else:

                    if hex['road_dist'][i]==0:
                        if demand_center_list[d][3] == '500 bar':
                            h2_costs_to_demand.append(h2_prod_costs[i]+cheapest_dist_option_pipeline('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[0])
                            
                            transport_state = cheapest_dist_option_pipeline('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[1]
                            transport_type.append(transport_state)


                            if  transport_state != 'NH3' and transport_state != 'LOHC' and transport_state != "Small Pipeline" and transport_state != "Medium Pipeline" and transport_state != "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                            elif transport_state == "Small Pipeline" or transport_state == "Medium Pipeline" or transport_state == "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i])
                            else:
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                        if demand_center_list[d][3] == 'LH2':
                            h2_costs_to_demand.append(h2_prod_costs[i]+cheapest_dist_option_pipeline('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[0])
                            
                            transport_state = cheapest_dist_option_pipeline('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[1]
                            transport_type.append(transport_state)

                            if  transport_state != 'NH3' and transport_state != 'LOHC' and transport_state != "Small Pipeline" and transport_state != "Medium Pipeline" and transport_state != "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                            elif transport_state == "Small Pipeline" or transport_state == "Medium Pipeline" or transport_state == "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i])
                            else:
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                        
                        if demand_center_list[d][3] == 'NH3':
                            h2_costs_to_demand.append(h2_prod_costs[i]+cheapest_dist_option_pipeline('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[0])
                            
                            transport_state = cheapest_dist_option_pipeline('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand, min(cheapest_elec_cost_grid)/1000,days_of_storage)[1]
                            transport_type.append(transport_state)

                            if  transport_state != 'NH3' and transport_state != 'LOHC' and transport_state != "Small Pipeline" and transport_state != "Medium Pipeline" and transport_state != "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                            elif transport_state == "Small Pipeline" or transport_state == "Medium Pipeline" or transport_state == "Large Pipeline":
                                h2_costs_incl_conversion.append(h2_prod_costs[i])     
                            else:
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                    else:
                        if demand_center_list[d][3] == '500 bar':
                            pipeline_dist = storage_costs('500 bar',hydrogen_quantity,days_of_storage,interest) + transport_pipeline(distance_to_demand[i],hydrogen_quantity,min(cheapest_elec_cost_grid)/1000,interest)[0] + h2_conversion_stand('500 bar', hydrogen_quantity, elec_costs_at_demand, 0.03, interest)[2]
                        
                            
                            transport_type.append(transport_pipeline(distance_to_demand[i],hydrogen_quantity,min(cheapest_elec_cost_grid)/1000,interest)[1])


                            h2_costs_to_demand.append(h2_prod_costs[i]+(pipeline_dist/hydrogen_quantity))
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand('500 bar', hydrogen_quantity, elec_costs_at_demand, 0.03, interest)[2]/hydrogen_quantity)

                        if demand_center_list[d][3] == 'LH2':
                            pipeline_dist = storage_costs('LH2',hydrogen_quantity,days_of_storage,interest) + transport_pipeline(distance_to_demand[i],hydrogen_quantity,min(cheapest_elec_cost_grid)/1000,interest)[0] + h2_conversion_stand('LH2', hydrogen_quantity, elec_costs_at_demand, 0.03, interest)[2]
                            
                            transport_type.append(transport_pipeline(distance_to_demand[i],hydrogen_quantity,min(cheapest_elec_cost_grid)/1000,interest)[1])

                            h2_costs_to_demand.append(h2_prod_costs[i]+(pipeline_dist/hydrogen_quantity))
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand('LH2', hydrogen_quantity, elec_costs_at_demand, 0.03, interest)[2]/hydrogen_quantity)
                            
                        if demand_center_list[d][3] == 'NH3':
                            pipeline_dist = storage_costs('NH3',hydrogen_quantity,days_of_storage,interest) + transport_pipeline(distance_to_demand[i],hydrogen_quantity,min(cheapest_elec_cost_grid)/1000,interest)[0] + h2_conversion_stand('NH3_load', hydrogen_quantity, elec_costs_at_demand, 0.03, interest)[2]
                            
                            transport_type.append(transport_pipeline(distance_to_demand[i],hydrogen_quantity,min(cheapest_elec_cost_grid)/1000,interest)[1])

                            h2_costs_to_demand.append(h2_prod_costs[i]+(pipeline_dist/hydrogen_quantity))
                            h2_costs_incl_conversion.append(h2_prod_costs[i])#+(h2_conversion_stand('NH3_load', hydrogen_quantity, elec_costs_at_demand, 0.03, interest)[2]/hydrogen_quantity))
                        
                        #h2_costs_to_demand.append((nan))
                        #h2_costs_incl_conversion.append(nan)

    else:
        if values3['Road construction'] == True:
            for i in range(len(hex)):

                if i == demand_fid:

                    h2_costs_incl_conversion.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                    h2_costs_to_demand.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)

                else:
                
                    if demand_center_list[d][3] == '500 bar':
                        h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)+h2_prod_costs[i]+cheapest_dist_option('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[0])

                        transport_state = cheapest_dist_option('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[1]
                        transport_type.append(transport_state)


                        if  transport_state != 'NH3' and transport_state != 'LOHC':
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                

                    if demand_center_list[d][3] == 'LH2':
                        h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)+h2_prod_costs[i]+cheapest_dist_option('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[0])
                    
                        transport_state = cheapest_dist_option('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[1]
                        transport_type.append(transport_state)


                        if  transport_state != 'NH3' and transport_state != 'LOHC':
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                    if demand_center_list[d][3] == 'NH3':
                        h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)+h2_prod_costs[i]+cheapest_dist_option('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[0])
                        
                        transport_state = cheapest_dist_option('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[1]
                        transport_type.append(transport_state)

                        if  transport_state != 'NH3' and transport_state != 'LOHC':
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
            

        else: 
            for i in range(len(hex)):

                if i == demand_fid:

                    h2_costs_incl_conversion.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                    h2_costs_to_demand.append(hex['h2_prod_costs'][i]+h2_conversion_stand(demand_center_list[d][3],hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)

                else:
    
                    if hex['road_dist'][i]==0:
                        if demand_center_list[d][3] == '500 bar':
                            h2_costs_to_demand.append(h2_prod_costs[i]+cheapest_dist_option('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[0])
                            
                            transport_state = cheapest_dist_option('500 bar', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[1]
                            transport_type.append(transport_state)

                            if  transport_state != 'NH3' and transport_state != 'LOHC':
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                            else:
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                        if demand_center_list[d][3] == 'LH2':
                            h2_costs_to_demand.append(h2_prod_costs[i]+cheapest_dist_option('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[0])
                        
                            transport_state = cheapest_dist_option('LH2', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[1]
                            transport_type.append(transport_state)

                            if  transport_state != 'NH3' and transport_state != 'LOHC':
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                            else:
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                        
                        if demand_center_list[d][3] == 'NH3':
                            h2_costs_to_demand.append(h2_prod_costs[i]+cheapest_dist_option('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[0])
                
                            transport_state = cheapest_dist_option('NH3', hydrogen_quantity, distance_to_demand[i], cheapest_elec_cost[i]/1000, 0.03, interest, elec_costs_at_demand,days_of_storage)[1]
                            transport_type.append(transport_state)

                            if  transport_state != 'NH3' and transport_state != 'LOHC':
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state,hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                            else:
                                h2_costs_incl_conversion.append(h2_prod_costs[i]+h2_conversion_stand(transport_state+'_load',hydrogen_quantity,cheapest_elec_cost[i]/1000,0.03,interest)[2]/hydrogen_quantity)
                
                    else:
                        h2_costs_to_demand.append((nan))
                        h2_costs_incl_conversion.append(nan)

    #h2_costs_to_demand = [round(num, 1) for num in h2_costs_to_demand]

    hex['h2_costs_to_demand' + str(d)] = h2_costs_to_demand
    hex['h2_costs_incl_conv' + str(d)] = h2_costs_incl_conversion
    hex['transport_type' + str(d)] = transport_type


#print(distance_to_demand)
NH3_costs_to_demand = []

if values3['ammonia_map'] == True:
    for i in range(len(hex)):
        # !!! add weight of h2 in NH3 as a varible, maybe in excel?
        h2_weight_NH3 = 0.178       
        NH3_costs_to_demand.append((h2_costs_to_demand[i]+(h2_conversion_stand('NH3_unload',hydrogen_quantity,elec_costs_at_demand/1000,0.03,interest)[2]/hydrogen_quantity))*h2_weight_NH3)


    hex['NH3_costs_to_demand'] = NH3_costs_to_demand


#See results and prepare plots
#print(hex)
hex_sort = {}
for d in range(len(demand_center_list)):

    hex_sort[d] = hex.sort_values(by=['h2_costs_to_demand' + str(d)]).head(100)

#print(hex_sort)
#gpd.hex_sort.to_excel()
#print(values[6])

h2_costs_to_demand = []

for i in range(len(hex)):
    options = []

    for d in range(len(demand_center_list)):
        options.append((hex['h2_costs_to_demand' + str(d)][i]))
    
    h2_costs_to_demand.append(min(options))
        

hex['h2_costs_to_demand'] = h2_costs_to_demand
# create excel writer object
#minimal_cost = pd.ExcelWriter('output.xlsx')
# write dataframe to excel
#hex_sort.to_excel(minimal_cost)

#minimal_cost.save()

kenya_shp = gpd.read_file('Data/kenyan-counties/County.shp')


#plotting

if values3['prod_loc'] == True:
    fig, ax = plt.subplots(figsize=(10, 8))
    for d in range(len(demand_center_list)):
        ax.plot(float(demand_center_list[d][1]), float(demand_center_list[d][0]), color= "black", markersize= 10, marker = "1", label='Demand center',linestyle="None")
    kenya_shp.plot(ax=ax, color='whitesmoke', legend = True, edgecolor='dimgrey', linewidth = .5)

    colour = ['Blues_r','Greens_r','Greys_r','Oranges_r','Purples_r','Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']


    for d in range(len(demand_center_list)):
        hex_sort[d].plot(column='h2_costs_to_demand' + str(d),legend = True, edgecolor='dimgrey', linewidth = .5 ,ax=ax, cmap= colour[d])


list_keys = []
for d in range(len(demand_center_list)):
    list_keys.append(d)


supply_demand = dict.fromkeys(list_keys)

supply_loc = {}

hex['h2_cap_after_supply'] = hex['h2_potential']

for d in range(len(demand_center_list)):
    supply = 0
    nr = 0
    supply_list = []
    index_list = list(hex_sort[d]['h2_costs_to_demand' +str(d)].index.values)

    while supply < float(demand_center_list[d][2]):
        if hex['h2_potential'][index_list[nr]]*1000000 > 0:
            if hex['h2_potential'][index_list[nr]]*1000000 >= (float(demand_center_list[d][2]) - supply):
                coverage = (float(demand_center_list[d][2]) - supply)
                supply = supply + coverage
                hex['h2_potential'][index_list[nr]] = hex['h2_potential'][index_list[nr]] - (coverage/1000000)
                supply_list.append(index_list[nr])
                nr = nr + 1
            else:
                coverage = hex['h2_potential'][index_list[nr]]*1000000
                supply = supply + hex['h2_potential'][index_list[nr]]*1000000
                hex['h2_potential'][index_list[nr]] = 0 
                supply_list.append(index_list[nr])
                nr = nr + 1
        else:
            nr = nr + 1

    supply_loc[d] = supply_list

print(supply_loc)


if values3['demand_cov'] == True:
    fig, ax = plt.subplots(figsize=(10, 8))
    for d in range(len(demand_center_list)):
        
        
        ax.plot(float(demand_center_list[d][1]), float(demand_center_list[d][0]), color= "black", markersize= 10, marker = "1", label='Demand center',linestyle="None")
        
        for s in supply_loc[d]:
            poly = shapely.wkt.loads(str(hex['geometry'][s]))
            center = poly.centroid
            ax.plot(center.x, center.y, markersize= 7, marker = "h", label='Demand center',linestyle="None")
            plt.plot([float(demand_center_list[d][1]),center.x], [float(demand_center_list[d][0]),center.y])

    kenya_shp.plot(ax=ax, color='whitesmoke', legend = True, edgecolor='dimgrey', linewidth = .5)


                
map_list = ['h2_prod_costs', 'h2_potential', 'cheapest_elec_cost', 'power_potential']
map_titel = {'h2_prod_costs':'LCOH [€/kg]', 'h2_potential':'H2 Capacity [kt/a]', 'cheapest_elec_cost':'LCOE [€/MWh]', 'power_potential': 'Power generation potential [TWh/a]'}
colour_dict = {'h2_prod_costs':'Greens_r', 'h2_potential':'Oranges', 'cheapest_elec_cost':'Blues_r', 'power_potential': 'Purples'}

map_count = 0
for i in map_list:
    if values3[i] == True:
        map_count = map_count + 1



for i in range(len(map_list)):
    if values3[map_list[i]] == True: 
        figure, ax = plt.subplots(1, 1)
        hex.plot(column= map_list[i], legend = True, cmap= colour_dict[map_list[i]], ax=ax)

        plt.title(map_titel[map_list[i]])


result_list = ['Supplied by hexagons:', 'Amount per hexagon:', 'Total cost for H2 from hexagons', 'Average price:', 'Average price breakup:']

result_dict = {}
for d in range(len(demand_center_list)):

    result_dict[d] = []
    result_dict[d].append(str(supply_loc[d]))




# create excel writer object
#output = pd.ExcelWriter('model_output.xlsx')

# write dataframe to excel
#for d in range(len(demand_center_list)):
    #hex = hex.sort_values(by=['h2_costs_to_demand' + str(d)])
    #hex['h2_prod_costs'].to_excel(output ,sheet_name = str(d), startcol=0)
    #hex['transport_type'+ str(d)].to_excel(output, sheet_name = str(d), startcol=2, index=False)
    #hex['h2_costs_incl_conv' + str(d)].to_excel(output, sheet_name = str(d), startcol=3, index=False)
    #hex['h2_costs_to_demand' + str(d)].to_excel(output, sheet_name = str(d), startcol=4, index=False)

#    result_dict[d].to_excel(output ,sheet_name = str(d), startcol=0)


#output.save()

workbook = Workbook('Ecl.xlsx')
for d in range(len(demand_center_list)):   
    Report_Sheet = workbook.add_worksheet()

# Write the column headers if required.
    #Report_Sheet.write(0, 0, 'Column1')
#Report_Sheet.write(0, 1, 'Column2')

# Write the column data.
    Report_Sheet.write_column(1, 0, result_list)
    Report_Sheet.write_column(1, 1, result_dict[d])
#Report_Sheet.write_column(1, 1, Column2)

workbook.close()



plt.show()


