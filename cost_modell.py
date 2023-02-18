'''
Originally written by Leander Müller at RWTH Aachen University
Updated by Claire Halloran at University of Oxford
Supported by the Climate Compatible Growth programme

Calculates the cost of hydrogen production in different regions
to meet different end demands

'''

# from calendar import c, prcal
# from turtle import color, distance
import geopandas as gpd
# from matplotlib import markers
from numpy import nanmin
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

#%% Data Input
# Hexagon file
hexagon = gpd.read_file('Data/hex_final.geojson')

# load shapefile of Kenyan counties
kenya_shp = gpd.read_file('Data/kenyan-counties/County.shp')

# Excel file with technology parameters
technology_parameters = "Data/technology_parameters.xlsx"

#%% load data from technology parameters Excel file
# 2D data is a dataframe, 1D data is a series
elec_tech_data = pd.read_excel(technology_parameters,
                               sheet_name= 'Electricity',
                               index_col='Technology')

ely_tech_data = pd.read_excel(technology_parameters,
                              sheet_name= 'Electrolysis',
                              index_col='Parameter'
                              ).squeeze("columns")

wind_tech_data = pd.read_excel(technology_parameters,
                               sheet_name='Wind turbine',
                               index_col='Parameter'
                               ).squeeze("columns")

infra_data = pd.read_excel(technology_parameters,
                           sheet_name='Infra',
                           index_col='Infrastructure')

global_data = pd.read_excel(technology_parameters,
                            sheet_name='Global',
                            index_col='Parameter'
                            ).squeeze("columns")

water_data = pd.read_excel(technology_parameters,
                            sheet_name='Water',
                            index_col='Parameter'
                            ).squeeze("columns")
#%% Fixed parameter declaration from excel file values

# water_spec_cost = 1.2                       # €/m3
# h2_en_den = 33.33                           #kWh/kgh2
# days_of_storage = 3
# interest = 0.08

water_spec_cost = global_data['Water specific cost (euros/m3)'] # €/m3
h2_en_den = global_data['H2 energy density (kWh/ kg H2)']     #kWh/kgh2
days_of_storage = global_data['Storage duration (days)']
interest = global_data['Interest rate']

pv_lifetime = elec_tech_data.at['PV','Lifetime']
pv_opex = elec_tech_data.at['PV','OPEX']                #€/kWp*a
pv_capex = elec_tech_data.at['PV','CAPEX']               #€/kWp

wind_lifetime = elec_tech_data.at['Wind','Lifetime']
wind_opex = elec_tech_data.at['Wind','OPEX']            #€/kWp*a
wind_capex = elec_tech_data.at['Wind','CAPEX']           #€/kW

ely_capex = ely_tech_data['CAPEX']               # €/kW
ely_stack_replacement = ely_tech_data['Stack']   # €/kW
ely_opex = ely_tech_data['OPEX']                 # % CAPEX/a
ely_lt = ely_tech_data['Lifetime']               # a
ely_eff = ely_tech_data['Efficiency']                                              
ely_water = ely_tech_data['Water_cons']          #liter/kg
flh_pv = int(ely_tech_data['Fullload_pv'])            #h
flh_wind = int(ely_tech_data['Fullload_wind'])        #h



grid_capex = infra_data.at['Grid','CAPEX']               
elec_trans_costs = infra_data.at['Grid','OPEX']                #€/MWh PLatzhalter see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8661478/pdf/main.pdf
grid_lifetime = infra_data.at['Grid','Lifetime']                  #years PLATZHALTER

road_capex_long = infra_data.at['Long road','CAPEX']            #€/km from John Hine, converted to Euro (Assumed earth to paved road)
road_capex_short = infra_data.at['Short road','CAPEX']         #€/km for raods < 10 km, from John Hine, converted to Euro (Assumed earth to paved road)
road_opex = infra_data.at['Short road','OPEX']                 #€/km/year from John Hine, converted to Euro (Assumed earth to paved road)
road_lifetime = infra_data.at['Short road','Lifetime']             #years, assumption

cp = wind_tech_data['cp']                            #Coefficient of performance wind turbine 
den_air = wind_tech_data['air density']              #Air density in kg/m3
d_rot = wind_tech_data['rotor diameter']             #Diameter of rotor in kg/m3
field_eff = wind_tech_data['field efficiency'] 
availability = wind_tech_data['availability'] 
p_turb = wind_tech_data['Power turbine'] 
start_speed = wind_tech_data['Start']                #[m/s]
switchoff_speed = wind_tech_data['Switch off']       #[m/s]

#%% prompt user for number of demand scenarios (what are those scenarios?)
sg.theme('Reddit')
#!!! rename values1-3 to something more descriptive
#button_menu = ['File',['Hydrogen state at destination', ['500 bar', 'Liquid H2']]]


# How many demand scenarios (specify location of demand, unsure which kind of demand)
# for each scenario, user can specify demand lat-lon, demand for hydrogen in kg
# !!! change demand to tonnes? kt? try to keep consistent units throughout
# and hydrogen state at destination
layout1 = [[sg.Text('Amount of demand scenarios'), sg.InputText()],
           [sg.Submit(), sg.Exit()]]  

window1 = sg.Window('Demand scenarios', layout1, size=(500,100),font=("Helvetica", 15))

event, values1 = window1.read(close=True)   

# Where are demand locations?
amount_list = []

for i in range(int(values1[0])):
    demand_center = str(i+1)
    amount_list.append([sg.Text('Demand location '+demand_center +":",font=("Helvetica", 12)),
                        sg.Text('Lat',font=("Helvetica", 12)), 
                        sg.InputText(default_text='-1.286', size=(10,10),font=("Helvetica", 12)),
                        sg.Text('Lon',font=("Helvetica", 12)), 
                        sg.InputText(default_text='36.817',size=(10,10),font=("Helvetica", 12))])
    amount_list.append([sg.Text('Hydrogen demand '+demand_center +" [kg]:",font=("Helvetica", 12)),
                        sg.InputText(default_text='10000000',font=("Helvetica", 12))])
    amount_list.append([sg.Text("Hydrogen state at destination",
                                font=("Helvetica", 12)), 
                                sg.Combo(['500 bar','LH2','NH3'], 
                                         default_value='500 bar',
                                         font=("Helvetica", 12))])

    amount_list.append([sg.Text()])


layout2 = [
           [sg.Column(amount_list, size=(500,500), scrollable=True, vertical_scroll_only=True)],
           [sg.Submit(), sg.Exit()]
           ]

window2 = sg.Window('Demand locations').Layout(layout2)

event, values2 = window2.read(close=True)   

# adjust whether different construction is allowed, whether to show cost of ammonia
# adjust maps produced
# adjust demand coverage and cheapest production locations (unsure what this means)

layout3 = [
                 [sg.Text("Allow grid construction"), 
                  sg.Radio('Yes', 'group 3', key='Grid construction', default= True), 
                  sg.Radio('No', 'group 3')],
                 [sg.Text("Allow road construction"), 
                  sg.Radio('Yes', 'group 4', key='Road construction', default= True), 
                  sg.Radio('No', 'group 4')],
                 [sg.Text("Allow pipeline construction"), 
                  sg.Radio('Yes', 'group 5', default= True, key='Pipeline construction'), 
                  sg.Radio('No', 'group 5')],
                 [sg.Text("Show costs of ammonia at destination (only chose yes, if state of destination is NH3)"), 
                  sg.Radio('Yes', 'group 2', key='ammonia_map'), 
                  sg.Radio('No', 'group 2',default= True)],
                 
                 [sg.Text('Maps to produce')],
                 [sg.Text('Locally production')],
                 [sg.Checkbox('H2 Production costs [€/kg]', key = 'h2_prod_costs'), 
                  sg.Checkbox('LCOE [€/MWh]', key = 'cheapest_elec_cost')],
                 [sg.Checkbox('H2 capacity [kt/a]', key = 'h2_potential'), 
                  sg.Checkbox('Power generation capacity [TWh/a]', key = 'power_potential')],

                 [sg.Text('Demand coverage')],        
                 [sg.Checkbox('Cheapest production locations', key = 'prod_loc'), 
                  sg.Checkbox('Demand coverage', key = 'demand_cov')],
       


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
    
    demand_center_list.append([float(values2[(4*i)]), 
                               float(values2[(4*i+1)]), 
                               float(values2[(4*i+2)]), 
                               values2[(4*i+3)]])

#%% calculate LCOE for PV and wind in each hexagon
#PV electricty costs
pv_elec_cost = []
pv_ely_ratio = []

for i in range(len(hexagon)):
    pv_hourly_output = []
    input_ely_pv = []

    pv_hourly_output = []
    input_ely_pv = []

    pv_elec_cost_hex = (((pv_capex/RBF(interest,pv_lifetime))+pv_opex)/hexagon['pv'][i]/365) * 1000
    pv_elec_cost.append(pv_elec_cost_hex)
    #!!! replace this hardcoded number (where does it come from?)
    pv_ely_ratio.append(0.5916057743717814)


hexagon['pv_ely_ratio'] = pv_ely_ratio
hexagon['pv_elec_cost'] = pv_elec_cost

#Wind electricty costs
wind_elec_cost = []
wind_yearly_output = []
wind_ely_ratio = []
wind_ely_size = []

#%% calculate wind output in each hexagon based on assumed distribution of wind speed and 
# also calculate LCOE of wind electricity
for i in range(len(hexagon)):

    wind_hourly_output = []
    input_ely_wind = []
    
    v_m = hexagon['wind'][i]                       # [m/s] mean wind speed
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

hexagon['wind_ely_ratio'] = wind_ely_ratio
hexagon['wind_output'] = wind_yearly_output
hexagon['wind_elec_cost'] = wind_elec_cost

#%% identify lowest-cost electricity souce (wind vs solar)
#!!! this is where I will need to input optimized cost based on wind-solar mix
#Cheapest electricity costs
cheapest_elec_cost = []
cheapest_elec_tech = []

for i in range(len(hexagon)):
    cheapest_elec_cost.append(nanmin([hexagon['wind_elec_cost'][i], hexagon['pv_elec_cost'][i]]))
    # cheapest_elec_cost.append(cheapest(hexagon['wind_elec_cost'][i], hexagon['pv_elec_cost'][i]))
    if hexagon['wind_elec_cost'][i] == nan:
        cheapest_elec_tech.append('PV')
    elif hexagon['wind_elec_cost'][i]<hexagon['pv_elec_cost'][i]:
        cheapest_elec_tech.append('Wind')
    else:
        cheapest_elec_tech.append('PV')

hexagon['cheapest_elec_tech'] = cheapest_elec_tech
hexagon['cheapest_elec_cost'] = cheapest_elec_cost

#%% calculating power potential from cheapest electricity technology
elec_power_pot = []
# what units are used here? seem to be doing conversions
for i in range(len(hexagon)):
    if hexagon['cheapest_elec_tech'][i] == 'PV' and hexagon['theo_pv'][i] >= 1:
        elec_power_pot.append((hexagon['theo_pv'][i]*hexagon['pv'][i]*1000)/1000000) # kWh to MWh?
        
    elif hexagon['cheapest_elec_tech'][i] == 'Wind' and hexagon['theo_turbines'][i] >= 1:
        elec_power_pot.append((hexagon['theo_turbines'][i]*hexagon['wind_output'][i])/1000000) # to MWh?

    else:
        elec_power_pot.append(nan)

hexagon['power_potential'] = elec_power_pot
#%% calculating hydrogen potential in each hexagon
h2_potential = []
# !!! seems like energy density parameter is hardcoded here to convert from energy to H2 mass    
for i in range(len(hexagon)):
    if hexagon['cheapest_elec_tech'][i] == 'PV' and hexagon['theo_pv'][i] > 0:
        h2_potential.append(hexagon['power_potential'][i]*hexagon['pv_ely_ratio'][i]*ely_eff/33.33)
    
    elif hexagon['cheapest_elec_tech'][i] == 'Wind' and hexagon['theo_turbines'][i] > 0:
        h2_potential.append(hexagon['power_potential'][i]*hexagon['wind_ely_ratio'][i]*ely_eff/33.33)

    else:
        h2_potential.append(nan)

hexagon['h2_potential'] = h2_potential
        
#%% calculate cost of expanding grid to each hexagon
    
#Transmission line: Possibility to trasnmit electricity over grid in other hexagons 

if values3['Grid construction'] == True:

    elec_cost_to_connect = []
    for i in range(len(hexagon)):
        if hexagon['grid_dist'][i] != 0:
            #!!! where are these numbers coming from? possibly efficiency of grid transmission
            elec_cost_to_connect.append(hexagon['cheapest_elec_cost'][i]
                                        +((hexagon['grid_dist'][i]*grid_capex/RBF(interest,grid_lifetime))
                                          /(2000*8760*0.95*0.9))) # unsure what conversion these numbers are doing
        else:
            elec_cost_to_connect.append(hexagon['cheapest_elec_cost'][i])

    
    elec_cost_at_grid = []
    for i in range(len(hexagon)):
        if hexagon['grid_dist'][i] == 0:
            elec_cost_at_grid.append(hexagon['cheapest_elec_cost'][i])
    
    
    cheapest_elec_cost_grid = []
    if min(elec_cost_at_grid) == min(elec_cost_to_connect):

        for i in range(len(hexagon)):
            if hexagon['grid_dist'][i] == 0:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], (min(elec_cost_at_grid)+elec_trans_costs)))
                if (min(elec_cost_at_grid)+elec_trans_costs) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'
            else:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], 
                                                   (min(elec_cost_at_grid)
                                                    +elec_trans_costs
                                                    +((hexagon['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime)))))
                #cheapest_elec_cost_grid.append(cheapest_elec_cost[i])
                if (min(elec_cost_at_grid)
                    +elec_trans_costs
                    +((hexagon['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime))) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'
    else:
        for i in range(len(hexagon)):
            if hexagon['grid_dist'][i] == 0:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i],min(elec_cost_to_connect)+elec_trans_costs))
                if (min(elec_cost_to_connect)+elec_trans_costs) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'

            elif i == min(range(len(elec_cost_to_connect)), key=elec_cost_to_connect.__getitem__):
                cheapest_elec_cost_grid.append(cheapest_elec_cost[i])
            else:
                cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], 
                                                   (min(elec_cost_to_connect)
                                                    +elec_trans_costs
                                                    +((hexagon['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime)))))
                if (min(elec_cost_to_connect)
                    +elec_trans_costs+((hexagon['grid_dist'][i]*grid_capex)/RBF(interest,grid_lifetime))) < cheapest_elec_cost[i]:
                    cheapest_elec_tech[i] = 'Grid'
else:
    cheapest_elec_cost_grid = []
    elec_cost_at_grid = []

    for i in range(len(hexagon)):
        if hexagon['grid_dist'][i] == 0:
            elec_cost_at_grid.append(hexagon['cheapest_elec_cost'][i])

    for i in range(len(hexagon)):
        if hexagon['grid_dist'][i] == 0:
            cheapest_elec_cost_grid.append(min(cheapest_elec_cost[i], (min(elec_cost_at_grid)+elec_trans_costs)))
            if (min(elec_cost_at_grid)+elec_trans_costs) < cheapest_elec_cost[i]:
                cheapest_elec_tech[i] = 'Grid'
        else:
            cheapest_elec_cost_grid.append(cheapest_elec_cost[i])
            #cheapest_elec_cost_grid.append(cheapest_elec_cost[i])


hexagon['cheapest_elec_cost'] = cheapest_elec_cost_grid


#for i in range(len(hexagon)):

#    if cheapest_elec_cost_grid[i] != cheapest_elec_cost[i]:
#        cheapest_elec_tech[i] = 'Grid'


hexagon['cheapest_elec_tech'] = cheapest_elec_tech
#%% water cost for each hexagon

h2o_costs_dom_water_bodies = []
h2o_costs_ocean = []
h2o_costs = []

# electricity_demand_h2o_treatment = 0.4                                  #kWh/          https://betterbuildingssolutioncenter.energy.gov/sites/default/files/Primer%20on%20energy%20efficiency%20in%20water%20and%20wastewater%20plants_0.pdf
# electricity_demand_h2o_ocean_treatment = 3.7                            #kWh/m3     //https://www.pnas.org/doi/epdf/10.1073/pnas.1902335116, https://essay.utwente.nl/78100/1/Antonyan%2C%20M.%201817078%20_openbaar.pdf, 
# water_transport_costs = 0.1                                             #€/100km/3
# water_spec_cost = 1.25                                                  #€/m3

electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
electricity_demand_h2o_ocean_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
water_transport_costs = water_data['Water transport cost (euros/100 km/m3)']
water_spec_cost = water_data['Water specific cost (euros/m3)']

for i in range(len(hexagon)):
    h2o_costs_dom_water_bodies.append(((water_spec_cost 
                                        + (water_transport_costs/100)*min(hexagon['waterbody_dist'][i],
                                                                          hexagon['waterway_dist'][i]) 
                                        + electricity_demand_h2o_treatment*(hexagon['cheapest_elec_cost'][i]/1000)*ely_water
                                        )/1000))
# for i in range(len(hexagon)):
    h2o_costs_ocean.append(((water_spec_cost 
                             + (water_transport_costs/100)*hexagon['ocean_dist'][i] 
                             + electricity_demand_h2o_ocean_treatment*(hexagon['cheapest_elec_cost'][i]/1000)*ely_water
                             )/1000))
# for i in range(len(hexagon)):
    h2o_costs.append(min(h2o_costs_dom_water_bodies[i],h2o_costs_ocean[i]))

#%% hydrogen production cost for each hexagon

h2_prod_costs = []

for i in range(len(hexagon)):
    if hexagon['cheapest_elec_tech'][i] == 'PV' and hexagon['theo_pv'][i] >= 1:
        ely_costs = ((((ely_capex+ely_stack_replacement)/RBF(interest,ely_lt))/(flh_pv))*(h2_en_den/ely_eff))*(1 + ely_opex)
        h2_prod_costs.append(((hexagon['cheapest_elec_cost'][i]/1000)* (h2_en_den/ely_eff)) + ely_costs + h2o_costs[i]*ely_water)
    
    elif hexagon['cheapest_elec_tech'][i] == 'Wind' and hexagon['theo_turbines'][i] >= 1:
        ely_costs = ((((ely_capex+ely_stack_replacement)/RBF(interest,ely_lt))/(flh_wind))*(h2_en_den/ely_eff))*(1 + ely_opex)
        h2_prod_costs.append(((hexagon['cheapest_elec_cost'][i]/1000)* (h2_en_den/ely_eff)) + ely_costs + h2o_costs[i]*ely_water)
    
    elif hexagon['cheapest_elec_tech'][i] == 'Grid':
        ely_costs = ((((ely_capex+ely_stack_replacement)/RBF(interest,ely_lt))/(flh_wind))*(h2_en_den/ely_eff))*(1 + ely_opex)
        h2_prod_costs.append(((hexagon['cheapest_elec_cost'][i]/1000)* (h2_en_den/ely_eff)) + ely_costs + h2o_costs[i]*ely_water)        

    else: 
        h2_prod_costs.append(nan)

hexagon['h2_prod_costs'] = h2_prod_costs

#print(hexagon)
#print(max(h2_prod_costs))
#print(min(h2o_costs_dom_water_bodies))

#%% calculate cost of hydrogen state conversion and transportation for demand
#distance to demand: Mombassa coordinates:  39.66007221012109, -4.039286378400124 info: flipped from google maps

distance_dict = {}
total_demand = 0

#%% loop through all demand centers
for d in range(len(demand_center_list)):

    total_demand = total_demand + demand_center_list[d][2]

    lat = demand_center_list[d][0]
    lon = demand_center_list[d][1]
    demand_location = Point(float(lon),float(lat))

    distance_to_demand = []
    
    
    hydrogen_quantity = demand_center_list[d][2]
    h2_costs_incl_conversion = []
    h2_costs_to_demand = []
    road_construction_costs = []
    transport_type = []
    demand_fid = 0
    #%% loop through all hexagons
    for i in range(len(hexagon)):
        # %% calculate distance to demand for each hexagon
        
        #Method takes to long so far
        #supply centre
        #poly = shapely.wkt.loads(str(hexagon['geometry'][i]))
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

        # calculate distance between each hexagon and demand center
        poly = shapely.wkt.loads(str(hexagon['geometry'][i]))
        center = poly.centroid
        #coords_1 = (-4.039286378400124, 39.66007221012109)
        coords_1 = (lat,lon)
        coords_2 = (center.y, center.x)
        dist = geopy.distance.geodesic(coords_1, coords_2).km

        #print(dist)

        distance_to_demand.append(dist)

        #%% label demand location under consideration
        if hexagon['geometry'][i].contains(demand_location) == True:
            demand_fid = i
        # determine elec_cost at demand to determine potential energy costs
        elec_costs_at_demand = float(hexagon['cheapest_elec_cost'][demand_fid])/1000
        #%% calculate cost of constructing a road to each hexagon
        if hexagon['road_dist'][i]==0:
            road_construction_costs.append(0)
        elif hexagon['road_dist'][i]!=0 and hexagon['road_dist'][i]<10:
            road_construction_costs.append(((hexagon['road_dist'][i]*road_capex_short)
                                            /(RBF(interest,road_lifetime)))+(hexagon['road_dist'][i]*road_opex))
        else:
            road_construction_costs.append(((hexagon['road_dist'][i]*road_capex_long)
                                            /(RBF(interest,road_lifetime)))+(hexagon['road_dist'][i]*road_opex))
        #%% calculate cost of transportation and conversion for all hexagons
    for i in range(len(hexagon)):
        #%% calculate cost of meeting hydrogen local demand within the same hexagon
        if i == demand_fid:
            # calculate cost of converting hydrogen to ammonia for local demand (i.e. no transport)
            if demand_center_list[d][3] == 'NH3':
            # !!! where are the 0.03 values coming from? it's the cost of heat in unknown units
                local_conversion_cost = h2_conversion_stand(demand_center_list[d][3]+'_load',
                                                            hydrogen_quantity,
                                                            cheapest_elec_cost[i]/1000,
                                                            0.03,
                                                            interest
                                                            )[2]/hydrogen_quantity
                h2_costs_incl_conversion.append(hexagon['h2_prod_costs'][i]
                                                +local_conversion_cost)
                h2_costs_to_demand.append(hexagon['h2_prod_costs'][i]
                                          +local_conversion_cost)
            else:
                local_conversion_cost = h2_conversion_stand(demand_center_list[d][3],
                                     hydrogen_quantity,
                                     cheapest_elec_cost[i]/1000,
                                     0.03,
                                     interest
                                     )[2]/hydrogen_quantity
                h2_costs_incl_conversion.append(hexagon['h2_prod_costs'][i] + local_conversion_cost)
                h2_costs_to_demand.append(hexagon['h2_prod_costs'][i] + local_conversion_cost)
            transport_type.append('None')
        #%% calculate total cost of hydrogen (production plus transport costs) if pipeline construction is allowed
        # calculate cost of hydrogen production plus transportation if road construction is allowed
        # !!! now that pipeline is a parameter in cheapest_transport_strategy, can condense code further
        # NON-LOCAL DEMAND
        elif (values3['Pipeline construction'] == True
            and values3['Road construction'] == True): 
                demand_state = demand_center_list[d][3]
                if demand_state in ['500 bar','LH2','NH3']:
                    # determine lowest cost transport state in the pipeline
                    transport_cost, transport_state = \
                        cheapest_transport_strategy(demand_state,
                                                       hydrogen_quantity,
                                                       distance_to_demand[i],
                                                       cheapest_elec_cost[i]/1000,
                                                       0.03, # heat costs?
                                                       interest, 
                                                       elec_costs_at_demand, 
                                                       min(cheapest_elec_cost_grid)/1000,
                                                       days_of_storage,
                                                       pipeline = True
                                                       )
                    h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)
                                              +h2_prod_costs[i]
                                              +transport_cost)
                    # determine lowest cost transport state in the pipeline
                    transport_type.append(transport_state)
                    # very long if not statement-- what transport states is this supposed to capture?
                    if transport_state not in ['NH3', 'LOHC', "Small Pipeline",
                                                "Medium Pipeline", "Large Pipeline"]:
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +h2_conversion_stand(transport_state,
                                                                             hydrogen_quantity,
                                                                             cheapest_elec_cost[i]/1000,
                                                                             0.03,interest
                                                                             )[2]/hydrogen_quantity)
                    # cost of conversion if transporting in pipeline
                    elif transport_state in ["Small Pipeline", "Medium Pipeline","Large Pipeline"]:
                            h2_costs_incl_conversion.append(h2_prod_costs[i])
                    else:
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +h2_conversion_stand(transport_state+'_load',
                                                                             hydrogen_quantity,
                                                                             cheapest_elec_cost[i]/1000,
                                                                             0.03,
                                                                             interest
                                                                             )[2]/hydrogen_quantity)
                else:
                    print('{} demand not supported.'.format(demand_state))
        #%% calculate transport cost if road construction isn't allowed but pipeline construction is allowed
        elif (values3['Pipeline construction'] == True
              and values3['Road construction'] != True):
                # if hydrogen is produced roadside
                if hexagon['road_dist'][i]==0:
                    demand_state = demand_center_list[d][3]
                    if demand_state in ['500 bar','LH2','NH3']:
                        distribution_cost, transport_state = \
                            cheapest_transport_strategy(demand_state, 
                                                          hydrogen_quantity, 
                                                          distance_to_demand[i], 
                                                          cheapest_elec_cost[i]/1000, 
                                                          0.03, 
                                                          interest, 
                                                          elec_costs_at_demand, 
                                                          min(cheapest_elec_cost_grid)/1000,
                                                          days_of_storage,
                                                          pipeline=True
                                                          )
                        h2_costs_to_demand.append(h2_prod_costs[i] + distribution_cost)
                        transport_type.append(transport_state)

                        if transport_state not in ['NH3', 'LOHC', "Small Pipeline",
                                                   "Medium Pipeline","Large Pipeline"]:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                            +h2_conversion_stand(transport_state,
                                                                                 hydrogen_quantity,
                                                                                 cheapest_elec_cost[i]/1000,
                                                                                 0.03,interest
                                                                                 )[2]/hydrogen_quantity)
                        elif transport_state in ["Small Pipeline", "Medium Pipeline","Large Pipeline"]:
                            h2_costs_incl_conversion.append(h2_prod_costs[i])
                        else:
                            h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                            +h2_conversion_stand(transport_state+'_load',
                                                                                 hydrogen_quantity,
                                                                                 cheapest_elec_cost[i]/1000,
                                                                                 0.03,interest
                                                                                 )[2]/hydrogen_quantity)
                # if hydrogen is produced off-road
                elif hexagon['road_dist'][i]>0: 
                    demand_state = demand_center_list[d][3]
                    if demand_state in ['500 bar','LH2','NH3']:
                        conversion_cost = h2_conversion_stand(demand_state, 
                                              hydrogen_quantity, 
                                              elec_costs_at_demand, 
                                              0.03, 
                                              interest
                                              )[2]
                        pipeline_dist = sum(storage_costs(demand_state,
                                                      hydrogen_quantity,
                                                      days_of_storage,
                                                      interest), 
                                            pipeline_costs(distance_to_demand[i],
                                             hydrogen_quantity,
                                             min(cheapest_elec_cost_grid)/1000,
                                             interest
                                             )[0],
                                            conversion_cost)
                    
                        
                        transport_type.append(pipeline_costs(distance_to_demand[i],
                                                                 hydrogen_quantity,
                                                                 min(cheapest_elec_cost_grid)/1000,
                                                                 interest)[1])
                        h2_costs_to_demand.append(h2_prod_costs[i]
                                                  +(pipeline_dist/hydrogen_quantity))
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +conversion_cost/hydrogen_quantity)

                    else:
                        print('Negative or non-numeric distance from road.')

                        #h2_costs_to_demand.append((nan))
                        #h2_costs_incl_conversion.append(nan)
        #%% pipeline construction not allowed but road construction allowed
        elif (values3['Road construction'] == True
              and values3['Pipeline construction'] != True):
                demand_state = demand_center_list[d][3]
                if demand_state in ['500 bar','LH2','NH3']:
                    transport_cost, transport_state =\
                        cheapest_transport_strategy(demand_state,
                                              hydrogen_quantity,
                                              distance_to_demand[i],
                                              cheapest_elec_cost[i]/1000,
                                              0.03,
                                              interest,
                                              elec_costs_at_demand,
                                              days_of_storage,
                                              pipeline = False
                                              )
                    h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)
                                              +h2_prod_costs[i]
                                              +transport_cost)
                    transport_type.append(transport_state)
                    # cost when transport state isn't ammonia or LOHC
                    if  transport_state != 'NH3' and transport_state != 'LOHC':
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +h2_conversion_stand(transport_state,
                                                                             hydrogen_quantity,
                                                                             cheapest_elec_cost[i]/1000,
                                                                             0.03,interest
                                                                             )[2]/hydrogen_quantity)
                    elif transport_state in ['NH3','LOHC']: # cost when transport state is ammonia or LOHC
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +h2_conversion_stand(transport_state+'_load',
                                                                             hydrogen_quantity,
                                                                             cheapest_elec_cost[i]/1000,
                                                                             0.03,interest
                                                                             )[2]/hydrogen_quantity)

       #%% pipeline construction and road construction not allowed         
        elif (values3['Pipeline construction'] != True
              and values3['Road construction'] != True):
            if hexagon['road_dist'][i]==0.: # hydrogen produced roadside
                demand_state = demand_center_list[d][3]
                if demand_state in ['500 bar', 'LH2', 'NH3']:
                    transport_cost, transport_state = cheapest_transport_strategy(demand_state,
                                                                           hydrogen_quantity,
                                                                           distance_to_demand[i],
                                                                           cheapest_elec_cost[i]/1000,
                                                                           0.03,
                                                                           interest,
                                                                           elec_costs_at_demand,
                                                                           days_of_storage,
                                                                           pipeline=False
                                                                           )
                    h2_costs_to_demand.append(h2_prod_costs[i]
                                              +transport_cost)
                    transport_type.append(transport_state)
                    if  transport_state != 'NH3' and transport_state != 'LOHC':
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +h2_conversion_stand(transport_state,
                                                                             hydrogen_quantity,
                                                                             cheapest_elec_cost[i]/1000,
                                                                             0.03,interest
                                                                             )[2]/hydrogen_quantity)
                    elif transport_state in ['NH3','LOHC']: # cost when transport state is ammonia or LOHC
                        h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                        +h2_conversion_stand(transport_state+'_load',
                                                                             hydrogen_quantity,
                                                                             cheapest_elec_cost[i]/1000,
                                                                             0.03,
                                                                             interest
                                                                             )[2]/hydrogen_quantity)
            elif hexagon['road_dist'][i]>0: 
                h2_costs_to_demand.append((nan))
                h2_costs_incl_conversion.append(nan)
                transport_type.append(nan)
            else:
                print('Negative or non-numeric distance from road.')
    #%% variables to save for each demand scenario
    #h2_costs_to_demand = [round(num, 1) for num in h2_costs_to_demand]

    # save costs for meeting demand for each demand center
    hexagon['h2_costs_to_demand' + str(d)] = h2_costs_to_demand
    hexagon['h2_costs_incl_conv' + str(d)] = h2_costs_incl_conversion
    hexagon['transport_type' + str(d)] = transport_type


#print(distance_to_demand)

#%%
NH3_costs_to_demand = []

if values3['ammonia_map'] == True:
    for i in range(len(hexagon)):
        # !!! add weight of h2 in NH3 as a varible, maybe in excel?
        h2_weight_NH3 = 0.178       
        NH3_costs_to_demand.append((h2_costs_to_demand[i]+(h2_conversion_stand('NH3_unload',hydrogen_quantity,elec_costs_at_demand/1000,0.03,interest)[2]/hydrogen_quantity))*h2_weight_NH3)

    hexagon['NH3_costs_to_demand'] = NH3_costs_to_demand

#See results and prepare plots
#print(hexagon)
hex_sort = {}
for d in range(len(demand_center_list)):

    hex_sort[d] = hexagon.sort_values(by=['h2_costs_to_demand' + str(d)]).head(100)

#print(hex_sort)
#gpd.hex_sort.to_excel()
#print(values[6])

h2_costs_to_demand = []

for i in range(len(hexagon)):
    options = []

    for d in range(len(demand_center_list)):
        options.append((hexagon['h2_costs_to_demand' + str(d)][i]))
    
    h2_costs_to_demand.append(min(options))
        

hexagon['h2_costs_to_demand'] = h2_costs_to_demand
# create excel writer object
#minimal_cost = pd.ExcelWriter('output.xlsx')
# write dataframe to excel
#hex_sort.to_excel(minimal_cost)

#minimal_cost.save()

#%% plotting

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

hexagon['h2_cap_after_supply'] = hexagon['h2_potential']

for d in range(len(demand_center_list)):
    supply = 0
    nr = 0
    supply_list = []
    index_list = list(hex_sort[d]['h2_costs_to_demand' +str(d)].index.values)

    while supply < float(demand_center_list[d][2]):
        if hexagon['h2_potential'][index_list[nr]]*1000000 > 0:
            if hexagon['h2_potential'][index_list[nr]]*1000000 >= (float(demand_center_list[d][2]) - supply):
                coverage = (float(demand_center_list[d][2]) - supply)
                supply = supply + coverage
                hexagon['h2_potential'][index_list[nr]] = hexagon['h2_potential'][index_list[nr]] - (coverage/1000000)
                supply_list.append(index_list[nr])
                nr = nr + 1
            else:
                coverage = hexagon['h2_potential'][index_list[nr]]*1000000
                supply = supply + hexagon['h2_potential'][index_list[nr]]*1000000
                hexagon['h2_potential'][index_list[nr]] = 0 
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
            poly = shapely.wkt.loads(str(hexagon['geometry'][s]))
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
        hexagon.plot(column= map_list[i], legend = True, cmap= colour_dict[map_list[i]], ax=ax)

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
    #hexagon = hexagon.sort_values(by=['h2_costs_to_demand' + str(d)])
    #hexagon['h2_prod_costs'].to_excel(output ,sheet_name = str(d), startcol=0)
    #hexagon['transport_type'+ str(d)].to_excel(output, sheet_name = str(d), startcol=2, index=False)
    #hexagon['h2_costs_incl_conv' + str(d)].to_excel(output, sheet_name = str(d), startcol=3, index=False)
    #hexagon['h2_costs_to_demand' + str(d)].to_excel(output, sheet_name = str(d), startcol=4, index=False)

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
