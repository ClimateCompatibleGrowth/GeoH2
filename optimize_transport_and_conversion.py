#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 16:52:59 2023

@author: Claire Halloran, University of Oxford
Contains code originally written by Leander MÃ¼ller, RWTH Aachen University

Calculates the cost-optimal hydrogen transportation strategy to the nearest demand center.
Also includes the optimal transportation strategy for both pipelines and trucking
to input as a demand profile into hydrogen plant optimization


Calculate cost of road transport and demand profile based on optimal schedule of shipments
Calculate cost of pipeline transport and demand profile based on optimal size



"""

# from calendar import c, prcal
# from turtle import color, distance
import geopandas as gpd
# from matplotlib import markers
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt                 #see: https://geopandas.org/en/stable/docs/user_guide/mapping.html for plotting
from functions import NPV, cheapest_transport_strategy, h2_conversion_stand, pipeline_costs
# from cmath import nan, pi
from shapely.geometry import Point
import shapely.geometry
import shapely.wkt
import geopy.distance
# import PySimpleGUI as sg 
# import math
# from xlsxwriter import Workbook

#%% Data Input
# Hexagon file
hexagon = gpd.read_file('Data/hex_final.geojson')

# Excel file with technology parameters
technology_parameters = "Parameters/technology_parameters.xlsx"
demand_parameters = 'Parameters/demand_parameters.xlsx'
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
demand_center_list = pd.read_excel(demand_parameters,
                                   sheet_name='Demand centers',
                                   index_col='Demand center',
                                   )

interest = global_data['Interest rate']
# !!! do we want to include spatial electricity costs? for now just using this number picked at random
cheapest_elec_cost = global_data['Electricity cost (euros/kWh)']
cheapest_elec_cost_grid = global_data['Electricity cost (euros/kWh)']
pipeline_construction = global_data['Pipeline construction allowed']
road_construction = global_data['Road construction allowed']


grid_capex = infra_data.at['Grid','CAPEX']               
elec_trans_costs = infra_data.at['Grid','OPEX']                #â¬/MWh PLatzhalter see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8661478/pdf/main.pdf
grid_lifetime = infra_data.at['Grid','Lifetime']                  #years PLATZHALTER

road_capex_long = infra_data.at['Long road','CAPEX']            #â¬/km from John Hine, converted to Euro (Assumed earth to paved road)
road_capex_short = infra_data.at['Short road','CAPEX']         #â¬/km for raods < 10 km, from John Hine, converted to Euro (Assumed earth to paved road)
road_opex = infra_data.at['Short road','OPEX']                 #â¬/km/year from John Hine, converted to Euro (Assumed earth to paved road)
road_lifetime = infra_data.at['Short road','Lifetime']             #years, assumption
days_of_storage = 0

#%% calculate cost of hydrogen state conversion and transportation for demand
# loop through all demand centers-- maybe each hexagon should only calculate this for the nearest demand 
# center, because on the continential scale this is a lot of calculations for distant demand centers
for d in demand_center_list.index:
    # !!! should these empty lists be numpy arrays of the length of the hexagon dataframe?
    demand_location = Point(demand_center_list.loc[d,'Lat [deg]'], demand_center_list.loc[d,'Lon [deg]'])
    distance_to_demand = np.zeros(len(hexagon))

    hydrogen_quantity = demand_center_list.loc[d,'Annual demand [kg/a]']
    h2_costs_incl_conversion = np.zeros(len(hexagon))
    h2_costs_to_demand = np.zeros(len(hexagon))
    road_construction_costs = np.zeros(len(hexagon))
    transport_type = np.zeros(len(hexagon))
    demand_fid = 0
    demand_state = demand_center_list.loc[d,'Demand state']

    #%% loop through all hexagons
    for i in range(len(hexagon)):
        # %% calculate distance to demand for each hexagon
        # calculate distance between each hexagon and demand center
        poly = shapely.wkt.loads(str(hexagon['geometry'][i]))
        center = poly.centroid
        demand_coords = (demand_center_list.loc[d,'Lat [deg]'], demand_center_list.loc[d,'Lon [deg]'])
        hexagon_coords = (center.y, center.x)
        dist = geopy.distance.geodesic(demand_coords, hexagon_coords).km
        
        distance_to_demand[i] = dist

        #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
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
                                            /(NPV(interest,road_lifetime)))+(hexagon['road_dist'][i]*road_opex))
        else:
            road_construction_costs.append(((hexagon['road_dist'][i]*road_capex_long)
                                            /(NPV(interest,road_lifetime)))+(hexagon['road_dist'][i]*road_opex))
        #%% calculate cost of transportation and conversion for all hexagons
    for i in range(len(hexagon)):
        #%% calculate cost of meeting hydrogen local demand within the same hexagon
        if i == demand_fid:
            # calculate cost of converting hydrogen to ammonia for local demand (i.e. no transport)
            if demand_state == 'NH3':
            # !!! where are the 0.03 values coming from? it's the cost of heat in unknown units
                local_conversion_cost = h2_conversion_stand(demand_state+'_load',
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
                local_conversion_cost = h2_conversion_stand(demand_state,
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
        # when road construction allowed
        # elif (values3['Pipeline construction'] == True
        #     and values3['Road construction'] == True): 
        elif road_construction == True:
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
                                                   pipeline = pipeline_construction
                                                   )

                # determine lowest cost transport state in the pipeline
                transport_type.append(transport_state)
                if transport_state not in ["Small Pipeline", "Medium Pipeline","Large Pipeline"]:
                    #!!! should this only be included if road transport is cheapest?
                    h2_costs_to_demand.append((road_construction_costs[i]/hydrogen_quantity)
                                              +h2_prod_costs[i]
                                              +transport_cost)
                else:
                    h2_costs_to_demand.append(h2_prod_costs[i] + transport_cost)
                if transport_state in ['500 bar','LH2']:
                    h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                    +h2_conversion_stand(transport_state,
                                                                         hydrogen_quantity,
                                                                         cheapest_elec_cost[i]/1000,
                                                                         0.03,
                                                                         interest
                                                                         )[2]/hydrogen_quantity)
                # cost of conversion if transporting in pipeline
                elif transport_state in ["Small Pipeline", "Medium Pipeline","Large Pipeline"]:
                    # no converstion costs apply for pipelines
                    h2_costs_incl_conversion.append(h2_prod_costs[i])
                else:
                    # need to specify loading transport state for conversion costs
                    h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                    +h2_conversion_stand(transport_state+'_load',
                                                                         hydrogen_quantity,
                                                                         cheapest_elec_cost[i]/1000,
                                                                         0.03,
                                                                         interest
                                                                         )[2]/hydrogen_quantity)
            else:
                raise NotImplementedError('{} demand not supported.'.format(demand_state))
        #%% calculate transport cost if road construction isn't allowed but pipeline construction is allowed
        # road construction not allowed
        elif road_construction == False:
        # elif (values3['Pipeline construction'] == True
        #       and values3['Road construction'] != True):
                # if hydrogen is produced roadside
                if hexagon['road_dist'][i]==0:
                    if demand_state in ['500 bar','LH2','NH3']: # why ensuring demand state is implemented
                        transport_cost, transport_state = \
                            cheapest_transport_strategy(demand_state, 
                                                          hydrogen_quantity, 
                                                          distance_to_demand[i], 
                                                          cheapest_elec_cost[i]/1000, 
                                                          0.03, 
                                                          interest, 
                                                          elec_costs_at_demand, 
                                                          min(cheapest_elec_cost_grid)/1000,
                                                          days_of_storage,
                                                          pipeline = pipeline_construction
                                                          )
                        h2_costs_to_demand.append(h2_prod_costs[i] + transport_cost)
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
                # if hydrogen is produced off-road, only pipeline options are available
                elif hexagon['road_dist'][i]>0: 
                    if pipeline_construction== True:
                        if demand_state in ['500 bar','LH2','NH3']:
                            conversion_cost = h2_conversion_stand(demand_state, 
                                                  hydrogen_quantity, 
                                                  elec_costs_at_demand, 
                                                  0.03, 
                                                  interest
                                                  )[2]
                            pipeline_dist = sum([pipeline_costs(distance_to_demand[i],
                                                 hydrogen_quantity,
                                                 min(cheapest_elec_cost_grid)/1000,
                                                 interest
                                                 )[0],
                                                conversion_cost])
                        
                            
                            transport_type.append(pipeline_costs(distance_to_demand[i],
                                                                     hydrogen_quantity,
                                                                     min(cheapest_elec_cost_grid)/1000,
                                                                     interest)[1])
                            h2_costs_to_demand.append(h2_prod_costs[i]
                                                      +(pipeline_dist/hydrogen_quantity))
                            h2_costs_incl_conversion.append(h2_prod_costs[i]
                                                            +conversion_cost/hydrogen_quantity)
                    else:
                        h2_costs_to_demand.append((np.nan))
                        h2_costs_incl_conversion.append(np.nan)
                        transport_type.append(np.nan)
                else:
                    print('Negative or non-numeric distance from road.')

                        #h2_costs_to_demand.append((nan))
                        #h2_costs_incl_conversion.append(nan)

    #%% variables to save for each demand scenario
    #h2_costs_to_demand = [round(num, 1) for num in h2_costs_to_demand]

    # save costs for meeting demand for each demand center
    hexagon['h2_costs_to_demand' + str(d)] = h2_costs_to_demand
    hexagon['h2_costs_incl_conv' + str(d)] = h2_costs_incl_conversion
    hexagon['transport_type' + str(d)] = transport_type

