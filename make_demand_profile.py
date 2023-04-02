#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 10:06:44 2023

@author: Claire Halloran, University of Oxford

Supported by the Climate Compatible Growth Programme

This script calculates an hourly demand schedule based on quantity of hydrogen demanded
for different transport methods.

"""
import pandas as pd

def demand_schedule(quantity, transport_state, transport_excel_path,
                             weather_excel_path):
    '''
    calculates hourly hydrogen demand for truck shipment.

    Parameters
    ----------
    quantity : float
        annual amount of hydrogen to transport.
    transport_state : string
        state hydrogen is transported in, one of '500 bar', 'LH2', 'LOHC', or 'NH3'.
    transport_excel_path : string
        path to transport_parameters.xlsx file
    weather_excel_path : string
        path to transport_parameters.xlsx file
            
    Returns
    -------
    trucking_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for hydrogen trucking.
    '''
    transport_parameters = pd.read_excel(transport_excel_path,
                                         sheet_name = transport_state,
                                         index_col = 'Parameter'
                                         ).squeeze('columns')
    weather_parameters = pd.read_excel(weather_excel_path,
                                       index_col = 'Parameters',
                                       ).squeeze('columns')
    truck_capacity = transport_parameters['Net capacity (kg H2)']
    start_date = weather_parameters['Start date']
    end_date = weather_parameters['End date (not inclusive)']
    
    # schedule for trucking
    annual_deliveries = quantity/truck_capacity
    quantity_per_delivery = quantity/annual_deliveries
    index = pd.date_range(start_date, end_date, periods=annual_deliveries)
    trucking_demand_schedule = pd.Series(quantity_per_delivery, index=index)
    trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)
    
    # schedule for pipeline
    index = pd.date_range(start_date, end_date, freq = 'H')
    pipeline_hourly_quantity = quantity/index.size
    pipeline_hourly_demand_schedule = pd.Series(pipeline_hourly_quantity, index=index)

    return trucking_hourly_demand_schedule, pipeline_hourly_demand_schedule

transport_excel_path = "Parameters/transport_parameters.xlsx"
weather_excel_path = "Parameters/weather_parameters.xlsx"
demand_excel_path = 'Parameters/demand_parameters.xlsx'

demand_parameters = pd.read_excel(demand_excel_path,
                                  index_col='Demand center',
                                  )
# calculate for each demand center
demand_centers = demand_parameters.index

for location in demand_centers:
    trucking_demand, pipeline_demand = demand_schedule(demand_parameters.loc[location,'Annual demand [kg/a]'],
                                               demand_parameters.loc[location,'Demand state'],
                                               transport_excel_path, weather_excel_path)
    trucking_demand.to_excel(f'Resources/{location} trucking demand.xlsx')
    pipeline_demand.to_excel(f'Resources/{location} pipeline demand.xlsx')
