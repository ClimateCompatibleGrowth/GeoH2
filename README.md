# GEOH2
Geospatial analysis of hydrogen production costs

GEOH2 helps to identify regional hydrogen production costs. In the code provided, the specific use case of Kenya is investiagte. As the code is written in a generalized way, it is possible to analyse all sorts of regions. For this, adjustments to the input data file, which can be found in the folder "Data", have to be made. Other techno-economic parameter changes can be done in the 'technology_parameters.xlsx' file or directly in the code.

# Analysis of desired region

To analyse a different area of interest, the input hexagon file needs to be changed, but needs to follow the logic of the one provided. An explanaition how a H3-Hexagon file can be created can be found in the following repo:

https://github.com/carderne/ccg-spider

The hexagon file needs to filled with the following attributes:

  - waterbody_dist: Distance to selected waterbodies in area of interest
  - waterway_dist: Distance to selected waterways in area of interest
  - ocean_dist: Distance to ocean coastline 
  
  - grid_dist: Distance to transmissin network
  
  - road_dist: Distance to road network
  
  - wind: Average windspeed at 100 meter [m/s]
  - pv: PV output potential, daily average [kWh/kWp/day]
  
  - theo_pv: Theoretical Potential of standarized PV plants       --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes
  - theo_wind: Theoretical Potential of standarized PV plants     --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes
  
  
  
  
  
