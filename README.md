# GEOH2
Geospatial analysis of hydrogen production costs

GEOH2 calculates the locational cost of green hydrogen production, storage, transport, and conversion to meet demand in a specified location. These costs can be compared to current or projected prices for energy and chemical feedstocks in the region to assess the competitiveness of green hydrogen. Currently, different end-uses, such as fertilizer production, export shipping, and steel production, are not modeled.

Required input parameters include the spatial area of interest, total annual demand for hydrogen, and prices and cost of capital for infrastructure investments. These values can be either current values or projected values for a single snapshot in time. The parameter values for running the model can be specified in a set of Excel files in the Parameters folder.

- **Demand parameters:** `demand_parameters.xlsx` includes a list of demand centers. For each demand center, its lat-lon location, annual demand, and hydrogen state for that demand must be specified. If multiple forms are hydrogen are demanded in one location, differentiate the demand center name (e.g. Nairobi LH2 and Nairobi NH3) to avoid problems from duplicate demand center names.

- ** Investment parameters:** `investment_parameters.xlsx` includes country- and technology-specific interest rates and asset lifetimes.
    - Interest rates should be expressed as a decimal, e.g. 5% as 0.05.
    - Asset lifetimes should be in years.
    
- **Basic H2 plant** in this folder, there are several csv files containing the global parameters for optimizing the plant design. All power units are MW and all energy units are MWh.

- **Weather parameters:** `weather_parameters.xlsx` includes the locational and time range to download historical weather data from the ERA5 reanalysis dataset to calculate wind and solar generation potential. At least a year is recommended for the weather time range to capture seasonal variation in renewable potential, but longer time ranges will increase the computation time for optimizing the design of a hydrogen plant.
    - Note that the end date is not inclusive, so if you want to download weather data for all of 2022, you should input 01-01-2023 as the end date.

- **Technology parameters:** ` technology_parameters.xlsx` includes the price of water, interest rate (i.e. cost of capital), and road infrastructure price. *Note: interest rates will soon be moved to their own sheet for country-specific and potentially technology-specific values.*
    - Sheets to ignore (to be deprecated): Electricity, Electrolysis, Wind Turbine, and interest rate parameter on Global sheet.
    
- **Pipeline parameters:** `pipeline_parameters.xlsx` includes the price, capacity, and lifetime data for different sizes of hydrogen pipeline.

- **Storage parameters:** `storage_parameters.xlsx` includes the price and lifetime of hydrogen storage. To be deprecated.

- **Transport parameters:** `transport_parameters.xlsx` includes the parameters related to road transport of hydrogen, including truck speed, cost, lifetime, and capacity.


The model outputs the levelized cost of hydrogen (LCOH) at the demand location including production, storage, transport, and conversion costs. 

In the code provided, the specific use case of Kenya is investigated. As the code is written in a generalized way, it is possible to analyse all sorts of regions.


# Installation instructions
## Clone the repository
First, clone the GEOH2 repository using `git`. 

`/some/other/path % cd /some/path/without/spaces`
`/some/path/without/spaces % git clone https://github.com/leandermue/GEOH2.git`

## Install Python dependencies
The python package requirements are in the `environment.yaml` file. You can install these requirements in a new environment using `conda` package and environment manager (available for installation [here](https://docs.conda.io/en/latest/miniconda.html)): 

` .../GEOH2 % conda env create -f envs/environment.yaml`

Then activate this new environment using

`.../GEOH2 % conda activate geoh2`

# Analysis of desired region
To analyse a different area of interest, the input hexagon file needs to be changed, but needs to follow the logic of the one provided. An explanation how a H3-Hexagon file can be created can be found in the following repo:

https://github.com/carderne/ccg-spider

The hexagon file needs to filled with the following attributes:

  - waterbody_dist: Distance to selected waterbodies in area of interest
  - waterway_dist: Distance to selected waterways in area of interest
  - ocean_dist: Distance to ocean coastline 
  
  - grid_dist: Distance to transmission network
  
  - road_dist: Distance to road network
  
  - wind: Average windspeed at 100 meter [m/s]
  - pv: PV output potential, daily average [kWh/kWp/day]
  
  - theo_pv: Theoretical Potential of standarized PV plants       --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes
  - theo_wind: Theoretical Potential of standarized PV plants     --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes
  
  
  
  
  
