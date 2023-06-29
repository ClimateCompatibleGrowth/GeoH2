# GEOH2
Geospatial analysis of hydrogen production costs

GEOH2 calculates the locational cost of green hydrogen production, storage, transport, and conversion to meet demand in a specified location. These costs can be compared to current or projected prices for energy and chemical feedstocks in the region to assess the competitiveness of green hydrogen. Currently, different end-uses, such as fertilizer production, export shipping, and steel production, are not modeled.

The model outputs the levelized cost of hydrogen (LCOH) at the demand location including production, storage, transport, and conversion costs. 

In the code provided, the specific use case of Kenya is investigated. As the code is written in a generalized way, it is possible to analyse all sorts of regions.


# Installation instructions
## Clone the repository
First, clone the GEOH2 repository using `git`. 

`/some/other/path % cd /some/path/without/spaces`

`/some/path/without/spaces % git clone https://github.com/ClimateCompatibleGrowth/GeoH2.git`

## Install Python dependencies
The python package requirements are in the `environment.yaml` file. You can install these requirements in a new environment using `conda` package and environment manager (available for installation [here](https://docs.conda.io/en/latest/miniconda.html)): 

` .../GEOH2 % conda env create -f environment.yaml`

Then activate this new environment using

`.../GEOH2 % conda activate geoh2`

# Preparing input data

## Hexagons 
To analyse a different area of interest, the input hexagon file needs to be changed, but needs to follow the logic of the one provided. An explanation of how to create a H3-Hexagon file can be found in the following repo:

https://github.com/carderne/ccg-spider

The hexagon file needs to filled with the following attributes:

  - waterbody_dist: Distance to selected waterbodies in area of interest
  - waterway_dist: Distance to selected waterways in area of interest
  - ocean_dist: Distance to ocean coastline
  - grid_dist: Distance to transmission network
  - road_dist: Distance to road network
  - theo_pv: Theoretical potential of standarized PV plants       --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes
  - theo_wind: Theoretical potential of standarized wind turbines     --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes
  
Once you have created a hexagon file with these features, save it in the `Data` folder.
  
## Input parameter Excel files

Required input parameters include the spatial area of interest, total annual demand for hydrogen, and prices and cost of capital for infrastructure investments. These values can be either current values or projected values for a single snapshot in time. The parameter values for running the model can be specified in a set of Excel files in the Parameters folder.

- **Demand parameters:** `demand_parameters.xlsx` includes a list of demand centers. For each demand center, its lat-lon location, annual demand, and hydrogen state for that demand must be specified. If multiple forms are hydrogen are demanded in one location, differentiate the demand center name (e.g. Nairobi LH2 and Nairobi NH3) to avoid problems from duplicate demand center names.

- **Country parameters:** `country_parameters.xlsx` includes country- and technology-specific interest rates, heat and electricity costs, and asset lifetimes.
    - Interest rates should be expressed as a decimal, e.g. 5% as 0.05.
    - Asset lifetimes should be in years.
    
- **Basic H2 plant** in this folder, there are several csv files containing the global parameters for optimizing the plant design. All power units are MW and all energy units are MWh. For more information on these parameters, refer to the [PyPSA documentation](https://pypsa.readthedocs.io/en/latest/components.html).

- **Weather parameters:** `weather_parameters.xlsx` includes the locational and time range to download historical weather data from the ERA5 reanalysis dataset to calculate wind and solar generation potential. At least a year is recommended for the weather time range to capture seasonal variation in renewable potential, but longer time ranges will increase the computation time for optimizing the design of a hydrogen plant.
    - Note: The end date is not inclusive, so if you want to download weather data for all of 2022, you should input 2023/01/01 as the end date. 
    - Note: The dates need to be presented in the following format: YYYY/MM/DD. If this is auto-corrected to a different format by excel, be sure to manually re-adjust it to be this format, or the download script will not recognise the date properly.

- **Technology parameters:** ` technology_parameters.xlsx` includes water parameters, road infrastructure parameters, and whether road and hydrogen pipeline construction is allowed.

- **Pipeline parameters:** `pipeline_parameters.xlsx` includes the price, capacity, and lifetime data for different sizes of hydrogen pipeline.

- **Transport parameters:** `transport_parameters.xlsx` includes the parameters related to road transport of hydrogen, including truck speed, cost, lifetime, and capacity.

- **Conversion parameters**: `conversion_parameters.xlsx` includes parameters related to converting between states of hydrogen.

# Assign countries
Assigns country-specific interest rates, technology lifetimes, and heat and electricity prices from `country_parameters.xlsx` to different hexagons based on their country.

You can run this script by entering the following command in your terminal:

`.../GEOH2 % python assign_country.py`

This script should run in a few seconds or less.

# Download weather data

Download the weather data to calculate wind and solar potential for the geographical area and time period of interest using the `get_weather_data.py` script. Ensure that you have specified your weather parameters in `Parameters/weather_parameters.xlsx` 
  
This script downloads the relevant historical weather data from the ERA-5 reanalysis dataset using [Atlite](https://atlite.readthedocs.io/en/latest/) to create a cutout. For this process to work, you need to register and set up your CDS API key as described on the [Climate Data Store website](https://cds.climate.copernicus.eu/api-how-to).

You can run this script by entering the following command in your terminal:

`.../GEOH2 % python get_weather_data.py`

Depending the length of your time period and size of your area as well as your internet connection, this script could take from a few minutes to several hours to run.

The downloaded data is saved in the Cutouts folder. Ensure that you have space on your computer to store the data, which can be several GB depending on the length of your time period and size of your area.

# Optimize conversion and transportation

Calculates the cost of the optimal hydrogen transportation and conversion strategy from each hexagon to each demand center using both pipelines and road transport using parameters from `technology_parameters.xlsx`, `demand_parameters.xlsx`, and `country_parameters.xlsx`.

You can run this script by entering the following command in your terminal: 

`.../GEOH2 % python optimize_transport_and_conversion.py`

This script should take a few minutes to run, depending on the size of the geographic area you are analyzing.

# Optimize green hydrogen plant design

Design green hydrogen plant to meet the hydrogen demand profile for each demand center for each transportation method to each demand center using the `optimize_hydrogen_plant.py` script. Ensure that you have specified your hydrogen plant parameters in the CSV files in the `Parameters/Basic_H2_plant` folder, your investment parameters in `Parameters/investment_parameters.xlsx`, and your demand centers in `Parameters/demand_parameters.xlsx`.

## Installing a solver
To use this script, you will need a solver installed on your computer. You can use any solver that works with [PyPSA](https://pypsa.readthedocs.io/en/latest/installation.html), such as [Cbc](https://github.com/coin-or/Cbc), a free, open-source solver, or [Gurobi](https://www.gurobi.com/), a commerical solver with free academic licenses available. Install your solver of choice following the instructions for use with Python and your operating system in the solver's documentation.

## Running script

You can run this script by entering the following command in your terminal:

`.../GEOH2 % python optimize_hydrogen_plant.py`

This script will take several minutes to several hours to run, depending on the size of the geographic area you are analyzing, the length of your weather data time period, and the number of demand centers.

# Water costs

Calculate water costs from the ocean and freshwater bodies for hydrogen production in each hexagon using `Parameters/technology_parameters.xlsx` and `Parameters/country_parameters.xlsx`.

You can run this script by entering the following command in your terminal:

`.../GEOH2 % python water_cost.py`

This script will take a few seconds to run.

# Total hydrogen costs
Combines the results `optimize_transport_and_conversion.py`,`optimize_hydrogen_plant.py`, and `water_cost.py` to find the lowest-cost method of producing, transporting, and converting hydrogen for each demand center.

You can run this script by entering the following command in your terminal:

`.../GEOH2 % python total_hydrogen_cost.py`

This script will take a few seconds to run.

# Visualizing results

Visualizes the spatial variation in different costs per kilogram of hydrogen.

You can run this script by entering the following command in your terminal: 

`.../GEOH2 % python map_costs.py`

This script will take a few seconds to run. 

If you encounter the following error in visualisation, you may have a package conflict with geopandas, matplotlib and cartopy: 
`AttributeError: 'GeoAxesSubplot' object has no attribute '_autoscaleXon'.` 
To get around this, try creating a new environment using the `environment-plot.yml` file using the same command presented above. This will create an environment named `geoh2-plot` which you can activate and use to run the visualisation script.

# Limitations

This model considers only greenfield wind and solar plants for hydrogen production. Therefore it does not consider using grid electricity or existing generation for hydrogen production.

While the design of the green hydrogen plant is convex and therefore guarenteed to find the global optimum solution if it exists, the selection of the trucking strategy is greedy to avoid the long computation times and potential computational intractability associated with a mixed-integer optimization problem.

Currently, only land transport is considered in the model. To calculate the cost of hydrogen production for export, any additional costs for conversion and transport via ship or undersea pipeline must be added in post-processing.

The availability of water for electrolysis is not limited in regions that could potentially face drought, and a single prices for freshwater and ocean water are used throughout the modeled area.

