# GEOH2
Geospatial analysis of hydrogen production costs

GEOH2 helps to identify regional hydrogen production costs. In the code provided, the specific use case of Kenya is investiagte. As the code is written in a generalized way, it is possible to analyse all sorts of regions. For this, adjustments to the input data file, which can be found in the folder "Data", have to be made. Other techno-economic parameter changes can be done in the 'technology_parameters.xlsx' file or directly in the code.

# Installation instructions
## Clone the repository
First, clone the GEOH2 repository using `git`. 

`/some/other/path % cd /some/path/without/spaces
`/some/path/without/spaces % git clone https://github.com/leandermue/GEOH2.git

## Install Python dependencies
The python package requirements are in the `environment.yaml` file. You can install these requirements in a new environment using `conda` package and environment manager: 
` .../GEOH2 % conda env create -f envs/environment.yaml`
Then activate this new environment using
`.../GEOH2 % conda activate geoh2`


Note: `gurobipy` package is only available from PyPI through pip, so create the GEOH2 environment and then run 
`.../GEOH2 pip install gurobipy`

# Analysis of desired region
To analyse a different area of interest, the input hexagon file needs to be changed, but needs to follow the logic of the one provided. An explanation how a H3-Hexagon file can be created can be found in the following repo:

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
  
  
  
  
  
