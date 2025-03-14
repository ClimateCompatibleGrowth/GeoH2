>[!NOTE]
>We now have a new repository for you to use! The new [GeoX](https://github.com/ClimateCompatibleGrowth/GeoX.git) repository contains both this repository and the [GeoNH3](https://github.com/ClimateCompatibleGrowth/GeoNH3.git) repository.


# GEOH2
**Geospatial analysis of hydrogen production costs**

GEOH2 calculates the locational cost of green hydrogen production, storage, transport, and conversion to meet demand in a specified location. These costs can be compared to current or projected prices for energy and chemical feedstocks in the region to assess the competitiveness of green hydrogen. Currently, different end-uses, such as fertilizer production, export shipping, and steel production, are not modeled.

The model outputs the levelized cost of hydrogen (LCOH) at the demand location including production, storage, transport, and conversion costs. 

In the code provided, the specific use case of Namibia is investigated. 
Parameter references for this case are attached.
However, as the code is written in a generalized way, it is possible to analyse all sorts of regions.

GeoH2 builds upon a preliminary code iteration produced by Leander Müller, available under a CC-BY-4.0 licence: [https://github.com/leandermue/GEOH2](https://github.com/leandermue/GEOH2).
It also integrates code produced by Nick Salmon under an MIT licence: 
[https://github.com/nsalmon11/LCOH_Optimisation](https://github.com/nsalmon11/LCOH_Optimisation)
___

# Setup instructions

## Clone the repository
First, clone the GeoH2 repository using `git`. 

`... % git clone https://github.com/ClimateCompatibleGrowth/GeoH2.git`

## Environment setup
The python package requirements are in the `environment.yaml` file. You can install these requirements in a new environment using `mamba` package and environment manager (installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)): 

` .../GEOH2 % mamba env create -f environment.yaml`

Then activate this new environment using

`.../GEOH2 % mamba activate geoh2`

## CDS API setup
The `get_weather_data` rule downloads the relevant historical weather data from the ERA-5 reanalysis dataset using [Atlite](https://atlite.readthedocs.io/en/latest/) to create a cutout. For this process to work, you need to register and set up your CDS API key as described on the [Climate Data Store website](https://cds.climate.copernicus.eu/api-how-to).

**Note:** Ensure the API key and URL are affiliated with CDS-Beta.

## Solver setup
For the `optimize_hydrogen_plant` rule to work, you will need a solver installed on your computer. You can use any solver that works with [PyPSA](https://pypsa.readthedocs.io/en/latest/installation.html), such as [Cbc](https://github.com/coin-or/Cbc), a free, open-source solver, or [Gurobi](https://www.gurobi.com/), a commerical solver with free academic licenses available. Install your solver of choice following the instructions for use with Python and your operating system in the solver's documentation. 

In `Scripts/optimize_hydrogen_plant.py` line 160, the solver is set to `gurobi`. This must be changed if you choose to use a different solver.

**Note**: Snakemake uses Cbc, which will be installed upon environment setup. To check, activate your environment and enter `mamba list` in your terminal for the environment's list of packages.
___

# Preparing input data

## Hexagons 
To analyse a different area of interest, the input hexagon file needs to be changed, but needs to follow the logic of the one provided. 

A full walkthrough on all the tools to create these hexagons are in the [GeoH2-data-prep](https://github.com/ClimateCompatibleGrowth/GeoH2-data-prep) repo.

An explanation of how to create a H3-Hexagon file can be found in the following repo:

https://github.com/carderne/ccg-spider

The hexagon file needs to filled with the following attributes:

  - waterbody_dist: Distance to selected waterbodies in area of interest
  - waterway_dist: Distance to selected waterways in area of interest
  - ocean_dist: Distance to ocean coastline
  - grid_dist: Distance to transmission network
  - road_dist: Distance to road network
  - theo_pv: Theoretical PV potential --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes. Note that this value should be in MW.
  - theo_wind: Theoretical wind turbine potential --> Possible to investigate with: https://github.com/FZJ-IEK3-VSA/glaes. Note that this value should be in MW.
  
Once you have created a hexagon file with these features, save it in the `Data` folder as `hex_final_[COUNTRY ISO CODE].geojson`. 

**Note:** `COUNTRY ISO CODE` is the country's ISO standard 2-letter abbreviation.
  
## Input parameter Excel files

Required input parameters include the spatial area of interest, total annual demand for hydrogen, and prices and cost of capital for infrastructure investments. These values can be either current values or projected values for a single snapshot in time. The parameter values for running the model can be specified in a set of Excel files in the Parameters folder.

- **Basic H2 plant:** In this folder, there are several csv files containing the global parameters for optimizing the plant design. All power units are MW and all energy units are MWh. For more information on these parameters, refer to the [PyPSA documentation](https://pypsa.readthedocs.io/en/latest/components.html).

**Note:** The excel files must be kept in a folder with the title matching the Country ISO Code. From the use case of Namibia, we have them in a folder titled "NA"

- **Conversion parameters:** `conversion_parameters.xlsx` includes parameters related to converting between states of hydrogen.

- **Country parameters:** `country_parameters.xlsx` includes country- and technology-specific interest rates, heat and electricity costs, and asset lifetimes.
    - Interest rates should be expressed as a decimal, e.g. 5% as 0.05.
    - Asset lifetimes should be in years.

- **Demand parameters:** `demand_parameters.xlsx` includes a list of demand centers. For each demand center, its lat-lon location, annual demand, and hydrogen state for that demand must be specified. If multiple forms of hydrogen are demanded in one location, differentiate the demand center name (e.g. Nairobi LH2 and Nairobi NH3) to avoid problems from duplicate demand center names.

- **Pipeline parameters:** `pipeline_parameters.xlsx` includes the price, capacity, and lifetime data for different sizes of hydrogen pipeline.

- **Technology parameters:** ` technology_parameters.xlsx` includes water parameters, road infrastructure parameters, and whether road and hydrogen pipeline construction is allowed.

- **Transport parameters:** `transport_parameters.xlsx` includes the parameters related to road transport of hydrogen, including truck speed, cost, lifetime, and capacity.
___

# Snakemake

This repository uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to automate its workflow (for a gentle introduction to Snakemake, see [Getting Started with Snakemake](https://carpentries-incubator.github.io/workflows-snakemake/) on The Carpentries Incubator).

## Wildcards
Wildcards specify the data used in the workflow. This workflow uses two wildcards: `country` (an ISO standard 2-letter abbreviation) and `weather_year` (a 4-digit year between 1940 and 2023 included in the ERA5 dataset).

## Config file

High-level workflow settings are controlled in the config file: `config.yaml`. 

Multiple wildcard values are specified in the `scenario` section. These can be changed to match the `country` and `weather_year` you are analysing.

Renewable generators considered for hydrogen plant construction are included in the `generators` section.

In the `transport` section, `pipeline_construction` and `road_construction` can be switched from `True` to `False`, as needed.

 **Note:** `country` and `weather_year` can be a list of more than one, depending on how many countries and years you are analysing. You must ensure all other files that need for each country run are where they should be.

## Rules

Rules can be run multiple ways using Snakemake. Below, you will be able to run rules by entering the rule name or their output in the terminal. Snakemake will run all necessary rules and their corresponding scripts to create an output. While all rules are discussed here for completeness, **you do not need to enter each rule one-by-one and can simply enter the output you're interested in or one of the run all rules.** Rules are defined in the `Snakefile`.

Snakemake requires a specification of the `number of cores to be used`; this can be up to 4.

### Run time

The `get_weather_data` rule, depending on country size and your internet connection, could take from a few minutes to several hours to run. Ensure that you have space on your computer to store the data, which can be several GB.

The `optimize_hydrogen_plant` rule, depending on country size and the number of demand centers, could take from several minutes to several hours to run.

The `optimize_transport_and_conversion` rule, depending on country size, should take a few minutes to run.

All other rules take a few seconds to run.

### Rule to remove all files

**Note:** This rule does not work on Windows, as of yet. Please manually remove the files you need to.

This rule is important to know first, as it will remove all the files that the below rules will create as well as the file you initially saved into the `Data` folder as `hex_final_[COUNTRY ISO CODE].geojson`.

This is to allow for a quicker transition to analyse more data and to clear up space. Make sure you save the created files that you need elsewhere before running the following rule into the terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] clean
```

### Run optimisation and mapping rules

This section can be used to run most rules, without having to run exact output files. If any files are changed after a completed run, the same command can be used again and Snakemake will only run the necessary scripts to ensure the results are up to date.

The total hydrogen cost for all scenarios can be run by entering the following rule into the terminal (make sure you have the necessary weather file in the Cutouts folder before running, you might have to run the get_weather_data rule):
```
snakemake -j [NUMBER OF CORES TO BE USED] optimise_all
```
Similarly, you can map hydrogen costs for all scenarios with the following rule:
```
snakemake -j [NUMBER OF CORES TO BE USED] map_all
```

### `assign_country` rule

Assign country-specific interest rates, technology lifetimes, and heat and electricity prices from `country_parameters.xlsx` to different hexagons based on their country.

You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Data/hexagons_with_country_[COUNTRY ISO CODE].geojson
```

### `get_weather_data` rule

**Note:** This rule will also create the assign_country rule's output, as it uses that file.
You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Cutouts/[COUNTRY ISO CODE]_[WEATHER YEAR].nc
```

### `optimize_transport_and_conversion` rule

Calculate the cost of the optimal hydrogen transportation and conversion strategy from each hexagon to each demand center, using both pipelines and road transport, using parameters from `technology_parameters.xlsx`, `demand_parameters.xlsx`, and `country_parameters.xlsx`.

You can run this rule by entering the following command in your terminal: 
```
snakemake -j [NUMBER OF CORES TO BE USED] Resources/hex_transport_[COUNTRY ISO CODE].geojson
```

### `calculate_water_costs` rule

Calculate water costs from the ocean and freshwater bodies for hydrogen production in each hexagon using `Parameters/technology_parameters.xlsx` and `Parameters/country_parameters.xlsx`.

You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Resources/hex_water_[COUNTRY ISO CODE].geojson
```

### `optimize_hydrogen_plant` rule

Design green hydrogen plant to meet the hydrogen demand profile for each demand center for each transportation method to each demand center using the `optimize_hydrogen_plant.py` script. Ensure that you have specified your hydrogen plant parameters in the CSV files in the `Parameters/Basic_H2_plant` folder, your investment parameters in `Parameters/investment_parameters.xlsx`, and your demand centers in `Parameters/demand_parameters.xlsx`.

You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Resources/hex_lcoh_[COUNTRY ISO CODE]_[WEATHER YEAR].geojson
```

### `calculate_total_hydrogen_cost` rule
Combine results to find the lowest-cost method of producing, transporting, and converting hydrogen for each demand center.

You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Results/hex_total_cost_[COUNTRY ISO CODE]_[WEATHER YEAR].geojson
```

### `calculate_cost_components` rule

Calculate the cost for each type of equipment in each polygon. 

You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Results/hex_cost_components_[COUNTRY ISO CODE]_[WEATHER YEAR].geojson
```

### `map_costs` rule

Visualize the spatial variation in different costs per kilogram of hydrogen.

You can run this rule by entering the following command in your terminal:
```
snakemake -j [NUMBER OF CORES TO BE USED] Plots/[COUNTRY ISO CODE]_[WEATHER YEAR]
```
___

# Limitations

This model considers only greenfield wind and solar plants for hydrogen production. Therefore it does not consider using grid electricity or existing generation for hydrogen production. The model further assumes that all excess electricity is curtailed.

While the design of the green hydrogen plant is convex and therefore guarenteed to find the global optimum solution if it exists, the selection of the trucking strategy is greedy to avoid the long computation times and potential computational intractability associated with a mixed-integer optimization problem.

Currently, only land transport is considered in the model. To calculate the cost of hydrogen production for export, any additional costs for conversion and transport via ship or undersea pipeline must be added in post-processing.

Transport costs are calculated from the center of the hexagon to the demand center. When using large hexagon sizes, this assumption may over- or underestimate transportation costs significantly. Additionally, only path length is considered when calculating the cost of road and pipeline construction. Additional costs due to terrain are not considered.

The availability of water for electrolysis is not limited in regions that could potentially face drought, and a single prices for freshwater and ocean water are used throughout the modeled area.
___

# Citation

If you decide to use GeoH2, please kindly cite us using the following: 

*Halloran, C., Leonard, A., Salmon, N., Müller, L., & Hirmer, S. (2024). 
GeoH2 model: Geospatial cost optimization of green hydrogen production including storage and transportation. 
MethodsX, 12, 102660. https://doi.org/10.1016/j.mex.2024.102660.*

```commandline
@article{Halloran_GeoH2_model_Geospatial_2024,
author = {Halloran, Claire and Leonard, Alycia and Salmon, Nicholas and Müller, Leander and Hirmer, Stephanie},
doi = {10.1016/j.mex.2024.102660},
journal = {MethodsX},
month = jun,
pages = {102660},
title = {{GeoH2 model: Geospatial cost optimization of green hydrogen production including storage and transportation}},
volume = {12},
year = {2024}
}
```
___

# Case study parameters

This repository includes sample parameters for a hydrogen production case in Namibia.
References for these parameters are included in the tables below for reference.
For the results of this case, please refer to the model MethodsX article: https://doi.org/10.1016/j.mex.2024.102660. 

**Green hydrogen plant parameters:**

| Hardware                   | Parameter             | Value     | Units                | Ref.                                                                                                                            |
|----------------------------|-----------------------|-----------|----------------------|---------------------------------------------------------------------------------------------------------------------------------|
| Solar photovoltaic         | Capex                 | 1,470,000 | €/MW                 | [Allington et al., 2021](https://doi.org/10.1016/j.dib.2022.108021)                                                             |
| Wind turbines              | Capex                 | 1,580,000 | €/MW                 | [Allington et al., 2021](https://doi.org/10.1016/j.dib.2022.108021)                                                             |
| Hydrogen electrolysis      | Capex                 | 1,250,000 | €/MW                 | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                           |
| Hydrogen electrolysis      | Efficiency            | 0.59      | MWh H2/MWh el        | [Taibi et al., 2020](https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2020/Dec/IRENA_Green_hydrogen_cost_2020.pdf)  |
| Hydrogen compression       | Isentropic efficiency | 0.051     | MWh el/MWh H2        | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                           |
| Hydrogen storage unloading | Efficiency            | 1         | MWh H2/MWh H2-stored | Assumption                                                                                                                      |
| Battery                    | Capex                 | 95,000    | €/MW                 | [BloombergNEF, 2022](https://about.bnef.com/blog/lithium-ion-battery-pack-prices-rise-for-first-time-to-an-average-of-151-kwh/) |
| Hydrogen storage           | Capex                 | 21,700    | €/MWh                | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                           |

**Conversion parameters:**

| Process                | Parameter                    | Value       | Units             | Ref.                                                                                                                                             |
|------------------------|------------------------------|-------------|-------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| 500 bar compression    | Heat capacity                | 0.0039444   | kWh/kg/K          | Kurzweil and Dietlmeier, 2016                                                                                                                    |
| 500 bar compression    | Input temperature            | 298.15      | K                 | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| 500 bar compression    | Input pressure               | 25          | bar               | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| 500 bar compression    | Isentropic exponent          | 1.402       |                   | Kurzweil and Dietlmeier, 2016                                                                                                                    |
| 500 bar compression    | Isentropic efficiency        | 0.8         |                   | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| 500 bar compression    | Compressor lifetime          | 15          | years             | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar compression    | Compressor capex coefficient | 40,035      | €/kg H2/day       | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar compression    | Compressor opex              | 4           | % capex/year      | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| Hydrogen liquification | Electricity demand           | 9.93        | kWh/kg H2         | [Ausfelder and Dura](https://dechema.de/dechema_media/Downloads/Positionspapiere/2021_DEC_P2X_III_Technischer_Anhang.pdf)                        |
| Hydrogen liquification | Capex quadratic coefficient  | -0.0002     | €/(kg H2)^2       | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Hydrogen liquification | Capex linear coefficient     | 1,781.9     | €/kg  H2          | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Hydrogen liquification | Capex constant               | 300,000,000 | €                 | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Hydrogen liquification | Opex                         | 8           | % capex/year      | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| Hydrogen liquification | Plant lifetime               | 20          | years             | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| LOHC hydrogenation     | Electricity demand           | 0.35        | kWh/kg H2         | [Andersson and Grönkvist, 2019](https://doi.org/10.1016/j.ijhydene.2019.03.063)                                                                  |
| LOHC hydrogenation     | Heat demand                  | -9          | kWh/kg H2         | [Hydrogenious, 2022](https://hydrogenious.net/how/)                                                                                              |
| LOHC hydrogenation     | Capex coefficient            | 0.84        | kWh/kg H2/year    | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC hydrogenation     | Opex                         | 4           | % capex/year      | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC hydrogenation     | Plant lifetime               | 25          | years             | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC hydrogenation     | Carrier costs                | 2           | €/kg carrier      | [Clark, 2020](https://hydrogenious.net/lohc-global-hydrogen-opportunity/)                                                                        |
| LOHC hydrogenation     | Carrier ratio                | 16.1        | kg carrier/kg  H2 | [Arlt and Obermeier, 2017](https://www.encn.de/fileadmin/user_upload/Studie_Wasserstoff_und_Schwerlastverkehr_WEB.pdf)                           |
| LOHC dehydrogenation   | Electricity demand           | 0.35        | kWh/kg H2         | [Andersson and Grönkvist, 2019](https://doi.org/10.1016/j.ijhydene.2019.03.063)                                                                  |
| LOHC dehydrogenation   | Heat demand                  | 12          | kWh/kg H2         | [Hydrogenious, 2022](https://hydrogenious.net/how/)                                                                                              |
| LOHC dehydrogenation   | Capex coefficient            | 2.46        | kWh/kg H2         | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC dehydrogenation   | Opex                         | 4           | % capex/year      | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC dehydrogenation   | Plant lifetime               | 25          | years             | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia synthesis      | Electricity demand           | 2.809       | kWh/kg H2         | [IEA, 2021](https://www.iea.org/reports/ammonia-technology-roadmap)                                                                              |
| Ammonia synthesis      | Capex coefficient            | 0.75717     | kWh/g H2/year     | [IEA, 2021](https://www.iea.org/reports/ammonia-technology-roadmap)                                                                              |
| Ammonia synthesis      | Opex                         | 1.5         | % capex/year      | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia synthesis      | Plant lifetime               | 25          | years             | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia cracking       | Heat demand                  | 4.2         | kWh/kg H2         | [Andersson and Grönkvist, 2019](https://doi.org/10.1016/j.ijhydene.2019.03.063)                                                                  |
| Ammonia cracking       | Capex coefficient            | 17,262,450  | kWh/g H2/hour     | [Cesaro et al., 2021](https://doi.org/10.1016/j.apenergy.2020.116009)                                                                            |
| Ammonia cracking       | Opex                         | 2           | % capex/year      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Ammonia cracking       | Plant lifetime               | 25          | years             | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |

**Trucking parameters:**

| Hardware                 | Parameter                  | Value   | Units        | Ref.                                                                                                                                             |
|--------------------------|----------------------------|---------|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| All trucks               | Average truck speed        | 70      | km/h         | Assumption                                                                                                                                       |
| All trucks               | Working hours              | 24      | h/day        | Assumption                                                                                                                                       |
| All trucks               | Diesel price               | 1.5     | €/L          | Assumption                                                                                                                                       |
| All trucks               | Driver wage                | 2.85    | €/h          | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| All trucks               | Working days               | 365     | days/year    | Assumption                                                                                                                                       |
| All trucks               | Max driving distance       | 160,000 | km/year      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| All trucks               | Truck capex                | 160,000 | €            | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Truck Opex                 | 12      | % capex/year | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Diesel consumption         | 35      | L/100 km     | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Truck lifetime             | 8       | years        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Trailer lifetime           | 12      | years        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| 500 bar hydrogen trailer | Trailer capex              | 660,000 | €            | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar hydrogen trailer | Trailer opex               | 2       | % capex/year | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar hydrogen trailer | Trailer capacity           | 1,100   | kg	H2        | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar hydrogen trailer | Loading and unloading time | 1.5     | hours        | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| Liquid hydrogen trailer  | Trailer capex              | 860,000 | €            | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Liquid hydrogen trailer  | Trailer opex               | 2       | % capex/year | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Liquid hydrogen trailer  | Trailer capacity           | 4,300   | kg H2        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Liquid hydrogen trailer  | Loading and unloading time | 3       | hours        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| LOHC trailer             | Trailer capex              | 660,000 | €            | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC trailer             | Trailer opex               | 2       | % capex/year | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| LOHC trailer             | Trailer capacity           | 1,800   | kg H2        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| LOHC trailer             | Loading and unloading time | 1.5     | hours        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Ammonia trailer          | Trailer capex              | 210,000 | €            | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia trailer          | Trailer opex               | 2       | % capex/year | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia trailer          | Trailer capacity           | 2,600   | kg H2        | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia trailer          | Loading and unloading time | 1.5     | hours        | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |

**Road parameters:**

| Road length         | Parameter | Value      | Units     | Ref.                                                                  |
|---------------------|-----------|------------|-----------|-----------------------------------------------------------------------|
| Short road (<10 km) | Capex     | 626,478.45 | €/km      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219) |
| Long road (>10 km)  | Capex     | 481,866.6  | €/km      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219) |
| All roads           | Opex      | 7,149.7    | €/km/year | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219) |

**Pipeline parameters:**

| Pipeline size   | Parameter           | Value     | Units        | Ref.                                                                                             |
|-----------------|---------------------|-----------|--------------|--------------------------------------------------------------------------------------------------|
| All pipelines   | Opex                | 1.25      | % capex/year | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| All pipelines   | Availability        | 95        | %            | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                            |
| All pipelines   | Pipeline lifetime   | 42.5      | years        | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| All pipelines   | Compressor lifetime | 24        | years        | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| All pipelines   | Electricity demand  | 0.000614  | kWh/kg H2/km | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Large pipeline  | Maximum capacity    | 13        | GW           | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Large pipeline  | Pipeline capex      | 2,800,000 | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Large pipeline  | Compressor capex    | 620,000   | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Medium pipeline | Maximum capacity    | 4.7       | GW           | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Medium pipeline | Pipeline capex      | 2,200,000 | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Medium pipeline | Compressor capex    | 310,000   | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Small pipeline  | Maximum capacity    | 1.2       | GW           | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Small pipeline  | Pipeline capex      | 90,000    | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Small pipeline  | Compressor capex    | 90,000    | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |

**Water parameters:**

| Type        | Parameter                    | Value | Units              | Ref.                                                                                                                                                                              |
|-------------|------------------------------|-------|--------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Freshwater  | Treatment electricity demand | 0.4   | kWh/m^3 water      | [US Dept. of Energy, 2016](https://betterbuildingssolutioncenter.energy.gov/sites/default/files/Primer%20on%20energy%20efficiency%20in%20water%20and%20wastewater%20plants_0.pdf) |
| Ocean water | Treatment electricity demand | 3.7   | kWh/m^3 water      | [Patterson et al., 2019](https://doi.org/10.1073/pnas.1902335116)                                                                                                                 |
| All water   | Transport cost               | 0.1   | €/100 km/m^3 water | [Zhou and Tol, 2005](https://doi.org/10.1029/2004WR003749)                                                                                                                        |
| All water   | Water specific cost          | 1.25  | €/m^3 water        | [Wasreb, 2019](https://wasreb.go.ke/wasrebsystems/tariffs/about-us.html)                                                                                                          |
| All water   | Water demand                 | 21    | L water/kg H2      | [Taibi et al., 2020](https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2020/Dec/IRENA_Green_hydrogen_cost_2020.pdf)                                                    |

**Country-specific parameters:**

| Country | Parameter                    | Value   | Units | Ref.                                                                               |
|---------|------------------------------|---------|-------|------------------------------------------------------------------------------------|
| Namibia | Electricity price            | 0.10465 | €/kWh | [GlobalPetrolPrices.com]({https://www.globalpetrolprices.com/electricity_prices/}) |
| Namibia | Heat price                   | 0.02    | €/kWh | Assumption                                                                         |
| Namibia | Solar interest rate          | 6       | %     | Assumption                                                                         |
| Namibia | Solar lifetime               | 20      | years | Assumption                                                                         |
| Namibia | Wind interest rate           | 6       | %     | Assumption                                                                         |
| Namibia | Wind lifetime                | 20      | years | Assumption                                                                         |
| Namibia | Plant interest rate          | 6       | %     | Assumption                                                                         |
| Namibia | Plant lifetime               | 20      | years | Assumption                                                                         |
| Namibia | Infrastructure interest rate | 6       | %     | Assumption                                                                         |
| Namibia | Infrastructure lifetime      | 50      | years | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)              |
