# GeoX
**Geospatial analysis of commodity production costs**

GeoX calculates the locational cost of production, storage, transport, and conversion of different commodities to meet demand in a specified location.
These costs can be compared to current or projected prices in the region under study to assess competitiveness. 
Currently, green hydrogen and green ammonia are implemented as commodities in GeoX. 
However, the code is written in a modular way to allow further commodities to be implemented.

The model outputs the levelised cost (LC) of the specified commodity at the specified demand location including production, storage, transport, and conversion costs from each production location. 

In the code provided, the use case of hydrogen and ammonia production in Namibia is provided as an example.
Parameter references for the hydrogen case are attached.
However, as the code is written in a generalized way, it is possible to analyse all sorts of regions.

GeoX builds upon a preliminary code iteration produced by Leander Müller, available under a CC-BY-4.0 licence: [https://github.com/leandermue/GEOH2](https://github.com/leandermue/GEOH2).
It also integrates code produced by Nick Salmon under an MIT licence: 
[https://github.com/nsalmon11/LCOH_Optimisation](https://github.com/nsalmon11/LCOH_Optimisation)

___

# Setup instructions

Follow the steps below before running GeoX. 
Note that GeoX is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/).
This lets you automate code execution using different "rules".
When we refer to rules below, we're just talking about running different parts of the codebase that execute different functionalities.

## 1) Clone the repository
First, clone the GeoX repository using `git`. 

`... % git clone https://github.com/ClimateCompatibleGrowth/GeoH2.git`

## 2) Set up the environment
The python package requirements are in the `environment.yaml` file. 
You can install these requirements in a new environment using the `mamba` package and environment manager (installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)): 

` .../GEOX % mamba env create -f environment.yaml`

Then activate this new environment using:

`.../GEOX % mamba activate geox`

## 3) Get a Climate Data Store API key
To simulate commodities which use renewable energy for production, historical climate data is needed as an input.
The `get_weather_data` rule can be used to download the relevant data from the ERA-5 reanalysis dataset using [Atlite](https://atlite.readthedocs.io/en/latest/) to create a cutout. 
For this process to work, you need to register and set up your CDS API key as described on the [Climate Data Store website](https://cds.climate.copernicus.eu/api-how-to).

## 4) Install a solver
For the `plant_optimization` rule to work, you will need a solver installed on your computer. 
You can use any solver that works with [PyPSA](https://pypsa.readthedocs.io/en/latest/installation.html), such as [Cbc](https://github.com/coin-or/Cbc), a free, open-source solver, or [Gurobi](https://www.gurobi.com/), a commerical solver with free academic licenses available. 
Install your solver of choice following the instructions for use with Python and your operating system in the solver's documentation. 

> [!NOTE]
> Snakemake uses Cbc as a solver by default. 
> This will be installed upon environment setup. 
> To check that this is installed, activate your environment and enter `mamba list` in your terminal for the environment's list of packages.
> If you choose to use Cbc, no additional solver needs to be installed.
___

# Preparing input data

Next, you'll need to prepare input data. 
GeoX has two types of input data: (1) spatial input data as a hexagon shapefile; and (2) techno-economic parameters in spreadsheets.

## 1) Prepare input hexagons 
First, prepare spatial input data as a set of hexagons using the [H3 standard](https://h3geo.org/).
These hexagons are provided in the repository for the illustrative case study of Namibia.
To analyse a different area of interest, the hexagon file needs to be changed, but needs to follow the logic of the one provided. 

A full walkthrough on making these hexagons, including tools to create them, are available in the [GeoH2-data-prep](https://github.com/ClimateCompatibleGrowth/GeoH2-data-prep) repository.
The hexagon file needs to have the following attributes:

  - waterbody_dist: Distance to selected waterbodies in area of interest
  - waterway_dist: Distance to selected waterways in area of interest
  - ocean_dist: Distance to ocean coastline
  - grid_dist: Distance to transmission network
  - road_dist: Distance to road network
  - theo_pv: Theoretical potential of standardized PV plants (can be determined with [GLAES](https://github.com/FZJ-IEK3-VSA/glaes))
  - theo_wind: Theoretical potential of standardized wind turbines (can be determined with [GLAES](https://github.com/FZJ-IEK3-VSA/glaes))

The [ccg-spider](https://github.com/carderne/ccg-spider) repository can be used to combine different spatial data into standard H3 hexagons.
Once you have created a hexagon file with these features, save it in the `data` folder as `hex_final_[COUNTRY ISO CODE].geojson`. 

> [!IMPORTANT]
> `COUNTRY ISO CODE` is the country's ISO standard 2-letter abbreviation.
  
## 2) Prepare input parameter Excel files
Next, prepare the techno-economic input parameter spreadsheets.
Required input parameters include the spatial area of interest, total annual demand for the commodity, and prices and cost of capital for infrastructure investments. 
These values can be either current values or projected values for a single snapshot in time. 
The parameter values for running the model can be specified in a set of Excel files in the `parameters` folder.

- **Basic plant:** The `basic_[COMMODITY ABBREVIATION]_plant` folder contains several csv files containing the global parameters for optimizing the plant design. 
All power units are MW and all energy units are MWh. 
For more information on these parameters, refer to the [PyPSA documentation](https://pypsa.readthedocs.io/en/latest/components.html).

> [!IMPORTANT]
> `COMMODITY ABBREVIATION` can be 'h2' or 'nh3' for current commodities. 
> The excel files must be kept in a folder with the commodity name within another folder with the title matching the Country ISO Code. 
> As currently implemented, the commodity must be either "hydrogen" or "ammonia".
> For the illustrative case of Namibia, we have them in a folder titled "NA" with two sub-folders "hydrogen" and "ammonia".

- **Conversion parameters:** `conversion_parameters.xlsx` includes parameters related to converting between states of the commodity.

- **Country parameters:** `country_parameters.xlsx` includes country- and technology-specific interest rates, heat and electricity costs, and asset lifetimes.
    - Interest rates should be expressed as a decimal, e.g. 5% as 0.05.
    - Asset lifetimes should be in years.

- **Demand parameters:** `demand_parameters.xlsx` includes a list of demand centers. 
For each demand center, its lat-lon location, annual demand, and commodity state for that demand must be specified. 
If multiple forms of the commodity are demanded in one location, differentiate the demand center name (e.g. Nairobi LH2 and Nairobi NH3) to avoid problems from duplicate demand center names.

- **Pipeline parameters:** `pipeline_parameters.xlsx` includes the price, capacity, and lifetime data for different sizes of pipeline.

- **Technology parameters:** ` technology_parameters.xlsx` includes water parameters, road infrastructure parameters, and whether road and pipeline construction is allowed.

- **Transport parameters:** `transport_parameters.xlsx` includes the parameters related to road transport of the commodity, including truck speed, cost, lifetime, and capacity.
___

# Running GeoX with Snakemake
Once you've done the repository setup and prepared your input data, you're ready to run GeoX!

This repository uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to automate its workflow (for a gentle introduction to Snakemake, see [Getting Started with Snakemake](https://carpentries-incubator.github.io/workflows-snakemake/) on The Carpentries Incubator).

## 1) Modify the config file

High-level workflow settings are controlled in the config file: `config.yaml`.

**Wildcards**: These are used in the `scenario` section of the config file to specify the data used in the workflow. 
They can be changed to match case you're analysing. They are: 
- `country`: an ISO standard 2-letter abbreviation
- `weather_year`: a 4-digit year between 1940 and 2023 included in the ERA5 dataset
- `plant_type`: the commodity to be produced (currently either "hydrogen" or "ammonia")

**Weather data**:
The amount of years you want to download weather data for should be added into `years_to_check`.
For instance, if you want your weather data to start in 2015 and span five years, `weather_year` should be 2015 and `years_to_check` should be 5.
You can set the frequency of data to be used in optimisation using `freq` (i.e., "1H" for hourly, "3H" for three-hourly, etc.)

**Generators**:
The types of renewable generators considered for plant construction are included in the `generators_dict` section. Currently, only Solar and Wind can be considered.
Dependent on which generators you are using, you can change the `panel` value for Solar and the `turbine` value for Wind.
In the `gen_capacity` section, you will find both `solar` and `wind`, which can be changed to match values that you are analysing.

**Other**:
You will have to set the `solver` to the solver name that you are going to be using. 
You will also have to set whether a `water_limit` is `True` or `False` (i.e., whether you want to consider water scarcity in your process).
In the `transport` section, `pipeline_construction` and `road_construction` can be switched from `True` to `False`, as needed.

> [!NOTE]
> `country` and `weather_year` can be a list of more than one, depending on how many countries and years you are analysing. 
> You must ensure all other input files that are needed for each country are where they should be as previously described.

## 2) Run GeoX rules

GeoX assumes that you will run it using Snakemake by entering rule names in the terminal. 
When you do this, Snakemake will run all the necessary rules to achieve the rule you entered and create the desired output.
Rules are defined in the `Snakefile`.

Snakemake requires a specification of the `NUMBER OF CORES TO BE USED` when you run a rule. This can be up to 4.

Each rule has a different expected run time:
- The `get_weather_data` rule, depending on country size and your internet connection, could take from a few minutes to several hours to run. 
Ensure that you have space on your computer to store the data, which can be several GB.
- The `optimize_plant` rule, depending on country size and the number of demand centers, could take from several minutes to several hours to run.
- The `optimize_transport` rule, depending on country size, should take a few minutes to run. 
- All other rules take a few seconds to run.

### 1) Get weather data (if needed)

If you already have weather data for your area of interest, just use step 2!
If not, you need to run two rules.
First, run the prep rule to: (1) assign country-specific interest rates, technology lifetimes, and heat and electricity prices from `country_parameters.xlsx` to different hexagons based on their country; and (2) drop any duplicated hexagons that do not belong to the country that is be run for.
```
snakemake -j [NUMBER OF CORES TO BE USED] run_prep
```
Then, run the weather data which will download the data you need from the CDS API: 

```
snakemake -j [NUMBER OF CORES TO BE USED] run_weather
```

### 2) Run optimisation or mapping rules

Finally, run the GeoX optimisation process. 
These rules are used to run the whole process without having go through each rule step-by-step.
If any input files are changed after a completed run, the same command can be used again and Snakemake will only run the necessary scripts to ensure the results are up to date.

The total commodity cost for all scenarios can be run by entering the following rule into the terminal.
Make sure you have the necessary weather file(s) in the cutouts folder before running. 
```
snakemake -j [NUMBER OF CORES TO BE USED] optimise_all
```
Similarly, you can map commodity costs for all scenarios with the following rule:
```
snakemake -j [NUMBER OF CORES TO BE USED] map_all
```
___
# Definitions of all rules

While all rules are discussed here for completeness, **you do not need to enter each rule one-by-one**. 
You can simply run one of the optimization or mapping rules, and Snakemake will ensure that all required rules are completed to achieve this. 
Please refer to the steps in the previous section for most usages.

### `prep_main` rule
This will assign country-specific interest rates, technology lifetimes, and heat and electricity prices from `country_parameters.xlsx` to different hexagons based on their country. As well as, drop any duplicated hexagons that do not belong to the country that is to be run for.
```
snakemake -j [NUMBER OF CORES TO BE USED] run_prep
```

### `get_weather_data` rule
This will download weather data from your area of interest for generation optimisation.
Make sure you have run the `prep_main` rule before running this rule as it needs the outputted file as an input.
```
snakemake -j [NUMBER OF CORES TO BE USED] run_weather
```

### `optimize_transport` rule
This will calculate the cost of the optimal commodity transportation and conversion strategy from each hexagon to each demand center, using both pipelines and road transport, using parameters from `technology_parameters.xlsx`, `demand_parameters.xlsx`, and `country_parameters.xlsx`.
```
snakemake -j [NUMBER OF CORES TO BE USED] resources/hex_transport_[COUNTRY ISO CODE].geojson
```

### `calculate_water_costs` rule
This will calculate water costs from the ocean and freshwater bodies for commodity production in each hexagon using `parameters/technology_parameters.xlsx` and `parameters/country_parameters.xlsx`.
```
snakemake -j [NUMBER OF CORES TO BE USED] resources/hex_water_[COUNTRY ISO CODE].geojson
```

### `optimize_plant` rule
This will design a plant to meet the commodity demand profile for each demand center for each transportation method to each demand center. 
Ensure that you have specified your plant parameters in the `parameters/basic_[COMMODITY ABBREVIATION]_plant` folder, your investment parameters in `parameters/investment_parameters.xlsx`, and your demand centers in `parameters/demand_parameters.xlsx`.
```
snakemake -j [NUMBER OF CORES TO BE USED] resources/hex_lcoh_[COUNTRY ISO CODE]_[WEATHER YEAR].geojson
```

### `calculate_total_cost` rule
This will combine results to find the lowest-cost method of producing, transporting, and converting the commodity for each demand center.
```
snakemake -j [NUMBER OF CORES TO BE USED] results/hex_total_cost_[COUNTRY ISO CODE]_[WEATHER YEAR].geojson
```

### `calculate_cost_components` rule
This will calculate the cost for each type of equipment in each polygon.
```
snakemake -j [NUMBER OF CORES TO BE USED] results/hex_cost_components_[COUNTRY ISO CODE]_[WEATHER YEAR].geojson
```

### `calculate_map_costs` rule
This will visualize the spatial variation in different costs per kilogram of hydrogen.
```
snakemake -j [NUMBER OF CORES TO BE USED] plots/[COUNTRY ISO CODE]_[WEATHER YEAR]
```
___

# Limitations

This model considers only greenfield generation. 
Therefore it does not consider using grid electricity or existing generation. 
The model further assumes that all excess electricity is curtailed.

While the design of the commodity plant is convex and therefore guarenteed to find the global optimum solution if it exists, the selection of the trucking strategy is greedy to avoid the long computation times and potential computational intractability associated with a mixed-integer optimization problem.

Currently, only land transport is considered in the model. 
To calculate the cost of commodity production for export, any additional costs for conversion and transport via ship or undersea pipeline must be added in post-processing.

Transport costs are calculated from the center of the hexagon to the demand center. 
When using large hexagon sizes, this assumption may over- or underestimate transportation costs significantly. 
Additionally, only path length is considered when calculating the cost of road and pipeline construction. 
Additional costs due to terrain are not considered.
___

# Citation

If you decide to use GeoX, please kindly cite us using the following: 

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
For the results of this case, please refer to the GeoH2 model MethodsX article: https://doi.org/10.1016/j.mex.2024.102660. 

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
