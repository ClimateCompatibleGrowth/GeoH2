configfile: "config.yaml"

wildcard_constraints:
    # ISO alpha-2 country code
    country="[A-Z]{2}",
    # ERA5 weather year (1940-2023)
    weather_year="(19[4-9]\d|20[0-1]\d|202[0-3])",

# rule to delete all necessary files to allow reruns
rule clean:
    shell: 'rm -r Cutouts/*.nc Data/*.geojson Resources/*.geojson Results/*.geojson temp/*.nc Results/*.csv Plots/'
    
# bulk run rule to run all countries and years listed in config file
rule optimise_all:
   input:
        expand('Results/hex_total_cost_{country}_{weather_year}.geojson',
        **config["scenario"]
        ),
        
# bulk run rule to map all countries and years listed in config file
rule map_all:
   input:
       expand('Plots/{country}_{weather_year}',
        **config["scenario"]
        ),

rule assign_country:
    input:
        "Data/hex_final_{country}.geojson",
    output:
        "Data/hexagons_with_country_{country}.geojson",
    script:
        "Scripts/assign_country.py"
        
rule get_weather_data:
    input:
        hexagons = "Data/hexagons_with_country_{country}.geojson",
    output:
        "Cutouts/{country}_{weather_year}.nc",
    script:
        'Scripts/get_weather_data.py'

rule optimize_transport_and_conversion:
    input:
        hexagons = 'Data/hexagons_with_country_{country}.geojson',
        technology_parameters = "Parameters/{country}/technology_parameters.xlsx",
        demand_parameters = 'Parameters/{country}/demand_parameters.xlsx',
        country_parameters = 'Parameters/{country}/country_parameters.xlsx',
        conversion_parameters = "Parameters/{country}/conversion_parameters.xlsx",
        transport_parameters = "Parameters/{country}/transport_parameters.xlsx",
        pipeline_parameters = "Parameters/{country}/pipeline_parameters.xlsx"
    output:
        'Resources/hex_transport_{country}.geojson'
    script:
        'Scripts/optimize_transport_and_conversion.py'

rule calculate_water_costs:
    input:
        technology_parameters = "Parameters/{country}/technology_parameters.xlsx",
        country_parameters = 'Parameters/{country}/country_parameters.xlsx',
        hexagons = 'Resources/hex_transport_{country}.geojson'
    output:
        'Resources/hex_water_{country}.geojson'
    script:
        'Scripts/water_cost.py'
        

rule optimize_hydrogen_plant:
    input:
        transport_parameters = "Parameters/{country}/transport_parameters.xlsx",
        country_parameters = 'Parameters/{country}/country_parameters.xlsx',
        demand_parameters = 'Parameters/{country}/demand_parameters.xlsx',
        # cutout = "Cutouts/{country}_{weather_year}.nc",
        hexagons = 'Resources/hex_water_{country}.geojson'
    output:
        'Resources/hex_lcoh_{country}_{weather_year}.geojson'
    script:
        'Scripts/optimize_hydrogen_plant.py'

rule calculate_total_hydrogen_cost:
    input:
        hexagons = 'Resources/hex_lcoh_{country}_{weather_year}.geojson',
        demand_parameters = 'Parameters/{country}/demand_parameters.xlsx'
    output:
        'Results/hex_total_cost_{country}_{weather_year}.geojson'
    script:
        'Scripts/total_hydrogen_cost.py'

rule calculate_cost_components:
    input:
        hexagons = 'Results/hex_total_cost_{country}_{weather_year}.geojson',
        demand_parameters = 'Parameters/{country}/demand_parameters.xlsx',
        country_parameters = 'Parameters/{country}/country_parameters.xlsx',
        stores_parameters = 'Parameters/Basic_H2_plant/stores.csv',
        storage_parameters = 'Parameters/Basic_H2_plant/storage_units.csv',
        links_parameters = 'Parameters/Basic_H2_plant/links.csv',
        generators_parameters = 'Parameters/Basic_H2_plant/generators.csv'
    output:
        'Results/hex_cost_components_{country}_{weather_year}.geojson',
        'Results/hex_cost_components_{country}_{weather_year}.csv'
    script:
        'Scripts/costs_by_component.py'

rule map_costs:
    input:
        hexagons = 'Results/hex_cost_components_{country}_{weather_year}.geojson',
        demand_parameters = 'Parameters/{country}/demand_parameters.xlsx'
    output:
        directory('Plots/{country}_{weather_year}')
    script:
        'Scripts/map_costs.py'

