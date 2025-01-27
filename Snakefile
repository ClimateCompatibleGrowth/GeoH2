configfile: 'config.yaml'

wildcard_constraints:
    # ISO alpha-2 country code
    country = '[A-Z]{2}',
    # ERA5 weather year (1940-2023)
    weather_year = '(19[4-9]\d|20[0-1]\d|202[0-3])',

# rule to delete all necessary files to allow reruns
rule clean:
    shell: 'rm -r cutouts/*.nc data/*.geojson resources/*.geojson results/*.geojson temp/*.nc results/*.csv plots/'
    
# bulk run rule to run all countries and years listed in config file
rule optimise_all:
   input:
        expand('results/hex_total_cost_{country}_{weather_year}.geojson',
        **config['scenario']
        ),
        
# bulk run rule to map all countries and years listed in config file
rule map_all:
   input:
       expand('plots/{country}_{weather_year}',
        **config['scenario']
        ),

rule prep_main:
    input:
        hexagons = 'data/hex_final_{country}.geojson',
        country_parameters = expand('parameters/{country}/{plant_type}/country_parameters.xlsx', 
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"])
    output:
        'data/hexagons_with_country_{country}.geojson',
    script:
        'src/prep/main.py'
        
rule get_weather_data:
    input:
        hexagons = 'data/hexagons_with_country_{country}.geojson',
    output:
        'cutouts/{country}_{weather_year}.nc',
    script:
        'src/prep/get_weather_data.py'

rule transport_optimization:
    input:
        hexagons = 'data/hexagons_with_country_{country}.geojson',
        technology_parameters = expand('parameters/{country}/{plant_type}/technology_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        demand_parameters = expand('parameters/{country}/{plant_type}/demand_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        country_parameters = expand('parameters/{country}/{plant_type}/country_parameters.xlsx', 
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        transport_parameters = expand('parameters/{country}/{plant_type}/transport_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        pipeline_parameters = expand('parameters/{country}/{plant_type}/pipeline_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"])
    output:
        'resources/hex_transport_{country}.geojson'
    script:
        'src/main/transport_optimization.py'

rule calculate_water_costs:
    input:
        technology_parameters = expand('parameters/{country}/{plant_type}/technology_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        country_parameters = expand('parameters/{country}/{plant_type}/country_parameters.xlsx', 
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        hexagons = 'resources/hex_transport_{country}.geojson'
    output:
        'resources/hex_water_{country}.geojson'
    script:
        'src/main/water_cost.py'

rule plant_optimization:
    input:
        transport_parameters = expand('parameters/{country}/{plant_type}/transport_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        country_parameters = expand('parameters/{country}/{plant_type}/country_parameters.xlsx', 
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        demand_parameters = expand('parameters/{country}/{plant_type}/demand_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        hexagons = 'resources/hex_water_{country}.geojson'
    output:
        'resources/hex_lc_{country}_{weather_year}.geojson'
    script:
        'src/main/plant_optimization.py'

rule calculate_total_costs:
    input:
        hexagons = 'resources/hex_lc_{country}_{weather_year}.geojson',
        demand_parameters = expand('parameters/{country}/{plant_type}/demand_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"])
    output:
        'results/hex_total_cost_{country}_{weather_year}.geojson'
    script:
        'src/main/total_costs.py'

rule calculate_cost_components:
    input:
        hexagons = 'results/hex_total_cost_{country}_{weather_year}.geojson',
        demand_parameters = expand('parameters/{country}/{plant_type}/demand_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"]),
        country_parameters = expand('parameters/{country}/{plant_type}/country_parameters.xlsx', 
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"])
    output:
        'results/hex_cost_components_{country}_{weather_year}.geojson',
        'results/hex_cost_components_{country}_{weather_year}.csv'
    script:
        'src/main/costs_by_component.py'

rule map_costs:
    input:
        hexagons = 'results/hex_cost_components_{country}_{weather_year}.geojson',
        demand_parameters = expand('parameters/{country}/{plant_type}/demand_parameters.xlsx',
                                        plant_type=config["plant_type"].lower(), 
                                        country=config["scenario"]["country"])
    output:
        directory('plots/{country}_{weather_year}')
    script:
        'src/main/map_costs.py'

