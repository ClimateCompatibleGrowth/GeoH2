configfile: 'config.yaml'

wildcard_constraints:
    # ISO alpha-2 country code
    country = '[A-Z]{2}',
    # ERA5 weather year (1940-2023)
    weather_year = '(19[4-9]\d|20[0-1]\d|202[0-3])',
    plant_type = str()

# bulk run rule to run all countries and years listed in config file
rule optimise_all:
   input:
        expand('results/hex_total_cost_{country}_{weather_year}_{plant_type}.geojson',
        **config['scenario']
        ),
        
# bulk run rule to map all countries and years listed in config file
rule map_all:
   input:
       expand('plots/{country}_{weather_year}_{plant_type}',
        **config['scenario']
        ),
        
rule run_prep:
    input:
        expand('data/hexagons_with_country_{country}_{plant_type}.geojson',
        **config['scenario']
        ),

rule run_weather:
    input:
        expand('cutouts/{country}_{weather_year}.nc',
        **config['scenario']
        ),

rule prep_main:
    input:
        hexagons = 'data/hex_final_{country}.geojson',
        country_parameters = 'parameters/{country}/{plant_type}/country_parameters.xlsx'
    output:
        'data/hexagons_with_country_{country}_{plant_type}.geojson',
    script:
        'src/prep/main.py'

rule get_weather_data:
    output:
        'cutouts/{country}_{weather_year}.nc',
    script:
        'src/prep/get_weather_data.py'

rule optimize_transport:
    input:
        hexagons = 'data/hexagons_with_country_{country}_{plant_type}.geojson',
        technology_parameters = 'parameters/{country}/{plant_type}/technology_parameters.xlsx',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx',
        country_parameters = 'parameters/{country}/{plant_type}/country_parameters.xlsx',
        transport_parameters = 'parameters/{country}/{plant_type}/transport_parameters.xlsx',
        pipeline_parameters = 'parameters/{country}/{plant_type}/pipeline_parameters.xlsx'
    output:
        'resources/hex_transport_{country}_{plant_type}.geojson'
    script:
        'src/main/transport_optimization.py'

rule calculate_water_costs:
    input:
        technology_parameters = 'parameters/{country}/{plant_type}/technology_parameters.xlsx',
        country_parameters = 'parameters/{country}/{plant_type}/country_parameters.xlsx',
        hexagons = 'resources/hex_transport_{country}_{plant_type}.geojson'
    output:
        'resources/hex_water_{country}_{plant_type}.geojson'
    script:
        'src/main/water_cost.py'

rule optimize_plant:
    input:
        transport_parameters = 'parameters/{country}/{plant_type}/transport_parameters.xlsx',
        country_parameters = 'parameters/{country}/{plant_type}/country_parameters.xlsx',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx',
        hexagons = 'resources/hex_water_{country}_{plant_type}.geojson'
    output:
        'resources/hex_lc_{country}_{weather_year}_{plant_type}.geojson'
    script:
        'src/main/plant_optimization.py'

rule calculate_total_costs:
    input:
        hexagons = 'resources/hex_lc_{country}_{weather_year}_{plant_type}.geojson',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx'
    output:
        'results/hex_total_cost_{country}_{weather_year}_{plant_type}.geojson'
    script:
        'src/main/total_costs.py'

rule calculate_cost_components:
    input:
        hexagons = 'results/hex_total_cost_{country}_{weather_year}_{plant_type}.geojson',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx',
        country_parameters = 'parameters/{country}/{plant_type}/country_parameters.xlsx'
    output:
        'results/hex_cost_components_{country}_{weather_year}_{plant_type}.geojson',
        'results/hex_cost_components_{country}_{weather_year}_{plant_type}.csv'
    script:
        'src/main/costs_by_component.py'

rule calculate_map_costs:
    input:
        hexagons = 'results/hex_cost_components_{country}_{weather_year}_{plant_type}.geojson',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx'
    output:
        directory('plots/{country}_{weather_year}_{plant_type}')
    script:
        'src/main/map_costs.py'

