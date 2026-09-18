CHELSA_PARQ = 'results/chelsa_climate.parquet'
CHELSA_CSV = 'results/chelsa_climate.csv'
SOIL_CSV = 'results/soilgrids_sampled.csv'
SOIL_PARQ = 'results/soilgrids_sampled.parquet'
OTHER_CSV = 'results/other_data.csv'
OTHER_PARQ = 'results/other_data.parquet'

rule all:
    input:
        ALL_DATA_PARQ, ALL_DATA_CSV

rule download_chelsa:
    input:
        config['points_csv']
    output:
        CHELSA_PARQ, CHELSA_CSV
    conda:
        "envs/conda_env_r_chelsa.yml"
    shell:
        """
        Rscript scripts/download_chelsa.R
        """

rule download_soilgrids:
    input:
        config['points_csv']
    output:
        SOIL_CSV, SOIL_PARQ
    conda:
        "envs/conda_env_soilgrids.yml"
    shell:
        """
        python scripts/download_soilgrids.py
        """

rule download_other_data:
    input:
        config['points_csv']
    output:
        OTHER_CSV, OTHER_PARQ
    conda:
        "envs/conda_env_others.yml"
    shell:
        """
        python scripts/download_others.py
        """

rule download_earth_data:
    input:
        CHELSA_PARQ, CHELSA_CSV, SOIL_CSV, SOIL_PARQ, OTHER_CSV, OTHER_PARQ
    output:
        ALL_DATA_PARQ, ALL_DATA_CSV
    conda:
        "envs/conda_env_merge.yml"
    shell:
        """
        Rscript scripts/merge_extract.R
        """
