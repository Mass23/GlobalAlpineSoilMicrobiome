CHELSA_PARQ = 'results/download_data/chelsa_climate.parquet'
CHELSA_CSV = 'results/download_data/chelsa_climate.csv'
SOIL_CSV = 'results/download_data/soilgrids_sampled.csv'
SOIL_PARQ = 'results/download_data/soilgrids_sampled.parquet'
OTHER_CSV = 'results/download_data/other_data.csv'
OTHER_PARQ = 'results/download_data/other_data.parquet'

ALL_DATA_PARQ = 'results/all_sampled_data.parquet'
ALL_DATA_CSV  = 'results/all_sampled_data.csv'

DATA_POINTS = 'data/points.csv'

rule all:
    input:
        ALL_DATA_PARQ, ALL_DATA_CSV

rule download_chelsa:
    input:
        DATA_POINTS
    output:
        CHELSA_PARQ, CHELSA_CSV
    conda:
        "envs/conda_env_chelsa.yml"
    shell:
        """
        Rscript scripts/download_chelsa.R
        """

rule download_soilgrids:
    input:
        DATA_POINTS
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
        DATA_POINTS
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
