CHELSA_CSV = '../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/chelsa_climate.csv'
SOIL_CSV = '../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/oilgrids_sampled.csv'
OTHER_CSV = '../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/other_data.csv'

ALL_DATA_CSV  = '../epfl-altshuler/GlobalAlpine/results/all_sampled_data.csv'

DATA_POINTS = 'data/points.csv'

rule all:
    input:
        ALL_DATA_PARQ, ALL_DATA_CSV

rule download_chelsa:
    input:
        DATA_POINTS
    output:
        CHELSA_CSV
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
        SOIL_CSV
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
        OTHER_CSV
    conda:
        "envs/conda_env_others.yml"
    shell:
        """
        python scripts/download_others.py
        """

rule download_earth_data:
    input:
        CHELSA_CSV, SOIL_CSV, OTHER_CSV
    output:
        ALL_DATA_CSV
    conda:
        "envs/conda_env_merge.yml"
    shell:
        """
        Rscript scripts/merge_extract.R
        """
