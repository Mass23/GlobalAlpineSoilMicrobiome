configfile: "config/extract_config.yaml"

# Simplified pipeline: separate rules and one final merge. Each rule points at its env YAML in envs/.
FINAL_PARQ = config["output"]["path"]
FINAL_CSV = config["output"].get("r_path", "results/sampled_soilgrids_r.csv")

CONDA_CHELSA = "envs/conda_env_r_chelsa.yml"
CONDA_SOIL = "envs/conda_env_soilgrids.yml"
CONDA_OTHER = "envs/conda_env_other.yml"
CONDA_MERGE = "envs/conda_env_r_merge.yml"

CHELSA_PARQ = 'results/chelsa_climate.parquet'
CHELSA_CSV = 'results/chelsa_climate.csv'
SOIL_CSV = 'results/soilgrids_sampled.csv'
SOIL_PARQ = 'results/soilgrids_sampled.parquet'
OTHER_CSV = 'results/other_data.csv'
OTHER_PARQ = 'results/other_data.parquet'

rule all:
    input:
        FINAL_PARQ, FINAL_CSV

rule download_chelsa:
    input:
        config['points_csv']
    output:
        CHELSA_PARQ, CHELSA_CSV
    conda:
        CONDA_CHELSA
    shell:
        """
        Rscript scripts/chelsa_extract.R
        """

rule download_soilgrids:
    input:
        config['points_csv']
    output:
        SOIL_CSV, SOIL_PARQ
    conda:
        CONDA_SOIL
    shell:
        """
        python scripts/soilgrids_extract.py
        """

rule download_other_data:
    input:
        config['points_csv']
    output:
        OTHER_CSV, OTHER_PARQ
    conda:
        CONDA_OTHER
    shell:
        """
        python scripts/other_extract.py
        """

rule download_earth_data:
    input:
        CHELSA_PARQ, CHELSA_CSV, SOIL_CSV, SOIL_PARQ, OTHER_CSV, OTHER_PARQ
    output:
        FINAL_PARQ, FINAL_CSV
    conda:
        CONDA_MERGE
    shell:
        """
        Rscript scripts/merge_extract.R
        """
