configfile: "config/extract_config.yaml"

# Minimal workflow: one rule doing all Earth-data sampling inside R.
OUTPUT_R = config["output"].get("r_path", "results/sampled_soilgrids_r.csv")
CONDA_ENV = "envs/conda_env_r_geodata.yml"

rule all:
    output:
        OUTPUT_R

rule download_earth_data:
    output:
        OUTPUT_R
    conda:
        CONDA_ENV
    shell:
        "Rscript scripts/extract_r.R {output[0]}"
