configfile: "config/extract_config.yaml"

OUTPUT = config["output"]["path"]

rule all:
    output:
        OUTPUT

rule extract:
    output:
        OUTPUT
    conda:
        "envs/conda_env.yml"
    shell:
        "python scripts/run_extract.py config/extract_config.yaml {output[0]}"
