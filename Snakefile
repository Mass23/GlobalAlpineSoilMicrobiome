configfile: "config/extract_config.yaml"

# Minimal workflow: one rule doing all Earth-data sampling inside R.
OUTPUT = config["output"]["path"]
OUTPUT_R = config["output"].get("r_path", "results/sampled_soilgrids_r.csv")
CONDA_ENV = "envs/conda_env_r_geodata-download.yml"

rule all:
    input:
        OUTPUT, OUTPUT_R

rule download_earth_data:
    input:
        config['points_csv']
    output:
        [OUTPUT, OUTPUT_R]
    conda:
        CONDA_ENV
    shell:
        """
        set -euo pipefail
        mkdir -p logs
        mkdir -p $(dirname {output[0]})
        mkdir -p $(dirname {output[1]})
        echo 'Running R extractor; logs -> logs/extract_r.out and logs/extract_r.err'
        Rscript scripts/extract.R
        if [ ! -f {output[0]} ] || [ ! -f {output[1]} ]; then
            echo 'Expected outputs not created:' {output[0]} {output[1]} >&2
            ls -la $(dirname {output[0]}) >&2 || true
            ls -la $(dirname {output[1]}) >&2 || true
            exit 1
        fi
        """