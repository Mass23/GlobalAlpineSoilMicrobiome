configfile: "config/extract_config.yaml"

# Minimal workflow: one rule doing all Earth-data sampling inside R.
OUTPUT = config["output"]["path"]
OUTPUT_R = config["output"].get("r_path", "results/sampled_soilgrids_r.csv")
CONDA_ENV = "envs/conda_env_r_geodata-download.yml"

rule all:
    output:
        [OUTPUT, OUTPUT_R]

R_PKG_DONE = 'envs/.r_packages_installed'

rule install_r_packages:
    output:
        R_PKG_DONE
    conda:
        CONDA_ENV
    run:
        import subprocess, os
        print('Installing CRAN/GitLab R packages inside conda env...')
        subprocess.check_call(['bash', 'scripts/install_r_packages.sh'])
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        open(output[0], 'w').close()

rule download_earth_data:
    input:
        R_PKG_DONE, config['points_csv']
    output:
        [OUTPUT, OUTPUT_R]
    conda:
        CONDA_ENV
    run:
        import os, subprocess
        # ensure log and result dirs exist so Snakemake can observe outputs
        os.makedirs('logs', exist_ok=True)
        os.makedirs(os.path.dirname(str(output[0])), exist_ok=True)
        os.makedirs(os.path.dirname(str(output[1])), exist_ok=True)
        cmd = ['Rscript', 'scripts/extract_r.R']
        with open('logs/extract_r.out', 'wb') as out, open('logs/extract_r.err', 'wb') as err:
            print('Running R extractor; logs -> logs/extract_r.{out,err}')
            subprocess.check_call(cmd, stdout=out, stderr=err)
        # verify outputs exist
        for f in output:
            if not os.path.exists(str(f)):
                raise Exception('Expected output not created: %s' % f)

