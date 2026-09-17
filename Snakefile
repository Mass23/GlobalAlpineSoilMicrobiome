configfile: "config/extract_config.yaml"

OUTPUT = config["output"]["path"]
OUTPUT_R = config["output"].get("r_path", "results/sampled_soilgrids_r.csv")

# Collect remote URLs from config and map to local tile paths (data/tiles/{basename}).
# Downloads are marked as temporary so Snakemake will remove them after the workflow.

def _collect_urls(cfg):
    urls = []
    targets = []
    import os
    def add(p):
        if isinstance(p, str) and (p.startswith('http://') or p.startswith('https://')):
            fn = os.path.basename(p.split('?')[0])
            tgt = os.path.join('data', 'tiles', fn)
            urls.append(p)
            targets.append(tgt)
    # soil_vars nested
    for var_map in cfg.get('soil_vars', {}).values():
        for p in var_map.values():
            add(p)
    for p in cfg.get('chelsa', {}).values():
        add(p)
    dem = cfg.get('dem')
    if isinstance(dem, str):
        add(dem)
    for p in cfg.get('landcover', {}).values():
        add(p)
    for p in cfg.get('ndvi', {}).values():
        add(p)
    for p in cfg.get('snow', {}).values():
        add(p)
    # deduplicate preserving order
    seen = set()
    n_urls = []
    n_targets = []
    for u,t in zip(urls, targets):
        if u not in seen:
            seen.add(u)
            n_urls.append(u)
            n_targets.append(t)
    return n_urls, n_targets

DOWNLOAD_URLS, DOWNLOAD_TARGETS = _collect_urls(config)

# mark downloaded tiles as temporary
DOWNLOAD_TARGETS_TEMP = temp(DOWNLOAD_TARGETS) if DOWNLOAD_TARGETS else []

rule all:
    output:
        OUTPUT_R

# Marker file to indicate CRAN-only R packages have been installed into the R env
R_PKG_DONE = 'envs/.r_packages_installed'

rule download_tiles:
    # outputs are temporary tiles
    output:
        DOWNLOAD_TARGETS_TEMP
    run:
        import os, subprocess
        os.makedirs('data/tiles', exist_ok=True)
        for url, out in zip(DOWNLOAD_URLS, output):
            if os.path.exists(out):
                continue
            print(f"Downloading {url} -> {out}")
            subprocess.check_call(["curl", "-fSL", "--retry", "3", "-o", out, url])

rule install_r_packages:
    output:
        R_PKG_DONE
    conda:
        "envs/conda_env_r_geodata.yml"
    run:
        import subprocess, os
        # Run the installer script in the repo root; it installs CRAN-only R packages into the env created by Snakemake.
        print('Installing CRAN-only R packages inside conda env...')
        subprocess.check_call(['bash', 'scripts/install_r_packages.sh'])
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        open(output[0], 'w').close()

rule extract_r:
    input:
        DOWNLOAD_TARGETS,
        R_PKG_DONE
    output:
        OUTPUT_R
    conda:
        "envs/conda_env_r_geodata.yml"
    shell:
        "Rscript scripts/extract_r.R {output[0]}"
