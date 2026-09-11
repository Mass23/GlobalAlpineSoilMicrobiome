configfile: "config/extract_config.yaml"

OUTPUT = config["output"]["path"]

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
        OUTPUT

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

rule extract:
    input:
        DOWNLOAD_TARGETS
    output:
        OUTPUT
    conda:
        "envs/conda_env.yml"
    shell:
        "python scripts/run_extract.py config/extract_config.yaml {output[0]}"
