import os

configfile: "smk_config.yaml"

# Where the data and the results live is set in smk_config.yaml (e.g. under epfl-altshuler):
#   data_folder:    raw MicrobeAtlas inputs AND derived data (filtered tables, points, CHELSA / SoilGrids / other csv)
#   results_folder: all_sampled_data.csv and the figures (optional; default: a "results" folder next to data_folder)
DATA = config["data_folder"].rstrip("/")
RESULTS = config.get("results_folder", os.path.join(os.path.dirname(DATA), "results")).rstrip("/")

CHELSA_CSV = f"{DATA}/chelsa_data.csv"
SOIL_CSV = f"{DATA}/soilgrids_data.csv"
OTHER_CSV = f"{DATA}/other_data.csv"

SAMPLE_DATA_FILTERED = f"{DATA}/sample_data_filtered.csv"
SAMPLE_ACCESSIONS = f"{DATA}/sample_accessions.txt"
DATA_POINTS = f"{DATA}/points.csv"

# raw BIOM and the filtered one
RAW_BIOM = f"{DATA}/microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease.biom"
FILTERED_BIOM = f"{DATA}/microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom"

ALL_DATA_CSV = f"{RESULTS}/all_sampled_data.csv"

DONE_1 = f"{RESULTS}/1_alpha_overview/fig1_panels.rds"
DONE_2 = f"{RESULTS}/2_alpha_beta_elevation_latitude/.done"

rule all:
    input:
        DONE_1,
        DONE_2


rule analysis_step_1:
    input:
        SAMPLE_DATA_FILTERED
    output:
        DONE_1
    conda:
        "envs/r_analyses/environment.yml"
    params:
        data_folder=DATA,
        results_folder=RESULTS
    shell:
        """
        Rscript scripts/1_alpha_overview.R {params.data_folder} {params.results_folder}
        """

rule analysis_step_2:
    input:
        SAMPLE_DATA_FILTERED,
        ALL_DATA_CSV,
        FILTERED_BIOM
    output:
        DONE_2
    conda:
        "envs/r_analyses/environment.yml"
    params:
        data_folder=DATA,
        results_folder=RESULTS
    shell:
        """
        Rscript scripts/2_alpha_beta_elevation_latitude.R {params.data_folder} {params.results_folder}
        """

# The R script builds the sample table. Its command line is unchanged, and the accession list is NOT an
# output of this rule on purpose: Snakemake deletes a rule's existing outputs before it reruns it, so
# declaring both here would delete an existing CSV before the script could reuse it.
rule prepare_microbeatlas_data:
    output:
        SAMPLE_DATA_FILTERED
    conda:
        "envs/microbeatlas/environment.yml"
    params:
        data_folder=DATA
    shell:
        """
        Rscript scripts/microbeatlas_filter_data.R {params.data_folder}
        """

# One accession per line, taken from the sample table (the script only reads the CSV when it exists)
rule sample_accessions:
    input:
        SAMPLE_DATA_FILTERED
    output:
        SAMPLE_ACCESSIONS
    conda:
        "envs/microbeatlas/environment.yml"
    params:
        data_folder=DATA
    shell:
        """
        Rscript scripts/microbeatlas_filter_data.R {params.data_folder}
        """

# Subset the BIOM to the filtered samples. -s takes a FILE with one ID per line (not a comma-separated
# string); this also avoids the 128 KB limit on a single command-line argument.
rule filter_biom:
    input:
        biom = RAW_BIOM,
        samples = SAMPLE_ACCESSIONS
    output:
        FILTERED_BIOM
    conda:
        "envs/biom_filter/environment.yml"
    shell:
        """
        biom subset-table -i {input.biom} -a sample -s {input.samples} -o {output}
        """

rule prepare_points:
    input:
        SAMPLE_DATA_FILTERED
    output:
        DATA_POINTS
    run:
        import pandas as pd

        # needs a QUOTED csv: run microbeatlas_filter_data.R once with REWRITE_QUOTED <- TRUE if it is not
        df = pd.read_csv(input[0])

        points = df[["Sample_acc", "Latitude", "Longitude", "Depth_m", "Date_collected"]].copy()
        points.columns = ["site", "lat", "lon", "depth_m", "sample_date"]
        points["depth_cm"] = points["depth_m"] * 100
        points = points.drop(columns=["depth_m"])      # was "depth_M": the column is lower case
        points = points.drop_duplicates()              # several read sets can share one sample

        points.to_csv(output[0], index=False)

rule download_chelsa:
    input:
        DATA_POINTS
    output:
        CHELSA_CSV
    conda:
        "envs/chelsa/environment.yml"
    shell:
        """
        python3 scripts/download_chelsa.py
        """

rule download_soilgrids:
    input:
        DATA_POINTS
    output:
        SOIL_CSV
    conda:
        "envs/soilgrids/environment.yml"
    shell:
        """
        python3 scripts/download_soilgrids.py
        """

rule download_other_data:
    input:
        DATA_POINTS
    output:
        OTHER_CSV
    conda:
        "envs/others/environment.yml"
    shell:
        """
        python3 scripts/download_others.py
        """

rule download_earth_data:
    input:
        CHELSA_CSV, SOIL_CSV, OTHER_CSV
    output:
        ALL_DATA_CSV
    run:
        import pandas as pd
        from pathlib import Path

        chelsa = pd.read_csv(input[0])
        soil = pd.read_csv(input[1])
        other = pd.read_csv(input[2])

        key_cols = [
            col for col in chelsa.columns
            if col in soil.columns and col in other.columns
        ]

        if not key_cols:
            raise ValueError(
                "The CHELSA/SoilGrids/Other tables share no columns to join on."
            )

        print("Joining on:", ", ".join(key_cols))

        merged = chelsa.merge(
            soil,
            how="outer",
            on=key_cols,
            suffixes=("", ".soil"),
        )

        merged = merged.merge(
            other,
            how="outer",
            on=key_cols,
            suffixes=("", ".other"),
        )

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output[0], index=False)

        print(
            f"Wrote {output[0]}: "
            f"{len(merged)} rows, {len(merged.columns)} columns"
        )
