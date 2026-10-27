import csv
import os

configfile: "smk_config.yaml"

# Where the data and the results live is set in smk_config.yaml (e.g. under epfl-altshuler):
#   data_folder:    raw MicrobeAtlas inputs AND derived data (filtered table, points, CHELSA / SoilGrids / other csv)
#   results_folder: all_sampled_data.csv and the figures (optional; default: a "results" folder next to data_folder)
DATA = config["data_folder"].rstrip("/")
RESULTS = config.get("results_folder", os.path.join(os.path.dirname(DATA), "results")).rstrip("/")

CHELSA_CSV = f"{DATA}/chelsa_data.csv"
SOIL_CSV = f"{DATA}/soilgrids_data.csv"
OTHER_CSV = f"{DATA}/other_data.csv"

# ALL FILE NAMES ARE DEFINED HERE ONLY: the rules hand the exact paths to the scripts (as arguments), so a name
# can be changed in this block without touching any script.
# The sample table is a TSV (its text fields contain commas): read it with quote = "" / QUOTE_NONE
SAMPLE_DATA_FILTERED = f"{DATA}/sample_data_filtered.tsv"
SAMPLE_ACCESSIONS = f"{DATA}/sample_accessions.txt"
DATA_POINTS = f"{DATA}/points.csv"
LOWLAND_SETS = f"{DATA}/lowland_sets.csv"          # one row per sample: group + the lowland replicates (rep_01 ... rep_20)

# raw BIOM and the filtered one
RAW_BIOM = f"{DATA}/microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease.biom"
FILTERED_BIOM = f"{DATA}/microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom"

ALL_DATA_CSV = f"{RESULTS}/all_sampled_data.csv"

# the flag files below are written by the analysis scripts themselves, inside the folder they are given
DONE_1 = f"{RESULTS}/1_alpha_overview/fig1_panels.rds"
DONE_2 = f"{RESULTS}/2_alpha_beta_elevation_latitude/.done"
DONE_3 = f"{RESULTS}/3_distance_decay/.done"
DONE_4 = f"{RESULTS}/4_variance_partitioning/.done"
DONE_5 = f"{RESULTS}/5_family_differential_abundance/.done"
DONE_6 = f"{RESULTS}/6_alpine_specialist_otus/.done"

# The first rule is the default target of `snakemake`: everything
rule all:
    input:
        DONE_1,
        DONE_2,
        DONE_3,
        DONE_4,
        DONE_5,
        DONE_6

# Script 3: distance decay of community dissimilarity in alpine / peri-alpine / lowland samples. Needs perialpine_alpine and cluster_id (written
# by prepare_microbeatlas_data) and the filtered BIOM; no climate / soil / terrain tables
rule analysis_step_3:
    input:
        table=SAMPLE_DATA_FILTERED,
        biom=FILTERED_BIOM
    output:
        DONE_3
    conda:
        "envs/r_analyses/environment.yml"
    params:
        out_folder=os.path.dirname(DONE_3)
    shell:
        """
        Rscript scripts/3_distance_decay.R {input.table} {input.biom} {params.out_folder}
        """

# Script 4: variation partitioning of the alpine microbiome (composition and alpha diversity) between the environment (climate + soil + terrain,
# all_sampled_data.csv) and geography; needs the soil download through all_sampled_data.csv
rule analysis_step_4:
    input:
        table=SAMPLE_DATA_FILTERED,
        env=ALL_DATA_CSV,
        biom=FILTERED_BIOM
    output:
        DONE_4
    conda:
        "envs/r_analyses/environment.yml"
    params:
        out_folder=os.path.dirname(DONE_4)
    shell:
        """
        Rscript scripts/4_variance_partitioning.R {input.table} {input.env} {input.biom} {params.out_folder}
        """


# Script 5: taxonomy (family level) and differential abundance alpine vs lowland (GAMs on CLR abundances); needs the taxonomy file in data/microbeatlas
rule analysis_step_5:
    input:
        table=SAMPLE_DATA_FILTERED,
        env=ALL_DATA_CSV,
        biom=FILTERED_BIOM,
        step4=DONE_4                                   # panel C reads the envfit results of step 4
    output:
        DONE_5
    conda:
        "envs/r_analyses/environment.yml"
    threads: 4
    params:
        out_folder=os.path.dirname(DONE_5),
        envfit=os.path.join(os.path.dirname(DONE_4), "tables", "envfit_summary.csv")
    shell:
        """
        Rscript scripts/5_family_differential_abundance.R {input.table} {input.env} {input.biom} {params.out_folder} {params.envfit}
        """


# Script 6: alpine specialist OTUs (niche breadth against a permutation null + alpine affinity): taxonomy (A), space (B), season and soil depth (C). Needs the taxonomy file in data/microbeatlas.
rule analysis_step_6:
    input:
        table=SAMPLE_DATA_FILTERED,
        env=ALL_DATA_CSV,
        biom=FILTERED_BIOM
    output:
        DONE_6
    conda:
        "envs/r_analyses/environment.yml"
    threads: 4
    params:
        out_folder=os.path.dirname(DONE_6)
    shell:
        """
        Rscript scripts/6_alpine_specialist_otus.R {input.table} {input.env} {input.biom} {params.out_folder}
        """


rule analysis_step_1:
    input:
        table=SAMPLE_DATA_FILTERED
    output:
        DONE_1
    conda:
        "envs/r_analyses/environment.yml"
    params:
        out_folder=os.path.dirname(DONE_1)
    shell:
        """
        Rscript scripts/1_alpha_overview.R {input.table} {params.out_folder}
        """

rule analysis_step_2:
    input:
        table=SAMPLE_DATA_FILTERED,
        env=ALL_DATA_CSV,
        biom=FILTERED_BIOM
    output:
        DONE_2
    conda:
        "envs/r_analyses/environment.yml"
    params:
        out_folder=os.path.dirname(DONE_2)
    shell:
        """
        Rscript scripts/2_alpha_beta_elevation_latitude.R {input.table} {input.env} {input.biom} {params.out_folder}
        """

# The R script writes both outputs: the sample table (TSV) and the list of accessions to keep
rule prepare_microbeatlas_data:
    output:
        table=SAMPLE_DATA_FILTERED,
        accessions=SAMPLE_ACCESSIONS
    conda:
        "envs/microbeatlas/environment.yml"
    params:
        data_folder=DATA
    shell:
        """
        Rscript scripts/microbeatlas_filter_data.R {params.data_folder} {output.table} {output.accessions}
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

# points.csv: one row per sample x location x depth x date (site, lat, lon, sample_date, depth_cm)
rule prepare_points:
    input:
        SAMPLE_DATA_FILTERED
    output:
        DATA_POINTS
    run:
        import pandas as pd

        # tab-separated, no quoting: text fields can contain commas and double quotes
        df = pd.read_csv(input[0], sep="\t", quoting=csv.QUOTE_NONE, low_memory=False)

        points = df[["Sample_acc", "Latitude", "Longitude", "Depth_m", "Date_collected"]].copy()
        points.columns = ["site", "lat", "lon", "depth_m", "sample_date"]
        points["depth_cm"] = points["depth_m"] * 100
        points = points.drop(columns=["depth_m"])      # was "depth_M" in an earlier version: lower case
        points = points.drop_duplicates()              # several read sets can share one sample

        points.to_csv(output[0], index=False)

# Lowland replicates (needs the DEM elevation, so it runs after download_other_data): n_rep random draws of lowland samples that are
# EVENLY SPACED over latitude / longitude (5-degree cells, one sample per cell in turn, at most one per thin_m metres), each with as
# many lowland samples as there are alpine samples (seed-controlled, reproducible). Columns rep_01 ... rep_<n_rep> (yes / no);
# all alpine samples are "yes" in every replicate; subset_any = alpine or in at least one replicate (what the soil download
# fetches: the UNION of the replicates; see lowland_sets_summary.csv for its size). Optional keys in smk_config.yaml (defaults shown):
#   lowland_sets_seed: 23        lowland_sets_ratio: 1.0 (lowland samples per alpine sample)
#   lowland_sets_cell_deg: 5     lowland_sets_n_rep: 20       lowland_sets_thin_m: 100
rule define_lowland_sets:
    input:
        table=SAMPLE_DATA_FILTERED,
        other=OTHER_CSV
    output:
        LOWLAND_SETS
    conda:
        "envs/others/environment.yml"
    params:
        seed=config.get("lowland_sets_seed", 23),
        ratio=config.get("lowland_sets_ratio", 1.0),
        cell=config.get("lowland_sets_cell_deg", 5),
        n_rep=config.get("lowland_sets_n_rep", 20),
        thin=config.get("lowland_sets_thin_m", 100)
    shell:
        """
        python3 scripts/define_lowland_sets.py {input.table} {input.other} {output} \
            --seed {params.seed} --ratio {params.ratio} --cell-deg {params.cell} --n-rep {params.n_rep} --thin-m {params.thin}
        """

rule download_chelsa:
    input:
        DATA_POINTS
    output:
        CHELSA_CSV
    conda:
        "envs/chelsa/environment.yml"
    shell:
        """
        python3 scripts/download_chelsa.py {input} {output}
        """

# Soil is fetched for the samples of lowland_sets.csv with subset_any = yes only (the others stay NA).
# soil_properties (optional, in smk_config.yaml): the SoilGrids properties to fetch, e.g.
#   soil_properties: [phh2o, soc, nitrogen, clay, sand, silt, cec, bdod]
# Each property costs 3 requests per coordinate, so fewer properties = a much shorter run. Empty / absent = all 13.
rule download_soilgrids:
    input:
        points=DATA_POINTS,
        sets=LOWLAND_SETS
    output:
        SOIL_CSV
    conda:
        "envs/soilgrids/environment.yml"
    params:
        properties=",".join(config.get("soil_properties", []))
    shell:
        """
        python3 scripts/download_soilgrids.py {input.points} {output} "{params.properties}" {input.sets}
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
        python3 scripts/download_others.py {input} {output}
        """

rule download_earth_data:
    input:
        chelsa=CHELSA_CSV,
        soil=SOIL_CSV,
        other=OTHER_CSV,
        sets=LOWLAND_SETS
    output:
        ALL_DATA_CSV
    run:
        import pandas as pd
        from pathlib import Path

        chelsa = pd.read_csv(input.chelsa)
        soil = pd.read_csv(input.soil)
        other = pd.read_csv(input.other)

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

        # which samples are alpine, and in which lowland replicates (rep_01 ...) they are drawn (one row per sample, joined by site)
        sets = pd.read_csv(input.sets)
        sets = sets[[c for c in sets.columns if c in ("site", "group", "n_rep_selected", "subset_any") or c.startswith("rep_")]]
        merged = merged.merge(sets, on="site", how="left")

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output[0], index=False)

        print(
            f"Wrote {output[0]}: "
            f"{len(merged)} rows, {len(merged.columns)} columns"
        )
