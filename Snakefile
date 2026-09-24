CHELSA_CSV = '../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/chelsa_climate.csv'
SOIL_CSV = '../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/oilgrids_sampled.csv'
OTHER_CSV = '../epfl-altshuler/GlobalAlpineSoilMicrobiome/data/other_data.csv'

ALL_DATA_CSV  = '../epfl-altshuler/GlobalAlpine/results/all_sampled_data.csv'

DATA_POINTS = 'data/points.csv'

rule all:
    input:
        ALL_DATA_CSV

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

