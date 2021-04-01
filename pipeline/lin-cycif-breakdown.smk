

lin_metadata = pd.read_csv("data-raw/metadata/lin-cycif-metadata.tsv", sep = "\t")

cancer_types = lin_metadata["Anatomic site"].unique()

lin_cycif_breakdown = {
    'sces': expand(output_dir + "sces/Lin-cycif-{cancer}.rds", cancer = cancer_types) 
}

rule create_sces:
    input:
        expand(output_path + "lin-cycif_processed/{core}.rds", core = lin_metadata[lin_metadata["Anatomic site"] == {cancer}]["core"].to_numpy())
    output:
        output_dir + "sces/Lin-cycif-{cancer}.rds"
    shell:
        "Rscript -e"