## All analysis for the Basel dataset

basel_metadata = pd.read_csv(os.path.join(config['basel']['base_dir'], config['basel']['metadata_file']))
basel_cores = list(basel_metadata.core)
tmp_basel_output = expand(output_path + "basel_processed/{core}.rds", core=basel_cores)


# Basel 15k reduced set for benchmarking
basel_15k_cells = pd.read_csv("data-raw/cellsets/cell_ids_15k.csv")
cores_15k = list(set(list(basel_15k_cells.core)))
basel_15k_csvs = expand(output_path + "basel_15k_subset/{core}.csv", core=cores_15k)

basel_output = {
    'csv_rds': tmp_basel_output,
    'loom': output_path + "looms/basel.loom"
}


rule read_jackson_2020_basel:
    input:
        scdat=config['basel']['base_dir'] + "SC_dat.csv",
        scloc=config['basel']['base_dir'] + "Basel_SC_locations.csv",
    output:
        expand(output_path + "basel_processed/{core}.rds", core=basel_cores),
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
    shell:
        "Rscript pipeline/jackson-2020/jackson-raw-to-sce-server.R "
        "--input_sc {input.scdat} "
        "--input_loc {input.scloc} "
        "--output_dir {output_path}/basel_processed "


rule basel_to_loom:
    input:
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
    output:
        output_path + "looms/basel.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/basel_processed "
        "{output} "


rule subset_basel_cells:
    input:
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
        cell_labels="data-raw/cellsets/cell_ids_15k.csv"
    output:
        basel_15k_csvs
    shell:
        "mkdir -p {output} && "
        "python pipeline/subset-expression-data.py "
        "--input_dir {output_path}/basel_processed "
        "--cell_labels {input.cell_labels} "
        "--output_dir {output_path}/basel_15k_subset "


rule basel_to_15k:
    input:
        astir_output=output_path + "astir_assignments/basel_astir_assignments.csv",
        cell_list="data-raw/15k_cells.csv",
    output:
        output_path + "summarized_assigments/basel_15k_all.csv"
    shell:
        "python pipeline/subset-to-15k.py "
        "{input.astir_output} "
        "{input.cell_list} "
        "{output}"