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
    'loom': output_path + "looms/basel.loom",
    'subset': output_path + "basel_subset/basel_subset_expression.csv"
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


rule basel_subset_cells:
    input:
        csvs=expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
        assignments=output_path + "astir_assignments/basel_astir_assignments.csv",
    output:
        assignments=output_path + "basel_subset/basel_subset_assignments.csv",
        expression=output_path + "basel_subset/basel_subset_expression.csv"
    run:
        import pandas as pd

        dfs = [pd.read_csv(f,index_col=0) for f in input.csvs]
        df = pd.concat(dfs)

        df = df.sample(n=config['n_subsample'],
        replace=False,random_state=1234)

        subset_cell_ids = list(df.index)

        df.to_csv(output.expression)

        ## Subsample assignments
        assignments = pd.read_csv(input.assignments, index_col=0)
        assignments_subset = assignments.loc[subset_cell_ids]
        assignments_subset.to_csv(output.assignments)


