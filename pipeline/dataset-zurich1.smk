## Analysis for zurich1 dataset


zurich1_metadata = pd.read_csv(os.path.join(config['zurich1']['base_dir'], config['zurich1']['metadata_file']))
zurich1_cores = list(zurich1_metadata.core)
tmp_zurich1_output = expand(output_path + "zurich1_processed/{core}.rds", core=zurich1_cores)
zurch1_subset_cell_ids = cell_ids=output_path + "zurich1_subset/zurich1_subset_cells.csv"

zurich1_output = {
    'rds_csv': tmp_zurich1_output,
    'loom': output_path + "looms/zurich1.loom",
    'subset': output_path + "zurich1_subset/zurich1_subset_expression.csv",
    'subset_sce': output_path + "zurich1_subset/zurich1_subset_sce.rds"
}


rule read_jackson_2020_zurich1:
    input:
        scdat=config['zurich1']['base_dir'] + "SC_dat.csv",
        scloc=config['zurich1']['base_dir'] + "Zurich_SC_locations.csv",
    output:
        expand(output_path + "zurich1_processed/{core}.rds", core=zurich1_cores),
        expand(output_path + "zurich1_processed/{core}.csv", core=zurich1_cores),
    shell:
        "Rscript pipeline/jackson-2020/jackson-raw-to-sce-server.R "
        "--input_sc {input.scdat} "
        "--input_loc {input.scloc} "
        "--output_dir {output_path}/zurich1_processed "


rule zurich1_to_loom:
    input:
        expand(output_path + "zurich1_processed/{core}.csv", core=zurich1_cores),
    output:
        output_path + "looms/zurich1.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/zurich1_processed "
        "{output} "



rule zurich1_subset_cells:
    input:
        csvs=expand(output_path + "zurich1_processed/{core}.csv", core=zurich1_cores),
        assignments_type=output_path + "astir_assignments/zurich1_astir_assignments.csv",
        assignments_state=output_path + "astir_assignments/zurich1_astir_assignments_state.csv",
    output:
        assignments_type=output_path + "zurich1_subset/zurich1_subset_assignments_type.csv",
        assignments_state=output_path + "zurich1_subset/zurich1_subset_assignments_state.csv",
        expression=output_path + "zurich1_subset/zurich1_subset_expression.csv"
    run:
        import pandas as pd

        dfs = [pd.read_csv(f,index_col=0) for f in input.csvs]
        df = pd.concat(dfs)

        df = df.sample(n=config['n_subsample'],
        replace=False,random_state=1234)

        subset_cell_ids = list(df.index)

        df.to_csv(output.expression)

        ## Subsample type assignments
        assignments = pd.read_csv(input.assignments_type, index_col=0)
        assignments_subset = assignments.loc[subset_cell_ids]
        assignments_subset.to_csv(output.assignments_type)

        ## Subsample astatessignments
        assignments = pd.read_csv(input.assignments_state, index_col=0)
        assignments_subset = assignments.loc[subset_cell_ids]
        assignments_subset.to_csv(output.assignments_state)

rule subset_zurich1_sce:
    params:
        input_dir = output_path + "zurich1_processed"
    input:
        expand(output_path + "zurich1_processed/{core}.rds", core=zurich1_cores),
        csv=output_path + "zurich1_subset/zurich1_subset_expression.csv",
    output:
        output_path + "zurich1_subset/zurich1_subset_sce.rds"
    shell:
        "Rscript pipeline/subset-sces.R "
        "--input_dir {params.input_dir} "
        "--input_csv {input.csv} "
        "--output_rds {output} "


