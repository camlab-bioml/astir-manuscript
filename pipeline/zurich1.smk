## Analysis for zurich1 dataset


zurich1_metadata = pd.read_csv(os.path.join(config['zurich1']['base_dir'], config['zurich1']['metadata_file']))
zurich1_cores = list(zurich1_metadata.core)
tmp_zurich1_output = expand(output_path + "zurich1_processed/{core}.rds", core=zurich1_cores)
zurch1_subset_cell_ids = cell_ids=output_path + "zurich1_subset/zurich1_subset_cells.csv"

zurich1_output = {
    'rds_csv': tmp_zurich1_output,
    'loom': output_path + "looms/zurich1.loom"
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




rule subset_zurich1:
    input:
        loom=output_path + "looms/zurich1.loom",
        sc_dat="data-raw/jackson-2020/SingleCell_and_Metadata/ZurichTMA/" + "SC_dat.csv",
        assignments=output_path + "astir_assignments/zurich1_astir_assignments.csv",
    output:
        cell_ids=output_path + "zurich1_subset/zurich1_subset_cells.csv",
        assignments=output_path + "zurich1_subset/zurich1_subset_assignments.csv",
    run:
        import numpy as np
        import pandas as pd
        sc_dat = pd.read_csv(input.sc_dat)

        ## Cell IDs we're going to subsample
        cell_ids = sc_dat.id.unique()
        np.random.seed(1234)
        cell_ids = list(np.random.choice(cell_ids, config['zurich1']['n_subsample'], replace=False))
        cell_df = pd.DataFrame({'id': cell_ids})
        cell_df.to_csv(output.cell_ids)

        ## Subsample assignments
        assignments = pd.read_csv(input.assignments, index_col=0)
        assignments_subset = assignments.loc[cell_ids]
        assignments_subset.to_csv(output.assignments)

