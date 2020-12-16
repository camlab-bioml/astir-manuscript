## All analysis for the schapiro dataset

schapiro_samples = pd.read_csv(config['schapiro']['sample_txt'], header=None)
schapiro_samples = list(schapiro_samples[schapiro_samples.columns[0]])


## Alternative quantifications
schapiro_alt_mask_samples = ['Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35']
schapiro_users = ['Catena', 'Jackson', 'Schulz']


schapiro_output = {
    'csv_rds': expand(output_path + "schapiro_processed/{sample}.csv", sample = schapiro_samples),
    'loom': output_path + "looms/schapiro.loom",
    'alt_mask_csvs': expand(output_path + "schapiro_processed_alt_mask/{alt_s}_{user}.csv",
    alt_s=schapiro_alt_mask_samples,user=schapiro_users)
}


rule read_schapiro_2017:
    input:
        csv=config['schapiro']['data_dir'] + "{sample}.csv"
    output:
        sce=output_path + "schapiro_processed/{sample}.rds",
        csv=output_path + "schapiro_processed/{sample}.csv"
    shell:
        "Rscript pipeline/schapiro-2017/schapiro-to-rds.R "
        "--input_csv {input.csv} "
        "--output_sce {output.sce} "
        "--output_csv {output.csv} "

rule read_schapiro_alt_masks:
    input:
        csv = config['schapiro']['data_dir'] + "alt-masks/{alt_s}_{user}.csv"
    output:
        sce=output_path + "schapiro_processed_alt_mask/{alt_s}_{user}.rds",
        csv=output_path + "schapiro_processed_alt_mask/{alt_s}_{user}.csv",
    shell:
        "Rscript pipeline/schapiro-2017/schapiro-to-rds.R "
        "--input_csv {input.csv} "
        "--output_sce {output.sce} "
        "--output_csv {output.csv} "      


rule schapiro_to_loom:
    input:
        expand(output_path + "schapiro_processed/{sample}.csv", sample=schapiro_samples),
    output:
        output_path + "looms/schapiro.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/schapiro_processed "
        "{output} "


rule schapiro_subset_cells:
    input:
        csvs=expand(output_path + "schapiro_processed/{core}.csv", core=schapiro_samples[0:30]),
        assignments_type=output_path + "astir_assignments/schapiro_astir_assignments.csv",
        assignments_state=output_path + "astir_assignments/schapiro_astir_assignments_state.csv",
    output:
        assignments_type=output_path + "schapiro_subset/schapiro_subset_assignments_type.csv",
        assignments_state=output_path + "schapiro_subset/schapiro_subset_assignments_state.csv",
        expression=output_path + "schapiro_subset/schapiro_subset_expression.csv"
    run:
        import pandas as pd

        core_dict = dict(zip(schapiro_samples, input.csvs))

        dfs = {c: pd.read_csv(f,index_col=0) for (c,f) in core_dict.items()}

        df = pd.concat(dfs.values())

        # df = df.sample(n=config['n_subsample'],
        # replace=False,random_state=1234)

        subset_cell_ids = list(df.index)

        df.to_csv(output.expression)

        # ## Output each to csv
        # for c in dfs.keys():
        #     df_core = dfs[c]
        #     cells_select = [ x for x in list(df.index) if x in list(df_core.index)]
        #     df_core = df_core.loc[ cells_select ]
        #     output_csv = output_path + f"schapiro_subset_separate_csvs/{c}.csv"
        #     df_core.to_csv(output_csv)

        ## Subsample type assignments
        assignments = pd.read_csv(input.assignments_type, index_col=0)
        assignments_subset = assignments.loc[subset_cell_ids]
        assignments_subset.to_csv(output.assignments_type)
        
        ## Subsample state assignments
        assignments = pd.read_csv(input.assignments_state, index_col=0)
        assignments_subset = assignments.loc[subset_cell_ids]
        assignments_subset.to_csv(output.assignments_state)

rule subset_schapiro_sce:
    params:
        input_dir = output_path + "schapiro_processed"
    input:
        expand(output_path + "schapiro_processed/{sample}.rds", sample=schapiro_samples),
        csv=output_path + "schapiro_subset/schapiro_subset_expression.csv",
    output:
        output_path + "schapiro_subset/schapiro_subset_sce.rds"
    shell:
        "Rscript pipeline/subset-sces.R "
        "--input_dir {params.input_dir} "
        "--input_csv {input.csv} "
        "--output_rds {output} "

