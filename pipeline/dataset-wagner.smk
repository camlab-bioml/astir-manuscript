## All analysis for the Wagner et al. dataset


# Wagner 2019 metadata
wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])


tmp_wagner_output = expand(output_path + 
        "wagner-2019_processed/{sample}.{ext}",
        sample=wagner_samples, ext=['csv','rds'])

wagner_output = {
    'rds_csv': tmp_wagner_output,
    'loom': output_path + "looms/wagner.loom",
    'subset': output_path + "wagner_subset/wagner_subset_expression.csv"
}


rule read_wagner_2019:
    input:
        config['wagner']['fcs_dir'] + "/{sample}.fcs"
    output:
        rds = output_path + "wagner-2019_processed/{sample}.rds",
        csv = output_path + "wagner-2019_processed/{sample}.csv",
    shell:
        "Rscript pipeline/wagner-2019/wagner-fcs-to-csv.R "
        "--input_fcs {input} "
        "--output_rds {output.rds} "
        "--output_csv {output.csv} "

rule wagner_to_loom:
    input:
        expand(output_path + "wagner-2019_processed/{core}.csv", core=wagner_samples),
    output:
        output_path + "looms/wagner.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/wagner-2019_processed "
        "{output} "



rule subset_wagner:
    input:
        csvs=expand(output_path + "wagner-2019_processed/{core}.csv", core=wagner_samples),
        assignments_type=output_path + "astir_assignments/wagner_astir_assignments.csv",
        assignments_state=output_path + "astir_assignments/wagner_astir_assignments_state.csv",
    output:
        assignments_type=output_path + "wagner_subset/wagner_subset_assignments_type.csv",
        assignments_state=output_path + "wagner_subset/wagner_subset_assignments_state.csv",
        expression=output_path + "wagner_subset/wagner_subset_expression.csv"
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

        ## Subsample state assignments
        assignments = pd.read_csv(input.assignments_state, index_col=0)
        assignments_subset = assignments.loc[subset_cell_ids]
        assignments_subset.to_csv(output.assignments_state)