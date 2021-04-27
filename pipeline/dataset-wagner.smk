## All analysis for the Wagner et al. dataset


# Wagner 2019 metadata
wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])


tmp_wagner_output = expand(output_path + 
        "wagner_processed/{sample}.{ext}",
        sample=wagner_samples, ext=['csv','rds'])

wagner_output = {
    'rds_csv': tmp_wagner_output,
    'loom': output_path + "looms/wagner.loom",
    # 'clusters': output_path + "wagner_clusters/wagner_clusters.tsv"
    # 'adata': output_path + "anndata/wagner.h5ad"
    # 'subset': output_path + "wagner_subset/wagner_subset_expression.csv"# ,
    # 'subset_sce': output_path + "wagner_subset/wagner_subset_sce.rds"
}


rule read_wagner_2019:
    input:
        fcs=config['wagner']['fcs_dir'] + "/{sample}.fcs",
        # cluster_identities = output_path + "wagner_clusters/wagner_clusters.tsv",
    output:
        rds = output_path + "wagner_processed/{sample}.rds",
        csv = output_path + "wagner_processed/{sample}.csv",
    shell:
        "Rscript pipeline/wagner/wagner-fcs-to-csv.R "
        "--input_fcs {input.fcs} "
        # "--cluster_identities {input.cluster_identities} "
        "--output_rds {output.rds} "
        "--output_csv {output.csv} "

rule wagner_to_loom:
    input:
        expand(output_path + "wagner_processed/{core}.csv", core=wagner_samples),
    output:
        output_path + "looms/wagner.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/wagner_processed "
        "{output} "

# rule wagner_to_ad:
#     input:
#         wagner_output['rds_csv'],
#     output:
#         wagner_output['adata'],
#     shell:
#         "R pipeline/cla/dir-of-csvs-to-scanpy.py "
#         "{output_path}/wagner_processed "
#         "{output} "

# rule get_wagner_clusters:
#     output:
#         output_path + "wagner_clusters/wagner_clusters.tsv"
#     shell:
#         "Rscript pipeline/wagner/get-wagner-clusters.R "
#         "--input_yaml imc-config.yml "
#         "--output_tsv {output} "


# rule subset_wagner:
#     input:
#         csvs=expand(output_path + "wagner_processed/{core}.csv", core=wagner_samples),
#         assignments_type=output_path + "astir_assignments/wagner_astir_assignments.csv",
#         assignments_state=output_path + "astir_assignments/wagner_astir_assignments_state.csv",
#     output:
#         assignments_type=output_path + "wagner_subset/wagner_subset_assignments_type.csv",
#         assignments_state=output_path + "wagner_subset/wagner_subset_assignments_state.csv",
#         expression=output_path + "wagner_subset/wagner_subset_expression.csv"
#     run:
#         import pandas as pd

#         dfs = [pd.read_csv(f,index_col=0) for f in input.csvs]
#         df = pd.concat(dfs)

#         df = df.sample(n=config['n_subsample'],
#         replace=False,random_state=1234)

#         subset_cell_ids = list(df.index)

#         df.to_csv(output.expression)

#         ## Subsample type assignments
#         assignments = pd.read_csv(input.assignments_type, index_col=0)
#         assignments_subset = assignments.loc[subset_cell_ids]
#         assignments_subset.to_csv(output.assignments_type)

#         ## Subsample state assignments
#         assignments = pd.read_csv(input.assignments_state, index_col=0)
#         assignments_subset = assignments.loc[subset_cell_ids]
#         assignments_subset.to_csv(output.assignments_state)

rule subset_wagner_sce:
    params:
        input_dir = output_path + "wagner_processed"
    input:
        expand(output_path + "wagner_processed/{sample}.rds", sample=wagner_samples),
        csv=output_path + "wagner_subset/wagner_subset_expression.csv",
    output:
        output_path + "wagner_subset/wagner_subset_sce.rds"
    shell:
        "Rscript pipeline/subset-sces.R "
        "--input_dir {params.input_dir} "
        "--input_csv {input.csv} "
        "--output_rds {output} "