
# ------------------ Cell level accuracy (CLA) quantifaction ----------------- #

def get_core_list(cohort):
    if cohort == "basel":
        return basel_cores
    elif cohort == "zurich1":
        return zurich1_cores
    elif cohort == "wagner":
        return wagner_samples
    elif cohort == "schapiro":
        return schapiro_samples
    elif cohort == "lin-cycif":
        return cycif_samples
    else:
        return None

# def get_core_list_acdc(cohort):
#     if cohort == "basel":
#         return basel_cores
#     elif cohort == "zurich1":
#         return zurich1_cores
#     elif cohort == "wagner":
#         return wagner_samples
#     elif cohort == "Lin-CyCif":
#         return cycif_samples
#     elif cohort == "schapiro":
#         return schapiro_samples
#     else:
#         return None

cla_methods = ['cytofLDA', 'acdc-absent', 'acdc-no-consider']
cla_cohorts = config['cla'].keys()

annotation_results = []



for cohort in cla_cohorts:
    for annotator in config['cla'][cohort]['annotators'].keys():
        annotation_results.append(
            output_path + f"cla/annotations_cytofLDA_{cohort}_{annotator}.tsv"
        )

for method in [ 'acdc-absent', 'acdc-no-consider']:
    for cohort in cla_cohorts:
        if cohort != "lin-cycif":
            annotation_results.append(
                output_path + f"cla/annotations_{method}_{cohort}.tsv"
            )

cla_outputs = {
    # 'jackson_annotations_fixed': [output_path + "cla/Basel_annotation_fixed.csv", output_path + "cla/Zurich_annotation_fixed.csv"],
    'anndata': expand(output_path + "anndata/{cohort}.h5ad", cohort=['basel','zurich1', 'wagner']),
    'annotation_results': annotation_results,
    'plots': expand(output_path + "cla/output_figs_tables/cla_{what}_{cohort}.png", cohort=['zurich1','basel','wagner', 'lin-cycif'], what=['annotation']),
    'wag_clear': output_path + "cla/wagner_clusters_by_cell_id.csv",
    'wag_clus': expand(output_path + "cla/output_figs_tables/wagner/cla_{what}_wagner.png",what=['cluster']),
    'plots2': expand(output_path + "cla/output_figs_tables/cla_{what}_{cohort}.png", cohort=['zurich1','basel'], what=['cluster']),
    'lin_subset': output_path + 'cla/lin-cycif-subset.h5ad',
    'lin_cluster_results': output_path + "cla/output_figs_tables/lin-cycif/cla_cluster_lin-cycif.tsv",
    'annotation_fig': output_path + "figures/cla/annotation.pdf",
    'clustering_fig': output_path + "figures/cla/clustering.pdf",
    'phenograph_fig': output_path + "figures/cla/phenograph.pdf",
    'lin_lung_clustered': output_path + 'cla/lin-lung-clustered.h5ad',
    'astir': output_path + "cla/astir_assignments/lin_cycif_lung_astir_assignments.csv",
    'acdc': expand(output_path + "cla/lin_acdc_assignments/lin_cycif_lung_acdc_assignments_{method}.csv",
    method=['absent', 'no-consider']),
    'cla_df': output_path + "figures/cla/cla_df.tsv",
    'cytof_lda_fig': output_path + "figures/cla/cytoflda.pdf"
}


rule fix_basel_cluster_annotations:
    input:
        clusters = "data-raw/jackson-clusters/{cohort}_metaclusters.csv",
        annotation = "data-raw/jackson-clusters/Metacluster_annotations.csv",
    output:
        output_path + "cla/{cohort}_annotation_fixed.csv",
    script:
        "fix_basel_cluster_annotations.R"

rule fix_wagner_cluster_annotations:
    input:
        rds=expand(output_path + "wagner_processed/{sample}.rds",sample=wagner_samples),
        cluster_map=config['cluster_mapping']['wagner']
    output:
        output_path + "cla/wagner_clusters_by_cell_id.csv"
    script:
        "fix_wagner_cluster_annotation.R"



# rule jackson_to_ad:
#     input:
#         lambda wildcards: expand(output_path + wildcards.cohort + "_processed/{core}.csv", core=get_core_list(wildcards.cohort)),
#     output:
#         output_path + "anndata/{cohort}.h5ad",
#     shell:
#         "python pipeline/cla/dir-of-csvs-to-scanpy.py "
#         "{output_path}/{wildcards.cohort}_processed "
#         "{output} "


rule run_cytofLDA:
    params:
        cytofLDA_path = config['cytofLDA_path'],
    input:
        h5ad=output_path + "anndata/{cohort}.h5ad",
        traintest = lambda wildcards: config['cla'][wildcards.cohort]['train_test'],
        labels = lambda wildcards: config['cla'][wildcards.cohort]['annotators'][wildcards.annotator],
    output:
        output_path + "cla/annotations_cytofLDA_{cohort}_{annotator}.tsv",
    shell:
        "Rscript pipeline/cla/run-cytofLDA.R "
        "--input_h5ad {input.h5ad} "
        "--input_traintest {input.traintest} "
        "--input_labels {input.labels} "
        "--cytofLDA_path {params.cytofLDA_path} "
        "--annotator {wildcards.annotator} "
        "--cohort {wildcards.cohort} "
        "--method cytofLDA "
        "--output_assignments {output} "

rule run_ACDC:
    input:
        h5ad=output_path + "anndata/{cohort}.h5ad",
        yaml=lambda wildcards: config[wildcards.cohort]['marker_file'],
    output:
        output_path + "cla/annotations_acdc-{method}_{cohort}.tsv"
    shell:
        "python pipeline/cla/run-acdc.py "
        "--input_h5ad {input.h5ad} "
        "--input_yaml {input.yaml} "
        "--output_assignments {output} "
        "--method {wildcards.method} "
        "--cohort {wildcards.cohort} "

rule copy_lin_assignments:
    input:
        output_path + "astir_assignments/lin_cycif_astir_assignments.csv",
    output:
        output_path + "astir_assignments/lin-cycif_astir_assignments.csv",
    shell:
        "cp {input} {output} "

rule graph_annotation_accuracy_lin_lung:
    params:
        cohort='lin-cycif',
        cytofLDA_path=output_path + "cla/",
        acdc_path=output_path + "cla/lin_acdc_assignments",
        other_workflow_path=output_path + "results/",
        taproom_path=config['taproom_path'],
    input:
        cla_outputs['annotation_results'],
        traintest = lambda wildcards: config['cla']['lin-cycif']['train_test'],
        annotations=lambda wildcards: config['cla']['lin-cycif']['annotators'].values(),
        astir_assignments=output_path+"cla/astir_assignments/lin_cycif_lung_astir_assignments.csv",
        # astir_assignments="astir_basel_assignments-withneg.csv",
    output:
        plot = output_path + "cla/output_figs_tables/cla_annotation_lin-cycif.png",
        tsv = output_path + "cla/output_figs_tables/cla_annotation_lin-cycif.tsv",
    script:
        "graph-accuracy-vs-annotated.R"

rule graph_annotation_accuracy:
    params:
        cohort='{cohort}',
        cytofLDA_path=output_path + "cla/",
        acdc_path=output_path + "cla/",
        other_workflow_path=output_path + "results/",
        taproom_path=config['taproom_path'],
    input:
        cla_outputs['annotation_results'],
        traintest = lambda wildcards: config['cla'][wildcards.cohort]['train_test'],
        annotations=lambda wildcards: config['cla'][wildcards.cohort]['annotators'].values(),
        astir_assignments=output_path+"astir_assignments/{cohort}_astir_assignments.csv"
        # astir_assignments="astir_basel_assignments-withneg.csv",
    output:
        plot = output_path + "cla/output_figs_tables/cla_annotation_{cohort}.png",
        tsv = output_path + "cla/output_figs_tables/cla_annotation_{cohort}.tsv",
    script:
        "graph-accuracy-vs-annotated.R"


rule graph_cluster_accuracy:
    params:
        cohort='{cohort}',
        cytofLDA_path=output_path + "cla/",
        acdc_path=output_path + "cla/",
        other_workflow_path=output_path + "results/",
        taproom_path=config['taproom_path'],
    input:
        coarse_fine_mapping=lambda wildcards: config['coarse_fine_mapping'][wildcards.cohort],
        cluster_mapping=lambda wildcards: config['cluster_mapping'][wildcards.cohort],
        jackson_clustering=output_path + "cla/{cohort}_annotation_fixed.csv",
        tmp=cla_outputs['annotation_results'],
        traintest = lambda wildcards: config['cla'][wildcards.cohort]['train_test'],
        annotations=lambda wildcards: config['cla'][wildcards.cohort]['annotators'].values(),
        astir_assignments=output_path+"astir_assignments/{cohort}_astir_assignments.csv"
    output:
        plot = output_path + "cla/output_figs_tables/cla_cluster_{cohort}.png",
        tsv = output_path + "cla/output_figs_tables/cla_cluster_{cohort}.tsv",
    script:
        "graph-accuracy-vs-clusters.R"

rule graph_cluster_accuracy_wagner:
    params:
        cohort='wagner',
        cytofLDA_path=output_path + "cla/",
        acdc_path=output_path + "cla/",
        other_workflow_path=output_path + "results/",
        taproom_path=config['taproom_path'],
    input:
        coarse_fine_mapping=config['coarse_fine_mapping']['wagner'],
        clustering=output_path + "cla/wagner_clusters_by_cell_id.csv",
        tmp=cla_outputs['annotation_results'],
        traintest = config['cla']['wagner']['train_test'],
        annotations= config['cla']['wagner']['annotators'].values(),
        astir_assignments=output_path+"astir_assignments/wagner_astir_assignments.csv"
    output:
        plot = output_path + "cla/output_figs_tables/wagner/cla_cluster_wagner.png",
        tsv = output_path + "cla/output_figs_tables/wagner/cla_cluster_wagner.tsv",
    script:
        "graph-accuracy-vs-wagner.R"

rule graph_cluster_accuracy_lin:
    params:
        cohort='lin-cycif',
        cytofLDA_path=output_path + "cla/",
        acdc_path=output_path + "cla/lin_acdc_assignments/",
        other_workflow_path=output_path + "results/",
        taproom_path=config['taproom_path'],
    input:
        clustering=config['cluster_mapping']['lin-cycif'],
        tmp=cla_outputs['annotation_results'],
        traintest = config['cla']['lin-cycif']['train_test'],
        annotations= config['cla']['lin-cycif']['annotators'].values(),
        astir_assignments=output_path+"cla/astir_assignments/lin_cycif_lung_astir_assignments.csv",
        acdc_assignments=cla_outputs['acdc'],
    output:
        plot = output_path + "cla/output_figs_tables/lin-cycif/cla_cluster_lin-cycif.png",
        tsv = output_path + "cla/output_figs_tables/lin-cycif/cla_cluster_lin-cycif.tsv",
    script:
        "graph-accuracy-vs-clusters-lin.R"

# rule subset_lin:
#     input:
#         output_path + 'anndata/lin-cycif.h5ad'
#     output:
#         output_path + 'cla/lin-cycif-subset.h5ad'
#     shell:
#         'python pipeline/cla/subset-lin.py '
#         '--input_h5ad {input} '
#         '--output_h5ad {output} '

rule cla_figure:
    params:
        output_fig_dir=output_path + "cla/output_figs_tables",
    input:
        cla_outputs['annotation_results'],
        cla_outputs['wag_clus'],
        cla_outputs['lin_cluster_results'],
        cla_outputs['plots'],
    output:
        annotation=cla_outputs['annotation_fig'],
        clustering=cla_outputs['clustering_fig'],
        phenograph=cla_outputs['phenograph_fig'],
        cla_df=cla_outputs['cla_df'],
    script:
        'cla-fig.R'

# rule run_ACDC_lin:
#     input:
#         h5ad=output_path + 'cla/lin-cycif-subset.h5ad',
#         yaml=ancient(lambda wildcards: config[wildcards.cohort]['marker_file']),
#     output:
#         output_path + "cla/annotations_acdc-{method}_{cohort}.tsv"
#     shell:
#         "python pipeline/cla/run-acdc.py "
#         "--input_h5ad {input.h5ad} "
#         "--input_yaml {input.yaml} "
#         "--output_assignments {output} "
#         "--method {wildcards.method} "
#         "--cohort {wildcards.cohort} "

rule run_astir_lung:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        anndata=output_path + 'anndata/lin_cycif-Lung.h5ad',
        markers=lambda wildcards: config['lin_cycif']['marker_file'],
    output:
        csv=output_path + "cla/astir_assignments/lin_cycif_lung_astir_assignments.csv",
    run:
         # fit
        from astir.data import from_anndata_yaml
        from datetime import datetime
        import scanpy as sc

        adata = sc.read_h5ad(input.anndata)

        print(adata.X[0:5,0:3])

        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_anndata_yaml(input.anndata, input.markers)
        print(f"{datetime.now()}\t Fitting model")

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)

        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate),
        n_init_epochs=int(params.n_initial_epochs))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.get_celltype_probabilities().to_csv(output.csv)

rule run_acdc_lung:
    input:
        h5ad=output_path + 'anndata/lin_cycif-Lung.h5ad',
        yaml=lambda wildcards: config['lin_cycif']['marker_file'],
    output:
        output_path + "cla/lin_acdc_assignments/lin_cycif_lung_acdc_assignments_{method}.csv",
    shell:
        "python pipeline/cla/run-acdc.py "
        "--input_h5ad {input.h5ad} "
        "--input_yaml {input.yaml} "
        "--output_assignments {output} "
        "--method {wildcards.method} "
        "--cohort lin_lung "

rule cluster_lin_lung:
    input:
        output_path + 'anndata/lin_cycif-Lung.h5ad'
    output:
        output_path + 'cla/lin-lung-clustered.h5ad'
    shell:
        'python pipeline/cla/assign-lin-lung-clusters.py '
        '--input_h5ad {input} '
        '--output_h5ad {output} '

rule cytoflda_fig:
    input:
        cla_outputs['cla_df'],
    output:
        cla_outputs['cytof_lda_fig'],
    script:
        'cytoflda-fig.R'
