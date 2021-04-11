
# ------------------ Cell level accuracy (CLA) quantifaction ----------------- #

cla_methods = ['cytofLDA']
cla_cohorts = config['cla'].keys()

annotation_results = []

for method in cla_methods:
    for cohort in cla_cohorts:
        for annotator in config['cla'][cohort]['annotators'].keys():
            annotation_results.append(
                output_path + f"cla/{method}/annotations_{method}_{cohort}_{annotator}.tsv"
            )

cla_outputs = {
    'jackson_annotations_fixed': output_path + "cla/jackson_annotation_fixed.csv",
    'basel_ad': output_path + "anndata/basel.h5ad",
    'annotation_results': annotation_results
}


rule fix_basel_cluster_annotations:
    input:
        clusters = "data-raw/jackson-clusters/Basel_metaclusters.csv",
        annotation = "data-raw/jackson-clusters/Metacluster_annotations.csv",
    output:
        cla_outputs['jackson_annotations_fixed'],
    script:
        "fix_basel_cluster_annotations.R"


rule basel_to_ad:
    input:
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
    output:
        cla_outputs['basel_ad'],
    shell:
        "python pipeline/cla/dir-of-csvs-to-scanpy.py "
        "{output_path}/basel_processed "
        "{output} "


rule run_cytofLDA:
    params:
        cytofLDA_path = config['cytofLDA_path'],
    input:
        h5ad=output_path + "anndata/{cohort}.h5ad",
        traintest = lambda wildcards: config['cla'][wildcards.cohort]['train_test'],
        labels = lambda wildcards: config['cla'][wildcards.cohort]['annotators'][wildcards.annotator],
    output:
        output_path + "cla/cytofLDA/annotations_cytofLDA_{cohort}_{annotator}.tsv",
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

