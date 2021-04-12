
# ------------------ Cell level accuracy (CLA) quantifaction ----------------- #

def get_core_list(cohort):
    if cohort == "basel":
        return basel_cores
    elif cohort == "zurich1":
        return zurich1_cores
    else:
        return None

cla_methods = ['cytofLDA', 'acdc-absent', 'acdc-no-consider']
cla_cohorts = config['cla'].keys()

annotation_results = []


for method in cla_methods:
    for cohort in cla_cohorts:
        for annotator in config['cla'][cohort]['annotators'].keys():
            annotation_results.append(
                output_path + f"cla/annotations_{method}_{cohort}_{annotator}.tsv"
            )

cla_outputs = {
    'jackson_annotations_fixed': [output_path + "cla/Basel_annotation_fixed.csv", output_path + "cla/Zurich_annotation_fixed.csv"],
    'anndata': expand(output_path + "anndata/{cohort}.h5ad", cohort=['basel','zurich1']),
    'annotation_results': annotation_results
}


rule fix_basel_cluster_annotations:
    input:
        clusters = "data-raw/jackson-clusters/{cohort}_metaclusters.csv",
        annotation = "data-raw/jackson-clusters/Metacluster_annotations.csv",
    output:
        output_path + "cla/{cohort}_annotation_fixed.csv",
    script:
        "fix_basel_cluster_annotations.R"


rule jackson_to_ad:
    input:
        lambda wildcards: expand(output_path + wildcards.cohort + "_processed/{core}.csv", core=get_core_list(wildcards.cohort)),
    output:
        output_path + "anndata/{cohort}.h5ad",
    shell:
        "python pipeline/cla/dir-of-csvs-to-scanpy.py "
        "{output_path}/{wildcards.cohort}_processed "
        "{output} "


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
        yaml=ancient(lambda wildcards: config[wildcards.cohort]['marker_file']),
    output:
        output_path + "cla/annotations_acdc-{method}_{cohort}_{annotator}.tsv"
    shell:
        "python pipeline/cla/run-acdc.py "
        "--input_h5ad {input.h5ad} "
        "--input_yaml {input.yaml} "
        "--output_assignments {output} "
        "--method {wildcards.method} "
        "--cohort {wildcards.cohort} "


rule graph_annotation_accuracy:
    params:
        cohort='{cohort}',
        cytofLDA_path=output_path + "cla/",
        other_workflow_path=output_path + "results/other-methods-cell-type-assignments/"
    input:
        cla_outputs['annotation_results'],
        traintest = lambda wildcards: config['cla'][wildcards.cohort]['train_test'],
        annotations=lambda wildcards: config['cla'][wildcards.cohort]['annotators'].values(),
        astir_assignment=output_path+"astir_assignments/{cohort}_astir_assignments.csv"
    output:
        plot = output_path + "cla/output_figs_tables/cla_annotation_{cohort}.png",
        tsv = output_path + "cla/output_figs_tables/cla_annotation_{cohort}.tsv",
    run:
        "pipeline/cla/graph-accuracy-vs-annotated.R"

