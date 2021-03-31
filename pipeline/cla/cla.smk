
# ------------------ Cell level accuracy (CLA) quantifaction ----------------- #

cla_outputs = {
    'jackson_annotations_fixed': output_path + "cla/jackson_annotation_fixed.csv"
}


rule fix_basel_cluster_annotations:
    input:
        clusters = "data-raw/jackson-clusters/Basel_metaclusters.csv",
        annotation = "data-raw/jackson-clusters/Metacluster_annotations.csv",
    output:
        cla_outputs['jackson_annotations_fixed'],
    script:
        "fix_basel_cluster_annotations.R"




