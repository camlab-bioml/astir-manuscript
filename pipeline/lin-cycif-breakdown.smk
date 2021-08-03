

lin_stratified_reports = output_dir_reports + "lin_stratified/"
FlowSOM_Cluster_sizes = ['4', '7', '8', '12', '15', '20']
Phenograph_Cluster_sizes = ['20', '30', '40', '50', '100', '500']

lin_GSVA_cell_type_assignments = []
lin_GSVA_assignments = [expand(output_dir_results + "lin_stratified/GSVA-assignment-{method}-Lin-CyCif-{cancer}-{markers}-{methodClusters}.csv", 
method = [m], cancer = cancer_types, markers = markers_spec, methodClusters = conditions[m]) for m in conditions.keys()]
for element in lin_GSVA_assignments:
    lin_GSVA_cell_type_assignments.extend(element)

acdc_lin_breakdown_results = []
for cancer in cancer_types:
    for options in [ 'absent', 'no-consider']:
        acdc_lin_breakdown_results.append(
            output_dir_results + f"lin_stratified/ACDC_clusters_Lin-CyCif_{cancer}_{options}.csv"
        )

lin_cycif_breakdown = {
    'sces': expand(output_path + "sces/Lin-cycif-{cancer}.rds", cancer = cancer_types),
    'Phenograph_analysis_report': expand(output_dir_reports + "lin_stratified/Rphenograph-analysis-Lin-CyCif-{cancer}-specified_markers_k{cluster}.html", cancer = cancer_types, cluster = Phenograph_Cluster_sizes),
    'Phenograph_clusters_output': expand(output_dir_results + "lin_stratified/Phenograph_clusters_Lin-CyCif_{cancer}_specified_markers_k{cluster}.csv", cancer = cancer_types, cluster = Phenograph_Cluster_sizes),
    'Phenograph_analysis_report_all': expand(output_dir_reports + "lin_stratified/Rphenograph-analysis-Lin-CyCif-{cancer}-all_markers_k{cluster}.html", cancer = cancer_types, cluster = Phenograph_Cluster_sizes),
    'Phenograph_clusters_output_all': expand(output_dir_results + "lin_stratified/Phenograph_clusters_Lin-CyCif_{cancer}_all_markers_k{cluster}.csv", cancer = cancer_types, cluster = Phenograph_Cluster_sizes),
    'FlowSOM_analysis_report': expand(output_dir_reports + "lin_stratified/FlowSOM-analysis-Lin-CyCif-{cancer}-specified_markers_k{cluster}.html", cancer = cancer_types, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_clusters_output': expand(output_dir_results + "lin_stratified/FlowSOM_clusters_Lin-CyCif_{cancer}_specified_markers_k{cluster}.csv", cancer = cancer_types, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_analysis_report_all': expand(output_dir_reports + "lin_stratified/FlowSOM-analysis-Lin-CyCif-{cancer}-all_markers_k{cluster}.html", cancer = cancer_types, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_clusters_output_all': expand(output_dir_results + "lin_stratified/FlowSOM_clusters_Lin-CyCif_{cancer}_all_markers_k{cluster}.csv", cancer = cancer_types, cluster = FlowSOM_Cluster_sizes),
    'ClusterX_analysis_report': expand(output_dir_reports + "lin_stratified/ClusterX-analysis-Lin-CyCif-{cancer}-specified_markers.html", cancer = cancer_types),
    'ClusterX_clusters_output': expand(output_dir_results + "lin_stratified/ClusterX_clusters_Lin-CyCif_{cancer}_specified_markers_default.csv", cancer = cancer_types),
    'ClusterX_analysis_report_all': expand(output_dir_reports + "lin_stratified/ClusterX-analysis-Lin-CyCif-{cancer}-all_markers.html", cancer = cancer_types),
    'ClusterX_clusters_output_all': expand(output_dir_results + "lin_stratified/ClusterX_clusters_Lin-CyCif_{cancer}_all_markers_default.csv", cancer = cancer_types),
    'acdc_cluster_output': acdc_lin_breakdown_results,
    
    'lin_GSVA_assignments': lin_GSVA_cell_type_assignments,

    'heatmap': output_path + 'figures/Final-heatmap-lin-tissue-types.pdf'
}

rule create_lin_cancerSubtype_sces:
    input:
        rds = lambda wildcards: expand(output_path + "lin-cycif_processed/{core}.rds", core = lin_metadata[lin_metadata["Anatomic site"] == wildcards.cancer]["core"].to_numpy())
    output:
        output_path + "sces/Lin-cycif-{cancer}.rds"
    container: "astir-manuscript.sif"
    shell:
        "Rscript pipeline/rds_to_sce_lin_cycif.R "
        "--rds_files '{input.rds}' "
        "--output_file '{output}' "


rule phenograph_analysis_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds",
        markers = config['lin_cycif']['marker_file']
    params:
        cohort = "Lin-CyCif",
        res_dir = output_dir_results + "lin_stratified/"
    resources: mem_mb=5000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "lin_stratified/Rphenograph-analysis-Lin-CyCif-{cancer}-specified_markers_k{cluster}.html",
        csvs = output_dir_results + "lin_stratified/Phenograph_clusters_Lin-CyCif_{cancer}_specified_markers_k{cluster}.csv"
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph-lin-breakdown.Rmd', output_file = '{output.html}', output_dir = '" + lin_stratified_reports + "', "
        "params = list(cells = '{input.cells}', cohort = '{params.cohort}', markers = 'specified_markers', markers_list = '{input.markers}', "
        "output_results = '{params.res_dir}', cluster_options = '{wildcards.cluster}', cancer = '{wildcards.cancer}'))\" "


rule phenograph_analysis_all_markers_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds"  
    params:
        cohort = "Lin-CyCif",
        res_dir = output_dir_results + "lin_stratified/"
    resources: mem_mb=5000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "lin_stratified/Rphenograph-analysis-Lin-CyCif-{cancer}-all_markers_k{cluster}.html",
        csvs = output_dir_results + "lin_stratified/Phenograph_clusters_Lin-CyCif_{cancer}_all_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph-lin-breakdown.Rmd', output_file = '{output.html}', output_dir = '" + lin_stratified_reports + "',"
        "params = list(cells = '{input.cells}', cohort = '{params.cohort}', markers = 'all_markers', output_results = '{params.res_dir}', "
        "cluster_options = '{wildcards.cluster}', cancer = '{wildcards.cancer}'))\" "


rule FlowSOM_analysis_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds",
        markers = config['lin_cycif']['marker_file']
    params:
        cohort = "Lin-CyCif",
        res_dir = output_dir_results + "lin_stratified/"
    resources: mem_mb=2000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "lin_stratified/FlowSOM-analysis-Lin-CyCif-{cancer}-specified_markers_k{cluster}.html",
        csvs = output_dir_results + "lin_stratified/FlowSOM_clusters_Lin-CyCif_{cancer}_specified_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM-lin-breakdown.Rmd', output_file = '{output.html}', output_dir = '" + lin_stratified_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', cohort = '{params.cohort}', markers_list = '{input.markers}', "
        "output_results = '{params.res_dir}', cluster_options = '{wildcards.cluster}', cancer = '{wildcards.cancer}'))\" "



rule FlowSOM_analysis_all_markers_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds"    
    params:
        cohort = "Lin-CyCif",
        res_dir = output_dir_results + "lin_stratified/"
    resources:
        mem_mb=2000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "lin_stratified/FlowSOM-analysis-Lin-CyCif-{cancer}-all_markers_k{cluster}.html",
        csvs = output_dir_results + "lin_stratified/FlowSOM_clusters_Lin-CyCif_{cancer}_all_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM-lin-breakdown.Rmd', output_file = '{output.html}', output_dir = '" + lin_stratified_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', cohort = '{params.cohort}', output_results = '{params.res_dir}', "
        "cluster_options = '{wildcards.cluster}', cancer = '{wildcards.cancer}'))\" "



rule ClusterX_analysis_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds",
        markers = config['lin_cycif']['marker_file'],
    params:
        cohort = "Lin-CyCif",
        res_dir = output_dir_results + "lin_stratified/"
    resources: mem_mb=4000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "lin_stratified/ClusterX-analysis-Lin-CyCif-{cancer}-specified_markers.html",
        csvs = output_dir_results + "lin_stratified/ClusterX_clusters_Lin-CyCif_{cancer}_specified_markers_default.csv",
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX-lin-breakdown.Rmd', output_file = '{output.html}', output_dir = '" + lin_stratified_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', cohort = '{params.cohort}', markers_list = '{input.markers}', "
        "output_results = '{params.res_dir}', cancer = '{wildcards.cancer}'))\" "

rule ClusterX_analysis_all_markers_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds",
    params:
        cohort = "Lin-CyCif",
        res_dir = output_dir_results + "lin_stratified/"
    resources: mem_mb=4000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "lin_stratified/ClusterX-analysis-Lin-CyCif-{cancer}-all_markers.html",
        csvs = output_dir_results + "lin_stratified/ClusterX_clusters_Lin-CyCif_{cancer}_all_markers_default.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX-lin-breakdown.Rmd', output_file = '{output.html}', output_dir = '" + lin_stratified_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', cohort = '{params.cohort}', output_results = '{params.res_dir}', "
        "cancer = '{wildcards.cancer}'))\" "

rule acdc_lin_breakdown:
    input:
        h5ad = output_path + "anndata/lin_cycif-{cancer}.h5ad",
        markers = lambda wildcards: config['lin_cycif']['marker_file']
    output:
        csv = output_dir_results + 'lin_stratified/ACDC_clusters_Lin-CyCif_{cancer}_{options}.csv'
    shell:
        "python pipeline/cla/run-acdc.py "
        "--input_h5ad {input.h5ad} "
        "--input_yaml {input.markers} "
        "--output_assignments {output.csv} "
        "--method {wildcards.options} "
        "--cohort {wildcards.cancer} "
        


rule assign_GSVA_cell_types_lin_breakdown:
    input:
        cells = output_path + "sces/Lin-cycif-{cancer}.rds",
        markers = config['lin_cycif']['marker_file'],
        clustering_result = output_dir_results + "lin_stratified/{method}_clusters_Lin-CyCif_{cancer}_{markers}_{methodClusters}.csv"
    output:
        #html = output_dir_reports + "lin_stratified/GSVA-assignment-{method}-Lin-CyCif-{cancer}-{markers}-{methodClusters}.html",
        csv = output_dir_results + "lin_stratified/GSVA-assignment-{method}-Lin-CyCif-{cancer}-{markers}-{methodClusters}.csv"
    container: "astir-manuscript.sif"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Assign-cell-type-to-GSVA-cluster.Rmd', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', markers = '{input.markers}', clustering_result = '{input.clustering_result}', output_file = '{output.csv}' ))\" "


lin_csvs = []
lin_csvs_tmp = [expand(output_dir_results + "lin_stratified/GSVA-assignment-{method}-Lin-CyCif-{tissue}-{markers}-{clusters}.csv", method = [m], markers = markers_spec, tissue = cancer_types, clusters = conditions[m]) for m in conditions.keys()]
for element in lin_csvs_tmp:
	lin_csvs.extend(element)
lin_csv_dict = {tissue: [f for f in lin_csvs if tissue in f] for tissue in cancer_types}\


rule lin_breakdown_comparison:
    input:
        bladder_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Bladder_astir_assignments.csv',
        breast_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Breast_astir_assignments.csv',
        GI_stomach_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_GI-Stomach_astir_assignments.csv',
        GI_colorectal_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_GI-Colonrectal_astir_assignments.csv',
        kidney_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Kidney_astir_assignments.csv',
        liver_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Liver_astir_assignments.csv',
        lung_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Lung_astir_assignments.csv',
        lymph_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Lymph-node_astir_assignments.csv',
        ovary_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Ovary_astir_assignments.csv',
        pancreas_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Pancreas_astir_assignments.csv',
        prostate_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Prostate_astir_assignments.csv',
        skin_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Skin_astir_assignments.csv',
        uterus_astir_assignments = output_path + 'lin_tissue_types_astir_assignments/lin_cycif_Uterus_astir_assignments.csv',

        bladder_sce = output_path + 'sces/Lin-cycif-Bladder.rds',
        breast_sce = output_path + 'sces/Lin-cycif-Breast.rds',
        GI_stomach_sce = output_path + 'sces/Lin-cycif-GI-Stomach.rds',
        GI_colorectal_sce = output_path + 'sces/Lin-cycif-GI-Colonrectal.rds',
        kidney_sce = output_path + 'sces/Lin-cycif-Kidney.rds',
        liver_sce = output_path + 'sces/Lin-cycif-Liver.rds',
        lung_sce = output_path + 'sces/Lin-cycif-Lung.rds',
        lymph_sce = output_path + 'sces/Lin-cycif-Lymph-node.rds',
        ovary_sce = output_path + 'sces/Lin-cycif-Ovary.rds',
        pancreas_sce = output_path + 'sces/Lin-cycif-Pancreas.rds',
        prostate_sce = output_path + 'sces/Lin-cycif-Prostate.rds',
        skin_sce = output_path + 'sces/Lin-cycif-Skin.rds',
        uterus_sce = output_path + 'sces/Lin-cycif-Uterus.rds',
        
        bladder_files = lin_csv_dict['Bladder'],
        breast_files = lin_csv_dict['Breast'],
        GI_stomach_files = lin_csv_dict['GI-Stomach'],
        GI_colorectal_files = lin_csv_dict['GI-Colonrectal'],
        kidney_files = lin_csv_dict['Kidney'],
        liver_files = lin_csv_dict['Liver'],
        lung_files = lin_csv_dict['Lung'],
        lymph_files = lin_csv_dict['Lymph-node'],
        ovary_files = lin_csv_dict['Ovary'],
        pancreas_files = lin_csv_dict['Pancreas'],
        prostate_files = lin_csv_dict['Prostate'],
        skin_files = lin_csv_dict['Skin'],
        uterus_files = lin_csv_dict['Uterus'],

        bladder_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Bladder_absent.csv",
        bladder_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Bladder_no-consider.csv",
        breast_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Breast_absent.csv",
        breast_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Breast_no-consider.csv",
        GI_stomach_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_GI-Stomach_absent.csv",
        GI_stomach_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_GI-Stomach_no-consider.csv",
        GI_colorectal_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_GI-Colonrectal_absent.csv",
        GI_colorectal_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_GI-Colonrectal_no-consider.csv",
        kidney_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Kidney_absent.csv",
        kidney_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Kidney_no-consider.csv",
        liver_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Liver_absent.csv",
        liver_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Liver_no-consider.csv",
        lung_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Lung_absent.csv",
        lung_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Lung_no-consider.csv",
        lymph_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Lymph-node_absent.csv",
        lymph_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Lymph-node_no-consider.csv",
        ovary_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Ovary_absent.csv",
        ovary_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Ovary_no-consider.csv",
        pancreas_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Pancreas_absent.csv",
        pancreas_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Pancreas_no-consider.csv",
        prostate_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Prostate_absent.csv",
        prostate_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Prostate_no-consider.csv",
        skin_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Skin_absent.csv",
        skin_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Skin_no-consider.csv",
        uterus_absent_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Uterus_absent.csv",
        uterus_no_consider_acdc = output_dir_results + "lin_stratified/ACDC_clusters_Lin-CyCif_Uterus_no-consider.csv",

        markers = config['lin_cycif']['marker_file']
    output:
        pdf = output_path + 'figures/Final-heatmap-lin-tissue-types.pdf'
    container: "astir-manuscript.sif"
    shell:
        "Rscript pipeline/clustering-comparison/Final-approaches-comparison-lin-stratified.R "
        "--bladder_astir_assignments {input.bladder_astir_assignments} "
        "--breast_astir_assignments {input.breast_astir_assignments} "
        "--GI_stomach_astir_assignments {input.GI_stomach_astir_assignments} "
        "--GI_colorectal_astir_assignments {input.GI_colorectal_astir_assignments} "
        "--kidney_astir_assignments {input.kidney_astir_assignments} "
        "--liver_astir_assignments {input.liver_astir_assignments} "
        "--lung_astir_assignments {input.lung_astir_assignments} "
        "--lymph_astir_assignments {input.lymph_astir_assignments} "
        "--ovary_astir_assignments {input.ovary_astir_assignments} "
        "--pancreas_astir_assignments {input.pancreas_astir_assignments} "
        "--prostate_astir_assignments {input.prostate_astir_assignments} "
        "--skin_astir_assignments {input.skin_astir_assignments} "
        "--uterus_astir_assignments {input.uterus_astir_assignments} "

        "--bladder_sce {input.bladder_sce} "
        "--breast_sce {input.breast_sce} "
        "--GI_stomach_sce {input.GI_stomach_sce} "
        "--GI_colorectal_sce {input.GI_colorectal_sce} "
        "--kidney_sce {input.kidney_sce} "
        "--liver_sce {input.liver_sce} "
        "--lung_sce {input.lung_sce} "
        "--lymph_sce {input.lymph_sce} "
        "--ovary_sce {input.ovary_sce} "
        "--pancreas_sce {input.pancreas_sce} "
        "--prostate_sce {input.prostate_sce} "
        "--skin_sce {input.skin_sce} "
        "--uterus_sce {input.uterus_sce} "

        "--bladder_files {input.bladder_files} "
        "--breast_files {input.breast_files} "
        "--GI_stomach_files {input.GI_stomach_files} "
        "--GI_colorectal_files {input.GI_colorectal_files} "
        "--kidney_files {input.kidney_files} "
        "--liver_files {input.liver_files} "
        "--lung_files {input.lung_files} "
        "--lymph_files {input.lymph_files} "
        "--ovary_files {input.ovary_files} "
        "--pancreas_files {input.pancreas_files} "
        "--prostate_files {input.prostate_files} "
        "--skin_files {input.skin_files} "
        "--uterus_files {input.uterus_files} "

        "--bladder_absent_acdc {input.bladder_absent_acdc} "
        "--bladder_no_consider_acdc {input.bladder_no_consider_acdc} "
        "--breast_absent_acdc {input.breast_absent_acdc} "
        "--breast_no_consider_acdc {input.breast_no_consider_acdc} "
        "--GI_stomach_absent_acdc {input.GI_stomach_absent_acdc} "
        "--GI_stomach_no_consider_acdc {input.GI_stomach_no_consider_acdc} "
        "--GI_colorectal_absent_acdc {input.GI_colorectal_absent_acdc} "
        "--GI_colorectal_no_consider_acdc {input.GI_colorectal_no_consider_acdc} "
        "--kidney_absent_acdc {input.kidney_absent_acdc} "
        "--kidney_no_consider_acdc {input.kidney_no_consider_acdc} "
        "--liver_absent_acdc {input.liver_absent_acdc} "
        "--liver_no_consider_acdc {input.liver_no_consider_acdc} "
        "--lung_absent_acdc {input.lung_absent_acdc} "
        "--lung_no_consider_acdc {input.lung_no_consider_acdc} "
        "--lymph_absent_acdc {input.lymph_absent_acdc} "
        "--lymph_no_consider_acdc {input.lymph_no_consider_acdc} "
        "--ovary_absent_acdc {input.ovary_absent_acdc} "
        "--ovary_no_consider_acdc {input.ovary_no_consider_acdc} "
        "--pancreas_absent_acdc {input.pancreas_absent_acdc} "
        "--pancreas_no_consider_acdc {input.pancreas_no_consider_acdc} "
        "--prostate_absent_acdc {input.prostate_absent_acdc} "
        "--prostate_no_consider_acdc {input.prostate_no_consider_acdc} "
        "--skin_absent_acdc {input.skin_absent_acdc} "
        "--skin_no_consider_acdc {input.skin_no_consider_acdc} "
        "--uterus_absent_acdc {input.uterus_absent_acdc} "
        "--uterus_no_consider_acdc {input.uterus_no_consider_acdc} "

        "--lin_markers {input.markers} "
        "--output_heatmap {output.pdf} "