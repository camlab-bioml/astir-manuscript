##### ALTERNATIVE CLUSTERING FOR ALTERNATIVE MASKS #####
alt_cores = ['Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35']
users = ['Catena', 'Jackson', 'Schulz']
output_dir_results_altMasks = output_path + "results/alternative_masks/"
output_dir_reports_altMasks = output_path + "reports/alternative_masks/"
output_dir_fig_altMasks = output_path + "figures/alternative_masks/"
alt_masks_shared_output_dir = "/home/ltri/campbell/share/projects/imc-2020/output/" + config['version'] + "/results/alternative_masks/"

markers_options = ['specified_markers', 'all_markers']

alt_methods_cell_assign = []
alt_methods_cell_assign_tmp = [expand(output_dir_reports_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.html", 
core = alt_cores, user = users, method = [m], clusters = conditions[m], markers = markers_options) for m in conditions.keys()]
for element in alt_methods_cell_assign_tmp:
	alt_methods_cell_assign.extend(element)

alt_methods_cell_assign_csv = []
alt_methods_cell_assign_csv_tmp = [expand(alt_masks_shared_output_dir + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.csv", 
core = alt_cores, user = users, method = [m], clusters = conditions[m], markers = markers_options) for m in conditions.keys()]
for element in alt_methods_cell_assign_csv_tmp:
	alt_methods_cell_assign_csv.extend(element)

alt_masks = {
    'sce': expand(output_path + "sces/schapiro_alt-{core}-{user}-sce.rds", core = alt_cores, user = users),
    'Phenograph_analysis_report': expand(output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-specified_markers.html", core = alt_cores, user = users),
    'Phenograph_clusters_output': expand(output_dir_results_altMasks + "Phenograph_clusters_{core}-{user}_specified_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = Phenograph_Cluster_sizes),
    'Phenograph_analysis_report_all': expand(output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-all_markers.html", core = alt_cores, user = users),
    'Phenograph_clusters_output_all': expand(output_dir_results_altMasks + "Phenograph_clusters_{core}-{user}_all_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = Phenograph_Cluster_sizes),
    'FlowSOM_analysis_report': expand(output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-specified_markers.html", core = alt_cores, user = users),
    'FlowSOM_clusters_output': expand(output_dir_results_altMasks + "FlowSOM_clusters_{core}-{user}_specified_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_analysis_report_all': expand(output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-all_markers.html", core = alt_cores, user = users),
    'FlowSOM_clusters_output_all': expand(output_dir_results_altMasks + "FlowSOM_clusters_{core}-{user}_all_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = FlowSOM_Cluster_sizes),
    'ClusterX_analysis_report': expand(output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-specified_markers.html", core = alt_cores, user = users),
    'ClusterX_clusters_output': expand(output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_specified_markers_default.csv", core = alt_cores, user = users),
    'ClusterX_analysis_report_all': expand(output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-all_markers.html", core = alt_cores, user = users),
    'ClusterX_clusters_output_all': expand(output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_all_markers_default.csv", core = alt_cores, user = users),
    
    'alt-mask-cell-type-identification-html': alt_methods_cell_assign,
    'alt-mask-cell-type-identification-csv': alt_methods_cell_assign_csv,
    
    # 'Other_approaches_report': expand(output_dir_reports_altMasks + "Approaches-comparison-{user}.html", user = users),
    # 'Other_approaches_heatmaps': expand(output_dir_results_altMasks + "Assessment-individual-heatmap-{user}.csv", user = users),
    # 'Other_approaches_summary_heatmap': expand(output_dir_results_altMasks + "Assessment-heatmap-{user}.csv", user = users),

    # 'final_benchmark_altMasks': output_dir_reports_altMasks + "Final-approaches-comparison.html",
    # 'final_benchmark_heatmap_altMasks': output_dir_fig_altMasks + "Final_benchmarking_heatmap.pdf"
    
}

rule create_sces_altMasks: 
    params:
        schapiro = output_path + "schapiro_processed_alt_mask/{core}_{user}.rds"
    output:
        schapiro = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"
    shell:
        "Rscript pipeline/rds_to_sce.R {params.schapiro} {output.schapiro}"


rule phenograph_analysis_altMasks:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        markers = lambda wildcards: config['schapiro']['marker_file'],
    params:
        res_dir = output_dir_results_altMasks
    output:
        html = output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-specified_markers.html",
        csvs = expand(output_dir_results_altMasks + "Phenograph_clusters_{{core}}-{{user}}_specified_markers_k{cluster}.csv", cluster = Phenograph_Cluster_sizes),
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "', "
        "params = list(cells = '{input.cells}', create_csv = TRUE, cohort = '{wildcards.core}-{wildcards.user}', markers = 'specified_markers', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule phenograph_analysis_all_markers_altMasks:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"
    params:
        res_dir = output_dir_results_altMasks
    output:
        html = output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-all_markers.html",
        csvs = expand(output_dir_results_altMasks + "Phenograph_clusters_{{core}}-{{user}}_all_markers_k{cluster}.csv", cluster = Phenograph_Cluster_sizes)
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', create_csv = TRUE, cohort = '{wildcards.core}-{wildcards.user}', markers = 'all_markers', output_results = '{params.res_dir}'))\" "


rule FlowSOM_analysis_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        markers = lambda wildcards: config['schapiro']['marker_file']
    params:
        res_dir = output_dir_results_altMasks    
    output:
        html = output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-specified_markers.html",
        csvs = expand(output_dir_results_altMasks + "FlowSOM_clusters_{{core}}-{{user}}_specified_markers_k{cluster}.csv", cluster=FlowSOM_Cluster_sizes),
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', "
        "create_csv = TRUE, cohort = '{wildcards.core}-{wildcards.user}',  markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "



rule FlowSOM_analysis_all_markers_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"    
    params:
        res_dir = output_dir_results_altMasks
    output:
        html = output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-all_markers.html",
        csvs = expand(output_dir_results_altMasks + "FlowSOM_clusters_{{core}}-{{user}}_all_markers_k{cluster}.csv", cluster=FlowSOM_Cluster_sizes)
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', create_csv = TRUE, cohort = '{wildcards.core}-{wildcards.user}', output_results = '{params.res_dir}'))\" "

    
rule ClusterX_analysis_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        markers = lambda wildcards: config['schapiro']['marker_file'],
    params:
        res_dir = output_dir_results_altMasks
    output:
        html = output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-specified_markers.html",
        csvs = output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_specified_markers_default.csv",
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', "
        "create_csv = TRUE, cohort = '{wildcards.core}-{wildcards.user}', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule ClusterX_analysis_all_markers_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
    params:
        res_dir = output_dir_results_altMasks
    output:
        html = output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-all_markers.html",
        csvs = output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_all_markers_default.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', create_csv = TRUE, cohort = '{wildcards.core}-{wildcards.user}', output_results = '{params.res_dir}'))\" "



# Create list of csvs
# csvs_altMasks = []
# csvs_tmp = [expand(output_dir_results_altMasks + "{method}_clusters-{core}-{user}-{markers}_{clusters}.csv", method = [m], markers = markers_spec, user = users, core = alt_cores, clusters = conditions[m]) for m in conditions.keys()]
# for element in csvs_tmp:
# 	csvs_altMasks.extend(element)

# print(csvs_altMasks)

# csv_dict_altMasks = {user: [f for f in csvs_altMasks if user in f] for user in users}



### Compare all the above approaches
rule get_altMasks_cell_types:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        csvs = output_dir_results_altMasks + "{method}_clusters_{core}-{user}_{markers}_{clusters}.csv",
        markers_yml = config['schapiro']['marker_file']
    output:
        html = output_dir_reports_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.html",
        csv = alt_masks_shared_output_dir + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/alt-masks-cluster-identification.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "', "
        "params = list(cells = '{input.cells}', csvs = '{input.csvs}', markers_list = '{input.markers_yml}', core = '{wildcards.core}', user = '{wildcards.user}', "
        "method = '{wildcards.method}', markers = '{wildcards.markers}', cluster = '{wildcards.clusters}', output_dir = '" + alt_masks_shared_output_dir + "' ))\" "


# rule final_approach_comparison_altMasks:
#     input:
#         Jackson = output_dir_results_altMasks + "Assessment-heatmap-Jackson.csv",
#         Catena = output_dir_results_altMasks + "Assessment-heatmap-Catena.csv",
#         Schulz = output_dir_results_altMasks + "Assessment-heatmap-Schulz.csv",
        
#         Jackson_indiv = output_dir_results_altMasks + "Assessment-individual-heatmap-Jackson.csv",
#         Catena_indiv = output_dir_results_altMasks + "Assessment-individual-heatmap-Catena.csv",
#         Schulz_indiv = output_dir_results_altMasks + "Assessment-individual-heatmap-Schulz.csv",

#     output:
#         html = output_dir_reports_altMasks + "Final-approaches-comparison.html",
#         heatmap = output_dir_fig_altMasks + "Final_benchmarking_heatmap.pdf"

#     shell:
#         "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
#         "rmarkdown::render('pipeline/clustering-comparison/Final-approaches-comparison_altMasks.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
#         "params = list(Jackson = '{input.Jackson}', Catena = '{input.Catena}', Schulz = '{input.Schulz}', Jackson_indiv = '{input.Jackson_indiv}', "
#         "Catena_indiv = '{input.Catena_indiv}', Schulz_indiv = '{input.Schulz_indiv}', output_dir = '" + output_dir_fig_altMasks + "'))\" "