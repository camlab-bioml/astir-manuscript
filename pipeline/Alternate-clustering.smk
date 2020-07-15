cohort_list = ['basel', 'zurich1', 'wagner']#, 'schapiro']
output_dir_reports = output_path + "reports/"
output_dir_results = output_path + "results/"

# Clustering parameters 
FlowSOM_Cluster_sizes = ['8', '12', '15', '20']
Phenograph_Cluster_sizes = ['20', '30', '40', '50']

# Generate list of files
basel_metadata = pd.read_csv(os.path.join(config["basel"]['base_dir'], config["basel"]['metadata_file']))
basel_cores = list(basel_metadata.core)

#schapiro_samples = pd.read_csv(config["schapiro"]['sample_txt'], header=None)
#schapiro_samples = list(schapiro_samples[schapiro_samples.columns[0]])

wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])

zurich1_metadata = pd.read_csv(os.path.join(config['zurich1']['base_dir'], config['zurich1']['metadata_file']))
zurich1_cores = list(zurich1_metadata.core)


alternate_approaches_output = {
    'Phenograph_analysis_report': expand(output_dir_reports + "Rphenograph-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
    'Phenograph_clusters_output': expand(output_dir_results + "Phenograph_clusters_{cohort}_specified_markers_k{cluster}.csv", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
    'Phenograph_analysis_report_all': expand(output_dir_reports + "Rphenograph-analysis-{cohort}-all_markers.html", cohort = cohort_list),
    'Phenograph_clusters_output_all': expand(output_dir_results + "Phenograph_clusters_{cohort}_all_markers_k{cluster}.csv", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
    'FlowSOM_analysis_report': expand(output_dir_reports + "FlowSOM-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
    'FlowSOM_clusters_output': expand(output_dir_results + "FlowSOM_clusters_{cohort}_specified_markers_k{cluster}.csv", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_analysis_report_all': expand(output_dir_reports + "FlowSOM-analysis-{cohort}-all_markers.html", cohort = cohort_list),
    'FlowSOM_clusters_output_all': expand(output_dir_results + "FlowSOM_clusters_{cohort}_all_markers_k{cluster}.csv", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),
    'ClusterX_analysis_report': expand(output_dir_reports + "ClusterX-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
    'ClusterX_clusters_output': expand(output_dir_results + "ClusterX_clusters_{cohort}_specified_markers_default.csv", cohort = cohort_list),
    'ClusterX_analysis_report_all': expand(output_dir_reports + "ClusterX-analysis-{cohort}-all_markers.html", cohort = cohort_list),
    'ClusterX_clusters_output_all': expand(output_dir_results + "ClusterX_clusters_{cohort}_all_markers_default.csv", cohort = cohort_list),   
}

#rule all:
    #input:
        #Rphenograph_analysis_report = expand(output_dir_reports + "Rphenograph-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
        #Rphenograph_clusters_output = expand(output_dir_results + "Phenograph_clusters_{cohort}_specified_markers_k{cluster}.csv", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
        #Rphenograph_analysis_report_all = expand(output_dir_reports + "Rphenograph-analysis-{cohort}-all_markers.html", cohort = cohort_list),
        #Rphenograph_clusters_output_all = expand(output_dir_results + "Phenograph_clusters_{cohort}_all_markers_k{cluster}.csv", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),

        #FlowSOM_analysis_report = expand(output_dir_reports + "FlowSOM-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
        #FlowSOM_clusters_output = expand(output_dir_results + "FlowSOM_clusters_{cohort}_specified_markers_k{cluster}.csv", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),
        #FlowSOM_analysis_report_all = expand(output_dir_reports + "FlowSOM-analysis-{cohort}-all_markers.html", cohort = cohort_list),
        #FlowSOM_clusters_output_all = expand(output_dir_results + "FlowSOM_clusters_{cohort}_all_markers_k{cluster}.csv", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),

        #ClusterX_analysis_report = expand(output_dir_reports + "ClusterX-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
        #ClusterX_clusters_output = expand(output_dir_results + "ClusterX_clusters_{cohort}_specified_markers_default.csv", cohort = cohort_list),
        #ClusterX_analysis_report_all = expand(output_dir_reports + "ClusterX-analysis-{cohort}-all_markers.html", cohort = cohort_list),
        #ClusterX_clusters_output_all = expand(output_dir_results + "ClusterX_clusters_{cohort}_all_markers_default.csv", cohort = cohort_list),

        #Other_approaches_report = expand(output_dir_reports + "Approaches-comparison-{cohort}.html", cohort = cohort_list),
        #Other_approaches_heatmaps = expand(output_dir_results + "Assessment-individual-heatmap-{cohort}.csv", cohort = cohort_list),
        #Other_approaches_summary_heatmap = expand(output_dir_results + "Assessment-heatmap-{cohort}.csv", cohort = cohort_list),

        #Final_approaches_comparison = expand(output_dir_reports + "Final-approaches-comparison.html")

#print(','.join(expand(output_path + "wagner_processed/{core}.rds", core = wagner_samples)))
rule create_sces: 
    params:
        basel = ','.join(expand(output_path + "basel_processed/{core}.rds", core = basel_cores)),
        zurich1 = ','.join(expand(output_path + "zurich1_processed/{core}.rds", core = zurich1_cores)),
        wagner = ','.join(expand(output_path + "wagner_processed/{core}.rds", core = wagner_samples))

    output:
        basel = output_path + "sces/basel_sce.rds",
        wagner = output_path + "sces/wagner_sce.rds",
        zurich1 = output_path + "sces/zurich1_sce.rds",

    shell:
        "Rscript pipeline/rds_to_sce.R {params.basel} {output.basel};"
        "Rscript pipeline/rds_to_sce.R {params.zurich1} {output.zurich1};"
        "Rscript pipeline/rds_to_sce.R {params.wagner} {output.wagner};"


rule phenograph_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
    
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results

    output:
        html = output_dir_reports + "Rphenograph-analysis-{cohort}-specified_markers.html",
        csvs = expand(output_dir_results + "Phenograph_clusters_{{cohort}}_specified_markers_k{cluster}.csv", cluster = Phenograph_Cluster_sizes),
        
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', "
        "create_csv = TRUE, cohort = '{params.cohort}', markers = 'specified_markers', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule phenograph_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
    
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results

    output:
        html = output_dir_reports + "Rphenograph-analysis-{cohort}-all_markers.html",
        csvs = expand(output_dir_results + "Phenograph_clusters_{{cohort}}_all_markers_k{cluster}.csv", cluster = Phenograph_Cluster_sizes)

    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}',"
        "create_csv = TRUE, cohort = '{params.cohort}', markers = 'all_markers', output_results = '{params.res_dir}'))\" "


rule FlowSOM_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']
    
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    
    output:
        html = output_dir_reports + "FlowSOM-analysis-{cohort}-specified_markers.html",
        csvs = expand(output_dir_results + "FlowSOM_clusters_{{cohort}}_specified_markers_k{cluster}.csv", cluster=FlowSOM_Cluster_sizes),
    
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', markers = 'specified_markers',"
        "create_csv = TRUE, cohort = '{params.cohort}',  markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "



rule FlowSOM_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
    
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    
    output:
        html = output_dir_reports + "FlowSOM-analysis-{cohort}-all_markers.html",
        csvs = expand(output_dir_results + "FlowSOM_clusters_{{cohort}}_all_markers_k{cluster}.csv", cluster=FlowSOM_Cluster_sizes)

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', markers = 'all_markers',"
        "create_csv = TRUE, cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "

    
rule ClusterX_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],

    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    
    output:
        html = output_dir_reports + "ClusterX-analysis-{cohort}-specified_markers.html",
        csvs = output_dir_results + "ClusterX_clusters_{cohort}_specified_markers_default.csv",

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', markers = 'specified_markers',"
        "create_csv = TRUE, cohort = '{params.cohort}', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule ClusterX_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",

    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    
    output:
        html = output_dir_reports + "ClusterX-analysis-{cohort}-all_markers.html",
        csvs = output_dir_results + "ClusterX_clusters_{cohort}_all_markers_default.csv"

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', markers = 'all_markers',"
        "create_csv = TRUE, cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "


# ### Compare all the above approaches

# rule compare_approaches:
#     input:
#         cells = output_path + "{cohort}_subset/{cohort}_subset_sce.rds",
#         celltypes = output_path + "{cohort}_subset/{cohort}_subset_assignments_type.csv"

#     params:
#         cohort = "{cohort}"

#     output:
#         html = output_dir_reports + "Approaches-comparison-{cohort}.html",
#         heatmap_individual = output_dir_results + "Assessment-individual-heatmap-{cohort}.csv",
#         heatmap_summary = output_dir_results + "Assessment-heatmap-{cohort}.csv"

#     shell:
#         "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
#         "rmarkdown::render('clustering-comparison/Approaches-comparison.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
#         "params = list(cells = '{input.cells}', celltypes = '{input.celltypes}', cohort = '{params.cohort}' ))\" "


# rule final_approach_comparison:
#     input:
#         basel = output_dir_results + "Assessment-heatmap-basel.csv",
#         wagner = output_dir_results + "Assessment-heatmap-wagner.csv",
#         zurich1 = output_dir_results + "Assessment-heatmap-zurich1.csv"

#     output:
#         html = output_dir_reports + "Final-approaches-comparison.html"

#     shell:
#         "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
#         "rmarkdown::render('clustering-comparison/Final-approaches-comparison.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
#         "params = list(basel = '{input.basel}', wagner = '{input.wagner}', zurich1 = '{input.zurich1}'))\" "
        
