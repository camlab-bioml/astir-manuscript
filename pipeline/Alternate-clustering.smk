cohort_list = ['basel', 'zurich1', 'wagner', 'schapiro']
reduced_cohort_list = ['basel', 'zurich1', 'wagner']
markers_spec = ['specified_markers', 'all_markers']

output_dir_reports = output_path + "reports/"
output_dir_results = output_path + "results/"
output_dir_fig = output_path + "figures/"

# Clustering parameters 
FlowSOM_Cluster_sizes = ['4', '7', '8', '12', '15', '20']
Phenograph_Cluster_sizes = ['20', '30', '40', '50']

conditions = { 
    'FlowSOM': ['k4', 'k7', 'k8', 'k12', 'k15', 'k20'],
    'Phenograph': ['k20', 'k30', 'k40', 'k50'],
    'ClusterX': ['default']
}

alluvialPDF = []
alluvial = [expand(output_dir_fig + 'Alluv_cohort_{cohort}_method_{method}_markers_{markers}_clusters_{clusters}.pdf', clusters=conditions[m], method = [m], cohort = cohort_list, markers = markers_spec) for m in conditions.keys()]
for element in alluvial:
    alluvialPDF.extend(element)

expressionPDF = []
expression = [expand(output_dir_fig + 'ExpressionHeatmap_cohort_{cohort}_method_{method}_markers_specified_markers_clusters_{clusters}.pdf', clusters=conditions[m], method = [m], cohort = cohort_list) for m in conditions.keys()]
for element in expression:
    expressionPDF.extend(element)

statesBoxPDF = []
statesBox = [expand(output_dir_fig + 'StatesBoxplot_cohort_{cohort}_method_{method}_clusters_{clusters}.png', clusters=conditions[m], method = [m], cohort = cohort_list) for m in conditions.keys()]
for element in statesBox:
    statesBoxPDF.extend(element)

# Generate list of files
basel_metadata = pd.read_csv(os.path.join(config["basel"]['base_dir'], config["basel"]['metadata_file']))
basel_cores = list(basel_metadata.core)

schapiro_samples = pd.read_csv(config["schapiro"]['sample_txt'], header=None)
schapiro_samples = list(schapiro_samples[schapiro_samples.columns[0]])

wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])

zurich1_metadata = pd.read_csv(os.path.join(config['zurich1']['base_dir'], config['zurich1']['metadata_file']))
zurich1_cores = list(zurich1_metadata.core)

# Create cohort metadata dictionary
cohort_metadata = {
    'basel': basel_cores,
    'zurich1': zurich1_cores,
    'wagner': wagner_samples,
    'schapiro': schapiro_samples
}

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

    # Analysis
    'Alluvial': alluvialPDF,
    'ExpressionHeatmaps': expressionPDF,
    'statesBoxPlots': statesBoxPDF,

    'Other_approaches_report': expand(output_dir_reports + "Approaches-comparison-{cohort}.html", cohort = cohort_list),
    'Other_approaches_heatmaps': expand(output_dir_results + "Assessment-individual-heatmap-{cohort}.csv", cohort = cohort_list),
    'Other_approaches_summary_heatmap': expand(output_dir_results + "Assessment-heatmap-{cohort}.csv", cohort = cohort_list),

    'final_benchmark': output_dir_reports + "Final-approaches-comparison.html",
    'GSVA': expand(output_dir_results + "Other_approaches_GSVA_{cohort}.csv", cohort = reduced_cohort_list)
}

rule create_sces: 
    params:
        basel = ','.join(expand(output_path + "basel_processed/{core}.rds", core = basel_cores)),
        zurich1 = ','.join(expand(output_path + "zurich1_processed/{core}.rds", core = zurich1_cores)),
        wagner = ','.join(expand(output_path + "wagner_processed/{core}.rds", core = wagner_samples)),
        schapiro = ','.join(expand(output_path + "schapiro_processed/{core}.rds", core = schapiro_samples))

    output:
        basel = output_path + "sces/basel_sce.rds",
        wagner = output_path + "sces/wagner_sce.rds",
        zurich1 = output_path + "sces/zurich1_sce.rds",
        schapiro = output_path + "sces/schapiro_sce.rds"

    shell:
        "Rscript pipeline/rds_to_sce.R {params.basel} {output.basel};"
        "Rscript pipeline/rds_to_sce.R {params.zurich1} {output.zurich1};"
        "Rscript pipeline/rds_to_sce.R {params.wagner} {output.wagner};"
        "Rscript pipeline/rds_to_sce.R {params.schapiro} {output.schapiro};"
        

rule phenograph_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
    
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results

    resources:
        mem_mb=2000

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

    resources:
        mem_mb=2000

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

    resources:
        mem_mb=2000
    
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

    resources:
        mem_mb=2000
    
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

    resources:
        mem_mb=3000
    
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

    resources:
        mem_mb=4000
    
    output:
        html = output_dir_reports + "ClusterX-analysis-{cohort}-all_markers.html",
        csvs = output_dir_results + "ClusterX_clusters_{cohort}_all_markers_default.csv"

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', markers = 'all_markers',"
        "create_csv = TRUE, cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "


# Create figures for all of the above methods
rule create_alluvials:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        clusters = output_path + "results/{method}_clusters_{cohort}_{markers}_{clusters}.csv"

    params:
        cohort = "{cohort}",
        method = "{method}",
        clustering_params = "{markers}_clusters_{clusters}"

    output:
        pdf = output_dir_fig + "Alluv_cohort_{cohort}_method_{method}_markers_{markers}_clusters_{clusters}.pdf"
 
    shell:
        "Rscript pipeline/clustering-comparison/alluvials.R {input.cells} {input.cellTypes} {input.cellStates} {input.clusters} "
        "{params.cohort} {params.method} {params.clustering_params} " + output_dir_fig +""


rule create_expression_heatmaps:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        clusters = output_path + "results/{method}_clusters_{cohort}_specified_markers_{clusters}.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],

    params:
        cohort = "{cohort}",
        method = "{method}",
        clustering_params = "markers_specified_markers_clusters_{clusters}"

    output:
        expression_pdf = output_dir_fig + "ExpressionHeatmap_cohort_{cohort}_method_{method}_markers_specified_markers_clusters_{clusters}.pdf"
 
    shell:
        "Rscript pipeline/clustering-comparison/expression-heatmaps.R {input.cells} {input.cellTypes} {input.cellStates} {input.clusters} "
        "{input.markers} {params.cohort} {params.method} {params.clustering_params} " + output_dir_fig +""


rule create_state_boxplots:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        clusters = output_path + "results/{method}_clusters_{cohort}_specified_markers_{clusters}.csv",

    params:
        cohort = "{cohort}",
        method = "{method}",
        clustering_params = "{clusters}"

    output:
        expression_pdf = output_dir_fig + "StatesBoxplot_cohort_{cohort}_method_{method}_clusters_{clusters}.png"
 
    shell:
        "Rscript pipeline/clustering-comparison/boxplot_states.R {input.cells} {input.cellTypes} {input.cellStates} {input.clusters} "
        "{params.cohort} {params.method} {params.clustering_params} " + output_dir_fig +""

# Create list of csvs
csvs = []
csvs_tmp = [expand(output_dir_results + "{method}_clusters_{cohort}_{markers}_{clusters}.csv", method = [m], markers = markers_spec, cohort = cohort_list, clusters = conditions[m]) for m in conditions.keys()]
for element in csvs_tmp:
	csvs.extend(element)

csv_dict = {cohort: [f for f in csvs if cohort in f] for cohort in cohort_list}


### Compare all the above approaches
rule compare_approaches:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
	    csvs = lambda wildcards: csv_dict[wildcards.cohort]

    params:
        cohort = "{cohort}"

    resources:
        mem_mb=5000

    threads: 15

    output:
        html = output_dir_reports + "Approaches-comparison-{cohort}.html",
        heatmap_individual = output_dir_results + "Assessment-individual-heatmap-{cohort}.csv",
        heatmap_summary = output_dir_results + "Assessment-heatmap-{cohort}.csv"

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Approaches-comparison.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cohort = '{params.cohort}', csvs = '{input.csvs}', "
        "markers_list = '{input.markers}', output_dir = '" + output_dir_results + "' ))\" "


rule final_approach_comparison:
    input:
        basel = output_dir_results + "Assessment-heatmap-basel.csv",
        schapiro = output_dir_results + "Assessment-heatmap-schapiro.csv",
        wagner = output_dir_results + "Assessment-heatmap-wagner.csv",
        zurich1 = output_dir_results + "Assessment-heatmap-zurich1.csv",

        basel_indiv = output_dir_results + "Assessment-individual-heatmap-basel.csv",
        schapiro_indiv = output_dir_results + "Assessment-individual-heatmap-schapiro.csv",
        wagner_indiv = output_dir_results + "Assessment-individual-heatmap-wagner.csv",
        zurich1_indiv = output_dir_results + "Assessment-individual-heatmap-zurich1.csv"


    output:
        html = output_dir_reports + "Final-approaches-comparison.html"

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Final-approaches-comparison.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(basel = '{input.basel}', schapiro = '{input.schapiro}', wagner = '{input.wagner}', zurich1 = '{input.zurich1}', "
        "basel_indiv = '{input.basel_indiv}', schapiro_indiv = '{input.schapiro_indiv}', wagner_indiv = '{input.wagner_indiv}', zurich1_indiv = '{input.zurich1_indiv}'))\" "

gsva_csvs = []
gsva_csvs_tmp = [expand(output_dir_results + "{method}_clusters_{cohort}_{markers}_{clusters}.csv", method = [m], markers = markers_spec, cohort = reduced_cohort_list, clusters = conditions[m]) for m in conditions.keys()]
for element in gsva_csvs_tmp:
	gsva_csvs.extend(element)

gsva_csv_dict = {cohort: [f for f in gsva_csvs if cohort in f] for cohort in reduced_cohort_list}      

rule GSVA:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
        csvs = lambda wildcards: gsva_csv_dict[wildcards.cohort]

    output:
        output_dir_results + "Other_approaches_GSVA_{cohort}.csv"

    shell:
        "Rscript pipeline/clustering-comparison/GSVA.R {input.cells} {input.markers} "
        "{input.csvs} {wildcards.cohort} " + output_dir_results

# Still need to define clusters here

# rule GSVA_alluvials:
#     input:
#         gsva_assignments = output_dir_results + "Other_approaches_GSVA_{cohort}.csv",
#         method = tmp_methods

#     output:
#         output_dir_fig + "robustness/other_approaches/Other_methods_robustness_{wildcards.method}_{cohort}.pdf"

#     shell:
#         "Rscript pipeline/clustering-comparison/Other_methods_robustness_alluvials.R {input.gsva_assignments} "
#         "{input.method} {wildcards.cohort} " + output_dir_fig + "robustness/other_approaches/"
