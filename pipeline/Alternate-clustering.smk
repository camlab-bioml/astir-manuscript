cohort_list = ['basel', 'zurich1', 'wagner', 'schapiro', 'lin_cycif']
reduced_cohort_list = ['basel', 'zurich1', 'wagner']
markers_spec = ['specified_markers', 'all_markers']

output_dir_reports = output_path + "reports/"
output_dir_results = output_path + "results/"
output_dir_fig = output_path + "figures/"

# Clustering parameters 
FlowSOM_Cluster_sizes = ['4', '7', '8', '12', '15', '20']
Phenograph_Cluster_sizes = ['10', '20', '30', '40', '50']

conditions = { 
    'FlowSOM': ['k4', 'k7', 'k8', 'k12', 'k15', 'k20'],
    'Phenograph': ['k10', 'k20', 'k30', 'k40', 'k50'],
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
    'schapiro': schapiro_samples,
    'lin_cycif': cycif_samples,
    'keren': keren_cores
}

# Create list containing markers to be removed from alternative approaches
markers_removal = ['None', 'Stromal', 'Stromal-Macrophage', 'Stromal-Macrophage-Endothelial']
markers_gsva = ['specified-markers', 'all-markers']

gsva_robustness_all = []
gsva_robustness_tmp_all = [expand(output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_all-markers_{markers_rem}.csv", 
cohort = reduced_cohort_list, method = [m], clusters = conditions[m], markers_rem = markers_removal) for m in conditions.keys()]
for element in gsva_robustness_tmp_all:
	gsva_robustness_all.extend(element)

gsva_assignments_reports_all = []
gsva_assignments_reports_tmp_all = [expand(output_dir_reports + "GSVA/GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_all-markers_{markers_rem}.html", 
cohort = reduced_cohort_list, method = [m], clusters = conditions[m], markers_rem = markers_removal) for m in conditions.keys()]
for element in gsva_assignments_reports_tmp_all:
	gsva_assignments_reports_all.extend(element)


gsva_robustness_spec = []
gsva_robustness_tmp_spec = [expand(output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_specified-markers_{markers_rem}.csv", 
cohort = reduced_cohort_list, method = [m], clusters = conditions[m], markers_rem = markers_removal) for m in conditions.keys()]
for element in gsva_robustness_tmp_spec:
	gsva_robustness_spec.extend(element)

gsva_assignments_reports_spec = []
gsva_assignments_reports_tmp_spec = [expand(output_dir_reports + "GSVA/GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_specified-markers_{markers_rem}.html", 
cohort = reduced_cohort_list, method = [m], clusters = conditions[m], markers_rem = markers_removal) for m in conditions.keys()]
for element in gsva_assignments_reports_tmp_spec:
	gsva_assignments_reports_spec.extend(element)

#out_dict = {cohort: [f for f in out_perm if cohort in f] for cohort in reduced_cohort_list}

robustness_alluvials = []
robustness_alluvials_tmp = [expand(output_dir_fig + "robustness/GSVA-robustness-alluv_{cohort}_{method}_{clusters}_{markers}.pdf", 
cohort = reduced_cohort_list, method = [m], clusters = conditions[m], markers = markers_gsva) for m in conditions.keys()]
#cohort = 'wagner', method = "FlowSOM", clusters = "k7", markers = "specified-markers")
for element in robustness_alluvials_tmp:
    robustness_alluvials.extend(element)


alternate_approaches_output = {
    'sces': expand(output_path + "sces/{cohort}_sce.rds", cohort = cohort_list),
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

    # # # Analysis
    'Alluvial': alluvialPDF,
    # # 'ExpressionHeatmaps': expressionPDF,
    # # 'statesBoxPlots': statesBoxPDF,

    'Other_approaches_report': expand(output_dir_reports + "Approaches-comparison-{cohort}.html", cohort = cohort_list),
    'Other_approaches_heatmaps': expand(output_dir_results + "Assessment-individual-heatmap-{cohort}.csv", cohort = cohort_list),
    'Other_approaches_summary_heatmap': expand(output_dir_results + "Assessment-heatmap-{cohort}.csv", cohort = cohort_list),

    'final_benchmark': output_dir_reports + "Final-approaches-comparison.html",
    'final_benchmark_heatmap': output_dir_fig + "Final_benchmarking_heatmap.pdf",
    ###'GSVA_csv_all': gsva_robustness_all,
    ###'GSVA_html_all': gsva_assignments_reports_all,
    #'GSVA_csv_spec': gsva_robustness_spec,
    #'GSVA_html_spec': gsva_assignments_reports_spec,
    #'GSVA_html': expand(output_dir_reports + "GSVA_{markers_rem}_{cohort}.html", cohort = reduced_cohort_list, markers_rem = markers_removal),
    #'GSVA_csv': expand(output_dir_results + "Other_approaches_GSVA_{markers_rem}_{cohort}.csv", cohort = reduced_cohort_list, markers_rem = markers_removal),

    #'GSVA_robust': robustness_alluvials,
    #'Astir-robustness': expand(output_dir_fig + "Astir-robustness-{cohort}.pdf", cohort = reduced_cohort_list)
}

rule create_sces: 
    params:
        basel = ','.join(expand(output_path + "basel_processed/{core}.rds", core = basel_cores)),
        zurich1 = ','.join(expand(output_path + "zurich1_processed/{core}.rds", core = zurich1_cores)),
        wagner = ','.join(expand(output_path + "wagner_processed/{core}.rds", core = wagner_samples)),
        schapiro = ','.join(expand(output_path + "schapiro_processed/{core}.rds", core = schapiro_samples)),
        lin = ','.join(expand(output_path + "lin-cycif_processed/{core}.rds", core = cycif_samples)),
        keren = ','.join(expand(output_path + "keren_processed/{core}.rds", core = keren_cores))
    output:
        basel = output_path + "sces/basel_sce.rds",
        wagner = output_path + "sces/wagner_sce.rds",
        zurich1 = output_path + "sces/zurich1_sce.rds",
        schapiro = output_path + "sces/schapiro_sce.rds",
        lin = output_path + "sces/lin_cycif_sce.rds",
        keren = output_path + "sces/keren_sce.rds"
    shell:
        "Rscript pipeline/rds_to_sce.R {params.basel} {output.basel};"
        "Rscript pipeline/rds_to_sce.R {params.zurich1} {output.zurich1};"
        "Rscript pipeline/rds_to_sce.R {params.wagner} {output.wagner};"
        "Rscript pipeline/rds_to_sce.R {params.schapiro} {output.schapiro};"
        "Rscript pipeline/rds_to_sce.R {params.lin} {output.lin};"
        "Rscript pipeline/rds_to_sce.R {params.keren} {output.keren};"
        

rule phenograph_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources:
        mem_mb=5000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "Rphenograph-analysis-{cohort}-specified_markers.html",
        csvs = expand(output_dir_results + "Phenograph_clusters_{{cohort}}_specified_markers_k{cluster}.csv", cluster = Phenograph_Cluster_sizes),  
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', cohort = '{params.cohort}', markers = 'specified_markers', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule phenograph_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds"  
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources:
        mem_mb=5000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "Rphenograph-analysis-{cohort}-all_markers.html",
        csvs = expand(output_dir_results + "Phenograph_clusters_{{cohort}}_all_markers_k{cluster}.csv", cluster = Phenograph_Cluster_sizes)
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', cohort = '{params.cohort}', markers = 'all_markers', output_results = '{params.res_dir}'))\" "


rule FlowSOM_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        #cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
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
        "params = list(cells = '{input.cells}', markers = 'specified_markers', "
        "create_csv = TRUE, cohort = '{params.cohort}',  markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "



rule FlowSOM_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        #cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
    
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
        "params = list(cells = '{input.cells}', markers = 'all_markers', "
        "create_csv = TRUE, cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "

    
rule ClusterX_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        #cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],

    params:
        cohort = "{cohort}",
        res_dir = output_dir_results

    resources:
        mem_mb=4000
    
    output:
        html = output_dir_reports + "ClusterX-analysis-{cohort}-specified_markers.html",
        csvs = output_dir_results + "ClusterX_clusters_{cohort}_specified_markers_default.csv",

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', "
        "create_csv = TRUE, cohort = '{params.cohort}', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule ClusterX_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        #cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",

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
        "params = list(cells = '{input.cells}', markers = 'all_markers', "
        "create_csv = TRUE, cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "


# Create figures for all of the above methods
rule create_alluvials:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        clusters = output_path + "results/{method}_clusters_{cohort}_{markers}_{clusters}.csv"

    params:
        cohort = "{cohort}",
        method = "{method}",
        clustering_params = "{markers}_clusters_{clusters}"

    output:
        pdf = output_dir_fig + "Alluv_cohort_{cohort}_method_{method}_markers_{markers}_clusters_{clusters}.pdf"
 
    shell:
        "Rscript pipeline/clustering-comparison/alluvials.R {input.cells} {input.cellTypes} {input.clusters} "
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
        lin = output_dir_results + "Assessment-heatmap-lin_cycif.csv",

        basel_indiv = output_dir_results + "Assessment-individual-heatmap-basel.csv",
        schapiro_indiv = output_dir_results + "Assessment-individual-heatmap-schapiro.csv",
        wagner_indiv = output_dir_results + "Assessment-individual-heatmap-wagner.csv",
        zurich1_indiv = output_dir_results + "Assessment-individual-heatmap-zurich1.csv",
        lin_indiv = output_dir_results + "Assessment-individual-heatmap-lin_cycif.csv",

    output:
        html = output_dir_reports + "Final-approaches-comparison.html",
        heatmap = output_dir_fig + "Final_benchmarking_heatmap.pdf"

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Final-approaches-comparison.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(basel = '{input.basel}', schapiro = '{input.schapiro}', wagner = '{input.wagner}', zurich1 = '{input.zurich1}', lin = '{input.lin}', "
        "basel_indiv = '{input.basel_indiv}', schapiro_indiv = '{input.schapiro_indiv}', wagner_indiv = '{input.wagner_indiv}', "
        "zurich1_indiv = '{input.zurich1_indiv}', lin_indiv = '{input.lin_indiv}', output_dir = '" + output_dir_fig + "'))\" "

gsva_csvs = []
gsva_csvs_tmp = [expand(output_dir_results + "{method}_clusters_{cohort}_{markers}_{clusters}.csv", method = [m], markers = markers_spec, cohort = reduced_cohort_list, clusters = conditions[m]) for m in conditions.keys()]
for element in gsva_csvs_tmp:
	gsva_csvs.extend(element)

gsva_csv_dict = {cohort: [f for f in gsva_csvs if cohort in f] for cohort in reduced_cohort_list}    

rule GSVA_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
        csvs = output_dir_results + "{method}_clusters_{cohort}_all_markers_{clusters}.csv"
        
    output:
        html = output_dir_reports + "GSVA/GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_all-markers_{markers_rem}.html",
        csv = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_all-markers_{markers_rem}.csv"

    resources:
        mem_mb=5000

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/GSVA.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "GSVA/', "
        "params = list(cells = '{input.cells}', markers = '{input.markers}', csvs = '{input.csvs}', cohort = '{wildcards.cohort}', "
        "markers_rem = '{wildcards.markers_rem}', output = '" + output_dir_results + "' ))\" "


rule GSVA_specified_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
        csvs = output_dir_results + "{method}_clusters_{cohort}_specified_markers_{clusters}.csv"
        
    output:
        html = output_dir_reports + "GSVA/GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_specified-markers_{markers_rem}.html",
        csv = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_specified-markers_{markers_rem}.csv"

    resources:
        mem_mb=5000

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/GSVA.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "GSVA/', "
        "params = list(cells = '{input.cells}', markers = '{input.markers}', csvs = '{input.csvs}', cohort = '{wildcards.cohort}', "
        "markers_rem = '{wildcards.markers_rem}', output = '" + output_dir_results + "' ))\" "


rule GSVA_create_alluvs:
    input:
        none_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_None.csv",
        stromal_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_Stromal.csv",
        strom_macr_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_Stromal-Macrophage.csv",
        strom_macr_end_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_Stromal-Macrophage-Endothelial.csv"

    output:
        output_dir_fig + "robustness/GSVA-robustness-alluv_{cohort}_{method}_{clusters}_{markers}.pdf"

    shell:
        "Rscript pipeline/clustering-comparison/Other_methods_robustness_alluvials.R {input.none_rem} {input.stromal_rem} {input.strom_macr_rem} {input.strom_macr_end_rem} "
        "{wildcards.cohort} {wildcards.method} {wildcards.clusters} {wildcards.markers} " + output_dir_fig + "robustness/" 