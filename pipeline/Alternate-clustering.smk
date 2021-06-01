cohort_list = ['basel', 'zurich1', 'wagner', 'schapiro', 'lin_cycif']
reduced_cohort_list = ['basel', 'zurich1', 'wagner']
markers_spec = ['specified_markers', 'all_markers']

output_dir_reports = output_path + "reports/"
output_dir_results = output_path + "results/"
output_dir_fig = output_path + "figures/"


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
    'lin_cycif': cycif_samples
}

# Create list containing markers to be removed from alternative approaches
markers_removal = ['None', 'Stromal', 'Stromal-Macrophage', 'Stromal-Macrophage-Endothelial']
markers_gsva = ['specified-markers', 'all-markers']


acdc_results = []
for cohort in ['basel', 'zurich1', 'wagner', 'schapiro', 'lin-cycif']:
    for options in [ 'absent', 'no-consider']:
        acdc_results.append(
            output_dir_results + f"ACDC_clusters_{cohort}_{options}.csv"
        )

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

robustness_alluvials = []
robustness_alluvials_tmp = [expand(output_dir_fig + "robustness/GSVA-robustness-alluv_{cohort}_{method}_{clusters}_{markers}.pdf", 
cohort = reduced_cohort_list, method = [m], clusters = conditions[m], markers = markers_gsva) for m in conditions.keys()]
for element in robustness_alluvials_tmp:
    robustness_alluvials.extend(element)

GSVA_cell_type_assignments = []
GSVA_assignments = [expand(output_dir_results + "GSVA-assignment-{method}-{cohort}-{markers}-{methodClusters}.csv", 
method = [m], cohort = cohort_list, markers = markers_spec, methodClusters = conditions[m]) for m in conditions.keys()]
for element in GSVA_assignments:
    GSVA_cell_type_assignments.extend(element)


alternate_approaches_output = {
    'sces': expand(output_path + "sces/{cohort}_sce.rds", cohort = cohort_list),
    'Phenograph_analysis_report': expand(output_dir_reports + "Rphenograph-analysis-{cohort}-specified_markers_k{cluster}.html", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
    'Phenograph_clusters_output': expand(output_dir_results + "Phenograph_clusters_{cohort}_specified_markers_k{cluster}.csv", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
    'Phenograph_analysis_report_all': expand(output_dir_reports + "Rphenograph-analysis-{cohort}-all_markers_k{cluster}.html", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
    'Phenograph_clusters_output_all': expand(output_dir_results + "Phenograph_clusters_{cohort}_all_markers_k{cluster}.csv", cohort = cohort_list, cluster = Phenograph_Cluster_sizes),
    'FlowSOM_analysis_report': expand(output_dir_reports + "FlowSOM-analysis-{cohort}-specified_markers_k{cluster}.html", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_clusters_output': expand(output_dir_results + "FlowSOM_clusters_{cohort}_specified_markers_k{cluster}.csv", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_analysis_report_all': expand(output_dir_reports + "FlowSOM-analysis-{cohort}-all_markers_k{cluster}.html", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),
    'FlowSOM_clusters_output_all': expand(output_dir_results + "FlowSOM_clusters_{cohort}_all_markers_k{cluster}.csv", cohort = cohort_list, cluster = FlowSOM_Cluster_sizes),

    'FlowSOM_100k': expand (output_dir_results + "Flowsom_k90/FlowSOM_clusters_basel_{markers}_k90.csv", markers = markers_spec),
    'FlowSOM_100k': expand(output_dir_results + "Flowsom_k90/GSVA-assignment-FlowSOM-basel-{markers}-k90.csv", markers = markers_spec),
    'FlowSOM_100k': expand(output_dir_fig + "Flowsom_k90/FlowSOM-overclustering-k90-heatmap-{markers}.pdf", markers = markers_spec),

    'ClusterX_analysis_report': expand(output_dir_reports + "ClusterX-analysis-{cohort}-specified_markers.html", cohort = cohort_list),
    'ClusterX_clusters_output': expand(output_dir_results + "ClusterX_clusters_{cohort}_specified_markers_default.csv", cohort = cohort_list),
    'ClusterX_analysis_report_all': expand(output_dir_reports + "ClusterX-analysis-{cohort}-all_markers.html", cohort = cohort_list),
    'ClusterX_clusters_output_all': expand(output_dir_results + "ClusterX_clusters_{cohort}_all_markers_default.csv", cohort = cohort_list),
    'acdc_clusters_output': acdc_results,

    'GSVA_cell_type_assignments': GSVA_cell_type_assignments,

    # # # Analysis
    'Alluvial': alluvialPDF,
    # # 'ExpressionHeatmaps': expressionPDF,
    # # 'statesBoxPlots': statesBoxPDF,

    'final_benchmark_heatmap': output_dir_fig + "Final_benchmarking_heatmap.pdf",
    ###'GSVA_csv_all': gsva_robustness_all,
    ###'GSVA_html_all': gsva_assignments_reports_all,
    'GSVA_csv_spec': gsva_robustness_spec,
    'GSVA_html_spec': gsva_assignments_reports_spec,
    #'GSVA_html': expand(output_dir_reports + "GSVA_{markers_rem}_{cohort}.html", cohort = reduced_cohort_list, markers_rem = markers_removal),
    #'GSVA_csv': expand(output_dir_results + "Other_approaches_GSVA_{markers_rem}_{cohort}.csv", cohort = reduced_cohort_list, markers_rem = markers_removal),

    'GSVA_robust': robustness_alluvials,
    #'Astir-robustness': expand(output_dir_fig + "Astir-robustness-{cohort}.pdf", cohort = reduced_cohort_list)
}

cores_dict = {
    'basel': expand(output_path + "basel_processed/{core}.rds", core = basel_cores),
    'zurich1': expand(output_path + "zurich1_processed/{core}.rds", core = zurich1_cores),
    'wagner': expand(output_path + "wagner_processed/{core}.rds", core = wagner_samples),
    'schapiro': expand(output_path + "schapiro_processed/{core}.rds", core = schapiro_samples),
    'lin_cycif': expand(output_path + "lin-cycif_processed/{core}.rds", core = cycif_samples)
}


rule create_sces: 
    input:
        lambda wildcards: cores_dict[wildcards.cohort]
    container: "astir-manuscript.sif"
    output:
        output_path + "sces/{cohort}_sce.rds"
    shell:
        "Rscript pipeline/rds_to_sce.R --rds {input} --output {output}"
        

rule phenograph_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources:
        mem_mb=40000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "Rphenograph-analysis-{cohort}-specified_markers_k{cluster}.html",
        csvs = output_dir_results + "Phenograph_clusters_{cohort}_specified_markers_k{cluster}.csv" 
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', cohort = '{params.cohort}', markers = 'specified_markers', "
        "markers_list = '{input.markers}', output_results = '{params.res_dir}', cluster_options = '{wildcards.cluster}'))\" "


rule phenograph_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds"  
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources:
        mem_mb=40000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "Rphenograph-analysis-{cohort}-all_markers_k{cluster}.html",
        csvs = output_dir_results + "Phenograph_clusters_{cohort}_all_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', cohort = '{params.cohort}', "
        "markers = 'all_markers', output_results = '{params.res_dir}', cluster_options = '{wildcards.cluster}' ))\" "


rule FlowSOM_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources: mem_mb=2000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "FlowSOM-analysis-{cohort}-specified_markers_k{cluster}.html",
        csvs = output_dir_results + "FlowSOM_clusters_{cohort}_specified_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', cluster_options = '{wildcards.cluster}', "
        "cohort = '{params.cohort}',  markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "



rule FlowSOM_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds"    
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources: mem_mb=2000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "FlowSOM-analysis-{cohort}-all_markers_k{cluster}.html",
        csvs = output_dir_results + "FlowSOM_clusters_{cohort}_all_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', cluster_options = '{wildcards.cluster}',"
        "cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "


rule FlowSOM_100k:
    input:
        cells = output_path +  "sces/basel_sce.rds",
        markers = config['basel']['marker_file']
    container: "astir-manuscript.sif"
    params:
        cohort = "basel",
        res_dir = output_dir_results + "Flowsom_k90/"
    output:
        output_dir_results + "Flowsom_k90/FlowSOM_clusters_basel_{markers}_k90.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', markers = '{wildcards.markers}', cluster_options = '90', "
        "cohort = 'basel',  markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule FlowSOM_k90_interpret:
    input:
        sce = output_path + "sces/basel_sce.rds",
        markers = config['basel']['marker_file'],
        FlowSOM = output_dir_results + "Flowsom_k90/FlowSOM_clusters_basel_{markers}_k90.csv"
    output:
        csv = output_dir_results + "Flowsom_k90/GSVA-assignment-FlowSOM-basel-{markers}-k90.csv",
        heatmap = output_dir_fig + "Flowsom_k90/FlowSOM-overclustering-k90-heatmap-{markers}.pdf"
    container: "astir-manuscript.sif"
    script:
        "clustering-comparison/FlowSOM_90k_overclustering.R"

    
rule ClusterX_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources: mem_mb=30000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "ClusterX-analysis-{cohort}-specified_markers.html",
        csvs = output_dir_results + "ClusterX_clusters_{cohort}_specified_markers_default.csv",
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', "
        "cohort = '{params.cohort}', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule ClusterX_analysis_all_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds"
    params:
        cohort = "{cohort}",
        res_dir = output_dir_results
    resources: mem_mb=30000
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports + "ClusterX-analysis-{cohort}-all_markers.html",
        csvs = output_dir_results + "ClusterX_clusters_{cohort}_all_markers_default.csv"

    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', "
        "cohort = '{params.cohort}', output_results = '{params.res_dir}'))\" "

        

rule create_anndata:
    input:
        lambda wildcards: expand(output_path + wildcards.cohort + "_processed/{core}.csv", core=get_core_list(wildcards.cohort)),
    output:
        output_path + "anndata/{cohort}.h5ad",
    shell:
        "python pipeline/cla/dir-of-csvs-to-scanpy.py "
        "{output_path}/{wildcards.cohort}_processed "
        "{output} "


rule acdc:
    input:
        h5ad = output_path + "anndata/{cohort}.h5ad",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']
    output:
        csv = output_dir_results + 'ACDC_clusters_{cohort}_{options}.csv'
    shell:
        "python pipeline/cla/run-acdc.py "
        "--input_h5ad {input.h5ad} "
        "--input_yaml {input.markers} "
        "--output_assignments {output.csv} "
        "--method {wildcards.options} "
        "--cohort {wildcards.cohort} "



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
    container: "astir-manuscript.sif"
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


rule assign_GSVA_cell_types:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
        clustering_result = output_dir_results + "{method}_clusters_{cohort}_{markers}_{methodClusters}.csv"
    output:
        html = output_dir_reports + "GSVA-assignment-{method}-{cohort}-{markers}-{methodClusters}.html",
        csv = output_dir_results + "GSVA-assignment-{method}-{cohort}-{markers}-{methodClusters}.csv"
    container: "astir-manuscript.sif"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Assign-cell-type-to-GSVA-cluster.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', markers = '{input.markers}', clustering_result = '{input.clustering_result}', output_file = '{output.csv}' ))\" "


# rule create_state_boxplots:
#     input:
#         cells = output_path + "sces/{cohort}_sce.rds",
#         cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
#         cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
#         clusters = output_path + "results/{method}_clusters_{cohort}_specified_markers_{clusters}.csv",
#     params:
#         cohort = "{cohort}",
#         method = "{method}",
#         clustering_params = "{clusters}"
#     output:
#         expression_pdf = output_dir_fig + "StatesBoxplot_cohort_{cohort}_method_{method}_clusters_{clusters}.png"
#     shell:
#         "Rscript pipeline/clustering-comparison/boxplot_states.R {input.cells} {input.cellTypes} {input.cellStates} {input.clusters} "
#         "{params.cohort} {params.method} {params.clustering_params} " + output_dir_fig +""

# Create list of csvs
csvs = []
csvs_tmp = [expand(output_dir_results + "GSVA-assignment-{method}-{cohort}-{markers}-{clusters}.csv", method = [m], markers = markers_spec, cohort = cohort_list, clusters = conditions[m]) for m in conditions.keys()]
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
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/Approaches-comparison.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cohort = '{params.cohort}', csvs = '{input.csvs}', "
        "markers_list = '{input.markers}', output_dir = '" + output_dir_results + "' ))\" "



rule final_approach_comparison_R:
    input:
        basel_astir_assignments = output_path + "astir_assignments/basel_astir_assignments.csv",
        schapiro_astir_assignments = output_path + "astir_assignments/schapiro_astir_assignments.csv",
        wagner_astir_assignments = output_path + "astir_assignments/wagner_astir_assignments.csv",
        zurich_astir_assignments = output_path + "astir_assignments/zurich1_astir_assignments.csv",
        lin_astir_assignments = output_path + "astir_assignments/lin_cycif_astir_assignments.csv",

        basel_sce = output_path + "sces/basel_sce.rds",
        schapiro_sce = output_path + "sces/schapiro_sce.rds",
        wagner_sce = output_path + "sces/wagner_sce.rds",
        zurich_sce = output_path + "sces/zurich1_sce.rds",
        lin_sce = output_path + "sces/lin_cycif_sce.rds",

        jackson_markers = config['basel']['marker_file'],
        schapiro_markers = config['schapiro']['marker_file'],
        wagner_markers = config['wagner']['marker_file'],
        lin_markers = config['lin_cycif']['marker_file'],

        basel_csvs = csv_dict['basel'],
        schapiro_csvs = csv_dict['schapiro'],
        wagner_csvs = csv_dict['wagner'],
        zurich_csvs = csv_dict['zurich1'],
        lin_csvs = csv_dict['lin_cycif'],
        
        acdc_absent_basel = output_dir_results + "ACDC_clusters_basel_absent.csv",
        acdc_no_consider_basel = output_dir_results + "ACDC_clusters_basel_no-consider.csv",
        acdc_absent_schapiro = output_dir_results + "ACDC_clusters_schapiro_absent.csv",
        acdc_no_consider_schapiro = output_dir_results + "ACDC_clusters_schapiro_no-consider.csv",
        acdc_absent_wagner = output_dir_results + "ACDC_clusters_wagner_absent.csv",
        acdc_no_consider_wagner = output_dir_results + "ACDC_clusters_wagner_no-consider.csv",
        acdc_absent_zurich = output_dir_results + "ACDC_clusters_zurich1_absent.csv",
        acdc_no_consider_zurich = output_dir_results + "ACDC_clusters_zurich1_no-consider.csv",
        acdc_absent_lin = output_dir_results + "ACDC_clusters_lin-cycif_absent.csv",
        acdc_no_consider_lin = output_dir_results + "ACDC_clusters_lin-cycif_no-consider.csv",
    output:
        heatmap = output_dir_fig + "Final_benchmarking_heatmap.pdf"
    resources:
        mem_mb=10000
    container: "astir-manuscript.sif"
    shell:
        "Rscript pipeline/clustering-comparison/Final-approaches-comparison.R "
        "   --basel_astir_assignments {input.basel_astir_assignments} "
        "   --schapiro_astir_assignments {input.schapiro_astir_assignments} "
        "   --wagner_astir_assignments {input.wagner_astir_assignments} "
        "   --zurich_astir_assignments {input.zurich_astir_assignments} "
        "   --lin_astir_assignments {input.lin_astir_assignments} "
        "   --basel_sce {input.basel_sce} "
        "   --schapiro_sce {input.schapiro_sce} "
        "   --wagner_sce {input.wagner_sce} "
        "   --zurich_sce {input.zurich_sce} "
        "   --lin_sce {input.lin_sce} "
        "   --jackson_markers {input.jackson_markers} "
        "   --schapiro_markers {input.schapiro_markers} "
        "   --wagner_markers {input.wagner_markers} "
        "   --lin_markers {input.lin_markers} "
        "   --basel_files {input.basel_csvs} "
        "   --schapiro_files {input.schapiro_csvs} "
        "   --wagner_files {input.wagner_csvs} "
        "   --lin_files {input.lin_csvs} "

        "   --acdc_absent_basel {input.acdc_absent_basel} "
        "   --acdc_no_consider_basel {input.acdc_no_consider_basel} "
        "   --acdc_absent_schapiro {input.acdc_absent_schapiro} "
        "   --acdc_no_consider_schapiro {input.acdc_no_consider_schapiro} "
        "   --acdc_absent_wagner {input.acdc_absent_wagner} "
        "   --acdc_no_consider_wagner {input.acdc_no_consider_wagner} "
        "   --acdc_absent_zurich {input.acdc_absent_zurich} "
        "   --acdc_no_consider_zurich {input.acdc_no_consider_zurich} "
        "   --acdc_absent_lin {input.acdc_absent_lin} "
        "   --acdc_no_consider_lin {input.acdc_no_consider_lin} "

        "   --zurich_files {input.zurich_csvs} "
        "   --output_heatmap {output.heatmap} "

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
    container: 'astir-manuscript.sif'
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/GSVA.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "GSVA/', "
        "params = list(cells = '{input.cells}', markers = '{input.markers}', csvs = '{input.csvs}', cohort = '{wildcards.cohort}', "
        "markers_rem = '{wildcards.markers_rem}', output = '{output.csv}' ))\" "


rule GSVA_specified_markers:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
        csvs = output_dir_results + "{method}_clusters_{cohort}_specified_markers_{clusters}.csv"
    output:
        html = output_dir_reports + "GSVA/GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_specified-markers_{markers_rem}.html",
        csv = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_specified-markers_{markers_rem}.csv"
    container: 'astir-manuscript.sif'
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/GSVA.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "GSVA/', "
        "params = list(cells = '{input.cells}', markers = '{input.markers}', csvs = '{input.csvs}', cohort = '{wildcards.cohort}', "
        "markers_rem = '{wildcards.markers_rem}', output = '{output.csv}' ))\" "


rule GSVA_create_alluvs:
    input:
        none_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_None.csv",
        stromal_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_Stromal.csv",
        strom_macr_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_Stromal-Macrophage.csv",
        strom_macr_end_rem = output_dir_results + "GSVA-Other-methods-robustness_{cohort}_{method}_{clusters}_{markers}_Stromal-Macrophage-Endothelial.csv"
    container: 'astir-manuscript.sif'
    output:
        output_dir_fig + "robustness/GSVA-robustness-alluv_{cohort}_{method}_{clusters}_{markers}.pdf"

    shell:
        "Rscript pipeline/clustering-comparison/Other_methods_robustness_alluvials.R {input.none_rem} {input.stromal_rem} {input.strom_macr_rem} {input.strom_macr_end_rem} "
        "{wildcards.cohort} {wildcards.method} {wildcards.clusters} {wildcards.markers} " + output_dir_fig + "robustness/" 