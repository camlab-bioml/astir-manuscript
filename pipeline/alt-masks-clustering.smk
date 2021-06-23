##### ALTERNATIVE CLUSTERING FOR ALTERNATIVE MASKS #####
alt_cores = ['Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35']
users = ['Catena', 'Jackson', 'Schulz']
acdc_options = [ 'absent', 'no-consider']
seeds = range(10)
output_dir_results_altMasks = output_path + "results/alternative_masks/"
output_dir_reports_altMasks = output_path + "reports/alternative_masks/"
output_dir_fig_altMasks = output_path + "figures/alternative_masks/"
alt_masks_shared_output_dir = "/home/campbell/share/projects/imc-2020/output/" + config['version'] + "/results/alternative_masks/"

markers_options = ['specified_markers', 'all_markers']

alt_methods_cell_assign = []
alt_methods_cell_assign_tmp = [expand(output_dir_reports_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.html", 
core = alt_cores, user = users, method = [m], clusters = conditions[m], markers = markers_options) for m in ['Phenograph', 'ClusterX']]
for element in alt_methods_cell_assign_tmp:
	alt_methods_cell_assign.extend(element)

alt_methods_cell_assign_csv = []
alt_methods_cell_assign_csv_tmp = [expand(output_dir_results_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.csv", 
core = alt_cores, user = users, method = [m], clusters = conditions[m], markers = markers_options) for m in ['Phenograph', 'ClusterX']]
for element in alt_methods_cell_assign_csv_tmp:
	alt_methods_cell_assign_csv.extend(element)


alt_methods_cell_assign_seed = []
alt_methods_cell_assign_seed_tmp = [expand(output_dir_reports_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-seed_{seed}-{markers}-{clusters}.html", 
core = alt_cores, user = users, method = [m], clusters = conditions[m], markers = markers_options, seed = seeds) for m in ['FlowSOM']]
for element in alt_methods_cell_assign_seed_tmp:
	alt_methods_cell_assign_seed.extend(element)

alt_methods_cell_assign_csv_seed = []
alt_methods_cell_assign_csv_seed_tmp = [expand(output_dir_results_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-seed_{seed}-{markers}-{clusters}.csv", 
core = alt_cores, user = users, method = [m], clusters = conditions[m], markers = markers_options, seed = seeds) for m in ['FlowSOM']]
for element in alt_methods_cell_assign_csv_seed_tmp:
	alt_methods_cell_assign_csv_seed.extend(element)


alt_masks = {
    'sce': expand(output_path + "sces/schapiro_alt-{core}-{user}-sce.rds", core = alt_cores, user = users),
    'Phenograph_analysis_report': expand(output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-specified_markers_k{cluster}.html", core = alt_cores, user = users, cluster = Phenograph_Cluster_sizes),
    'Phenograph_clusters_output': expand(output_dir_results_altMasks + "Phenograph_clusters_{core}-{user}_specified_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = Phenograph_Cluster_sizes),
    'Phenograph_analysis_report_all': expand(output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-all_markers_k{cluster}.html", core = alt_cores, user = users, cluster = Phenograph_Cluster_sizes),
    'Phenograph_clusters_output_all': expand(output_dir_results_altMasks + "Phenograph_clusters_{core}-{user}_all_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = Phenograph_Cluster_sizes),
    'FlowSOM_analysis_report': expand(output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-seed_{seed}-specified_markers_k{cluster}.html", core = alt_cores, user = users, cluster = FlowSOM_Cluster_sizes, seed = seeds),
    'FlowSOM_clusters_output': expand(output_dir_results_altMasks + "FlowSOM_clusters_{core}-{user}_seed-{seed}_specified_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = FlowSOM_Cluster_sizes, seed = seeds),
    'FlowSOM_analysis_report_all': expand(output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-seed_{seed}-all_markers_k{cluster}.html", core = alt_cores, user = users, cluster = FlowSOM_Cluster_sizes, seed = seeds),
    'FlowSOM_clusters_output_all': expand(output_dir_results_altMasks + "FlowSOM_clusters_{core}-{user}_seed-{seed}_all_markers_k{cluster}.csv", core = alt_cores, user = users, cluster = FlowSOM_Cluster_sizes, seed = seeds),
    'ClusterX_analysis_report': expand(output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-specified_markers.html", core = alt_cores, user = users, seed = seeds),
    'ClusterX_clusters_output': expand(output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_specified_markers_default.csv", core = alt_cores, user = users, seed = seeds),
    'ClusterX_analysis_report_all': expand(output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-all_markers.html", core = alt_cores, user = users, seed = seeds),
    'ClusterX_clusters_output_all': expand(output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_all_markers_default.csv", core = alt_cores, user = users, seed = seeds),
    
    'alt-mask-cell-type-identification-html': alt_methods_cell_assign,
    'alt-mask-cell-type-identification-csv': alt_methods_cell_assign_csv,
    'alt-mask-cell-type-identification-html-seed': alt_methods_cell_assign_seed,
    'alt-mask-cell-type-identification-csv-seed': alt_methods_cell_assign_csv_seed,

    'alt-mask-h5ads': expand(output_path + 'anndata/schapiro_alt-{core}-{user}.h5ad', core = alt_cores, user = users),
    'acdc': expand(output_dir_results_altMasks + 'ACDC_clusters_{core}-{user}_seed-{seed}_{options}.csv', core = alt_cores, user = users, options = acdc_options, seed = seeds),
    #'acdc_assignment': expand(output_dir_results_altMasks + "Alternative_masks-ACDC-cell_type_assignments-{core}-{user}-seed_{seed}.csv", core = alt_cores, user = users, seed = seeds)    
}

rule create_sces_altMasks: 
    params:
        schapiro = output_path + "schapiro_processed_alt_mask/{core}_{user}.rds"
    container: "astir-manuscript.sif"
    output:
        schapiro = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"
    shell:
        "Rscript pipeline/rds_to_sce.R --rds {params.schapiro} --output {output.schapiro}"


rule create_anndata_alt_masks:
    input:
        sce = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"
    output:
        h5ad = output_path + 'anndata/schapiro_alt-{core}-{user}.h5ad'
    script:
        "create_h5ad_from_sce.R"


rule phenograph_analysis_altMasks:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        markers = lambda wildcards: config['schapiro']['marker_file'],
    params:
        res_dir = output_dir_results_altMasks
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-specified_markers_k{cluster}.html",
        csvs = output_dir_results_altMasks + "Phenograph_clusters_{core}-{user}_specified_markers_k{cluster}.csv"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "', "
        "params = list(cells = '{input.cells}', cohort = '{wildcards.core}-{wildcards.user}', markers = 'specified_markers', "
        "markers_list = '{input.markers}', output_results = '{params.res_dir}', cluster_options = '{wildcards.cluster}'))\" "


rule phenograph_analysis_all_markers_altMasks:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"
    params:
        res_dir = output_dir_results_altMasks
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "Rphenograph-analysis-{core}-{user}-all_markers_k{cluster}.html",
        csvs = output_dir_results_altMasks + "Phenograph_clusters_{core}-{user}_all_markers_k{cluster}.csv"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/clustering-comparison/Rphenograph.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', cohort = '{wildcards.core}-{wildcards.user}', markers = 'all_markers', output_results = '{params.res_dir}', "
        "cluster_options = '{wildcards.cluster}' ))\" "


rule FlowSOM_analysis_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        markers = lambda wildcards: config['schapiro']['marker_file']
    params:
        res_dir = output_dir_results_altMasks  
    container: "astir-manuscript.sif"  
    output:
        html = output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-seed_{seed}-specified_markers_k{cluster}.html",
        csvs = output_dir_results_altMasks + "FlowSOM_clusters_{core}-{user}_seed-{seed}_specified_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', FlowSOM_seed = '{wildcards.seed}', cluster_options = '{wildcards.cluster}', "
        "cohort = '{wildcards.core}-{wildcards.user}',  markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "



rule FlowSOM_analysis_all_markers_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds"    
    params:
        res_dir = output_dir_results_altMasks
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "FlowSOM-analysis-{core}-{user}-seed_{seed}-all_markers_k{cluster}.html",
        csvs = output_dir_results_altMasks + "FlowSOM_clusters_{core}-{user}_seed-{seed}_all_markers_k{cluster}.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/FlowSOM.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', cohort = '{wildcards.core}-{wildcards.user}', output_results = '{params.res_dir}', "
        "cluster_options = '{wildcards.cluster}', FlowSOM_seed = '{wildcards.seed}'))\" "

    
rule ClusterX_analysis_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        markers = lambda wildcards: config['schapiro']['marker_file'],
    params:
        res_dir = output_dir_results_altMasks
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-specified_markers.html",
        csvs = output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_specified_markers_default.csv",
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'specified_markers', "
        "cohort = '{wildcards.core}-{wildcards.user}', markers_list = '{input.markers}', output_results = '{params.res_dir}'))\" "


rule ClusterX_analysis_all_markers_altMarkers:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
    params:
        res_dir = output_dir_results_altMasks
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "ClusterX-analysis-{core}-{user}-all_markers.html",
        csvs = output_dir_results_altMasks + "ClusterX_clusters_{core}-{user}_all_markers_default.csv"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/clustering-comparison/ClusterX.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "',"
        "params = list(cells = '{input.cells}', markers = 'all_markers', cohort = '{wildcards.core}-{wildcards.user}', output_results = '{params.res_dir}'))\" "


rule acdc_alt_masks:
    input:
        h5ad = output_path + 'anndata/schapiro_alt-{core}-{user}.h5ad',
        markers = config['schapiro']['marker_file']
    output:
        csv = output_dir_results_altMasks + 'ACDC_clusters_{core}-{user}_seed-{seed}_{options}.csv'
    shell:
        "python pipeline/epithelial-overclustering/run-acdc.py "
        "--input_h5ad {input.h5ad} "
        "--input_yaml {input.markers} "
        "--output_assignments {output.csv} "
        "--method {wildcards.options} "
        "--cohort schapiro "
        "--seed {wildcards.seed} "
        "--user {wildcards.user} "
        "--sample {wildcards.core} "

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
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.html",
        csv = output_dir_results_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-{markers}-{clusters}.csv"
    shell:
        "Rscript -e \" rmarkdown::render('pipeline/clustering-comparison/alt-masks-cluster-identification.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "', "
        "params = list(cells = '{input.cells}', csvs = '{input.csvs}', markers_list = '{input.markers_yml}', core = '{wildcards.core}', user = '{wildcards.user}', "
        "method = '{wildcards.method}', markers = '{wildcards.markers}', cluster = '{wildcards.clusters}', output_dir = '" + output_dir_results_altMasks + "' ))\" "


rule get_altMasks_cell_types_seed:
    input:
        cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
        csvs = output_dir_results_altMasks + "{method}_clusters_{core}-{user}_seed-{seed}_{markers}_{clusters}.csv",
        markers_yml = config['schapiro']['marker_file']
    container: "astir-manuscript.sif"
    output:
        html = output_dir_reports_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-seed_{seed}-{markers}-{clusters}.html",
        csv = output_dir_results_altMasks + "Alternative_masks-{method}-cell_type_assignments-{core}-{user}-seed_{seed}-{markers}-{clusters}.csv"
    shell:
        "Rscript -e \" rmarkdown::render('pipeline/clustering-comparison/alt-masks-cluster-identification.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports_altMasks + "', "
        "params = list(cells = '{input.cells}', csvs = '{input.csvs}', markers_list = '{input.markers_yml}', core = '{wildcards.core}', user = '{wildcards.user}', "
        "method = '{wildcards.method}', markers = '{wildcards.markers}', cluster = '{wildcards.clusters}', output_dir = '" + output_dir_results_altMasks + "', seed = '{wildcards.seed}' ))\" "



# rule get_altMasks_cell_types_acdc:
#     input:
#         cells = output_path + "sces/schapiro_alt-{core}-{user}-sce.rds",
#         csvs = output_dir_results_altMasks + "ACDC_clusters_{core}-{user}_seed-{seed}_{options}.csv",
#         markers_yml = config['schapiro']['marker_file']
#     container: "astir-manuscript.sif"
#     output:
#         csv = output_dir_results_altMasks + "Alternative_masks-ACDC-cell_type_assignments-{core}-{user}-seed_{seed}.csv"
#     shell:
#         "Rscript -e \" rmarkdown::render('pipeline/clustering-comparison/alt-masks-cluster-identification.Rmd', output_dir = '" + output_dir_reports_altMasks + "', "
#         "params = list(cells = '{input.cells}', csvs = '{input.csvs}', markers_list = '{input.markers_yml}', core = '{wildcards.core}', user = '{wildcards.user}', "
#         "method = 'ACDC', markers = '{wildcards.markers}', cluster = '{wildcards.clusters}', output_dir = '" + output_dir_results_altMasks + "', seed = '{wildcards.seed}' ))\" "
