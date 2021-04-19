epi_output_dir_reports = output_path + "reports/epithelial_overclustering/"
epi_output_dir_res = output_path + "results/epithelial_overclustering/"
epi_output_dir_fig = output_path + "figures/epithelial_overclustering/"
epi_shared_output_dir = "/home/campbell/share/projects/imc-2020/output/" + config['version'] + "/results/"
#epi_shared_output_dir = "/home/ltri/campbell/share/projects/imc-2020/output/squirrel/results/"

percent_luminal = [0.99, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25]
marker_options = ["all_markers", "specified_markers"]

UMAP_viz = []
UMAP_viz_tmp = [expand(epi_output_dir_fig + 'Epithelial_overclustering_viz-{method}-{marker_category}-{clusters}.pdf', clusters=conditions[m], method = [m], marker_category = marker_options) for m in conditions.keys()]
for element in UMAP_viz_tmp:
    UMAP_viz.extend(UMAP_viz_tmp)

epithelial_overclustering = {
    'basel_dataset_subsets': expand(epi_shared_output_dir + 'epithelial_overclustering/luminal-{perc}.csv', perc = percent_luminal),
    'ClusterX_assignments_html': expand(epi_output_dir_reports + 'Epithelial_overclustering_ClusterX_clusters-{marker_category}-default-{perc}.html', marker_category = marker_options, perc = percent_luminal),
    'ClusterX_assignments': expand(epi_shared_output_dir + 'epithelial_overclustering/Epithelial_overclustering_ClusterX_clusters-{marker_category}-default-{perc}.csv', marker_category = marker_options, perc = percent_luminal),
    'Phenograph_assignments_html': expand(epi_output_dir_reports + 'Epithelial_overclustering_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.html', marker_category = marker_options, clusters = Phenograph_Cluster_sizes, perc = percent_luminal),
    'Phenograph_assignments_csv': expand(epi_shared_output_dir + 'epithelial_overclustering/Epithelial_overclustering_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.csv', marker_category = marker_options, clusters = Phenograph_Cluster_sizes, perc = percent_luminal),

    #'Phenograph_cluster_enrichmen_assignments_html': expand(epi_output_dir_reports + 'Epithelial_overclustering_cluster_enrichment_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.html', marker_category = marker_options, clusters = Phenograph_Cluster_sizes, perc = percent_luminal),
    #'Phenograph_cluster_enrichment_assignments_csv': expand(epi_shared_output_dir + 'epithelial_overclustering/Epithelial_overclustering_cluster_enrichment_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.csv', marker_category = marker_options, clusters = Phenograph_Cluster_sizes, perc = percent_luminal),
    
    'FlowSOM_assignments_html': expand(epi_output_dir_reports + 'Epithelial_overclustering_FlowSOM_clusters-{marker_category}-k{clusters}-{perc}.html', marker_category = marker_options, clusters = FlowSOM_Cluster_sizes, perc = percent_luminal),
    'FlowSOM_assignments_csv': expand(epi_shared_output_dir + 'epithelial_overclustering/Epithelial_overclustering_FlowSOM_clusters-{marker_category}-k{clusters}-{perc}.csv', marker_category = marker_options, clusters = FlowSOM_Cluster_sizes, perc = percent_luminal),
    'UMAPs': expand(epi_shared_output_dir + "epithelial_overclustering/UMAPs/UMAP-luminal-{perc}-{marker_category}.csv", perc = percent_luminal, marker_category = marker_options),
    #'UMAP_viz': UMAP_viz,
}


rule create_cell_type_proportion_datasets:
    input:
        basel_types = output_path + "astir_assignments/basel_astir_assignments.csv"
    output:
        epi_shared_output_dir + "epithelial_overclustering/luminal-{perc}.csv"
    shell:
        "Rscript pipeline/epithelial-overclustering/create-cell-type-proportion-datasets.R {input.basel_types} {wildcards.perc} " + epi_shared_output_dir + "epithelial_overclustering/"


rule ClusterX_epithelial_overclustering:
    input: 
        cells = output_path + "sces/basel_sce.rds",
        cell_subset = output_path + "results/epithelial_overclustering/luminal-{perc}.csv",
        markers = lambda wildcards: config['basel']['marker_file']
    output:
        html = epi_output_dir_reports + "Epithelial_overclustering_ClusterX_clusters-{marker_category}-default-{perc}.html",
        csv = epi_shared_output_dir + "epithelial_overclustering/Epithelial_overclustering_ClusterX_clusters-{marker_category}-default-{perc}.csv"
    container:
        "astir-manuscript.sif"
    resources:
        mem_mb=5000
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/epithelial-overclustering/Epithelial-overclustering-ClusterX.Rmd', output_file = '{output.html}', "
        "output_dir = '" + epi_output_dir_reports + "', params = list(markers = '{wildcards.marker_category}', cells = '{input.cells}', "
        "cellSubset = '{input.cell_subset}', markers_list = '{input.markers}', percent_luminal = '{wildcards.perc}', output_results = '" + epi_shared_output_dir + "epithelial_overclustering/'))\" "


rule Phenograph_epithelial_overclustering:
    input:
        markers = lambda wildcards: config['basel']['marker_file'],
        cells = output_path + "sces/basel_sce.rds",
        cell_subset = output_path + "results/epithelial_overclustering/luminal-{perc}.csv"
    output:
        html = epi_output_dir_reports + "Epithelial_overclustering_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.html",
        csv = epi_shared_output_dir + "epithelial_overclustering/Epithelial_overclustering_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.csv"
    container: "astir-manuscript.sif"
    resources:
        mem_mb=5000
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/epithelial-overclustering/Epithelial-overclustering-Phenograph.Rmd', output_file = '{output.html}', "
        "output_dir = '" + epi_output_dir_reports + "', params = list(markers = '{wildcards.marker_category}', cells = '{input.cells}', "
        "cellSubset = '{input.cell_subset}', markers_list = '{input.markers}', percent_luminal = '{wildcards.perc}', "
        "output_results = '" + epi_shared_output_dir + "epithelial_overclustering/', cluster_options = '{wildcards.clusters}'))\" "

rule Phenograph_epithelial_overclustering_cluster_enrichment:
    input:
        markers = lambda wildcards: config['basel']['marker_file'],
        cells = output_path + "sces/basel_sce.rds",
        cell_subset = output_path + "results/epithelial_overclustering/luminal-{perc}.csv"
    output:
        html = epi_output_dir_reports + "Epithelial_overclustering_cluster_enrichment_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.html",
        csv = epi_shared_output_dir + "epithelial_overclustering/Epithelial_overclustering_cluster_enrichment_Phenograph_clusters-{marker_category}-k{clusters}-{perc}.csv"
    container: "astir-manuscript.sif"
    resources:
        mem_mb=5000
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/epithelial-overclustering/Epithelial-overclustering-cluster-enrichment-Phenograph.Rmd', output_file = '{output.html}', "
        "output_dir = '" + epi_output_dir_reports + "', params = list(markers = '{wildcards.marker_category}', cells = '{input.cells}', "
        "cellSubset = '{input.cell_subset}', markers_list = '{input.markers}', percent_luminal = '{wildcards.perc}', "
        "output_results = '" + epi_shared_output_dir + "epithelial_overclustering/', cluster_options = '{wildcards.clusters}'))\" "


rule FlowSOM_epithelial_overclustering:
    input:
        markers = lambda wildcards: config['basel']['marker_file'],
        cells = output_path + "sces/basel_sce.rds",
        cell_subset = output_path + "results/epithelial_overclustering/luminal-{perc}.csv"
    output:
        html = epi_output_dir_reports + "Epithelial_overclustering_FlowSOM_clusters-{marker_category}-k{clusters}-{perc}.html",
        csv = epi_shared_output_dir + "epithelial_overclustering/Epithelial_overclustering_FlowSOM_clusters-{marker_category}-k{clusters}-{perc}.csv"
    container: "astir-manuscript.sif"
    resources:
        mem_mb=5000
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/epithelial-overclustering/Epithelial-overclustering-FlowSOM.Rmd', output_file = '{output.html}', "
        "output_dir = '" + epi_output_dir_reports + "', params = list(markers = '{wildcards.marker_category}', cells = '{input.cells}', "
        "cellSubset = '{input.cell_subset}', markers_list = '{input.markers}', percent_luminal = '{wildcards.perc}', output_results = '" + epi_shared_output_dir + 
        "epithelial_overclustering/', cluster_options = '{wildcards.clusters}'))\" "



rule create_UMAP_embedding:
    input:
        cell_subset = epi_shared_output_dir + "epithelial_overclustering/luminal-{perc}.csv",
        cell_expression = output_path + "sces/basel_sce.rds",
        markers = lambda wildcards: config['basel']['marker_file'],
    output:
        epi_shared_output_dir + "epithelial_overclustering/UMAPs/UMAP-luminal-{perc}-{marker_category}.csv"
    shell:
        "Rscript pipeline/epithelial-overclustering/create-UMAPs.R {input.cell_subset} {input.cell_expression} " + epi_shared_output_dir + "epithelial_overclustering/ {wildcards.marker_category} {input.markers} {wildcards.perc}" 




umap = []
umapTMP = [expand(epi_output_dir_res + 'UMAPs/UMAP-luminal-{perc}-{{marker_category}}.csv', perc = [m]) for m in percent_luminal]
for element in umapTMP:
    umap.extend(element)

assignments_other_methods = []
assignments_other_methods_TMP = [expand(epi_output_dir_res + "Epithelial_overclustering_{{method}}_clusters-{{marker_category}}-{{clusters}}-{perc}.csv", perc = percent_luminal)]
for element in assignments_other_methods_TMP:
    assignments_other_methods.extend(element)


rule visualize_overclustering:
    input:
        umaps = umap,
        cluster_assignments = assignments_other_methods,
        astir_assignments = expand(epi_shared_output_dir + "imbalance/astir_imbalance_{perc}.csv", perc = percent_luminal)
    output:
        html = epi_output_dir_reports + "Epithelial_overclustering_viz-{method}-{marker_category}-{clusters}.html",
        pdf = epi_output_dir_fig + "Epithelial_overclustering_viz-{method}-{marker_category}-{clusters}.pdf"
    shell:
        "Rscript -e \" Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/epithelial-overclustering/Epithelial-overclustering-FlowSOM.Rmd', output_file = '{output.html}', "
        "output_dir = '" + epi_output_dir_reports + "', params = list(umaps = '{input.umaps}', assignments = '{input.cluster_assignments}', "
        "astir_assignments = '{input.astir_assignments}', output_results = '" + epi_output_dir_fig + "'))\""