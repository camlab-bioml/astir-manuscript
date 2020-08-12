


spatial_output = {
    'basel': expand(output_path + "spatial/{dataset}/{core}.csv", dataset=['basel'], core=basel_cores),
    'zurich': expand(output_path + "spatial/{dataset}/{core}.csv", dataset=['zurich1'], core=zurich1_cores),
    'reports': expand(output_path + "spatial/spatial_report_{dataset}.html", dataset=['basel','zurich1']),
    'kaplan_meier': output_path + "figures/spatial/km.pdf",
    'pathway_fig': output_path + "figures/spatial/pathway_fig.pdf",
    'spatial_heatmap': output_path + "figures/spatial/heatmap.pdf"
}

rule calc_dists:
    params:
    input:
        markers=config['basel']['marker_file'],
        assignments=output_path + "astir_assignments/{dataset}_astir_assignments.csv",
        locations="data-raw/jackson-2020/spatial-locations/{dataset}_SC_locations.csv"
    output:
        output_path + "spatial/{dataset}/{core}.csv"
    shell:
        "Rscript pipeline/spatial/spatial-celltypes.R "
        "--input_locations {input.locations} "
        "--input_assignments {input.assignments} "
        "--marker_yaml {input.markers} "
        "--core {wildcards.core} "
        "--output_csv {output}"

rule dist_report:
    params:
        curr_dir = os.getcwd(),
        output_path = output_path,
    input:
        expand(output_path + "spatial/basel/{core}.csv", core=basel_cores),
        expand(output_path + "spatial/zurich1/{core}.csv", core=zurich1_cores),
        metadata = lambda wildcards: config[wildcards.dataset]['base_dir']+ config[wildcards.dataset]['metadata_file'],
        celltype_assignments = output_path + "astir_assignments/{dataset}_astir_assignments.csv",
        state_assignments = output_path + "astir_assignments/{dataset}_astir_assignments_state.csv",
    output:
        html=output_path + "spatial/spatial_report_{dataset}.html",
        rds=output_path + "spatial/rds_output_{dataset}.rds",
        scree_plot=output_path + "figures/supplementary/scree-plot-spatial-{dataset}.png",
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin'); "
        "rmarkdown::render('pipeline/spatial/distance-clustering.Rmd',   "
        "output_file='{output.html}', "
        "params=list(input_dir='{params.output_path}/spatial/{wildcards.dataset}/', "
        "celltype_assignments='{input.celltype_assignments}', "
        "state_assignments='{input.state_assignments}', "  
        "output_rds='{output.rds}', "
        "scree_plot='{output.scree_plot}', "
        "dataset='{wildcards.dataset}', "      
        "metadata='{input.metadata}'), "
        "knit_root_dir='{params.curr_dir}', "
        "output_dir='{params.curr_dir}/{params.output_path}/spatial/'"
        ")\"  "

rule spatial_figures:
    params:
        curr_dir = os.getcwd(),
        output_path = output_path,
    input:
        basel=output_path + "spatial/rds_output_basel.rds",
        zurich1=output_path + "spatial/rds_output_zurich1.rds",
        basel_metadata=config['basel']['base_dir']+ config['basel']['metadata_file'],
    output:
        html=output_path + "spatial/spatial_figures.html",
        kaplan_meier = output_path + "figures/spatial/km.pdf",
        pathway_fig = output_path + "figures/spatial/pathway_fig.pdf",
        spatial_heatmap = output_path + "figures/spatial/heatmap.pdf",
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin'); "
        "rmarkdown::render('pipeline/spatial/spatial-figures.Rmd',   "
        "output_file='{output.html}', "
        "params=list(basel_output='{input.basel}', "
        "zurich1_output='{input.zurich1}', "
        "basel_metadata='{input.basel_metadata}', "  
        "kaplan_meier='{output.kaplan_meier}', "
        "pathway_fig='{output.pathway_fig}', "      
        "spatial_heatmap='{output.spatial_heatmap}'), "
        "knit_root_dir='{params.curr_dir}', "
        "output_dir='{params.curr_dir}/{params.output_path}/spatial/'"
        ")\"  "