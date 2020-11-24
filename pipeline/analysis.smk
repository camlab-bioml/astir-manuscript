output_dir_fig = output_path + "figures/"
output_dir_results = output_path + "results/"
output_dir_reports = output_path + "reports/"

cohort_list = ['basel', 'zurich1', 'wagner', 'schapiro', 'lin_cycif']
reduced_cohort_list = ['basel', 'zurich1', 'wagner', 'schapiro']
alpha_values = [1, 0.01, 0.005, 0.001, 0.000001]

analysis_output = {
    'tSNE_clustering_cellType': expand(output_dir_fig + 'tSNE_cellType-{alpha}.png', alpha =  0.000001),
    'tSNE_clustering_cohort': expand(output_dir_fig + 'tSNE_cohort-{alpha}.png', alpha =  0.000001),
    'tSNE_clustering_results': output_dir_results + "tsne_results.rds",
    'tSNE_unknown_cells': output_dir_fig + "tSNE_cellType_unknownCells.png",
    'tSNE_mTOR': expand(output_dir_fig + "tSNE_mTOR-{alpha}.png", alpha = 1),

    #'patient_analysis': expand(output_dir_reports + "Patient-level-analysis-{cohort}.html", cohort = reduced_cohort_list),
    'core_cell_identity': expand(output_dir_fig + "core_cell_identity_barchart_{cohort}.pdf", cohort = reduced_cohort_list),
    'core_cell_identity_allCohorts': output_dir_fig + "allCohort_core_cell_identity_barchart.pdf",
    'marker_correlation': expand(output_dir_fig + "expression_correlation_{cohort}.pdf", cohort = cohort_list),
    'unknownCells_score_hist': expand(output_dir_fig + "Unknown_celltype_mutual_exclusivity_score_hist_{cohort}.pdf", cohort = reduced_cohort_list),
    'unknownCells_score_qqplot': expand(output_dir_fig + "Unknown_celltype_mutual_exclusivity_score_qqplot_{cohort}.png", cohort = reduced_cohort_list),
    'probability_heatmap': expand(output_dir_fig + "celltype_probability_{cohort}.pdf", cohort = cohort_list),
    'astir_expression_heatmaps': expand(output_dir_fig + "Astir_expressionHeatmap_specifiedMarkers_{cohort}.pdf", cohort = cohort_list),
    'expression_heatmap': expand(output_dir_fig + "Unknown_celltype_expression_{cohort}.pdf", cohort = cohort_list), #created by celltypeProbability.R
    # 'EMT_expression_diagnostic': expand(output_dir_reports + "EMT_expression_{cohort}.html", cohort = ['basel', 'zurich1'])
    'Other-unknown-diff-expression': expand(output_dir_reports + "Unknown-other-differential-expression-{cohort}-{thresh}.html", cohort = cohort_list, thresh = [0.5]),
    'Other-unknown-diff-expression-pdf': expand(output_dir_fig + "Unknown-other-differential-expression-{cohort}-{thresh}.pdf", cohort = cohort_list, thresh = [0.5])
}

rule dim_reduct:
    input:
        basel_cells = output_path + "sces/basel_sce.rds",
        basel_types = output_path + "astir_assignments/basel_astir_assignments.csv",
        basel_states = output_path + "astir_assignments/basel_astir_assignments_state.csv",
        zurich1_cells = output_path + "sces/zurich1_sce.rds",
        zurich1_types = output_path + "astir_assignments/zurich1_astir_assignments.csv",
        zurich1_states = output_path + "astir_assignments/zurich1_astir_assignments_state.csv",
        wagner_cells = output_path + "sces/wagner_sce.rds",
        wagner_types = output_path + "astir_assignments/wagner_astir_assignments.csv",
        wagner_states = output_path + "astir_assignments/wagner_astir_assignments_state.csv",

    output:
        tsne_rds = output_dir_results + "tsne_results.rds"

    shell:
        "Rscript pipeline/analysis/dimReduct.R {input.basel_cells} {input.basel_types} {input.basel_states} "
        "{input.zurich1_cells} {input.zurich1_types} {input.zurich1_states} "
        "{input.wagner_cells} {input.wagner_types} {input.wagner_states} " + output_dir_results


rule dim_reduct_plot:
    input:
        tsne_rds = output_dir_results + "tsne_results.rds"

    output:
        cohort = output_dir_fig + "tSNE_cohort-{alpha}.png",
        cell_type = output_dir_fig + 'tSNE_cellType-{alpha}.png'

    shell:
        "Rscript pipeline/analysis/dimReductPlot.R {input.tsne_rds} {wildcards.alpha} " + output_dir_fig


rule dim_reduct_unknownCells:
    input:
        basel_cells = output_path + "sces/basel_sce.rds",
        basel_types = output_path + "astir_assignments/basel_astir_assignments.csv",
        basel_states = output_path + "astir_assignments/basel_astir_assignments_state.csv",
        zurich1_cells = output_path + "sces/zurich1_sce.rds",
        zurich1_types = output_path + "astir_assignments/zurich1_astir_assignments.csv",
        zurich1_states = output_path + "astir_assignments/zurich1_astir_assignments_state.csv",
        wagner_cells = output_path + "sces/wagner_sce.rds",
        wagner_types = output_path + "astir_assignments/wagner_astir_assignments.csv",
        wagner_states = output_path + "astir_assignments/wagner_astir_assignments_state.csv",

    output:
        cell_type = output_dir_fig + "tSNE_cellType_unknownCells.png"

    shell:
        "Rscript pipeline/analysis/dimReduct_unknownCells.R {input.basel_cells} {input.basel_types} {input.basel_states} "
        "{input.zurich1_cells} {input.zurich1_types} {input.zurich1_states} "
        "{input.wagner_cells} {input.wagner_types} {input.wagner_states} " + output_dir_fig


rule dim_reduct_mTOR:
    input:
        tsne_rds = output_dir_results + "tsne_results.rds",
        basel_states = output_path + "astir_assignments/basel_astir_assignments_state.csv",
        zurich1_states = output_path + "astir_assignments/zurich1_astir_assignments_state.csv",
        wagner_states = output_path + "astir_assignments/wagner_astir_assignments_state.csv",

    output:
        tsne = output_dir_fig + "tSNE_mTOR-{alpha}.png"
    
    resources:
        mem_mb=4000

    shell:
        "Rscript pipeline/analysis/mTOR_tSNE.R {input.tsne_rds} {input.basel_states} {input.zurich1_states} {input.wagner_states} "
        "{wildcards.alpha} " + output_dir_fig


rule patient_analysis:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        metadata = lambda wildcards: config[wildcards.cohort]['base_dir']+ config[wildcards.cohort]['metadata_file'],

    output:
        output_dir_reports + "Patient-level-analysis-{cohort}.html"

    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/analysis/Patient-level-analysis.Rmd', output_file = '{output}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', celltypes = '{input.cellTypes}', cellstates = '{input.cellStates}', "
        "metadata = '{input.metadata}', cohort = '{wildcards.cohort}'))\" "


rule cell_identity_by_core:
    input:
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv", 

    output:
        output_dir_fig + "core_cell_identity_barchart_{cohort}.pdf"

    shell:
        "Rscript pipeline/analysis/core_cell_identity_barchart.R {input.cellTypes} {input.cellStates} {wildcards.cohort} " + output_dir_fig


rule allCohort_core_barchart:
    input:
        basel = output_path + "astir_assignments/basel_astir_assignments.csv",
        schapiro = output_path + "astir_assignments/schapiro_astir_assignments.csv",
        wagner = output_path + "astir_assignments/wagner_astir_assignments.csv",
        zurich1 = output_path + "astir_assignments/zurich1_astir_assignments.csv",

    output:
        barchart = output_dir_fig + "allCohort_core_cell_identity_barchart.pdf"

    shell:
        "Rscript pipeline/analysis/core_cell_identity_barchart_allCohorts.R {input.basel} {input.schapiro} {input.wagner} "
        "{input.zurich1} " + output_dir_fig


rule astir_expression_heatmaps:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        #cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv", 
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],
    
    output:
        allMarkers = output_dir_fig + "Astir_expressionHeatmap_allMarkers_{cohort}.pdf",
        specifiedMarkers = output_dir_fig + "Astir_expressionHeatmap_specifiedMarkers_{cohort}.pdf"

    resources:
        mem_mb=2000

    shell:
        "Rscript pipeline/analysis/Astir_expressionHeatmaps.R {input.cells} {input.cellTypes} "# {input.cellStates} "
        "{input.markers} {wildcards.cohort} " + output_dir_fig


### Unknown celltype analysis
rule markers_correlation:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        #cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv", 
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],

    output:
        output_dir_fig + "expression_correlation_{cohort}.pdf"

    shell:
        "Rscript pipeline/analysis/marker_gene_correlations.R {input.cells} {input.cellTypes} "#{input.cellStates} "
        "{input.markers} {wildcards.cohort} " + output_dir_fig

rule EMT_expression:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv",
        metadata = lambda wildcards: config[wildcards.cohort]['base_dir']+ config[wildcards.cohort]['metadata_file'],
        markers = lambda wildcards: config[wildcards.cohort]['marker_file']

    output:
        output_dir_reports + "EMT_expression_{cohort}.html"

    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/analysis/EMT.Rmd', output_file = '{output}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', types = '{input.cellTypes}', states = '{input.cellStates}', "
        "metadata = '{input.metadata}', markers = '{input.markers}', cohort = '{wildcards.cohort}'))\" "



rule unknown_celltype_score:
    input:
        cells =  output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv", 
    
    output:
        hist = output_dir_fig + "Unknown_celltype_mutual_exclusivity_score_hist_{cohort}.pdf",
        qqplot = output_dir_fig + "Unknown_celltype_mutual_exclusivity_score_qqplot_{cohort}.png"
    
    shell:
        "Rscript pipeline/analysis/Unknown_celltype_score.R {input.cells} {input.cellTypes} "
        "{wildcards.cohort} " + output_dir_fig


rule cell_type_probability_heatmap:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv", 
        markers = lambda wildcards: config[wildcards.cohort]['marker_file'],

    output:
        probability_heatmap = output_dir_fig + "celltype_probability_{cohort}.pdf",
        expression_heatmap = output_dir_fig + "Unknown_celltype_expression_{cohort}.pdf",

    shell:
        "Rscript pipeline/analysis/celltypeProbability.R {input.cells} {input.cellTypes} "
        "{input.markers} {wildcards.cohort} " + output_dir_fig


rule other_unknown_expression:
    input:
        cells = output_path + "sces/{cohort}_sce.rds",
        cellTypes = output_path + "astir_assignments/{cohort}_astir_assignments.csv",
        #cellStates = output_path + "astir_assignments/{cohort}_astir_assignments_state.csv"

    output:
        html = output_dir_reports + "Unknown-other-differential-expression-{cohort}-{thresh}.html",
        pdf = output_dir_fig + "Unknown-other-differential-expression-{cohort}-{thresh}.pdf"

    resources:
        mem_mb=2000

    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/analysis/Unknown-other-differential-expression.Rmd', output_file = '{output.html}', output_dir = '" + output_dir_reports + "', "
        "params = list(cells = '{input.cells}', types = '{input.cellTypes}', cohort = '{wildcards.cohort}', thresh = '{wildcards.thresh}', output_dir = '" + output_dir_fig +"'))\" "