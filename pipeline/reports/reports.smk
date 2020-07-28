
import os

reports_output = {
    'expression_html': expand(output_path + "reports/expression_heatmap_{dataset}.html", dataset=datasets),
    'expression_figs': expand(output_path + "figures/supplementary/expression_heatmap_{dataset}.pdf", dataset=datasets)
}

rule subset_expression_reports:
    params:
        curr_dir = os.getcwd(),
    input:
        subset_exprs=output_path + "{dataset}_subset/{dataset}_subset_sce.rds",
        subset_assignments=output_path + "{dataset}_subset/{dataset}_subset_assignments_type.csv",
        marker_yml = lambda wildcards: config[wildcards.dataset]['marker_file'],
    output:
        html=output_path + "reports/expression_heatmap_{dataset}.html",
        figure=output_path + "figures/supplementary/expression_heatmap_{dataset}.pdf",
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin'); "
        "rmarkdown::render('pipeline/reports/subset-expression-figures.Rmd', "
        "output_file='{output.html}', "
        "params=list(marker_yml='{input.marker_yml}', "
        "subset_exprs='{input.subset_exprs}', "
        "expression_pdf='{output.figure}', "
        "subset_assignments='{input.subset_assignments}'), "
        "knit_root_dir='{params.curr_dir}', "
        "output_dir='{params.curr_dir}/{output_path}/reports/'"
        ")\"  "
