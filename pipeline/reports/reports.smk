
import os

reports_datasets = datasets + ["keren"]

reports_output = {
    'expression_html': expand(output_path + "reports/expression_heatmap_{dataset}.html", dataset=reports_datasets),
    'expression_figs': expand(output_path + "figures/supplementary/expression_heatmap_{dataset}.pdf", dataset=reports_datasets)
}

rule expression_reports:
    params:
        curr_dir = os.getcwd(),
        input_rds_dir=output_path + "{dataset}_processed/"
    input:
        assignments=output_path + "astir_assignments/{dataset}_astir_assignments.csv",
        marker_yml = lambda wildcards: config[wildcards.dataset]['marker_file'],
    output:
        html=output_path + "reports/expression_heatmap_{dataset}.html",
        figure=output_path + "figures/supplementary/expression_heatmap_{dataset}.pdf",
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin'); "
        "rmarkdown::render('pipeline/reports/subset-expression-figures.Rmd', "
        "output_file='{output.html}', "
        "params=list(marker_yml='{input.marker_yml}', "
        "input_rds_dir='{params.input_rds_dir}', "
        "expression_pdf='{output.figure}', "
        "assignments='{input.assignments}'), "
        "knit_root_dir='{params.curr_dir}', "
        "output_dir='{params.curr_dir}/{output_path}/reports/'"
        ")\"  "
