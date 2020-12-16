


## Snakefile for assessing the robustness of astir to undesirable implementations. Currently:
## (1)  create_reduced_sets: remove cell types from described markers and assess how well 
##      inference is performed
## (2)  

import yaml
import numpy as np
import string


def remove_punc(s):
    """Remove punctuation from a string"""
    s = s.replace(" ", "")
    return s.translate(str.maketrans('', '', string.punctuation))

cohorts = ['basel', 'zurich1', 'wagner', 'schapiro', 'keren']

celltypes_to_remove = {}

for cohort in cohorts:
    markers = None
    with open(config[cohort]['marker_file'], "r") as stream:
        markers = yaml.safe_load(stream)
    cell_types = list(markers['cell_types'].keys())
    ctr = {remove_punc(s): s for s in cell_types}
    ctr['None'] = []
    celltypes_to_remove[cohort] = ctr



reduced_assignments = {cohort: \
    expand(output_path + f"robustness/removed/{cohort}/{cohort}-assignments-{{removed}}.csv", 
    removed=celltypes_to_remove[cohort].keys()) for cohort in cohorts}

# print(list(config['fake_celltypes']['basel'].keys()))

added_assignments = {
    cohort: \
    expand(output_path + f"robustness/added/{cohort}/{cohort}-assignments-{{fake_celltype}}.csv", 
    fake_celltype=list(config['fake_celltypes'][cohort].keys())) \
    for cohort in cohorts
}



robustness_output = {
    'reduced_assignments': reduced_assignments.values(),
    'removed_dfs': expand(output_path + "robustness/removed/df_{cohort}.tsv", cohort=cohorts),
    'added_assignments': added_assignments.values(),
    'added_dfs': expand(output_path + "robustness/added/df_{cohort}.tsv", cohort=cohorts),
    'figure_added': output_path + "figures/robustness/added.pdf",
    'figure_removed': output_path + "figures/robustness/removed.pdf"
}



rule remove_celltypes:
    params:
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        input_dir = lambda wildcards: output_path + f"{wildcards.cohort}_processed"
    input:
        markers=lambda wildcards: config[wildcards.cohort]['marker_file'],
        csvs=asts_type.values(), # this just ensures we have necessary csvs
    output:
        markers=output_path + "robustness/removed/{cohort}/{removed}.yml",
        assignments = output_path + "robustness/removed/{cohort}/{cohort}-assignments-{removed}.csv",
        fig=output_path + "robustness/removed/{cohort}/{removed}_loss.png" ,       
        diagnostics=output_path + "robustness/removed/{cohort}/{removed}_diagnostics.csv"
    run:

        marker_dict = None
        with open(input.markers, "r") as stream:
            marker_dict = yaml.safe_load(stream)


        t = celltypes_to_remove[wildcards.cohort][wildcards.removed]
        if str(wildcards.removed) != "None":
            marker_dict['cell_types'].pop(t)

        
        with open(output.markers, 'w') as stream:
            yaml.dump(marker_dict, stream)
            
        from astir.data import from_csv_dir_yaml

        ast = from_csv_dir_yaml(params.input_dir, output.markers)

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)
        # batch_size = int(params.batch_size)

        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate))
        
        ast.get_celltype_probabilities().to_csv(output.assignments)

        ast.get_type_dataset().get_exprs_df().to_csv(
            os.path.join(
                output_path, "robustness",
                "tmp_exprs_" + wildcards.removed + ".csv"
            )
        )

        # run diagnostics
        ast.diagnostics_celltype().to_csv(output.diagnostics)

        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)

rule parse_removed_output_to_df:
    params:
        input_dir = output_path + "robustness/removed/{cohort}"
    input:
        markers=lambda wildcards: config[wildcards.cohort]['marker_file'],
        ass=lambda wildcards: reduced_assignments[wildcards.cohort],
    output:
        output_path + "robustness/removed/df_{cohort}.tsv",
    shell:
        "Rscript pipeline/robustness/parse-robustness-to-df.R "
        "--marker_yml {input.markers} "
        "--input_dir {params.input_dir} "
        "--cohort {wildcards.cohort} "
        "--output_file {output}"


rule add_celltypes:
    params:
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        input_dir = lambda wildcards: output_path + f"{wildcards.cohort}_processed"
    input:
        markers=lambda wildcards: config[wildcards.cohort]['marker_file'],
        csvs=asts_type.values(), # this just ensures we have necessary csvs
    output:
        markers=output_path + "robustness/added/{cohort}/{fake_celltype}.yml",
        assignments = output_path + "robustness/added/{cohort}/{cohort}-assignments-{fake_celltype}.csv",
        fig=output_path + "robustness/added/{cohort}/{fake_celltype}_loss.png" ,       
        diagnostics=output_path + "robustness/added/{cohort}/{fake_celltype}_diagnostics.csv"
    run:
        marker_dict = None
        with open(input.markers, "r") as stream:
            marker_dict = yaml.safe_load(stream)

        ct = wildcards.fake_celltype
        t = config['fake_celltypes'][wildcards.cohort][ct]

        marker_dict['cell_types'][ct] = t
        
        with open(output.markers, 'w') as stream:
            yaml.dump(marker_dict, stream)
            
        from astir.data import from_csv_dir_yaml

        ast = from_csv_dir_yaml(params.input_dir, output.markers)

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)
        # batch_size = int(params.batch_size)

        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate))


        ast.get_celltype_probabilities().to_csv(output.assignments)

        # run diagnostics
        ast.diagnostics_celltype().to_csv(output.diagnostics)

        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)

rule parse_added_output_to_df:
    params:
        input_dir = output_path + "robustness/added/{cohort}"
    input:
        markers=lambda wildcards: config[wildcards.cohort]['marker_file'],
        ass=lambda wildcards: added_assignments[wildcards.cohort],
    output:
        output_path + "robustness/added/df_{cohort}.tsv",
    shell:
        "Rscript pipeline/robustness/parse-fake-to-df.R "
        "--input_dir {params.input_dir} "
        "--cohort {wildcards.cohort} "
        "--output_file {output}"

rule make_figs:
    params:
        version=config['version'],
        curr_dir = os.getcwd(),
        input_dir_added = output_path + "robustness/added",
        input_dir_removed = output_path + "robustness/removed",
        output_path = output_path,
    input:
        robustness_output['added_dfs'],
        robustness_output['removed_dfs'],
    output:
        html=output_path + "robustness/robustness-figs.html",
        output_fig_removed=robustness_output['figure_removed'],
        output_fig_added=robustness_output['figure_added']
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin'); "
        "rmarkdown::render('pipeline/robustness/new-figs.Rmd',   "
        "output_file='{output.html}', "
        "params=list(input_dir_added='{params.input_dir_added}', "
        "input_dir_removed='{params.input_dir_removed}', "
        "output_fig_added='{output.output_fig_added}', "
        "output_fig_removed='{output.output_fig_removed}'), "
        "knit_root_dir='{params.curr_dir}', "
        "output_dir='{params.curr_dir}/{params.output_path}/robustness/'"
        ")\"  "