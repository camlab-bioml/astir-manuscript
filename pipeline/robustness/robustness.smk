


## Snakefile for assessing the robustness of astir to undesirable implementations. Currently:
## (1)  create_reduced_sets: remove cell types from described markers and assess how well 
##      inference is performed
## (2)  

import yaml
import numpy as np


## Removal of cell types
celltypes_to_remove = {
    'Stromal_only': ['Stromal'],
    'Stromal_Macrophage': ['Stromal', 'Macrophage'],
    'Stromal_Macrophage_Endothelial': ['Stromal', 'Macrophage', 'Endothelial'],
    'None': []
}

reduced_assignments = expand(output_path + "robustness/assignments-15k-removed-{removed}.csv", 
    removed=celltypes_to_remove.keys())

## Addition of cell types
fake_celltypes = {
    'CD3+Keratin+': [
        'pan Cytokeratin',
        'CD3'
    ],
    'CD20+CDH1+': [
        'CD20',
        'E-Cadherin'
    ],
    'CK5+VIM+': [
        'Cytokeratin 5',
        'Vimentin'
    ]
}

celltypes_to_add = {
    'One_additional': list(fake_celltypes.keys())[0:1],
    'Two_additional': list(fake_celltypes.keys())[0:2],
    'Three_additional': list(fake_celltypes.keys())[0:3]
}

added_assignments = expand(output_path + "robustness/assignments-15k-added-{added}.csv",
    added=celltypes_to_add.keys())

robustness_output = {
    'added_assignments': added_assignments,
    'reduced_assignments': reduced_assignments,
    'report': output_path + "robustness/robustness.html"
}

rule robustness_figures:
    params:
        version=config['version'],
        curr_dir = os.getcwd(),
        output_path = output_path,
    input:
        robustness_output['added_assignments'],
        robustness_output['reduced_assignments'],
    output:
        html=output_path + "robustness/robustness.html",
        remove_fig=output_path + "robustness/robustness_removed.pdf",
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin'); "
        "rmarkdown::render('pipeline/robustness/robustness-figures.Rmd',   "
        "output_file='{output.html}', "
        "params=list(version='{params.version}', "
        "input_dir='{params.output_path}robustness', "
        "output_fig_removed='{output.remove_fig}'), "
        "knit_root_dir='{params.curr_dir}', "
        "output_dir='{params.curr_dir}/{params.output_path}/robustness/'"
        ")\"  "

rule remove_celltypes:
    params:
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
    input:
        markers="markers/jackson-2020-markers.yml",
        csvs=basel_output['subset_csvs'],
    output:
        markers=output_path + "robustness/{removed}.yml",
        assignments = output_path + "robustness/assignments-15k-removed-{removed}.csv",
        fig=output_path + "robustness/{removed}_loss.png" ,       diagnostics=output_path + "robustness/{removed}_diagnostics.csv"
    run:

        marker_dict = None
        with open(input.markers, "r") as stream:
            marker_dict = yaml.safe_load(stream)
            for t in celltypes_to_remove[wildcards.removed]:
                marker_dict['cell_types'].pop(t)
        
        with open(output.markers, 'w') as stream:
            yaml.dump(marker_dict, stream)
            
        from astir.data import from_csv_dir_yaml

        ast = from_csv_dir_yaml(os.path.join(output_path, "basel_subset_separate_csvs"), output.markers)
        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = int(params.batch_size), 
        learning_rate = float(params.learning_rate))
        ast.type_to_csv(output.assignments)

        # run diagnostics
        ast.diagnostics_celltype().to_csv(output.diagnostics)

        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)


rule add_celltypes:
    params:
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
    input:
        markers="markers/jackson-2020-markers.yml",
        csvs=basel_output['subset_csvs'],
    output:
        markers=output_path + "robustness/{added}.yml",
        assignments = output_path + "robustness/assignments-15k-added-{added}.csv",
        fig=output_path + "robustness/{added}_loss.png",       diagnostics=output_path + "robustness/{added}_diagnostics.csv"
    run:
        marker_dict = None
        with open(input.markers, "r") as stream:
            marker_dict = yaml.safe_load(stream)
            for t in celltypes_to_add[wildcards.added]:
                marker_dict['cell_types'][t] = fake_celltypes[t]
        
        with open(output.markers, 'w') as stream:
            yaml.dump(marker_dict, stream)
            
        from astir.data import from_csv_dir_yaml

        ast = from_csv_dir_yaml(os.path.join(output_path, "basel_subset_separate_csvs"), output.markers)
        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = int(params.batch_size), 
        learning_rate = float(params.learning_rate))
        ast.type_to_csv(output.assignments)

        # run diagnostics
        ast.diagnostics_celltype().to_csv(output.diagnostics)

        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)
