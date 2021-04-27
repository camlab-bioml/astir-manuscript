
import pandas as pd
import numpy as np
import os
from anndata import AnnData


configfile: "imc-config.yml"
logdir = "logs/slurm"

output_path = "output/" + config['version'] + "/"

# Final assignments
datasets = ['basel', 'zurich1', 'wagner', 'schapiro']
asts_type = {d: output_path + f"astir_assignments/{d}_astir_assignments.csv" for d in datasets + ["lin_cycif"]}
asts_state = {d: output_path + f"astir_assignments/{d}_astir_assignments_state.csv" for d in datasets}

# Lin-cycif tissue type assignments
lin_metadata = pd.read_csv("data-raw/metadata/lin-cycif-metadata.tsv", sep = "\t")
cancer_types = lin_metadata["Anatomic site"].unique()
lin_tissue_types_anndata = expand(output_path + "anndata/lin_cycif-{anatomic_site}.h5ad", anatomic_site = cancer_types)
lin_tissue_types_asts_type = {cancer_types: output_path + f"lin_tissue_types_astir_assignments/lin_cycif_{cancer_types}_astir_assignments.csv" for cancer_types in cancer_types}


def create_annData_from_csv_lin(metadata_file, anatomic_site, input_file_path, output_file):
    # creates anndata files for lin_cycif stratified by tissue type
    metadata = pd.read_csv(metadata_file, sep = '\t')
    
    cores = list(metadata[metadata['Anatomic site'] == anatomic_site]['core'])
    
    csv_files = [f + '.csv' for f in cores]
    
    dfs = [pd.read_csv(os.path.join(input_file_path, f), index_col = 0) for f in csv_files]
    df = pd.concat(dfs, axis = 0)
    
    batch = list(np.concatenate([np.repeat(csv_files[i].replace(".csv", ""), dfs[i].shape[0]) for i in range(len(dfs))], axis=0))
    
    obs = pd.DataFrame({'cell_id': list(df.index), 'batch': batch})
    adata = AnnData(df.values, obs)
    adata.var_names = list(df.columns)
    adata.obs_names = list(df.index)
    
    adata.write(output_file)


## Dataset-specific pipelines
include: "pipeline/dataset-wagner.smk"
include: "pipeline/dataset-basel.smk"
include: "pipeline/dataset-zurich1.smk"
include: "pipeline/dataset-schapiro.smk"
include: "pipeline/dataset-lin-cycif.smk"
# include: "pipeline/dataset-keren.smk"

include: "analysis/analysis.smk"
include: "pipeline/benchmarking/benchmarking.smk"
include: "pipeline/Alternate-clustering.smk"
include: "pipeline/epithelial-overclustering.smk"    
include: "pipeline/alt-masks-clustering.smk" 
include: "pipeline/lin-cycif-breakdown.smk"

include: "pipeline/robustness/robustness.smk"
# include: "pipeline/reports/reports.smk"
include: "pipeline/spatial/spatial.smk"
include: "pipeline/imbalance/imbalance.smk"
include: "pipeline/segmentation/segmentation.smk"

include: "pipeline/cla/cla.smk"
#include: "pipeline/granularity/granularity.smk"


## Beginning of rules ----- 
rule all:
    input:        
        # asts_type.values(),
        # wagner_output.values(),
        # basel_output.values(),
        # zurich1_output.values(),
        # schapiro_output.values(),
        # keren_output.values(),
        # reports_output.values(),
        # spatial_output.values(),
        # lin_cycif_output.values(),
        lin_cycif_breakdown.values(),
        segmentation_output.values(),
        # imbalance_output.values(),
        # robustness_output.values(),
        # cla_outputs.values(),
        # benchmarking_output.values(),
        alternate_approaches_output.values(),
        # analysis_output.values(),

        #epithelial_overclustering.values(),
        alt_masks.values(),
        lin_tissue_types_anndata,
        lin_tissue_types_asts_type.values(),
        #granularity_outputs.values(),

rule create_lin_anndata_objects:
    input:
        metadata = 'data-raw/metadata/lin-cycif-metadata.tsv'
    params:
        input_csvs = output_path + 'lin-cycif_processed/'
    #container: 'astir-manuscript.sif'
    output:
        anndata = output_path + 'anndata/lin_cycif-{anatomic_site}.h5ad'
    run:
        create_annData_from_csv_lin(input.metadata, wildcards.anatomic_site, params.input_csvs, output.anndata)

# rule run_astir_type:
#     params:
#         op = output_path,
#         max_epochs = config['astir_opts']['max_epochs'],
#         batch_size = config['astir_opts']['batch_size'],
#         learning_rate = config['astir_opts']['learning_rate'],
#         n_initial_epochs = config['astir_opts']['n_initial_epochs']
#     input:
#         loom=ancient(output_path + "looms/{dataset}.loom"),
#         markers=ancient(lambda wildcards: config[wildcards.dataset]['marker_file']),
#     output:
#         csv=output_path + "astir_assignments/{dataset}_astir_assignments.csv",
#         fig=output_path + "astir_assignments/{dataset}_astir_loss.png",
#         diagnostics=output_path + "astir_assignments/{dataset}_diagnostics.csv"
#     run:
#         # fit
#         from astir.data import from_loompy_yaml
#         from datetime import datetime
#         print(f"{datetime.now()}\t Reading loom file ")
#         ast = from_loompy_yaml(input.loom, input.markers)
#         print(f"{datetime.now()}\t Fitting model")

#         N = ast.get_type_dataset().get_exprs_df().shape[0]
#         batch_size = int(N/100)


#         ast.fit_type(max_epochs = int(params.max_epochs), 
#         batch_size = batch_size, 
#         learning_rate = float(params.learning_rate),
#         n_init_epochs=int(params.n_initial_epochs))
#         print(f"{datetime.now()}\t Finished fitting model")
#         ast.get_celltype_probabilities().to_csv(output.csv)

#         # run diagnostics
#         ast.diagnostics_celltype().to_csv(output.diagnostics)

#         # plot loss
#         import matplotlib.pyplot as plt
#         plt.figure(figsize=(5,4))
#         plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
#         plt.ylabel("Loss")
#         plt.xlabel("Epoch")
#         plt.tight_layout()
#         plt.savefig(output.fig, dpi=300)

rule run_astir_type_lin_tissue_types:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        anndata=ancient(output_path + "anndata/lin_cycif-{cancer_types}.h5ad"),
        markers=ancient(lambda wildcards: config['lin_cycif']['marker_file']),
    output:
        csv=output_path + "lin_tissue_types_astir_assignments/lin_cycif_{cancer_types}_astir_assignments.csv",
        fig=output_path + "lin_tissue_types_astir_assignments/lin_cycif_{cancer_types}_astir_loss.png",
        diagnostics=output_path + "lin_tissue_types_astir_assignments/lin_cycif_{cancer_types}_diagnostics.csv"
    run:
         # fit
        from astir.data import from_anndata_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_anndata_yaml(input.anndata, input.markers)
        print(f"{datetime.now()}\t Fitting model")

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)

        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate),
        n_init_epochs=int(params.n_initial_epochs))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.get_celltype_probabilities().to_csv(output.csv)

        # run diagnostics
        ast.diagnostics_celltype().to_csv(output.diagnostics)

        # plot loss
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)
        


# rule run_astir_state:
#     params:
#         op = output_path,
#         max_epochs = config['astir_opts']['max_epochs'],
#         batch_size = config['astir_opts']['batch_size'],
#         learning_rate = config['astir_opts']['learning_rate'],
#     input:
#         loom=ancient(output_path + "looms/{dataset}.loom"),
#         markers=lambda wildcards: config[wildcards.dataset]['marker_file'],
#         type_diagnostic=output_path + "astir_assignments/{dataset}_diagnostics.csv", # make sure state doesn't run concurrently with type
#     output:
#         csv=output_path + "astir_assignments/{dataset}_astir_assignments_state.csv",
#         fig=output_path + "astir_assignments/{dataset}_astir_loss_state.png",
#         # diagnostics=output_path + "astir_assignments/{dataset}_diagnostics.csv"
#     run:
#         # fit
#         from astir.data import from_loompy_yaml
#         from datetime import datetime
#         print(f"{datetime.now()}\t Reading loom file ")
#         ast = from_loompy_yaml(input.loom, input.markers)
#         print(f"{datetime.now()}\t Fitting model")
#         ast.fit_state(max_epochs = int(params.max_epochs), 
#         batch_size = int(params.batch_size), 
#         learning_rate = float(params.learning_rate))
#         print(f"{datetime.now()}\t Finished fitting model")
#         ast.state_to_csv(output.csv)

#         # run diagnostics
#         # ast.diagnostics_celltype().to_csv(output.diagnostics)

#         # plot loss
#         import matplotlib.pyplot as plt
#         plt.figure(figsize=(5,4))
#         plt.plot(np.arange(len(ast.get_state_losses())), ast.get_state_losses())
#         plt.ylabel("Loss")
#         plt.xlabel("Epoch")
#         plt.tight_layout()
#         plt.savefig(output.fig, dpi=300)




