
import pandas as pd
import numpy as np
import os


configfile: "imc-config.yml"
logdir = "logs/slurm"

output_path = "output/" + config['version'] + "/"

# Final assignments
datasets = ['basel', 'zurich1', 'wagner']
asts_type = {d: output_path + f"astir_assignments/{d}_astir_assignments.csv" for d in datasets}
asts_state = {d: output_path + f"astir_assignments/{d}_astir_assignments_state.csv" for d in datasets}


## Dataset-specific pipelines
include: "pipeline/dataset-wagner.smk"
include: "pipeline/dataset-basel.smk"
include: "pipeline/dataset-zurich1.smk"

include: "analysis/analysis.smk"

# include: "pipeline/benchmarking/benchmarking.smk"
include: "pipeline/robustness/robustness.smk"
include: "pipeline/reports/reports.smk"
include: "pipeline/benchmarking/benchmarking.smk"

## Beginning of rules ----- 
rule all:
    input:
        asts_type.values(),
        # asts_state.values(),
        robustness_output.values(),
        analysis_output.values(),
        wagner_output.values(),
        basel_output.values(),
        zurich1_output.values(),
        # benchmarking_output.values(),
        # reports_output.values(),



rule run_astir_type:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        loom=ancient(output_path + "looms/{dataset}.loom"),
        markers=ancient(lambda wildcards: config[wildcards.dataset]['marker_file']),
    output:
        csv=output_path + "astir_assignments/{dataset}_astir_assignments.csv",
        fig=output_path + "astir_assignments/{dataset}_astir_loss.png",
        diagnostics=output_path + "astir_assignments/{dataset}_diagnostics.csv"
    run:
        # fit
        from astir.data_readers import from_loompy_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(input.loom, input.markers)
        print(f"{datetime.now()}\t Fitting model")
        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = int(params.batch_size), 
        learning_rate = float(params.learning_rate),
        n_init_epochs=int(params.n_initial_epochs))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.type_to_csv(output.csv)

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


rule run_astir_state:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
    input:
        loom=ancient(output_path + "looms/{dataset}.loom"),
        markers=lambda wildcards: config[wildcards.dataset]['marker_file'],
        type_diagnostic=output_path + "astir_assignments/{dataset}_diagnostics.csv", # make sure state doesn't run concurrently with type
    output:
        csv=output_path + "astir_assignments/{dataset}_astir_assignments_state.csv",
        fig=output_path + "astir_assignments/{dataset}_astir_loss_state.png",
        # diagnostics=output_path + "astir_assignments/{dataset}_diagnostics.csv"
    run:
        # fit
        from astir.data_readers import from_loompy_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(input.loom, input.markers)
        print(f"{datetime.now()}\t Fitting model")
        ast.fit_state(max_epochs = int(params.max_epochs), 
        batch_size = int(params.batch_size), 
        learning_rate = float(params.learning_rate))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.state_to_csv(output.csv)

        # run diagnostics
        # ast.diagnostics_celltype().to_csv(output.diagnostics)

        # plot loss
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_state_losses())), ast.get_state_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)



