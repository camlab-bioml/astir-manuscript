
import pandas as pd
import os


configfile: "imc-config.yml"
logdir = "logs/slurm"

output_path = "output/" + config['version'] + "/"

# Final assignments
datasets = ['basel', 'zurich1', 'wagner']
asts = {d: output_path + f"astir_assignments/{d}_astir_assignments.csv" for d in datasets}



## Dataset-specific pipelines
include: "pipeline/wagner.smk"
include: "pipeline/basel.smk"
include: "pipeline/zurich1.smk"

include: "analysis/analysis.smk"

# include: "pipeline/benchmarking/benchmarking.smk"
include: "pipeline/robustness/robustness.smk"

## Beginning of rules ----- 
rule all:
    input:
        asts.values(),
        robustness_output.values(),
        analysis_output.values(),
        wagner_output.values(),
        basel_output.values(),
        zurich1_output.values()



rule run_astir:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
    input:
        loom=ancient(output_path + "looms/{dataset}.loom"),
        markers=lambda wildcards: config[wildcards.dataset]['marker_file'],
    output:
        csv=output_path + "astir_assignments/{dataset}_astir_assignments.csv"
    run:
        from astir.data_readers import from_loompy_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(input.loom, input.markers, include_beta=False)
        print(f"{datetime.now()}\t Fitting model")
        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = int(params.batch_size), 
        learning_rate = float(params.learning_rate))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.type_to_csv(output.csv)

