
import pandas as pd
import os


configfile: "imc-config.yml"
logdir = "logs/slurm"

output_path = "output/" + config['version'] + "/"

# Final assignments
datasets = ['basel', 'zurich1']
asts = {d: output_path + f"astir_assignments/{d}_astir_assignments.csv" for d in datasets}


# include: "pipeline/benchmarking/benchmarking.smk"
include: "pipeline/robustness/robustness.smk"

## Dataset-specific pipelines
include: "pipeline/wagner.smk"
include: "pipeline/basel.smk"
include: "pipeline/zurich1.smk"

include: "analysis/analysis.smk"



## Beginning of rules ----- 
rule all:
    input:
        # tmp_zurich1_output,
        asts.values(),
        # expand(output_path + "looms/{dataset}.loom", dataset=datasets),
        # output_path + "summarized_assigments/basel_15k_all.csv",
        # reduced_assignments, added_assignments,
        zurch1_subset_cell_ids,
        analysis_deliverables,
        # geneset_files
        wagner_output.values()
        basel_output.values()



rule astir:
    params:
        op = output_path
    input:
        loom=ancient(output_path + "looms/{dataset}.loom"),
        markers="markers/jackson-2020-markers.yml"
    output:
        csv=output_path + "astir_assignments/{dataset}_astir_assignments.csv"
    run:
        from astir.data_readers import from_loompy_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(input.loom, input.markers, include_beta=False)
        print(f"{datetime.now()}\t Fitting model")
        ast.fit_type(max_epochs = 20, batch_size = 512, learning_rate = 1e-3)
        print(f"{datetime.now()}\t Finished fitting model")
        ast.type_to_csv(output.csv)

