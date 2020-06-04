
import pandas as pd
import os


configfile: "imc-config.yml"
logdir = "logs/slurm"

output_path = "output/" + config['version'] + "/"

# Jackson 2020 metadata
basel_metadata = pd.read_csv(os.path.join(config['basel']['base_dir'], config['basel']['metadata_file']))
basel_cores = list(basel_metadata.core)
tmp_basel_output = expand(output_path + "basel_processed/{core}.rds", core=basel_cores)

zurich1_metadata = pd.read_csv(os.path.join(config['zurich1']['base_dir'], config['zurich1']['metadata_file']))
zurich1_cores = list(zurich1_metadata.core)
tmp_zurich1_output = expand(output_path + "zurich1_processed/{core}.rds", core=zurich1_cores)
zurch1_subset_cell_ids = cell_ids=output_path + "zurich1_subset/zurich1_subset_cells.csv"

# Wagner 2019 metadata
wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])


tmp_wagner_output = expand(output_path + 
        "wagner-2019_processed/{sample}.{ext}",
        sample=wagner_samples, ext=['csv','rds'])

# Basel 15k reduced set for benchmarking
basel_15k_cells = pd.read_csv("data-raw/cellsets/cell_ids_15k.csv")
cores_15k = list(set(list(basel_15k_cells.core)))
basel_15k_csvs = expand(output_path + "basel_15k_subset/{core}.csv", core=cores_15k)

# Final assignments
datasets = ['basel', 'zurich1']
asts = {d: output_path + f"astir_assignments/{d}_astir_assignments.csv" for d in datasets}


# include: "pipeline/benchmarking/benchmarking.smk"
include: "pipeline/robustness/robustness.smk"
include: "analysis/analysis.smk"

## Beginning of rules ----- 
rule all:
    input:
        # tmp_basel_output,
        # tmp_zurich1_output,
        asts.values(),
        # expand(output_path + "looms/{dataset}.loom", dataset=datasets),
        # output_path + "summarized_assigments/basel_15k_all.csv",
        # reduced_assignments, added_assignments,
        zurch1_subset_cell_ids,
        analysis_deliverables,
        # geneset_files
        # tmp_wagner_output

rule read_wagner_2019:
    input:
        config['wagner']['fcs_dir'] + "/{sample}.fcs"
    output:
        rds = output_path + "wagner-2019_processed/{sample}.rds",
        csv = output_path + "wagner-2019_processed/{sample}.csv",
    shell:
        "Rscript pipeline/wagner-2019/wagner-fcs-to-csv.R "
        "--input_fcs {input} "
        "--output_rds {output.rds} "
        "--output_csv {output.csv} "


rule read_jackson_2020_basel:
    input:
        scdat=config['basel']['base_dir'] + "SC_dat.csv",
        scloc=config['basel']['base_dir'] + "Basel_SC_locations.csv",
    output:
        expand(output_path + "basel_processed/{core}.rds", core=basel_cores),
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
    shell:
        "Rscript pipeline/jackson-2020/jackson-raw-to-sce-server.R "
        "--input_sc {input.scdat} "
        "--input_loc {input.scloc} "
        "--output_dir {output_path}/basel_processed "

rule read_jackson_2020_zurich1:
    input:
        scdat=config['zurich1']['base_dir'] + "SC_dat.csv",
        scloc=config['zurich1']['base_dir'] + "Zurich_SC_locations.csv",
    output:
        expand(output_path + "zurich1_processed/{core}.rds", core=zurich1_cores),
        expand(output_path + "zurich1_processed/{core}.csv", core=zurich1_cores),
    shell:
        "Rscript pipeline/jackson-2020/jackson-raw-to-sce-server.R "
        "--input_sc {input.scdat} "
        "--input_loc {input.scloc} "
        "--output_dir {output_path}/zurich1_processed "

rule basel_to_loom:
    input:
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
    output:
        output_path + "looms/basel.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/basel_processed "
        "{output} "

rule zurich1_to_loom:
    input:
        expand(output_path + "zurich1_processed/{core}.csv", core=zurich1_cores),
    output:
        output_path + "looms/zurich1.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/zurich1_processed "
        "{output} "

rule subset_basel_cells:
    input:
        expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
        cell_labels="data-raw/cellsets/cell_ids_15k.csv"
    output:
        basel_15k_csvs
    shell:
        "mkdir -p {output} && "
        "python pipeline/subset-expression-data.py "
        "--input_dir {output_path}/basel_processed "
        "--cell_labels {input.cell_labels} "
        "--output_dir {output_path}/basel_15k_subset "

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

rule basel_to_15k:
    input:
        astir_output=output_path + "astir_assignments/basel_astir_assignments.csv",
        cell_list="data-raw/15k_cells.csv",
    output:
        output_path + "summarized_assigments/basel_15k_all.csv"
    shell:
        "python pipeline/subset-to-15k.py "
        "{input.astir_output} "
        "{input.cell_list} "
        "{output}"


rule subset_zurich1:
    input:
        loom=output_path + "looms/zurich1.loom",
        sc_dat="data-raw/jackson-2020/SingleCell_and_Metadata/ZurichTMA/" + "SC_dat.csv",
        assignments=output_path + "astir_assignments/zurich1_astir_assignments.csv",
    output:
        cell_ids=output_path + "zurich1_subset/zurich1_subset_cells.csv",
        assignments=output_path + "zurich1_subset/zurich1_subset_assignments.csv",
    run:
        import loompy
        sc_dat = pd.read_csv(input.sc_dat)

        ## Cell IDs we're going to subsample
        cell_ids = list(sc_dat.id.sample(config['zurich1']['n_subsample'], replace=False, random_state=1234))
        cell_df = pd.DataFrame({'id': cell_ids})
        cell_df.to_csv(output.cell_ids)

        ## Subsample assignments
        assignments = pd.read_csv(input.assignments, index_col=0)
        assignments_subset = assignments.loc[cell_ids]
        assignments_subset.to_csv(output.assignments)
