
## Snakefile for benchmarking state




repeats = 20
gsva_methods = ["plage", "astir"] # "gsva", "ssgsea", 

n_cells = [1000, 5000, 10000, 50000]

datasets = ['basel', 'zurich1', 'wagner', 'schapiro']

output_path_benchmarking = output_path + "benchmarking/geneset/"

geneset_files = expand(output_path_benchmarking + "geneset_benchmark_{d}_{m}_{n}_{r}.csv",
                d=datasets,
                m=gsva_methods, n=n_cells, r = range(repeats))

# geneset_files = geneset_files + \
#     expand(output_path_benchmarking + "geneset_benchmark_{d}_{m}_{n}_{r}.csv",
#             d=['basel', 'zurich1'],
#             m=['plage', 'astir'], n=[50000,100000,200000], r = range(repeats))

benchmarking_output = {
    'files': geneset_files
}

core_list = {
    'basel': basel_cores,
    'zurich1': zurich1_cores,
    'wagner': wagner_samples,
    'schapiro': schapiro_samples
}

rule run_gsva_benchmark:
    params:
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs'],
    input:
        csvs = lambda wildcards: ['{}{}_processed/'.format(output_path, wildcards.d) + f + '.csv' for f in core_list[wildcards.d]],
        markers = lambda wildcards: config[wildcards.d]['marker_file']
    output:
        output_path_benchmarking + "geneset_benchmark_{d}_{m}_{n}_{r}.csv"
    shell:
        "Rscript pipeline/benchmarking/benchmark-gsva-methods2.R "
        "--input_dir {output_path}/{wildcards.d}_processed/ "
        "--markers {input.markers} "
        "--n_cells {wildcards.n} "
        "--method {wildcards.m} "
        "--dataset {wildcards.d} "
        "--output_file {output} "
        "--max_epochs {params.max_epochs} "
        "--batch_size {params.batch_size} "
        "--learning_rate {params.learning_rate} "
        "--n_initial_epochs {params.n_initial_epochs} "

# Rscript pipeline/benchmarking/benchmark-gsva-methods.R --input_dir output/v1//basel_processed/ --markers data-raw/jackson-2020-markers.yml --n_cells 1000 --method astir --output_file test.csv