
## Snakefile for benchmarking

repeats = 20
gsva_methods = ["gsva", "ssgsea", "zscore", "plage", "astir"]

n_cells = [1000, 5000, 10000]

output_path_benchmarking = output_path + "benchmarking/geneset/"

geneset_files = expand(output_path_benchmarking + "geneset_benchmark_{m}_{n}_{r}.csv",
                m=gsva_methods, n=n_cells, r = range(repeats))

geneset_files = geneset_files + \
    expand(output_path_benchmarking + "geneset_benchmark_{m}_{n}_{r}.csv",
            m=['plage', 'astir'], n=[50000,100000,200000], r = range(repeats))



rule run_gsva_benchmark:
    input:
        all_csvs = expand(output_path + "basel_processed/{core}.csv", core=basel_cores),
        markers = "data-raw/jackson-2020-markers.yml"
    output:
        output_path_benchmarking + "geneset_benchmark_{m}_{n}_{r}.csv"
    shell:
        "Rscript pipeline/benchmarking/benchmark-gsva-methods.R "
        "--input_dir {output_path}/basel_processed/ "
        "--markers {input.markers} "
        "--n_cells {wildcards.n} "
        "--method {wildcards.m} "
        "--output_file {output} "

# Rscript pipeline/benchmarking/benchmark-gsva-methods.R --input_dir output/v1//basel_processed/ --markers data-raw/jackson-2020-markers.yml --n_cells 1000 --method astir --output_file test.csv