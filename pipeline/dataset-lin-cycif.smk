## lin-cycif

# cycif_samples = [
#     'TMA00n',
#     'TMA01n',
#     'TMA02n',
#     'TMA03n',
#     'TMA04n',
#     'TMA05n',
#     'TMA06n',
#     'TMA07n',
#     'TMA08n',
#     'TMA09n',
#     'TMA10n',
#     'TMA11n',
#     'TMA12n',
#     'TMA13n',
#     'TMA14n',
#     'TMA15n',
#     'TMA16n',
#     'TMA17n',
#     'TMA18n',
#     'TMA19n',
#     'TMA20n',
#     'TMA21n',
#     'TMA22n',
#     'TMA23n',
#     'TMA24n',
#     'TMA25n',
#     'TMA26n',
#     'TMA27n',
#     'TMA28n',
#     'TMA29n',
#     'TMA30n',
#     'TMA31n',
#     'TMA32n',
#     'TMA33n',
#     'TMA34n',
#     'TMA35n',
#     'TMA36n',
#     'TMA37n',
#     'TMA38n'
# ]


cycif_samples = [
    'TMA-0',
    'TMA-1',
    'TMA-2',
    'TMA-3',
    'TMA-4',
    'TMA-5',
    'TMA-6',
    'TMA-7',
    'TMA-8',
    'TMA-9',
    'TMA-10',
    'TMA-11',
    'TMA-12',
    'TMA-13',
    'TMA-14',
    'TMA-15',
    'TMA-16',
    'TMA-17',
    'TMA-18',
    'TMA-19',
    'TMA-20',
    'TMA-21',
    'TMA-22',
    'TMA-23',
    'TMA-24',
    'TMA-25',
    'TMA-26',
    'TMA-27',
    'TMA-28',
    'TMA-29',
    'TMA-30',
    'TMA-31',
    'TMA-32',
    'TMA-33',
    'TMA-34',
    'TMA-35',
    'TMA-36',
    'TMA-37',
    'TMA-38'
]



lin_cycif_output = {
    'sces': expand(output_path + "lin-cycif_processed/{sample_id}.rds", sample_id=cycif_samples),
    'loom': output_path + "looms/lin_cycif.loom",
}


rule read_lin_cycif:
    input:
        scdata="data-raw/lin-2019-deepcell/summarized_expression/expression_Composite-{sample_id}.tsv"
    output:
        csv = output_path + "lin-cycif_processed/{sample_id}.csv",
        rds = output_path + "lin-cycif_processed/{sample_id}.rds",
    shell:
        "Rscript pipeline/lin-cycif/parse-lin-cycif2.R "
        "--input_sc {input.scdata} "
        "--id {wildcards.sample_id} "
        "--output_csv {output.csv} "
        "--output_rds {output.rds} "


rule cycif_to_loom:
    input:
        expand(output_path + "lin-cycif_processed/{sample_id}.csv", sample_id=cycif_samples),
    output:
        output_path + "looms/lin_cycif.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/lin-cycif_processed "
        "{output} "

