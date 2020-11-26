## lin-cycif

cycif_samples = [
    'TMA00n',
    'TMA01n',
    'TMA02n',
    'TMA03n',
    'TMA04n',
    'TMA05n',
    'TMA06n',
    'TMA07n',
    'TMA08n',
    'TMA09n',
    'TMA10n',
    'TMA11n',
    'TMA12n',
    'TMA13n',
    'TMA14n',
    'TMA15n',
    'TMA16n',
    'TMA17n',
    'TMA18n',
    'TMA19n',
    'TMA20n',
    'TMA21n',
    'TMA22n',
    'TMA23n',
    'TMA24n',
    'TMA25n',
    'TMA26n',
    'TMA27n',
    'TMA28n',
    'TMA29n',
    'TMA30n',
    'TMA31n',
    'TMA32n',
    'TMA33n',
    'TMA34n',
    'TMA35n',
    'TMA36n',
    'TMA37n',
    'TMA38n'
]


lin_cycif_output = {
    'sces': expand(output_path + "lin-cycif_processed/{sample_id}.rds", sample_id=cycif_samples),
    'loom': output_path + "looms/lin_cycif.loom",
}


rule read_lin_cycif:
    input:
        scdata="data-raw/lin-2020/raw-data/TMApanels/{sample_id}.csv"
    output:
        csv = output_path + "lin-cycif_processed/{sample_id}.csv",
        rds = output_path + "lin-cycif_processed/{sample_id}.rds",
    shell:
        "Rscript pipeline/lin-cycif/parse-lin-cycif.R "
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

