## All analysis for the Keren dataset


keren_cores = pd.read_csv(config['keren']['samples_file'], header=None)
keren_cores = list(keren_cores[keren_cores.columns[0]])


keren_csvs = expand(output_path + "keren_processed/{core}.csv", core=keren_cores)


keren_output = {
    'csvs': keren_csvs,
    'loom': output_path + "looms/keren.loom"
}


rule parse_keren:
    input:
        "data-raw/keren-2018/{core}.csv"
    output:
        rds = output_path + "keren_processed/{core}.rds",
        csv = output_path + "keren_processed/{core}.csv",
    shell:
        "Rscript pipeline/keren-2018/keren-parse-csv.R "
        "--input_sc {input} "
        "--output_rds {output.rds} "
        "--output_csv {output.csv} "


rule keren_to_loom:
    input:
        keren_csvs
    output:
        output_path + "looms/keren.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/keren_processed "
        "{output} "
