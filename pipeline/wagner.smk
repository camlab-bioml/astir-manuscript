## All analysis for the Wagner et al. dataset


# Wagner 2019 metadata
wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])


tmp_wagner_output = expand(output_path + 
        "wagner-2019_processed/{sample}.{ext}",
        sample=wagner_samples, ext=['csv','rds'])

wagner_output = {
    'rds_csv': tmp_wagner_output,
    'loom': output_path + "looms/wagner.loom"
}


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

rule wagner_to_loom:
    input:
        expand(output_path + "wagner-2019_processed/{core}.csv", core=wagner_samples),
    output:
        output_path + "looms/wagner.loom"
    shell:
        "python pipeline/dir-of-csvs-to-loom.py "
        "{output_path}/wagner-2019_processed "
        "{output} "