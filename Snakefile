
import pandas as pd
import os

configfile: "imc-config.yml"
logdir = "logs/slurm"
output_path = "output/" + config['version'] + "/"

# Jackson 2020 metadata
basel_metadata = pd.read_csv(os.path.join(config['basel']['base_dir'], config['basel']['metadata_file']))
basel_cores = list(basel_metadata.core)

# Wagner 2019 metadata
wagner_sample_df = pd.read_csv(config['wagner']['sample_file'], header=None)
wagner_samples = list(wagner_sample_df[0])

tmp_basel_output = expand(output_path + "basel_processed/{core}.rds", core=basel_cores)

tmp_wagner_output = expand(output_path + 
        "wagner-2019_processed/{sample}.{ext}",
        sample=wagner_samples, ext=['csv','rds'])

## Beginning of rules ----- 
rule all:
    input:
        # tmp_basel_output
        tmp_wagner_output

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



