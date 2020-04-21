
import pandas as pd
import os

configfile: "imc-config.yml"

output_path = "output/" + config['version'] + "/"

basel_metadata = pd.read_csv(os.path.join(config['basel']['base_dir'], config['basel']['metadata_file']))

basel_cores = list(basel_metadata.core)


rule jackson2020basel:
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


tmp_basel_output = expand(output_path + "basel_processed/{core}.rds", core=basel_cores)

rule all:
    input:
        tmp_basel_output
