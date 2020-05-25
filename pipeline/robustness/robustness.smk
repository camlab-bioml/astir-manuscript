


## Snakefile for assessing the robustness of astir to undesirable implementations. Currently:
## (1)  create_reduced_sets: remove cell types from described markers and assess how well 
##      inference is performed
## (2)  

import yaml

celltypes_to_remove = {
    'Stromal_only': ['Stromal'],
    'Stromal_Macrophage': ['Stromal', 'Macrophage'],
    'Stromal_Macrophage_Monocyte': ['Stromal', 'Macrophage', 'Monocyte'],
    'None': []
}

reduced_assignments = expand(output_path + "robustness/assignments-15k-removed-{removed}.csv", 
    removed=celltypes_to_remove.keys())

rule create_reduced_sets:
    input:
        markers="markers/jackson-2020-markers.yml",
        csvs=basel_15k_csvs
    output:
        markers=output_path + "robustness/{removed}.yml",
        assignments = output_path + "robustness/assignments-15k-removed-{removed}.csv",
    run:
        marker_dict = None
        with open(input.markers, "r") as stream:
            marker_dict = yaml.safe_load(stream)
            for t in celltypes_to_remove[wildcards.removed]:
                marker_dict['cell_types'].pop(t)
        
        with open(output.markers, 'w') as stream:
            yaml.dump(marker_dict, stream)
            
        from astir.data_readers import from_csv_dir_yaml

        ast = from_csv_dir_yaml(os.path.join(output_path, "basel_15k_subset"), output.markers, include_beta=False)
        ast.fit_type(max_epochs = 20, batch_size = 512, learning_rate = 1e-3)
        ast.type_to_csv(output.assignments)



