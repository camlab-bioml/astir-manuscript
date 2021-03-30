
epithelial_props = [0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.99]

imbalance_output = {
    'astir_assignments': expand(output_path + "results/imbalance/astir_imbalance_{p}.csv",
    p = epithelial_props),
    'imbalance_fig': output_path + "figures/imbalance/imbalance.pdf"
}



rule imbalance_figures:
    params:
        input_dir_astir = output_path + "results/imbalance",
        input_dir_other = output_path + "results/epithelial_overclustering",
    input:
        imbalance_output['astir_assignments'],
    output:
        pdf = imbalance_output['imbalance_fig'],
    script:
        "imbalance-plots.R"


rule run_astir_imbalance:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        loom=ancient(output_path + "looms/basel.loom"),
        cells_to_sample=output_path + "results/epithelial_overclustering/luminal-{p}.csv",
        markers=config['basel']['markers_imbalance'],
    output:
        csv=output_path + "results/imbalance/astir_imbalance_{p}.csv",
        loom=output_path + "results/imbalance/looms/basel-{p}.loom"
    run:
        # fit
        from astir.data import from_loompy_yaml
        import pandas as pd
        import loompy
        import numpy as np
        from datetime import datetime

        cell_ids_to_sample = pd.read_csv(input.cells_to_sample)
        cell_ids_to_sample = list(cell_ids_to_sample.cell_id)

        ds = loompy.connect(input.loom)
        cell_names = list(ds.col_attrs['cell_name'])

        full_mat = ds[:,:]

        keep_cell_index = np.isin(cell_names, cell_ids_to_sample)
        full_mat = full_mat[:, keep_cell_index]

        new_batch = np.array(ds.col_attrs['batch'])[keep_cell_index]

        new_col_attr = {
            'batch': list(new_batch),
            'cell_name': cell_ids_to_sample
        }

        loompy.create(output.loom, full_mat, ds.row_attrs, new_col_attr)


        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(output.loom, input.markers)
        print(f"{datetime.now()}\t Fitting model")

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)


        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate),
        n_init_epochs=int(params.n_initial_epochs))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.get_celltype_probabilities().to_csv(output.csv)