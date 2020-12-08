


segmentation_output = {
    'astir_assignments': expand(
        output_path + "schapiro_astir_assignments_alt_mask/assignments_{alt_s}_{user}.csv",
        alt_s=schapiro_alt_mask_samples,user=schapiro_users
    ),
    'figure': output_path + "figures/segmentation/segmentation.pdf"
}





rule run_astir_type_alt_segmentation:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        csv=output_path + "schapiro_processed_alt_mask/{alt_s}_{user}.csv",
        markers= config['schapiro']['marker_file'],
    output:
        csv=output_path + "schapiro_astir_assignments_alt_mask/assignments_{alt_s}_{user}.csv",
        fig=output_path + "schapiro_astir_assignments_alt_mask/loss_{alt_s}_{user}.png",
        diagnostics=output_path + "schapiro_astir_assignments_alt_mask/diagnostics_{alt_s}_{user}.csv"
    run:
        # fit
        from astir.data import from_csv_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_csv_yaml(input.csv, input.markers)
        print(f"{datetime.now()}\t Fitting model")

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)


        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate),
        n_init_epochs=int(params.n_initial_epochs))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.get_celltype_probabilities().to_csv(output.csv)

        # run diagnostics
        ast.diagnostics_celltype().to_csv(output.diagnostics)

        # plot loss
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
        plt.ylabel("Loss")
        plt.xlabel("Epoch")
        plt.tight_layout()
        plt.savefig(output.fig, dpi=300)

rule segmentation_figures:
    params:
        dir_astir = output_path + "schapiro_astir_assignments_alt_mask/",
        dir_other = output_path + "results/alternative_masks/"
    input: 
        segmentation_output['astir_assignments'],
    output:
        pdf = segmentation_output['figure']
    script:
        "segmentation-consistency.R"

