
granularity_outputs = {
    'HR': expand(output_path + "granularity/HR-astir_assignments/{dataset}_astir_assignments.csv",dataset='basel'),
    'HR2': expand(output_path + "granularity/rare-HR-astir_assignments/rare-HR-{dataset}_astir_assignments.csv",dataset='basel')
}

rule run_astir_granularity:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        loom=ancient(output_path + "looms/{dataset}.loom"),
        markers="markers/jackson-2020-markers-HR.yml",
    output:
        csv=output_path + "granularity/HR-astir_assignments/{dataset}_astir_assignments.csv",
        fig=output_path + "granularity/HR-astir_assignments/{dataset}_astir_loss.png",
        diagnostics=output_path + "granularity/HR-astir_assignments/{dataset}_diagnostics.csv"
    run:
        # fit
        from astir.data import from_loompy_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(input.loom, input.markers)
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

rule run_astir_granularity2:
    params:
        op = output_path,
        max_epochs = config['astir_opts']['max_epochs'],
        batch_size = config['astir_opts']['batch_size'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        loom=ancient(output_path + "looms/{dataset}.loom"),
        markers="markers/jackson-2020-markers-rare.yml",
    output:
        csv=output_path + "granularity/rare-HR-astir_assignments/rare-HR-{dataset}_astir_assignments.csv",
        fig=output_path + "granularity/rare-HR-astir_assignments/rare-HR-{dataset}_astir_loss.png",
        diagnostics=output_path + "granularity/rare-HR-astir_assignments/rare-HR-{dataset}_diagnostics.csv"
    run:
        # fit
        from astir.data import from_loompy_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_loompy_yaml(input.loom, input.markers)
        print(f"{datetime.now()}\t Fitting model")

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/1000)

        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(5e-3),
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