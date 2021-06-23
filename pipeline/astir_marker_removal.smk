
markers_to_rem = {
    'Cytokeratin-7': ['Cytokeratin 7'],
    'Cytokeratin-7_Cytokeratin_19': ['Cytokeratin 7', 'Cytokeratin 19'],
    'all_luminal': ['Cytokeratin 7', 'Cytokeratin 19', 'Cytokeratin 8/18', 'pan Cytokeratin', 'E-Cadherin'],
    'Cytokeratin-19': ['Cytokeratin 19'],
    'Cytokeratin-8_18': ['Cytokeratin 8/18']
}

seeds = range(5)

astir_marker_removal = {
    'markers_remove': expand("markers/jackson-2020-markers-v4-removed-{remove_marker}.yml", remove_marker = list(markers_to_rem.keys())),
    'astir_rem': expand(output_path + "astir_iterative_marker_removal/basel-remove-{markers}-runSeed_{seed}.csv", markers = list(markers_to_rem.keys()), seed = seeds),
    'viz': expand(output_path + "figures/astir-marker-removal-runSeed_{seed}.pdf", seed = seeds),
    'viz2': expand(output_path + "figures/astir-marker-removal-{marker}-percent-{thresh}-runSeed_{seed}.pdf", marker = ['Cytokeratin-8_18', 'Cytokeratin-19', 'Cytokeratin-7'], thresh = [0.5, 0.9], seed = seeds),
    'heatmap': expand(output_path + "figures/astir-marker-removal-heatmap-{marker}-percent-{thresh}-runSeed_{seed}.pdf", marker = ['Cytokeratin-8_18', 'Cytokeratin-19', 'Cytokeratin-7'], thresh = [0.5, 0.9], seed = seeds)
}


rule create_markers:
    input:
        markers = "markers/jackson-2020-markers-v4.yml"
    output:
        m = "markers/jackson-2020-markers-v4-removed-{remove_marker}.yml"
    run:
        import yaml
        from pathlib import Path

        with open(input.markers) as file:
            markers = yaml.full_load(file)

        for x in markers_to_rem[wildcards.remove_marker]:
            markers['cell_types']['Epithelial (luminal)'].remove(x)

        with open(output.m, 'w') as file:
            yaml.dump(markers, file)



rule astir_marker_removal:
    params:
        max_epochs = config['astir_opts']['max_epochs'],
        learning_rate = config['astir_opts']['learning_rate'],
        n_initial_epochs = config['astir_opts']['n_initial_epochs']
    input:
        anndata=ancient(output_path + "anndata/basel.h5ad"),
        markers = "markers/jackson-2020-markers-v4-removed-{markers}.yml"
    output:
        csv = output_path + "astir_iterative_marker_removal/basel-remove-{markers}-runSeed_{seed}.csv",
        fig = output_path + "astir_iterative_marker_removal/basel-remove-{markers}_astir_loss_{seed}.png",
        diagnostics = output_path + "astir_iterative_marker_removal/basel-remove-{markers}-diagnostics_{seed}.csv",
        delta = output_path + 'astir_iterative_marker_removal/basel-remove-{markers}-delta-{seed}.tsv'
    run:
        from astir.data import from_anndata_yaml
        from datetime import datetime
        print(f"{datetime.now()}\t Reading loom file ")
        ast = from_anndata_yaml(input.anndata, input.markers, random_seed = int(wildcards.seed))
        print(f"{datetime.now()}\t Fitting model")

        N = ast.get_type_dataset().get_exprs_df().shape[0]
        batch_size = int(N/100)

        ast.fit_type(max_epochs = int(params.max_epochs), 
        batch_size = batch_size, 
        learning_rate = float(params.learning_rate),
        n_init_epochs=int(params.n_initial_epochs))
        print(f"{datetime.now()}\t Finished fitting model")
        ast.get_celltype_probabilities().to_csv(output.csv)

        # Save delta
        delta = ast._type_ast._variables['log_delta'].detach().numpy()
        delta = pd.DataFrame(delta)
        cell_types = ast._type_dset._classes

        # add column and row indexes
        delta.columns = cell_types + ['Other']
        delta.index = ast._type_dset._m_features
        delta.to_csv(output.delta, sep = '\t')

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
        

rule plot_astir_marker_removal:
    input:
        none = output_path + "astir_assignments/basel_astir_assignments.csv",
        cyto7 = output_path + "astir_iterative_marker_removal/basel-remove-Cytokeratin-7-runSeed_{seed}.csv",
        cyto7_19 = output_path + "astir_iterative_marker_removal/basel-remove-Cytokeratin-7_Cytokeratin_19-runSeed_{seed}.csv",
        all_luminal = output_path + "astir_iterative_marker_removal/basel-remove-all_luminal-runSeed_{seed}.csv"
    container: "astir-manuscript.sif"
    output:
        pdf = output_path + "figures/astir-marker-removal-runSeed_{seed}.pdf"
    script:
        "analysis/Astir-marker-removal.R"


rule plot_astir_marker_removal_pairwise:
    input:
        normal = output_path + "astir_assignments/basel_astir_assignments.csv",
        removed = output_path + "astir_iterative_marker_removal/basel-remove-{marker}-runSeed_{seed}.csv",
        sce = output_path + "sces/basel_sce.rds",
        markers = config['basel']['marker_file']
    container: "astir-manuscript.sif"
    output:
        pdf = output_path + "figures/astir-marker-removal-{marker}-percent-{thresh}-runSeed_{seed}.pdf",
        heatmap = output_path + "figures/astir-marker-removal-heatmap-{marker}-percent-{thresh}-runSeed_{seed}.pdf"
    script:
        "analysis/Astir-marker-removal-2.R"