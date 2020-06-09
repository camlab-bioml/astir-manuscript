
analysis_output = {
    'SP43_115_X4Y8_assignments': output_path + "assignments_subset_for_pathology/SP43_115_X4Y8.csv"
}

rule subset_SP45_115_X4Y8:
    input:
        output_path + "astir_assignments/basel_astir_assignments.csv"
    output:
        output_path + "assignments_subset_for_pathology/SP43_115_X4Y8.csv"
    run:
        df = pd.read_csv(input[0], index_col=0)
        print(f"df.shape: {df.shape}")
        df2 = df.loc[ [c for c in df.index if "SP43_115_X4Y8" in c] ]
        df2.to_csv(output[0])

