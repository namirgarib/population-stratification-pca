import sys
import numpy as np
import matplotlib.pyplot as plt

def load_eigenvalues(file_path):
    eigen_data = np.loadtxt(file_path, delimiter=",")
    components = eigen_data[:, 0]
    eigenvalues = eigen_data[:, 1]
    return components, eigenvalues

def plot_scree(components, eigenvalues, output_file):
    components = components[:4]
    eigenvalues = eigenvalues[:4]
    plt.figure(figsize=(6, 4))
    plt.bar(components, eigenvalues, color='blue', width=0.5)
    plt.title('Scree Plot')
    plt.xlabel('Principal Component')
    plt.xticks(components, [f'PC{i+1}' for i in range(len(components))])
    plt.ylabel('Eigenvalue')
    plt.savefig(output_file)
    plt.close()

def load_pca_results(file_path):
    pca_scores = np.loadtxt(file_path, delimiter=",")
    return pca_scores

def plot_pca(pca_scores, output_prefix):
    # Plot PC1 vs PC2
    plt.figure(figsize=(6, 6))
    plt.grid(color='gray', alpha=0.5, linestyle=':')
    #plt.xlim(-1, 1)
    #plt.ylim(-2, 2)
    plt.scatter(pca_scores[:, 0], pca_scores[:, 1], c='red', marker='o', zorder=2)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PC1 vs PC2')
    for i, txt in enumerate(['ind' + str(i+1) for i in range(pca_scores.shape[0])]):
        plt.annotate(txt, (pca_scores[i, 0], pca_scores[i, 1]), fontsize=8, ha='right')
    plt.savefig(f'{output_prefix}_pc1_pc2.png')
    plt.close()
    
    # Plot PC3 vs PC4 if available
    if pca_scores.shape[1] >= 4:
        plt.figure(figsize=(6, 6))
        #plt.xlim(-1, 1)
        #plt.ylim(-2, 2)
        plt.scatter(pca_scores[:, 2], pca_scores[:, 3], c='green', marker='o', zorder=2)
        plt.xlabel('PC3')
        plt.ylabel('PC4')
        plt.title('PC3 vs PC4')
        plt.grid(color='gray', alpha=0.5, linestyle=':')
        plt.savefig(f'{output_prefix}_pc3_pc4.png')

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot.py <timestamp_folder>")
        sys.exit(1)
    
    # Base path is hardcoded
    base_path = "/Users/namirgarib/Desktop/KANAZAWA UNIVERSITY/3Q3/research/program_c/results"
    sub_folder = sys.argv[1]  # e.g. "20250113133005_3k"
    
    # Construct the full path to the results folder
    results_folder = f"{base_path}/{sub_folder}"

    eigenvalues_path = f"{results_folder}/eigenvalues.csv"
    results_path = f"{results_folder}/results.csv"

    components, eigenvalues = load_eigenvalues(eigenvalues_path)
    plot_scree(components, eigenvalues, f'{results_folder}/scree_plot.png')

    pca_scores = load_pca_results(results_path)
    plot_pca(pca_scores, f'{results_folder}/pc')

if __name__ == "__main__":
    main()