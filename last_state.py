# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 17:21:35 2024

@author: kaders
"""

import fa2
import louvain
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as multi
from gprofiler import GProfiler


# importing data
adata = sc.read_h5ad("GSE243292_Microglia_GEO.h5ad")


# this code is for related to configuring settings for the Scanpy library in Python, which is commonly used for scRNA-seq data analysis.
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Filter nuclei based on expressed cells ;# Remove genes with 3 expression across all cells
sc.pp.filter_genes(adata, min_cells=3)

# Filter nuclei based on expressed genes
sc.pp.filter_cells(adata, min_genes=200)
# Filter out low-quality cells
sc.pp.filter_cells(adata, min_counts=1)  
# Filter out low-expressing genes
sc.pp.filter_genes(adata, min_counts=1)  


# Print the number of cells with zero counts
cells_with_zero_counts = (adata.X.sum(axis=1) == 0).sum()
print(f"Number of cells with zero counts: {cells_with_zero_counts}")

# I also run this code to take them unique as gene ids
adata.var_names_make_unique()

# Normalization
sc.pp.normalize_per_cell(adata)
# Log transformation
sc.pp.log1p(adata)
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pl.highly_variable_genes(adata)

#Reduce the dimensionality of the data by running PCA which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')

# to visualize them
sc.pl.violin(adata, ['n_genes', 'n_counts'],
             jitter=0.4, multi_panel=True)

# before clustering neighbors were determined.
sc.pp.neighbors(adata)

# determined cluster by using louvain method. Resolution was adjusted to 0.6 because I should find 7 clusters. 
sc.tl.louvain(adata, resolution= 0.6)
# also see them as umap
sc.tl.umap(adata)

# visualize them with differrent palettes and density levels
# Plot UMAP with adjusted resolution Louvain clusters
sc.pl.umap(adata, color='louvain', legend_loc='on data', palette='tab10', alpha=0.5, size=30, frameon=False, legend_fontsize=10)
sc.pl.umap(adata, color='louvain', legend_loc='on data', palette='viridis', alpha=0.7, size=30, frameon=False, legend_fontsize=10)

#most density
sc.pl.umap(adata, color='louvain', legend_loc='on data', palette='tab10', alpha=1.0, size=30, frameon=False, legend_fontsize=10)



# Create a list to store subcluster-specific AnnData objects
subcluster_adatas = []

# Iterate over each Louvain cluster and analyze separately
for cluster_label in ['0', '1', '2', '3','4', '5', '6']: # my clusters labeles 0 to 6
    # Filter AnnData for the current cluster
    adata_subcluster = adata[adata.obs['louvain'] == cluster_label, :]

    # Store the subcluster-specific AnnData object
    subcluster_adatas.append(adata_subcluster)


# Computing t-SNE to see them
sc.tl.tsne(adata)

# Our AnnData object with t-SNE coordinates and cell type annotations with louvain
tsne_coordinates = adata.obsm['X_tsne']

# Define a custom color palette for Louvain clusters
custom_palette = sns.color_palette("Set2", n_colors=len(adata.obs['louvain'].unique()))

# Creating a scatter plot with Seaborn using the custom color palette
plt.figure(figsize=(10, 6))
sns.scatterplot(x=tsne_coordinates[:, 0], y=tsne_coordinates[:, 1], hue=adata.obs['louvain'].astype(str), palette=custom_palette, legend='full')

# Adding labels, title, legend, etc.
plt.title('t-SNE of Microglia Subpopulations (Louvain Clusters)')
plt.xlabel('t-SNE Dimension 1')
plt.ylabel('t-SNE Dimension 2')
plt.legend(title='Louvain Cluster', loc='upper left',  bbox_to_anchor=(1, 1))

plt.show()



# Get unique louvain cluster labels
unique_clusters = adata.obs['louvain'].unique()

# Perform differential gene expression analysis for each cluster
for cluster_label in unique_clusters:
    sc.tl.rank_genes_groups(adata, groupby='louvain', groups=[cluster_label], method='wilcoxon', use_raw=False, corr_method='benjamini-hochberg')

    # Extract results
    results = adata.uns['rank_genes_groups']

    # Convert results to DataFrames
    df_names = pd.DataFrame(results['names'])
    df_scores = pd.DataFrame(results['scores'])
    df_pvals = pd.DataFrame(results['pvals'])
    df_pvals_adj = pd.DataFrame(results['pvals_adj'])
    df_logfoldchanges = pd.DataFrame(results['logfoldchanges'])

    # Get the cluster-specific DataFrame
    cluster_df = pd.DataFrame({
        'Gene': df_names[cluster_label],
        'Score': df_scores[cluster_label],
        'P-value': df_pvals[cluster_label],
        'Adjusted P-value': df_pvals_adj[cluster_label],
        'Log-fold Change': df_logfoldchanges[cluster_label]
    })

    # Filter for significant genes (you can adjust the threshold)
    significant_genes = cluster_df[cluster_df['Adjusted P-value'] < 0.05]

    # Sort by score or other criteria to get top genes
    top_markers = significant_genes.sort_values(by='Score', ascending=False).head(20)

    # Display the top markers for the current cluster
    print(f"Top 10 markers for cluster {cluster_label}:\n")
    print(top_markers)
    print("\n" + "="*50 + "\n")



# Create a dictionary to store top markers for each cluster
top_markers_by_cluster = {}

# Perform differential gene expression analysis for each cluster
for cluster_label in unique_clusters:
    sc.tl.rank_genes_groups(adata, groupby='louvain', groups=[cluster_label], method='wilcoxon', use_raw=False, corr_method='benjamini-hochberg')

    # Extract results
    results = adata.uns['rank_genes_groups']

    # Convert results to DataFrames
    df_names = pd.DataFrame(results['names'])
    df_scores = pd.DataFrame(results['scores'])
    df_pvals = pd.DataFrame(results['pvals'])
    df_pvals_adj = pd.DataFrame(results['pvals_adj'])
    df_logfoldchanges = pd.DataFrame(results['logfoldchanges'])

    # Get the cluster-specific DataFrame
    cluster_df = pd.DataFrame({
        'Gene': df_names[cluster_label],
        'Score': df_scores[cluster_label],
        'P-value': df_pvals[cluster_label],
        'Adjusted P-value': df_pvals_adj[cluster_label],
        'Log-fold Change': df_logfoldchanges[cluster_label]
    })

    # Filter for significant genes (you can adjust the threshold)
    significant_genes = cluster_df[cluster_df['Adjusted P-value'] < 0.05]

    # Sort by score or other criteria to get top genes
    top_markers = significant_genes.sort_values(by='Score', ascending=False).head(10)

    # Store in the dictionary
    top_markers_by_cluster[cluster_label] = top_markers


top_markers_by_cluster # cluster top dge genes

# cluster annotations were done according to marker genes and other article
cluster_to_group_mapping = {
    0: 'MicrogliaSub1', #ARM
    1: 'MicrogliaSub2', #Motile
    2: 'Oligo_Microglia',
    3: 'Mixed_Microglia',
    4: 'MicrogliaSub3',
    5: 'Astro_Microglia',
    6: 'MicrogliaSub4',
}

adata.obs['Subpopulation_Group'] = adata.obs['louvain'].astype(int).map(cluster_to_group_mapping)

sc.tl.rank_genes_groups(adata, groupby='louvain', groups=[0,1,2,3,4,5,6], method='wilcoxon', use_raw=False, corr_method='benjamini-hochberg')

sc.tl.rank_genes_groups(adata, groupby='louvain', groups=[0,1,2,3,4,5,6], method='wilcoxon', use_raw=False, corr_method='benjamini-hochberg')
sc.pl.rank_genes_groups_heatmap(
    adata,
    groupby= 'Subpopulation_Group',  # Assuming you are using Louvain clusters
    groups=['0','1','2','3','4','5','6'],
    n_genes= 10,  # Number of top genes to show in the heatmap
    cmap='viridis',  # Specify a colormap 
    show_gene_labels=True,  # Show gene labels on the heatmap
    dendrogram=True, 
    standard_scale='var',  # Standardize data by variance
    use_raw=False,
    )

# Combine all DataFrames into a single DataFrame
all_markers_df = pd.concat(top_markers_by_cluster.values(), keys=top_markers_by_cluster.keys(), names=['Cluster', 'Row'])
# Reset index for better indexing in the heatmap
# Replace 'cluster' with the correct column name containing your cluster information
all_markers_df['cluster'] = all_markers_df.groupby('Cluster').cumcount()
all_markers_df_reset = all_markers_df.reset_index()


# Removing Doublets

new_adata = adata[(adata.obs['Subpopulation_Group'] != 'Oligo_Microglia'), :]
new_adata = new_adata[(new_adata.obs['Subpopulation_Group'] != 'Mixed_Microglia'), :]
new_adata = new_adata[(new_adata.obs['Subpopulation_Group'] != 'Astro_Microglia'), :]

adata = new_adata
# Map categorical values to numeric values
category_mapping = {'R47H': 0, 'WT': 1}
adata.obs['trem2_numeric'] = adata.obs['trem2'].map(category_mapping)
adata.obs['trem2_numeric'] = pd.to_numeric(adata.obs['trem2_numeric'], errors='coerce')

plt.figure(figsize=(15, 6))
sns.countplot(x='Subpopulation_Group', hue='trem2_numeric', data=adata.obs)
legend_labels = {0: 'R47H', 1: 'WT'}
plt.legend(title='trem2_numeric', labels=[f'{key}: {value}' for key, value in legend_labels.items()])
plt.title('Distribution of trem2_numeric in Subgroups')
plt.show()

# Map categorical values to numeric values
apoe_mapping = {'E3/E3': 0, 'E3/E4': 1, 'E4/E4': 2}
adata.obs['apoe_numeric'] = adata.obs['apoe'].map(apoe_mapping)
adata.obs['apoe_numeric'] = pd.to_numeric(adata.obs['apoe_numeric'], errors='coerce')

plt.figure(figsize=(12, 8))
sns.countplot(x='Subpopulation_Group',hue='apoe_numeric', data=adata.obs)
# Create a legend using Matplotlib
legend_labels = {0: 'E3/E3', 1: 'E3/E4', 2: 'E4/E4'}
plt.legend(title='apoe_numeric', labels=[f'{key}: {value}' for key, value in legend_labels.items()])

plt.title('Average apoe_numeric in Subpopulation Groups')
plt.show()


# filtering ARM and Motile subgroups

filter_data = new_adata[(new_adata.obs['Subpopulation_Group'] != 'MicrogliaSub3'), :]

filter_data = filter_data[(filter_data.obs['Subpopulation_Group'] != 'MicrogliaSub4'), :]

# pseudotime trajectory analysis with PAGA
# adata is our current data without doublets

# Step 1: Compute Neighborhood Graph
sc.pp.neighbors(filter_data, n_neighbors=15, method='umap', metric='cosine')

# Step 2: Run PAGA
sc.tl.paga(filter_data)

# Step 3: Generate Graph and Plot PAGA
sc.pl.paga(filter_data)  # Adjust the threshold as needed

# Step 4: Run Diffusion Map
sc.tl.diffmap(filter_data)

# Step 5: Update Neighborhood Graph using Diffusion Map representation
sc.pp.neighbors(filter_data, n_neighbors=15, use_rep='X_diffmap')

# Step 6: Draw Graph
sc.tl.draw_graph(filter_data)

# Step 7: Plot Graph with color annotations
sc.pl.draw_graph(filter_data, color='Subpopulation_Group', legend_loc='on data')

# Step 8: Run PAGA with specific groups
sc.tl.paga(filter_data, groups='Subpopulation_Group')

# Step 9: Plot PAGA with color annotations
sc.pl.paga(filter_data, color=['Subpopulation_Group'], title=['Cell type'])

# Step 10: Draw Graph using PAGA initialization
sc.tl.draw_graph(filter_data, init_pos='paga')

# Step 11: Compare PAGA graphs
sc.pl.paga_compare(
    filter_data, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)

# Step 12: Set root for diffusion pseudotime
filter_data.uns['iroot'] = np.flatnonzero(filter_data.obs['Subpopulation_Group'] == 'MicrogliaSub1')[0]

# Step 13: Run diffusion pseudotime
sc.tl.dpt(filter_data)

# Step 14: Plot Graph with pseudotime and color annotations
sc.pl.draw_graph(filter_data, color=['Subpopulation_Group', 'dpt_pseudotime'], legend_loc='on data', size=20, alpha=0.8)

# Step 15: Plot PAGA with pseudotime annotations
sc.pl.paga(filter_data, color=['dpt_pseudotime'], threshold=0.1, node_size_scale=1.5, node_size_power=0.5)


# some different plots to see changes of pseudotime values btw ARM and Motile


sc.pl.draw_graph(filter_data, color=['Subpopulation_Group'], title=['Cluster'], legend_loc='on data',size=20, alpha=0.8,cmap='viridis', frameon=False)
sc.pl.draw_graph(filter_data, color=['trem2'], title=['TREM2'],size=20, alpha=0.8,cmap='viridis', frameon=False)
sc.pl.draw_graph(filter_data, color=['CD163', 'FGD4','CX3CR1', 'TREM2'], title=['CD163_ARM', 'FGD4_Motile','CX3CR1_Homeostatic', 'TREM2'],size=20, alpha=0.8,cmap='viridis', frameon=False)
sc.pl.draw_graph(filter_data, color=['trem2'], title=['TREM2'],size=20, alpha=0.8,cmap='viridis', frameon=False)
sc.pl.draw_graph(filter_data, color=['apoe'], title=['APOE'],size=20, alpha=0.8,cmap='viridis', frameon=False)
sc.pl.draw_graph(filter_data, color=['dpt_pseudotime'], title=['dpt_pseudotime'],size=20, alpha=0.8,cmap='viridis', frameon=False)

sc.pl.draw_graph(filter_data, color=['SPP1', 'APOE','CD74', 'CD44'], title=['SPP1', 'APOE','CD74', 'CD44'],size=20, alpha=0.8,cmap='viridis', frameon=False)
sc.pl.draw_graph(filter_data, color=['P2RY12', 'CCL3','CSF3R', 'TMEM119'], title=['P2RY12', 'CCL3','CSF3R', 'TMEM119'],size=20, alpha=0.8,cmap='viridis', frameon=False)


# linear mixed effect modeling for both group

ARM = filter_data[(filter_data.obs['Subpopulation_Group'] == 'MicrogliaSub1'), :]
Motile = filter_data[(filter_data.obs['Subpopulation_Group'] == 'MicrogliaSub2'), :]

# this linear mixed effect modeling are without seperation of APOE genotypes

# for ARM I prepared data to apply linear mixed effect modeling in R

# changing type of columns to as needed as they acceps as references
ARM.obs['amyloid_pathology'] = np.where(ARM.obs['atscore'].str.contains('A\+'), 1, 0)
ARM.obs['tau_pathology'] = np.where(ARM.obs['atscore'].str.contains('T\+'), 1, 0)
ARM.obs['APOE_genotype'] = np.where(ARM.obs['apoe'] == 'E3/E3', 0, 1)
ARM.obs['TREM2_genotype'] = np.where(ARM.obs['trem2'] == 'WT', 0, 1)
ARM.obs['sample_ID'] = ARM.obs['sampleID'].astype(str)

# for creating data 
pseudotime_values = ARM.obs['dpt_pseudotime']
amyloid_pathology = ARM.obs['amyloid_pathology']  
tau_pathology = ARM.obs['tau_pathology']      
APOE_genotype = ARM.obs['APOE_genotype']   
TREM2_genotype = ARM.obs['TREM2_genotype']  
sample_ID = ARM.obs['sample_ID']   

# Creating a dictionary with your data
data = {
    'pseudotime': pseudotime_values,
    'amyloid_pathology': amyloid_pathology,
    'tau_pathology': tau_pathology,
    'APOE_genotype': APOE_genotype,
    'TREM2_genotype': TREM2_genotype,
    'sample_ID': sample_ID
}

# Creating a pandas DataFrame from the dictionary
df = pd.DataFrame(data)
df['amyloid_pathology'] = df['amyloid_pathology'].astype(int)
df['tau_pathology'] = df['tau_pathology'].astype(int)
df['APOE_genotype'] = df['APOE_genotype'].astype(int)
df['TREM2_genotype'] = df['TREM2_genotype'].astype(int)
df['sample_ID'] = df['sample_ID'].astype(int)
df['pseudotime'] = df['pseudotime']

df.to_excel('ARM.xlsx', sheet_name='new_sheet_name')



# for Motile I prepared data to apply linear mixed effect modeling in R

# changing type of columns to as needed as they acceps as references
Motile.obs['amyloid_pathology'] = np.where(Motile.obs['atscore'].str.contains('A\+'), 1, 0)
Motile.obs['tau_pathology'] = np.where(Motile.obs['atscore'].str.contains('T\+'), 1, 0)
Motile.obs['APOE_genotype'] = np.where(Motile.obs['apoe'] == 'E3/E3', 0, 1)
Motile.obs['TREM2_genotype'] = np.where(Motile.obs['trem2'] == 'WT', 0, 1)
Motile.obs['sample_ID'] = Motile.obs['sampleID'].astype(str)

# for creating data 
pseudotime_values = Motile.obs['dpt_pseudotime']
amyloid_pathology = Motile.obs['amyloid_pathology']  
tau_pathology = Motile.obs['tau_pathology']      
APOE_genotype = Motile.obs['APOE_genotype']   
TREM2_genotype = Motile.obs['TREM2_genotype']  
sample_ID = Motile.obs['sample_ID']   

# Creating a dictionary with your data
data = {
    'pseudotime': pseudotime_values,
    'amyloid_pathology': amyloid_pathology,
    'tau_pathology': tau_pathology,
    'APOE_genotype': APOE_genotype,
    'TREM2_genotype': TREM2_genotype,
    'sample_ID': sample_ID
}

# Creating a pandas DataFrame from the dictionary
df = pd.DataFrame(data)
df['amyloid_pathology'] = df['amyloid_pathology'].astype(int)
df['tau_pathology'] = df['tau_pathology'].astype(int)
df['APOE_genotype'] = df['APOE_genotype'].astype(int)
df['TREM2_genotype'] = df['TREM2_genotype'].astype(int)
df['sample_ID'] = df['sample_ID'].astype(int)
df['pseudotime'] = df['pseudotime']

df.to_excel('Motile.xlsx', sheet_name='new_sheet_name')


# Now linear rmixed effect modeling for ARM with seperation of APOE genotypes

# for ARM I prepared data to apply linear mixed effect modeling in R

# changing type of columns to as needed as they acceps as references
ARM.obs['amyloid_pathology'] = np.where(ARM.obs['atscore'].str.contains('A\+'), 1, 0)
ARM.obs['tau_pathology'] = np.where(ARM.obs['atscore'].str.contains('T\+'), 1, 0)
ARM.obs['APOE_genotype[E3/E4]'] = np.where(ARM.obs['apoe'] == 'E3/E4', 1, 0)
ARM.obs['APOE_genotype[E4/E4]'] = np.where(ARM.obs['apoe'] == 'E4/E4', 1, 0)
ARM.obs['TREM2_genotype'] = np.where(ARM.obs['trem2'] == 'WT', 0, 1)
ARM.obs['sample_ID'] = ARM.obs['sampleID'].astype(str)

# for creating data 
pseudotime_values = ARM.obs['dpt_pseudotime']  
amyloid_pathology = ARM.obs['amyloid_pathology']  
tau_pathology = ARM.obs['tau_pathology']     
APOE_genotype1 = ARM.obs['APOE_genotype[E4/E4]']   
TREM2_genotype = ARM.obs['TREM2_genotype']
sample_ID = ARM.obs['sample_ID']
APOE_genotype2 = ARM.obs['APOE_genotype[E3/E4]']  

# Creating a dictionary with your data
data = {
    'pseudotime': pseudotime_values,
    'amyloid_pathology': amyloid_pathology,
    'tau_pathology': tau_pathology,
    'APOE_genotype[E4/E4]': APOE_genotype1,
    'APOE_genotype[E3/E4]': APOE_genotype2,
    'TREM2_genotype': TREM2_genotype,
    'sample_ID': sample_ID
}


# Creating a pandas DataFrame from the dictionary
df = pd.DataFrame(data)
df['amyloid_pathology'] = df['amyloid_pathology'].astype(int)
df['tau_pathology'] = df['tau_pathology'].astype(int)
df['APOE_genotype[E3/E4]'] = df['APOE_genotype[E3/E4]'].astype(int)
df['APOE_genotype[E4/E4]'] = df['APOE_genotype[E4/E4]'].astype(int)
df['TREM2_genotype'] = df['TREM2_genotype'].astype(int)
df['sample_ID'] = df['sample_ID'].astype(int)
df['pseudotime'] = df['pseudotime']

df.to_excel('ARM_with_APOEgenotype.xlsx', sheet_name='new_sheet_name')


# for Motile I prepared data to apply linear mixed effect modeling in R

# changing type of columns to as needed as they acceps as references
Motile.obs['amyloid_pathology'] = np.where(Motile.obs['atscore'].str.contains('A\+'), 1, 0)
Motile.obs['tau_pathology'] = np.where(Motile.obs['atscore'].str.contains('T\+'), 1, 0)
Motile.obs['APOE_genotype[E3/E4]'] = np.where(Motile.obs['apoe'] == 'E3/E4', 1, 0)
Motile.obs['APOE_genotype[E4/E4]'] = np.where(Motile.obs['apoe'] == 'E4/E4', 1, 0)
Motile.obs['TREM2_genotype'] = np.where(Motile.obs['trem2'] == 'WT', 0, 1)
Motile.obs['sample_ID'] = Motile.obs['sampleID'].astype(str)

# for creating data 
pseudotime_values = Motile.obs['dpt_pseudotime']  
amyloid_pathology = Motile.obs['amyloid_pathology']  
tau_pathology = Motile.obs['tau_pathology']     
APOE_genotype1 = Motile.obs['APOE_genotype[E4/E4]']   
TREM2_genotype = Motile.obs['TREM2_genotype']
sample_ID = Motile.obs['sample_ID']
APOE_genotype2 = Motile.obs['APOE_genotype[E3/E4]']  

# Creating a dictionary with your data
data = {
    'pseudotime': pseudotime_values,
    'amyloid_pathology': amyloid_pathology,
    'tau_pathology': tau_pathology,
    'APOE_genotype[E4/E4]': APOE_genotype1,
    'APOE_genotype[E3/E4]': APOE_genotype2,
    'TREM2_genotype': TREM2_genotype,
    'sample_ID': sample_ID
}


# Creating a pandas DataFrame from the dictionary
df = pd.DataFrame(data)
df['amyloid_pathology'] = df['amyloid_pathology'].astype(int)
df['tau_pathology'] = df['tau_pathology'].astype(int)
df['APOE_genotype[E3/E4]'] = df['APOE_genotype[E3/E4]'].astype(int)
df['APOE_genotype[E4/E4]'] = df['APOE_genotype[E4/E4]'].astype(int)
df['TREM2_genotype'] = df['TREM2_genotype'].astype(int)
df['sample_ID'] = df['sample_ID'].astype(int)
df['pseudotime'] = df['pseudotime']

df.to_excel('Motile_with_APOEgenotype.xlsx', sheet_name='new_sheet_name')

# I will use these dataframes for linear mixed effect modeling 

# Now linear regression analysis for ARM


# List to store regression results for each gene
regression_results = []

# 'dpt_pseudotime' is the pseudotime column in ARM.obs
pseudotime_values = ARM.obs['dpt_pseudotime'].values

# ARM.X is the matrix of gene expression
gene_expression_matrix = ARM.X.A # Convert to dense array

# Creating a DataFrame with gene expression and pseudotime
data_frame = pd.DataFrame(gene_expression_matrix, columns=ARM.var_names)
data_frame['pseudotime'] = pseudotime_values

# for each gene
for gene_column in data_frame.columns[:-1]:  # Exclude 'dpt_pseudotime'
    gene_expression_values = data_frame[gene_column]

    # Performing linear regression with stats.linregress
    slope, intercept, r_value, p_value, std_err = stats.linregress(pseudotime_values, gene_expression_values)

    # Apply Bonferroni correction to p-value
    p_value_adjusted = multi.multipletests([p_value], method='bonferroni')[1][0]

    # Calculate effect size which is beta and R^2
    effect_size = slope
    r_squared = r_value**2

    # Store results in a dictionary
    result_dict = {
        'gene': gene_column,
        'slope': slope,
        'intercept': intercept,
        'r_value': r_value,
        'p_value': p_value,
        'p_value_adjusted': p_value_adjusted,
        'std_err': std_err,
        'effect_size': effect_size,
        'r_squared': r_squared
    }

    regression_results.append(result_dict)

# Converting the list of dictionaries to a DataFrame
regression_results_df = pd.DataFrame(regression_results)

# Displaingy the regression results
print(regression_results_df)

# accorsding to threshold determining significantly changed genes above pseudotime
# very bad threshold for r_squared unfortunately

filtered_genes = regression_results_df[
    (abs(regression_results_df['slope']) > 0.1) &
    (regression_results_df['r_squared'] > 0.009) &
    (regression_results_df['p_value_adjusted'] < 0.01)
]

# Filter ypregulated genes (positive slope)
upregulated_genes = regression_results_df[
    (regression_results_df['slope'] > 0.1) &
    (regression_results_df['r_squared'] > 0.009) &
    (regression_results_df['p_value_adjusted'] < 0.01)
]

# Filter downregulated genes (negative slope)
downregulated_genes = regression_results_df[
    (regression_results_df['slope'] < -0.1) &
    (regression_results_df['r_squared'] > 0.009) &
    (regression_results_df['p_value_adjusted'] < 0.01)
]

# Display the results
print("Upregulated Genes:")
print(upregulated_genes)

print("\nDownregulated Genes:")
print(downregulated_genes)


# enrichment analysis with gprofiler
gene_names = filtered_genes['gene'].tolist()


from gprofiler import GProfiler

def perform_enrichment_analysis(gene_symbols):
    gp = GProfiler(return_dataframe=True)
    enrichment_results = gp.profile(organism='hsapiens', query=gene_symbols, sources=['GO:BP', 'GO:MF' , 'KEGG'], user_threshold=0.05)

    return enrichment_results

# gene symbols to anlaysis enrichment
gene_symbols = gene_names
enrichment_results = perform_enrichment_analysis(gene_symbols)
print(enrichment_results)


# enrichment results use it to visualize
plt.figure(figsize=(12, 8))
sns.heatmap(enrichment_results.pivot_table(index='name', columns='source', values='p_value'),
            cmap="YlGnBu", linewidths=.5, cbar_kws={'label': 'P-Value'})

# Adding labels and title
plt.title('Enrichment Results Heatmap')
plt.xlabel('Source')
plt.ylabel('Enriched Terms')

plt.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#finding top enrichment results significant
sorted_enrichments = enrichment_results.sort_values(by='p_value', ascending=True)
top_enrichments = sorted_enrichments.head(20)  # Select the top 10 entries, you can adjust this number

plt.figure(figsize=(12, 8))
sns.heatmap(top_enrichments.pivot_table(index='name', columns='source', values='p_value'),
            cmap="YlGnBu", linewidths=.5, cbar_kws={'label': 'P-Value'})

# Add labels and title
plt.title('Top Enrichment Results Heatmap')
plt.xlabel('Source')
plt.ylabel('Enriched Terms')

plt.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# other way to visualize
top_enrichments = enrichment_results.loc[enrichment_results['p_value'] <= 0.01]

plt.figure(figsize=(14, 10))
heatmap = sns.heatmap(top_enrichments.pivot_table(index='name', columns='source', values='p_value'),
                      cmap="coolwarm", linewidths=.5, cbar_kws={'label': 'P-Value'},
                      annot=False, fmt=".4g", square=True, cbar=True, vmin=0, vmax=0.01,
                      annot_kws={"size": 8}, robust=True)

# Adjust cell width
plt.setp(heatmap.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(heatmap.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor")

# Add labels and title
plt.title('Top Enrichment Results Heatmap')
plt.xlabel('Source')
plt.ylabel('Enriched Terms')

plt.show()


from wordcloud import WordCloud
import matplotlib.pyplot as plt


# it is for just visualize different way
# Combine the gene names into a single string for each category
upregulated_text = ' '.join(upregulated_genes['gene'])
downregulated_text = ' '.join(downregulated_genes['gene'])

# Create word clouds
up_wordcloud = WordCloud(width=800, height=400, background_color='white').generate(upregulated_text)
down_wordcloud = WordCloud(width=800, height=400, background_color='white').generate(downregulated_text)

# Plot the word clouds
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.imshow(up_wordcloud, interpolation='bilinear')
plt.title('Upregulated Genes Word Cloud')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(down_wordcloud, interpolation='bilinear')
plt.title('Downregulated Genes Word Cloud')
plt.axis('off')

plt.show()




import matplotlib.pyplot as plt

# this is for MA plot I tried to create 

# Set the threshold for distinguishing genes
threshold = 0.01 # to line but doent needed 

# Set the maximum limit for R-squared
max_r_squared = 0.04 # I dont want to see line therefore it is 0.04

# Create a figure with subplots
fig, ax = plt.subplots(figsize=(12, 8))

# Plot all genes with bulbs (markers)
scatter_all = ax.scatter(regression_results_df['r_squared'], regression_results_df['slope'],
                         color='grey', label='Other Genes', s=15, alpha=0.5)  # Decrease marker size

# Plot upregulated genes with gene names
scatter_upregulated = ax.scatter(upregulated_genes['r_squared'], upregulated_genes['slope'],
                                 color='green', label='_nolegend_', s=15, alpha=0.7)  # Decrease marker size
for i, txt in enumerate(upregulated_genes['gene']):
    ax.annotate(txt, (upregulated_genes.iloc[i]['r_squared'], upregulated_genes.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Highlight significant upregulated genes
significant_upregulated = upregulated_genes[upregulated_genes['p_value_adjusted'] < threshold]
scatter_significant_upregulated = ax.scatter(significant_upregulated['r_squared'], significant_upregulated['slope'],
                                             color='orange', marker='*', s=30, label='Significant Upregulated Genes',
                                             alpha=0.7)  # Decrease marker size
for i, txt in enumerate(significant_upregulated['gene']):
    ax.annotate(txt, (significant_upregulated.iloc[i]['r_squared'], significant_upregulated.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Plot downregulated genes with gene names
scatter_downregulated = ax.scatter(downregulated_genes['r_squared'], downregulated_genes['slope'],
                                   color='red', label='_nolegend_', s=15, alpha=0.7)  # Decrease marker size
for i, txt in enumerate(downregulated_genes['gene']):
    ax.annotate(txt, (downregulated_genes.iloc[i]['r_squared'], downregulated_genes.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Highlight significant downregulated genes
significant_downregulated = downregulated_genes[downregulated_genes['p_value_adjusted'] < threshold]
scatter_significant_downregulated = ax.scatter(significant_downregulated['r_squared'], significant_downregulated['slope'],
                                               color='blue', marker='*', s=30, label='Significant Downregulated Genes',
                                               alpha=0.7)  # Decrease marker size
for i, txt in enumerate(significant_downregulated['gene']):
    ax.annotate(txt, (significant_downregulated.iloc[i]['r_squared'], significant_downregulated.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Set y and x axis limits based on data excluding max values
y_min, y_max = -1, 2
x_min, x_max = 0, max_r_squared

# Add some padding to the limits
y_padding = (y_max - y_min) * 0.1
x_padding = (x_max - x_min) * 0.1
y_min -= y_padding
y_max += y_padding
x_min -= x_padding
x_max += x_padding

ax.set_ylim(y_min, y_max)
ax.set_xlim(x_min, x_max)

# Add lines for beta and R-squared thresholds
ax.axhline(y=0.1, color='black', linestyle='--', label='Beta Threshold (0.1)')
ax.axhline(y=-0.1, color='black', linestyle='--', label='Beta Threshold (-0.1)')
ax.axvline(x=0.09, color='red', linestyle='--', label='R-squared Threshold (0.009)')

# Add labels and legend
ax.set_ylabel('Beta (slope)')
ax.set_xlabel('R-squared')
ax.set_title('Gene Expression Analysis')
ax.legend()
ax.grid(True)

plt.tight_layout()  # Adjust layout to eliminate empty spaces
plt.show()




# in other way whihc means adata includes 4 microglia subtype, and I performed PAGA all of them


# Step 1: Compute Neighborhood Graph
sc.pp.neighbors(adata, n_neighbors=15, method='umap', metric='cosine')

# Step 2: Run PAGA
sc.tl.paga(adata)

# Step 3: Generate Graph and Plot PAGA
sc.pl.paga(adata)  # Adjust the threshold as needed

# Step 4: Run Diffusion Map
sc.tl.diffmap(adata)

# Step 5: Update Neighborhood Graph using Diffusion Map representation
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap')

# Step 6: Draw Graph
sc.tl.draw_graph(adata)

# Step 7: Plot Graph with color annotations
sc.pl.draw_graph(adata, color='Subpopulation_Group', legend_loc='on data')

# Step 8: Run PAGA with specific groups
sc.tl.paga(adata, groups='Subpopulation_Group')

# Step 9: Plot PAGA with color annotations
sc.pl.paga(adata, color=['Subpopulation_Group'], title=['Cell type'])

# Step 10: Draw Graph using PAGA initialization
sc.tl.draw_graph(adata, init_pos='paga')

# Step 11: Compare PAGA graphs
sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)

# Step 12: Set root for diffusion pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['Subpopulation_Group'] == 'MicrogliaSub1')[0]

# Step 13: Run diffusion pseudotime
sc.tl.dpt(adata)

# Step 14: Plot Graph with pseudotime and color annotations
sc.pl.draw_graph(adata, color=['Subpopulation_Group', 'dpt_pseudotime'], legend_loc='on data', size=20, alpha=0.8)

# Step 15: Plot PAGA with pseudotime annotations
sc.pl.paga(adata, color=['dpt_pseudotime'], threshold=0.1, node_size_scale=1.5, node_size_power=0.5)


# filtering ARM again to linear regression
ARM = adata[(adata.obs['Subpopulation_Group'] == 'MicrogliaSub1'), :]

# List to store regression results for each gene
regression_results = []

# 'dpt_pseudotime' is the pseudotime column in ARM.obs
pseudotime_values = ARM.obs['dpt_pseudotime'].values

# ARM.X is the matrix of gene expression
gene_expression_matrix = ARM.X.A # Convert to dense array

# Creating a DataFrame with gene expression and pseudotime
data_frame = pd.DataFrame(gene_expression_matrix, columns=ARM.var_names)
data_frame['pseudotime'] = pseudotime_values

# for each gene
for gene_column in data_frame.columns[:-1]:  # Exclude 'dpt_pseudotime'
    gene_expression_values = data_frame[gene_column]

    # Performing linear regression with stats.linregress
    slope, intercept, r_value, p_value, std_err = stats.linregress(pseudotime_values, gene_expression_values)

    # Apply Bonferroni correction to p-value
    p_value_adjusted = multi.multipletests([p_value], method='bonferroni')[1][0]

    # Calculate effect size which is beta and R^2
    effect_size = slope
    r_squared = r_value**2

    # Store results in a dictionary
    result_dict = {
        'gene': gene_column,
        'slope': slope,
        'intercept': intercept,
        'r_value': r_value,
        'p_value': p_value,
        'p_value_adjusted': p_value_adjusted,
        'std_err': std_err,
        'effect_size': effect_size,
        'r_squared': r_squared
    }

    regression_results.append(result_dict)

# Converting the list of dictionaries to a DataFrame
regression_results_df = pd.DataFrame(regression_results)

# Displaingy the regression results
print(regression_results_df)

# accorsding to threshold determining significantly changed genes above pseudotime
# very bad threshold for r_squared unfortunately

filtered_genes = regression_results_df[
    (abs(regression_results_df['slope']) > 0.1) &
    (regression_results_df['r_squared'] > 0.009) &
    (regression_results_df['p_value_adjusted'] < 0.01)
]

# Filter ypregulated genes (positive slope)
upregulated_genes = regression_results_df[
    (regression_results_df['slope'] > 0.1) &
    (regression_results_df['r_squared'] > 0.009) &
    (regression_results_df['p_value_adjusted'] < 0.01)
]

# Filter downregulated genes (negative slope)
downregulated_genes = regression_results_df[
    (regression_results_df['slope'] < -0.1) &
    (regression_results_df['r_squared'] > 0.009) &
    (regression_results_df['p_value_adjusted'] < 0.01)
]

# Display the results
print("Upregulated Genes:")
print(upregulated_genes)

print("\nDownregulated Genes:")
print(downregulated_genes)


# again enrichment analysis
gene_names = filtered_genes['gene'].tolist()


from gprofiler import GProfiler

def perform_enrichment_analysis(gene_symbols):
    gp = GProfiler(return_dataframe=True)
    enrichment_results = gp.profile(organism='hsapiens', query=gene_symbols, sources=['GO:BP', 'GO:MF' , 'KEGG'], user_threshold=0.05)

    return enrichment_results

# Example usage
gene_symbols = gene_names
enrichment_results = perform_enrichment_analysis(gene_symbols)

# Display the results
print(enrichment_results)

# visualizetion
plt.figure(figsize=(12, 8))
sns.heatmap(enrichment_results.pivot_table(index='name', columns='source', values='p_value'),
            cmap="YlGnBu", linewidths=.5, cbar_kws={'label': 'P-Value'})

# Add labels and title
plt.title('Enrichment Results Heatmap')
plt.xlabel('Source')
plt.ylabel('Enriched Terms')

plt.show()


sorted_enrichments = enrichment_results.sort_values(by='p_value', ascending=True)
top_enrichments = sorted_enrichments.head(30)  # Select the top 10 entries, you can adjust this number

plt.figure(figsize=(12, 8))
sns.heatmap(top_enrichments.pivot_table(index='name', columns='source', values='p_value'),
            cmap="YlGnBu", linewidths=.5, cbar_kws={'label': 'P-Value'})

# Add labels and title
plt.title('Top Enrichment Results Heatmap')
plt.xlabel('Source')
plt.ylabel('Enriched Terms')

plt.show()


# other visualization technique
top_enrichments = enrichment_results.loc[enrichment_results['p_value'] <= 0.01]

plt.figure(figsize=(20, 10))
heatmap = sns.heatmap(top_enrichments.pivot_table(index='name', columns='source', values='p_value'),
                      cmap="coolwarm", linewidths=.5, cbar_kws={'label': 'P-Value'},
                      annot=False, fmt=".4g", square=True, cbar=True, vmin=0, vmax=0.01,
                      annot_kws={"size": 8}, robust=True)

# Adjust cell width
plt.setp(heatmap.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
plt.setp(heatmap.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor", fontsize=10)

# Set explicit x-axis labels
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8)

# Add labels and title
plt.title('Top Enrichment Results Heatmap')
plt.xlabel('Source')
plt.ylabel('Enriched Terms')

plt.show()



from wordcloud import WordCloud
import matplotlib.pyplot as plt


# Combine the gene names into a single string for each category
upregulated_text = ' '.join(upregulated_genes['gene'])
downregulated_text = ' '.join(downregulated_genes['gene'])

# Create word clouds
up_wordcloud = WordCloud(width=800, height=400, background_color='white').generate(upregulated_text)
down_wordcloud = WordCloud(width=800, height=400, background_color='white').generate(downregulated_text)

# Plot the word clouds
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.imshow(up_wordcloud, interpolation='bilinear')
plt.title('Upregulated Genes Word Cloud')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(down_wordcloud, interpolation='bilinear')
plt.title('Downregulated Genes Word Cloud')
plt.axis('off')

plt.show()




import matplotlib.pyplot as plt

# Set the threshold for distinguishing genes
threshold = 0.1

# Set the maximum limit for R-squared
max_r_squared = 0.04

# Create a figure with subplots
fig, ax = plt.subplots(figsize=(12, 8))

# Plot all genes with bulbs (markers)
scatter_all = ax.scatter(regression_results_df['r_squared'], regression_results_df['slope'],
                         color='grey', label='Other Genes', s=15, alpha=0.5)  # Decrease marker size

# Plot upregulated genes with gene names
scatter_upregulated = ax.scatter(upregulated_genes['r_squared'], upregulated_genes['slope'],
                                 color='green', label='_nolegend_', s=15, alpha=0.7)  # Decrease marker size
for i, txt in enumerate(upregulated_genes['gene']):
    ax.annotate(txt, (upregulated_genes.iloc[i]['r_squared'], upregulated_genes.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Highlight significant upregulated genes
significant_upregulated = upregulated_genes[upregulated_genes['p_value_adjusted'] < threshold]
scatter_significant_upregulated = ax.scatter(significant_upregulated['r_squared'], significant_upregulated['slope'],
                                             color='orange', marker='*', s=30, label='Significant Upregulated Genes',
                                             alpha=0.7)  # Decrease marker size
for i, txt in enumerate(significant_upregulated['gene']):
    ax.annotate(txt, (significant_upregulated.iloc[i]['r_squared'], significant_upregulated.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Plot downregulated genes with gene names
scatter_downregulated = ax.scatter(downregulated_genes['r_squared'], downregulated_genes['slope'],
                                   color='red', label='_nolegend_', s=15, alpha=0.7)  # Decrease marker size
for i, txt in enumerate(downregulated_genes['gene']):
    ax.annotate(txt, (downregulated_genes.iloc[i]['r_squared'], downregulated_genes.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Highlight significant downregulated genes
significant_downregulated = downregulated_genes[downregulated_genes['p_value_adjusted'] < threshold]
scatter_significant_downregulated = ax.scatter(significant_downregulated['r_squared'], significant_downregulated['slope'],
                                               color='blue', marker='*', s=30, label='Significant Downregulated Genes',
                                               alpha=0.7)  # Decrease marker size
for i, txt in enumerate(significant_downregulated['gene']):
    ax.annotate(txt, (significant_downregulated.iloc[i]['r_squared'], significant_downregulated.iloc[i]['slope']),
                fontsize=8, alpha=0.7)

# Set y and x axis limits based on data excluding max values
y_min, y_max = -1, 2
x_min, x_max = 0, max_r_squared

# Add some padding to the limits
y_padding = (y_max - y_min) * 0.1
x_padding = (x_max - x_min) * 0.1
y_min -= y_padding
y_max += y_padding
x_min -= x_padding
x_max += x_padding

ax.set_ylim(y_min, y_max)
ax.set_xlim(x_min, x_max)

# Add lines for beta and R-squared thresholds
ax.axhline(y=0.1, color='black', linestyle='--', label='Beta Threshold (0.1)')
ax.axhline(y=-0.1, color='black', linestyle='--', label='Beta Threshold (-0.1)')
ax.axvline(x=0.09, color='red', linestyle='--', label='R-squared Threshold (0.009)')

# Add labels and legend
ax.set_ylabel('Beta (slope)')
ax.set_xlabel('R-squared')
ax.set_title('Gene Expression Analysis')
ax.legend()
ax.grid(True)

plt.tight_layout()  # Adjust layout to eliminate empty spaces
plt.show()

