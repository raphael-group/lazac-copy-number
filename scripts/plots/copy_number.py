import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
 
def plot_copy_number_heatmap(data_file, nodes_file):
    # Read the data from CSV
    data = pd.read_csv(data_file)

    # Read the nodes to plot from the text file
    with open(nodes_file, 'r') as file:
        nodes = [line.strip() for line in file.readlines()]

    # Filter the data for the specified nodes
    filtered_data = data[data['node'].isin(nodes)]

    # chromosome ordering
    chrom_ordering = {chrm : j for j, chrm in enumerate([str(x + 1) for x in range(22)] + ['X'])}

    # unique nodes and chromosomes
    chromosome_starts = filtered_data[['chrom', 'start']].drop_duplicates()
    chromosome_starts = chromosome_starts.sort_values(by='start').sort_values(by='chrom', key=lambda x: x.map(chrom_ordering), kind='stable')
    nodes = filtered_data['node'].unique()

    # convert to 2D matrix format
    copy_number_matrix = np.zeros((len(nodes), len(chromosome_starts)))
   
    chrom_start_to_index = {}
    index_j = 0
    for chrom, start in chromosome_starts.values:
        chrom_start_to_index[(chrom, start)] = index_j
        index_j += 1


    for i, node in enumerate(nodes):
        node_data = filtered_data[filtered_data['node'] == node]
        for j, row in node_data.iterrows():
            chrom = row['chrom']
            start = row['start']
            index_j = chrom_start_to_index[(chrom, start)]
            copy_number_matrix[i, index_j] = row['cn_a']

    copy_number_matrix = copy_number_matrix.astype(int)

    # Get the xticks
    chromosome_starts['index'] = chromosome_starts.apply(lambda row: chrom_start_to_index[(row['chrom'], row['start'])], axis=1)
    xtick_indices = chromosome_starts.groupby('chrom').median('index').reset_index()[['chrom', 'index']].sort_values(by='index')
    xticks = xtick_indices['index'].values
    xtick_labels = xtick_indices['chrom'].values

    # Create and plot figure
    fig, ax = plt.subplots(figsize=(15, 8))

    copy_number_states = np.unique(copy_number_matrix)
    n = len(copy_number_states)
    cmap_dict = dict(zip(copy_number_states, sns.color_palette("deep", n)))
    cmap = ListedColormap([cmap_dict[x] for x in copy_number_states])
    sns.heatmap(copy_number_matrix, ax=ax, cmap=cmap, yticklabels=nodes, vmin=copy_number_states.min() - 0.5, vmax=copy_number_states.max() + 0.5)

    chromosome_start_indices = chromosome_starts.groupby('chrom').min('index').reset_index()[['chrom', 'index']].sort_values(by='index')
    for _, row in chromosome_start_indices.iterrows():
        if row['chrom'] == '1':
            continue

        ax.axvline(row['index'], color='black', lw=1, ls='--')

    colorbar = ax.collections[0].colorbar
    r = colorbar.vmax - colorbar.vmin
    colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n) + r * i / (n) for i in range(n)])
    colorbar.set_ticklabels(copy_number_states)

    # Set xticks
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)

    fig.tight_layout()
    plt.show()

def plot_copy_number(data_file, nodes_file):
    # Read the data from CSV
    data = pd.read_csv(data_file)

    # Read the nodes to plot from the text file
    with open(nodes_file, 'r') as file:
        nodes = [line.strip() for line in file.readlines()]

    # Filter the data for the specified nodes
    filtered_data = data[data['node'].isin(nodes)]

    # Get unique nodes and chromosomes
    unique_nodes = filtered_data['node'].unique()
    unique_chromosomes = filtered_data['chrom'].unique()

    # Calculate the number of rows and columns for subplots
    nrows = len(unique_nodes)
    ncols = len(unique_chromosomes)


    fig = plt.figure(figsize=(22, 10))
    gs = gridspec.GridSpec(nrows, ncols, figure=fig, wspace=0, hspace=0.045)
    axes = []

    # Iterate over each node and chromosome combination
    for i, node in enumerate(unique_nodes):
        for j, chromosome in enumerate([str(x + 1) for x in range(21)] + ['X']):
            # Filter data for the specific node and chromosome
            node_chromosome_data = filtered_data[(filtered_data['node'] == node) & (filtered_data['chrom'] == chromosome)]

            # Create a subplot within the GridSpec
            if j == 0:
                ax = fig.add_subplot(gs[i, j])
                axes.append(ax)
            else:
                ax = fig.add_subplot(gs[i, j], sharey=axes[i])

            # Plot copy number as a function of x using seaborn
            sns.lineplot(data=node_chromosome_data, x='start', y='cn_a', ax=ax)

            # Remove the x-axis label for all but the bottom row
            if i < nrows - 1:
                ax.set_xticklabels([])

            # Remove the y-axis label for all but the leftmost column
            if j > 0:
                ax.tick_params(labelleft=False)
                # ax.set_yticklabels([])
                # ax.set_yticks([])

            if chromosome != 'X':
                ax.spines['right'].set_visible(False)

            ax.spines['left'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)

            # Set labels for each subplot
            ax.set_xticks([])
            ax.set_xlabel(chromosome)
            if j == 0:
                ax.set_ylabel("Copy Number")
                ax.set_ylabel("Copy Number")
            else:
                ax.set_ylabel(None)

            # Remove the unnecessary repeated information in the title
            if chromosome == 'X':
                node_name = f'{node[17:]}' 
                plt.text(1.1, 0.5, node_name, rotation=0, va='center', transform=ax.transAxes)

    # Show the plot
    plt.show()

if __name__ == '__main__':
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Plot copy number for each node and chromosome combination.')

    # Add the command-line arguments
    parser.add_argument('data_file', help='Path to the data file')
    parser.add_argument('nodes_file', help='Path to the text file containing the nodes')

    parser.add_argument("--type", help="Type of plot to create", choices=['heatmap', 'line'], default='line')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to plot the copy number
    if args.type == 'line':
        plot_copy_number(args.data_file, args.nodes_file)
    else:
        plot_copy_number_heatmap(args.data_file, args.nodes_file)
