import argparse
import logging
import os
import pandas as pd
from Bio import Phylo
import subprocess

def setup_logging(log_file_path):
    """Set up logging configuration to log both to console and file."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='a'),
            logging.StreamHandler()
        ]
    )

def print_pretty(header, message, border_char='-', border_length=60):
    """Print a message with a header and border."""
    border = border_char * border_length
    print(f"\n{border}")
    print(f"{header.center(border_length)}")
    print(f"{message}\n")

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Process a phylogeny in Newick format with metadata.')

    # Add arguments for input files and output path
    parser.add_argument('-n', '--newick_file', type=str, required=True,
                        help='Path to the Newick format phylogeny file.')

    parser.add_argument('-o', '--output_path', type=str, required=True,
                        help='Directory path where the output will be written.')

    parser.add_argument('-m', '--metadata_csv_path', type=str, required=True,
                        help='Path to the metadata CSV file.')

    parser.add_argument('-c', '--id_column_header', type=str, required=True,
                        help='The column header in the metadata CSV file that contains the identifiers.')

    # Add optional arguments for clustering
    parser.add_argument('--cluster', action='store_true',
                        help='Flag to indicate if clustering using TreeCluster should be executed.')

    parser.add_argument('--clustering_method', type=str, default='max_clade',
                        choices=['avg_clade', 'leaf_dist_avg', 'leaf_dist_max', 'leaf_dist_min', 'length', 
                                 'length_clade', 'max', 'max_clade', 'med_clade', 'root_dist', 'single_linkage', 
                                 'single_linkage_cut', 'single_linkage_union', 'sum_branch', 'sum_branch_clade'],
                        help='Clustering method to be used (default: max_clade).')

    parser.add_argument('-t', '--threshold', type=float, required='--cluster' in parser.parse_known_args()[0],
                        help='Threshold for the clustering method, indicating the maximum phylogenetic distance between two members of the same cluster.')

    parser.add_argument('-s', '--support', type=float, default=None,
                        help='Support for the clustering method, indicating the minimum branch support for a cluster to be considered.')

    # Add argument for existing clustering file
    parser.add_argument('--clustering_file', type=str, default=None,
                        help='Path to an existing clustering file to be used instead of performing clustering.')

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed arguments
    newick_file = args.newick_file
    output_path = args.output_path
    metadata_csv_path = args.metadata_csv_path
    id_column_header = args.id_column_header
    clustering_method = args.clustering_method
    threshold = args.threshold
    support = args.support
    clustering_file = args.clustering_file

    # Set up logging to file and console
    input_filename = os.path.basename(newick_file)
    file_root, _ = os.path.splitext(input_filename)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    log_file_path = os.path.join(output_path, f'{file_root}_clustering.log')
    setup_logging(log_file_path)

    # Log the input variables
    logging.info('Newick file: %s', newick_file)
    logging.info('Output path: %s', output_path)
    logging.info('Metadata CSV path: %s', metadata_csv_path)
    logging.info('ID Column Header: %s', id_column_header)
    logging.info('Clustering Method: %s', clustering_method)
    logging.info('Threshold: %f', threshold)
    if support is not None:
        logging.info('Support: %f', support)
    if clustering_file is not None:
        logging.info('Clustering file provided: %s', clustering_file)

    # Read the Newick file
    logging.info('Reading the Newick file...')
    tree = Phylo.read(newick_file, 'newick')

    # Read the metadata CSV
    logging.info('Reading the metadata CSV file...')
    metadata = pd.read_csv(metadata_csv_path)

    # Validate the presence of the ID column in the metadata
    if id_column_header not in metadata.columns:
        logging.error('The specified ID column header is not in the metadata CSV.')
        raise ValueError('The specified ID column header is not in the metadata CSV.')

    # Process the phylogeny and metadata (Example operation)
    matched_identifiers = metadata[id_column_header].tolist()
    logging.info('Matched identifiers count: %d', len(matched_identifiers))

    clustering_df = None  # Initialize the clustering_df variable

    # Perform clustering if the cluster flag is set
    if args.cluster:
        logging.info('Performing clustering with TreeCluster...')
        
        # Construct the output file name for clustering
        cluster_output_file = os.path.join(output_path, f"{file_root}_{clustering_method}_clusters.txt")

        # Construct the TreeCluster command
        treecluster_command = [
            "TreeCluster.py", 
            "-i", newick_file, 
            "-o", cluster_output_file, 
            "-m", clustering_method, 
            "-t", str(threshold)
        ]
        
        # Add support if provided
        if support is not None:
            treecluster_command.extend(["-s", str(support)])

        # Print and execute the TreeCluster command
        logging.info('TreeCluster command: %s', ' '.join(treecluster_command))
        print_pretty("TreeCluster Command", ' '.join(treecluster_command))
        
        try:
            subprocess.run(treecluster_command, check=True)
            logging.info('Clustering completed successfully. Output written to: %s', cluster_output_file)

            # Read the TreeCluster output into a Pandas DataFrame
            clustering_df = pd.read_csv(cluster_output_file, sep='\t')
            logging.info('TreeCluster output read into DataFrame successfully.')
        except subprocess.CalledProcessError as e:
            logging.error('Clustering failed: %s', e)
            print(f"Error: {e}")
    elif clustering_file:
        # If clustering is not performed, but a clustering file is provided
        logging.info('Reading the provided clustering file into DataFrame...')
        clustering_df = pd.read_csv(clustering_file, sep='\t')
        logging.info('Clustering file read into DataFrame successfully.')
    else:
        logging.info('Clustering not performed. Please ensure you use the --cluster flag to execute clustering or provide a clustering file.')
        print_pretty("Clustering Module Message", "Clustering not performed. Use the --cluster flag to execute clustering or provide a clustering file.")

    if clustering_df is not None:
        # Match metadata with clustering results
        def match_identifiers(identifier, clustering_df):
            for idx, row in clustering_df.iterrows():
                if identifier in row['SequenceName']:
                    return row['ClusterNumber']
            return None

        metadata['ClusterNumber'] = metadata[id_column_header].apply(lambda x: match_identifiers(x, clustering_df))

        # Export the combined DataFrame to a CSV file
        combined_output_file = os.path.join(output_path, f'{file_root}_combined_metadata_clustering.csv')
        metadata.to_csv(combined_output_file, index=False)
        logging.info('Combined metadata and clustering results written to: %s', combined_output_file)
        print_pretty("Export Successful", f"Combined metadata and clustering results exported to {combined_output_file}")

    # Output summary of the results
    output_summary = (
        f"Processed Newick file: {newick_file}\n"
        f"Metadata file used: {metadata_csv_path}\n"
        f"Output written to: {output_path}\n"
        f"Log file: {log_file_path}\n"
    )
    print_pretty("Clustering Module Summary", output_summary)

if __name__ == '__main__':
    main()
