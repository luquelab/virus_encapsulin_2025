import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from src.auxiliary_functions.VIBRANT_filtering_functions import (
    read_VIBRANT_custom_fasta_23,
    filter_VIBRANT_annotations_by_mcp,
    write_VIBRANT_contigs_w_mcp_fasta,
    filter_and_write_MCP_sequences
)

def setup_logging(log_file_path):
    """Set up logging configuration to log both to console and file."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='w'),
            logging.StreamHandler()
        ]
    )

def print_pretty(header, message, border_char='-', border_length=60):
    """Print a message with a header and border."""
    border = border_char * border_length
    print(f"\n{border}")
    print(f"{header.center(border_length)}")
    print(f"{border}")
    print(f"{message}\n")

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Process input and output paths for VIBRANT analysis.')

    # Add arguments for input paths
    parser.add_argument('-v', '--vibrant_output_path', type=str, default='/Users/abelardoaguilar/projects/github_repos/mini-devel/mini-devel/bin/Modules_refactoring/1_VIBRANT/results/VIBRANT_DTRs_20kb/',
                        help='Path to the VIBRANT output directory.')

    parser.add_argument('-m', '--metadata_csv_path', type=str, default=None,
                        help='Path to the metadata CSV file. If not provided, this will be None.')

    # Add argument for the column header if metadata is given
    parser.add_argument('-c', '--id_column_header', type=str, default=None,
                        help='The column header in the metadata CSV file that contains the ID to link to the VIBRANT output.')

    # Add arguments for parameters
    parser.add_argument('-t', '--mcp_terms_by_usr', nargs='+', default=None,
                        help='List of terms of proteins to be found. Defaults to a predefined list for major capsid proteins.')

    parser.add_argument('-f', '--false_terms_by_usr', nargs='+', default=None,
                        help='List of terms of protein annotations to be excluded. Defaults to a predefined list of false positives for major capsid proteins if None.')

    # Add argument for output path suffix
    parser.add_argument('-s', '--mcps_fasta_suffix', type=str, default='_MCPs.faa',
                        help='Suffix for extracted protein FASTA files. Output paths are created automatically within VIBRANT output path.')

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed arguments
    VIBRANT_output_path = args.vibrant_output_path
    metadata_csv_path = args.metadata_csv_path
    id_column_header = args.id_column_header
    mcp_terms_by_usr = args.mcp_terms_by_usr
    false_terms_by_usr = args.false_terms_by_usr
    mcps_fasta_suffix = args.mcps_fasta_suffix

    # Determine if metadata is provided
    metadata_provided = metadata_csv_path is not None

    # Validate the presence of the ID column header if metadata is provided
    if metadata_provided and not id_column_header:
        parser.error('If metadata is provided, the ID column header must also be specified with the -c flag.')

    # MAIN SCRIPT
    # create the output directory if it does not exist
    output_directory = os.path.join(VIBRANT_output_path, 'post_VIBRANT_filtered_output')

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # Set up logging to file and console
    log_file_path = os.path.join(output_directory, 'vibrant_post_analysis.log')
    setup_logging(log_file_path)

    logging.info('Output directory is ready: %s', output_directory)

    VIBRANT_output_path_elements = VIBRANT_output_path.split('/')
    VIBRANT_folder_name = next((element for element in reversed(VIBRANT_output_path_elements) if element), None)
    VIBRANT_input_name = '_'.join(VIBRANT_folder_name.split('_')[1:])

    # Read VIBRANT output key files and store them in pandas dataframes
    # output 38 is the VIBRANT annotation file
    output38_file_path = os.path.join(VIBRANT_output_path, f'VIBRANT_results_{VIBRANT_input_name}/VIBRANT_annotations_{VIBRANT_input_name}.tsv')
    VIBRANT_annotations = pd.read_csv(output38_file_path, sep='\t')
    # output 42 is the list of predicted genome quality and type
    output42_file_path = os.path.join(VIBRANT_output_path, f'VIBRANT_results_{VIBRANT_input_name}/VIBRANT_genome_quality_{VIBRANT_input_name}.tsv')
    VIBRANT_genome_quality = pd.read_csv(output42_file_path, sep='\t')
    # output 23 corresponds to all encoded proteins among identified phages
    output23_file_path = os.path.join(VIBRANT_output_path, f'VIBRANT_phages_{VIBRANT_input_name}/{VIBRANT_input_name}.phages_combined.faa')
    VIBRANT_phages_proteins = read_VIBRANT_custom_fasta_23(output23_file_path)
    count_VIBRANT_phages_proteins = len(VIBRANT_phages_proteins)
    # output 25 corresponds to all phages genomes
    output25_file_path = os.path.join(VIBRANT_output_path, f'VIBRANT_phages_{VIBRANT_input_name}/{VIBRANT_input_name}.phages_combined.fna')
    VIBRANT_phages_genomes = SeqIO.parse(output25_file_path, "fasta")
    VIBRANT_phages_genomes_as_list = list(SeqIO.parse(output25_file_path, "fasta"))
    count_VIBRANT_phages_genomes = sum(1 for _ in VIBRANT_phages_genomes)

    logging.info('Count of VIBRANT phages genomes: %d', count_VIBRANT_phages_genomes)
    logging.info('Count of VIBRANT phages proteins: %d', count_VIBRANT_phages_proteins)

    # Pretty print the input variables and summary
    input_summary = (
        f"VIBRANT output path: {VIBRANT_output_path}\n"
        f"Metadata CSV path: {metadata_csv_path}\n"
        f"ID Column Header: {id_column_header}\n"
        f"Metadata Provided: {metadata_provided}\n"
        f"MCP terms: {mcp_terms_by_usr}\n"
        f"False MCP terms: {false_terms_by_usr}\n"
        f"FASTA suffix: {mcps_fasta_suffix}\n"
        f"Number of VIBRANT phages genomes: {count_VIBRANT_phages_genomes}\n"
        f"Number of VIBRANT phages proteins: {count_VIBRANT_phages_proteins}\n"
    )
    print_pretty("VIBRANT Post Analysis Configuration", input_summary)

    # Operate using VIBRANT output
    # Filter the VIBRANT annotations DataFrame based on MCP-related terms
    filtered_VIBRANT_contigs_w_mcp = filter_VIBRANT_annotations_by_mcp(VIBRANT_annotations, mcp_terms=mcp_terms_by_usr, false_terms=false_terms_by_usr)
    filtered_VIBRANT_contigs_w_mcp_path = os.path.join(output_directory, f'{VIBRANT_input_name}_contigs_w_mcp.tsv')
    filtered_VIBRANT_contigs_w_mcp.to_csv(filtered_VIBRANT_contigs_w_mcp_path, sep='\t', index=False)
    logging.info('Contigs with MCP DataFrame written to: %s', filtered_VIBRANT_contigs_w_mcp_path)

    # Write filtered VIBRANT contigs to FASTA
    write_VIBRANT_contigs_w_mcp_fasta(filtered_VIBRANT_contigs_w_mcp, output25_file_path, output_directory, VIBRANT_input_name, VIBRANT_phages_proteins)

    # Filter the VIBRANT proteins and write the MCP sequences to a fasta file
    filter_and_write_MCP_sequences(filtered_VIBRANT_contigs_w_mcp, VIBRANT_phages_proteins, output_directory)

    # If metadata is provided, add the number of MCPs per contig and the quality of the contig to the metadata
    if metadata_provided:
        # Load the metadata
        metadata = pd.read_csv(metadata_csv_path)
        # Add columns to the metadata DataFrame to store the number of MCPs per contig
        if 'MCPs_count' in metadata.columns:
            metadata = metadata.drop(columns=['MCPs_count'])
        metadata['MCPs_count'] = metadata[id_column_header].apply(lambda x: filtered_VIBRANT_contigs_w_mcp[filtered_VIBRANT_contigs_w_mcp['scaffold'] == x].shape[0])

        # Add quality information
        if 'quality' in metadata.columns:
            metadata = metadata.drop(columns=['quality'])

        def get_quality(genome_id):
            quality_data = VIBRANT_genome_quality[VIBRANT_genome_quality['scaffold'] == genome_id]['Quality']
            return quality_data.values[0] if not quality_data.empty else None

        metadata['quality'] = metadata[id_column_header].apply(get_quality)

        # Save the updated metadata
        updated_metadata_path = os.path.join(output_directory, f'{VIBRANT_input_name}_updated_metadata.csv')
        metadata.to_csv(updated_metadata_path, index=False)
        logging.info('Updated metadata written to: %s', updated_metadata_path)

    # Output summary of the results
    output_summary = (
        f"Filtered contigs with MCPs TSV: {filtered_VIBRANT_contigs_w_mcp_path}\n"
        f"Filtered MCP FASTA files in: {output_directory}\n"
        f"Updated metadata CSV: {updated_metadata_path if metadata_provided else 'Not provided'}\n"
        f"Log file: {log_file_path}\n"
    )
    print_pretty("VIBRANT Post Analysis Output Summary", output_summary)

if __name__ == '__main__':
    main()
