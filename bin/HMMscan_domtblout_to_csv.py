import csv
import argparse

def parse_domtblout(input_file, output_file):
    """
    Parse a domtblout file from hmmscan and write the results to a CSV file.

    Args:
        input_file (str): Path to the input domtblout file.
        output_file (str): Path to the output CSV file.
    """
    # Open the input and output files
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        
        # Write the header to the CSV file
        writer.writerow([
            "target_name", "accession", "tlen", "query_name", "query_accession",
            "qlen", "E-value", "score", "bias", "#", "of", "c-Evalue",
            "i-Evalue", "score", "bias", "hmm_from", "hmm_to", "ali_from",
            "ali_to", "env_from", "env_to", "acc", "description"
        ])
        
        # Iterate through each line in the input file
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue  # Skip comments and empty lines
            
            # Split the line into columns based on whitespace
            columns = line.split()
            
            # Extract relevant columns from the domtblout file
            target_name = columns[0]
            accession = columns[1]
            tlen = columns[2]
            query_name = columns[3]
            query_accession = columns[4]
            qlen = columns[5]
            e_value = columns[6]
            score = columns[7]
            bias = columns[8]
            num = columns[9]
            of = columns[10]
            c_evalue = columns[11]
            i_evalue = columns[12]
            dom_score = columns[13]
            dom_bias = columns[14]
            hmm_from = columns[15]
            hmm_to = columns[16]
            ali_from = columns[17]
            ali_to = columns[18]
            env_from = columns[19]
            env_to = columns[20]
            acc = columns[21]
            description = ' '.join(columns[22:])  # Capture the rest as the description
            
            # Write the row to the CSV file
            writer.writerow([
                target_name, accession, tlen, query_name, query_accession,
                qlen, e_value, score, bias, num, of, c_evalue, i_evalue,
                dom_score, dom_bias, hmm_from, hmm_to, ali_from, ali_to,
                env_from, env_to, acc, description
            ])

def main():
    """
    Main function to handle argument parsing and call the parser function.
    """
    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description='Convert a domtblout file from hmmscan to a CSV file.'
    )
    parser.add_argument('-i', '--input', required=True, help='Path to the input domtblout file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function to parse the domtblout file
    parse_domtblout(args.input, args.output)

# Entry point of the script
if __name__ == '__main__':
    main()
