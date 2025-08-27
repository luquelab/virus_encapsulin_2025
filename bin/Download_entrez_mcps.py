from Bio import Entrez, SeqIO
import pandas as pd

# Set your email for NCBI Entrez
Entrez.email = "abelardoacm@ciencias.unam.mx"

# Path to save the results
save_path = "/Users/abelardoaguilar/projects/github_repos/mini-devel/mini-devel/bin/Modules_refactoring/X_isolates_verse/results"

# Search for major capsid proteins in Heunggongvirae within RefSeq
search_term = "Heunggongvirae[Organism] AND (major capsid protein[Title] OR mcp[Title]) AND srcdb_refseq[PROP]"

# Perform search in protein database
handle = Entrez.esearch(db="protein", term=search_term, retmax=5000)
record = Entrez.read(handle)
id_list = record["IdList"]

# Fetch the sequences and metadata
sequences = []
metadata = []
organisms_seen = set()

for protein_id in id_list:
    # Fetch the full record to get metadata and ensure uniqueness by organism
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
    gb_record = SeqIO.read(handle, "genbank")
    
    # Extract relevant metadata
    organism = gb_record.annotations.get('organism', 'Unknown')
    
    if organism in organisms_seen:
        continue  # Skip if we've already seen this organism
    organisms_seen.add(organism)
    
    lineage = gb_record.annotations.get('taxonomy', [])
    family = lineage[-1] if lineage else 'Unknown'
    host = 'Unknown'
    for feature in gb_record.features:
        if feature.type == "source":
            host = feature.qualifiers.get("host", ["Unknown"])[0]
    
    # Fetch the protein sequence in FASTA format
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    sequences.append(seq_record)
    
    # Add metadata to the list, including the protein accession
    metadata.append({
        "Protein ID": protein_id,
        "protein_accession": gb_record.annotations.get('accessions', ['Unknown'])[0],
        "Organism": organism,
        "Lineage": " > ".join(lineage),
        "Family": family,
        "Host": host,
        "Protein Length": len(seq_record.seq),
        "Viral Genome Length": gb_record.annotations.get('contig', 'Unknown')
    })

# Save the sequences to a fasta file
fasta_file_path = f"{save_path}/Heunggongvirae_major_capsid_proteins.fasta"
with open(fasta_file_path, "w") as file:
    SeqIO.write(sequences, file, "fasta")

# Convert metadata list to DataFrame
df = pd.DataFrame(metadata)

# Save DataFrame to a CSV file
csv_file_path = f"{save_path}/Heunggongvirae_major_capsid_proteins_metadata.csv"
df.to_csv(csv_file_path, index=False)

# Display the DataFrame
df.head()

# Extract unique virus names from the "Organism" column
unique_viruses = df['Organism'].unique()

# Function to fetch genome length from NCBI for a given virus name
def fetch_genome_length(virus_name):
    handle = Entrez.esearch(db="nucleotide", term=f"{virus_name}[Organism] AND complete genome[Title]", retmax=1)
    record = Entrez.read(handle)
    if record["IdList"]:
        nucleotide_id = record["IdList"][0]
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
        gb_record = SeqIO.read(handle, "genbank")
        return len(gb_record.seq)
    else:
        return "Unknown"

# Fetch genome lengths and update the DataFrame
genome_lengths = {}
for virus in unique_viruses:
    genome_lengths[virus] = fetch_genome_length(virus)

# Update the DataFrame with the fetched genome lengths
df['Viral Genome Length'] = df['Organism'].map(genome_lengths)

# Save the updated DataFrame to a CSV file
df.to_csv(csv_file_path, index=False)
