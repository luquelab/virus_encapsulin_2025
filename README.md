# virus_encapsulin_2025

This repository contains the code and analyses associated with the manuscript:

**Small viruses reveal multiple evolutionary transitions between HK97-fold viruses and encapsulins** *(under revision)*

Related preprint (bioRxiv, 2025):  
**Small uncultured viral-like entities redraw the origin of viruses and cellular compartments**  
https://www.biorxiv.org/content/10.1101/2025.06.18.659913v1

---

### What this repository contains

This repository provides small, self-contained scripts/apps used to generate and explore analyses shown in the manuscript, particularly those behind plots and summary outputs.

The repository is organized as follows:

- `bin/` contains analysis entrypoints and a utils folder with helper scripts
- `data/` contains input files used with those entrypoints.
- `results/` contains representative outputs produced from the example inputs (for full outputs please access the Figshare datasets indicated in the manuscript).

---

### File-to-file correspondence

Each main analysis has:
1) a script/app in `bin/`  
2) a matching input in `data/`  
3) a matching output folder in `results/`

| Analysis | Entrypoint (bin/) | Input (data/) | Example outputs (results/) |
|---|---|---|---|
| Domain co-occurrence network (Tkinter) | `Co_occurrence_domains_network.tkinter.py` | `Co_occurrence_domains_network_clan0373_data.csv` | `Co_occurrence_domains/` |
| Hexbin: encapsulin fraction vs genome length (Shiny) | `Hexbin_enc_fraction_vs_genome_length.app.R` | `Hexbin_enc_fraction_vs_genome_length_data.csv` | `Hexbin_enc_fraction_vs_genome_length/` |
| Violin plots: phage coat & encapsulin fraction (Shiny) | `Violin_plots_phage_coat_and_enc_fraction.app.R` | `Violin_plots_phage_coat_and_enc_fraction_data.csv` | `Violin_plots_phage_coat_and_enc_fraction/` |
| Statistical comparisons across groups (R script) | `Statistics_phage_coat_and_encapsulin_fraction.R` | `Statistics_phage_coat_and_encapsulin_fraction_data.csv` | `Statistics_phage_coat_and_encapsulin_fraction/` |

Additional included data:
- `data/DTRs_genomes.fna` is a FASTA file used in upstream parts of the broader project.

---

## Requirements

**Python (Tkinter app)**
- Python 3.x
- Tkinter

**R (Shiny apps + statistics)**
- R 4.x
- tidyverse
- RStudio is recommended for running Shiny apps interactively

---

## How to run the analyses

### 1) Co-occurrence network app (Python + Tkinter)

**Purpose**  
Builds a domain co-occurrence network from a binary matrix:

- rows = proteins  
- columns = Pfam domains  
- cell = 1 (domain present) / 0 (domain absent)

Nodes are domains; an edge connects two domains if they co-occur in at least one protein. The app also supports group assignment (e.g., “encapsulin” vs “viral”) and path evaluation with a simple cost function.

**Run**
```bash
python bin/Co_occurrence_domains_network.tkinter.py
```

**Using the app (4 tabs overview)**

1. **Input tab**
   - Select the binary matrix input file
   - Provide the row number where the matrix begins (if the CSV contains header/metadata rows)

2. **Domain grouping tab**
   - Assign domains to groups (e.g., `encapsulin` or `viral`)

3. **Network + export tab**
   - Build the co-occurrence graph
   - Export network

**Example input**
- `data/Co_occurrence_domains_network_clan0373_data.csv`

---

### 2) Hexbin app (R Shiny): encapsulin fraction vs genome length

**Purpose**  
Interactive hexbin plot where:
- x-axis = genome length (source genome of each protein)
- y-axis = encapsulin-fraction coverage (per protein)

Includes a smoothing spline and controls for hexbin/aesthetic parameters.

**Run from terminal**
```bash
R -e "shiny::runApp('bin/Hexbin_enc_fraction_vs_genome_length.app.R', launch.browser=TRUE)"
```

**Example input**
- `data/Hexbin_enc_fraction_vs_genome_length_data.csv`

---

### 3) Violin plots app (R Shiny): phage coat & encapsulin fraction

**Purpose**  
Interactive violin plots comparing groups of proteins using:
- phage coat coverage (%)
- encapsulin fraction coverage (%)

**Run from terminal**
```bash
R -e "shiny::runApp('bin/Violin_plots_phage_coat_and_enc_fraction.app.R', launch.browser=TRUE)"
```

**Example input**
- `data/Violin_plots_phage_coat_and_enc_fraction_data.csv`

---

### 4) Statistics script (R): group comparisons + post hoc tests

**Purpose**  
Runs a statistical testing workflow comparing:
- `phage_coat_coverage_percent`
- `encapsulin-fraction`

The script performs an omnibus test and then post hoc comparisons, producing pairwise p-value tables and summary outputs.

**Run:**
```bash
Rscript bin/Statistics_phage_coat_and_encapsulin_fraction.R data/Statistics_phage_coat_and_encapsulin_fraction_data.csv
```

---

## Helper utilities (bin/utils)

The `bin/utils/` folder contains helper scripts used in upstream or supporting steps (not figure apps themselves). Examples include:
- downloading reference Duplodnaviria MCPs,
- running VIBRANT in-site,
- converting `hmmscan` domtblout outputs to CSV,
- pulling TreeCluster outputs,
- wrappers to standardize running steps.

These scripts were used inside a broader pipeline context and may require environment- and path-specific configurations.

---

## License

See `LICENSE`.