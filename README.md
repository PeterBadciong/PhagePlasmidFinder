# PhagePlasmidFinder
A tool to take complete and partial genomes and searches for Phage Plasmid Hybrids using Genomad and HMMscan, and extracts them for further analysis.

Phage plasmids are hybridized mobile genetic elements found within bacterial genomes. This tool is used to discover novel phage plasmids from inputed genomes. By using a combined search of proviruses and plasmids obtained from genomad, and annotation of phage/plasmid like genes using HMMscan and a currated selection of HMMs of the opposite mobile genetic element. The scaffolds and contigs containing potenial phage plasmids are extracted for later analysis  

Usage

Installation:
Mamba installation: 

Download the PPF.py, PhagePlasmidFinder.py, and PlasmidBypass.py scripts

Execution

Running the PPF
  The PPF uses 3 scripts, PPF.py is the input wrapper script, while PlasmidBypass.py and PhagePlasmidFinder.py are the scripts that execute genomad and hmmscan, along with the parsing out the data and extracting the scaffolds
  To execute the PPF, required minimum for execution is below
  python3 PPF.py (Input.fasta) (Path/to/genomad_db/) -o (OutputFolder/) -j (Path/to/PhageProteins.hmm) -l (Path/to/PlasmidProteins.hmm) 

  The following inputs are optional commands for controlling the strictness of parameters
  -h, --help
    Opens the help menu
  -s, --splits
    Determines number of splits for genomad (default 8)
  -t, --output_dir
    Determines number of threads for genomad (default 10)
  -e, --evalue_cutoff
    Set E-value cutoff for hmmscan (default 1e-5)
  -g, --gene_min
    
  -p, --percent_min
  -m, --plasmid_threshold
  -c, --combined_threshold
  -x, --extract_toggle
  
  
