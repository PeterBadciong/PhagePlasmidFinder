# PhagePlasmidFinder
A tool to take complete and partial genomes and searches for Phage Plasmid Hybrids using Genomad and HMMscan, and extracts them for further analysis.

Phage plasmids are hybridized mobile genetic elements found within bacterial genomes. This tool is used to discover novel phage plasmids from inputed genomes. By using a combined search of proviruses and plasmids obtained from genomad, and annotation of phage/plasmid like genes using HMMscan and a currated selection of HMMs of the opposite mobile genetic element. The scaffolds and contigs containing potenial phage plasmids are extracted for later analysis  

Usage

Installation:

wget https://github.com/PeterBadciong/PhagePlasmidFinder/archive/refs/heads/PPF.zip

#Downloads the python scripts and database

unzip PPF.zip

#Unzips wget

cd PhagePlasmidFinder-PPF

#Changes Directory to where main script is present

unzip hmm_files/PhageProteins.hmm.zip -d hmm_files/

#Unzips compressed Phage hmms

rm hmm_files/PhageProteins.hmm.zip

#removes old zip file

hmmpress hmm_files/PlasmidProteins.hmm

hmmpress hmm_files/PhageProteins.hmm

#hmmpress the phage and plasmid hmms

conda env create -f PPF.yml

#Creates a conda enviroment for genomad dependancy

conda activate PPF

#Activates the conda enviroment

Execution

Running the PPF
  The PPF uses 3 scripts, PPF.py is the input wrapper script, while PlasmidBypass.py and PhagePlasmidFinder.py are the scripts that execute genomad and hmmscan, along with the parsing out the data and extracting the scaffolds
  To execute the PPF, required minimum for execution is below
  python3 PPF.py (Input.fasta) (Path/to/genomad_db/) -o (OutputFolder/) -j (Path/to/PhageProteins.hmm) -l (Path/to/PlasmidProteins.hmm) 

  The following inputs are optional commands for controlling the strictness of parameters
  
  -h, --help:
    Opens the help menu
    
  -s, --splits:
    Determines number of splits for genomad (default 8)
    
  -t, --output_dir:
    Determines number of threads for genomad (default 10)
  
  -e, --evalue_cutoff:
    Set E-value cutoff for hmmscan (default 1e-5)
  
  -g, --gene_min:
    Minimum amount of genes for a phage plasmid to be identified (default 15)
    
  -p, --percent_min:
    Minimum percent crossover of phages and plasmids for a phage plasmid to be identified (default 0.15)
    
  -m, --plasmid_threshold:
    Minimum plasmid_score needed to be have an HMMscan run (default 0.8)
  
  -c, --combined_threshold:
    Minimum plasmid_score + phage_score sum to have an HMMscan run (default 0.85)
  
  -x, --extract_toggle:
