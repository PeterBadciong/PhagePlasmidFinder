# PhagePlasmidFinder: Tool for Identification of Phage Plasmid Hybrids

## Introduction
Phage plasmids are hybridized mobile genetic elements found within bacterial genomes. This tool is used to discover novel phage plasmids from inputed genomes. By using a combined search of prophages and plasmids obtained from genomad, and annotation of phage/plasmid like genes using HMMscan and a currated selection of HMMs of the opposite mobile genetic element. The scaffolds and contigs containing potenial phage plasmids are extracted for further analysis. 

## Installation: 
### Conda + Setup
```
wget https://github.com/PeterBadciong/PhagePlasmidFinder/archive/refs/heads/PPF.zip
unzip PPF.zip
rm -r PPF.zip
cd PhagePlasmidFinder-PPF
unzip hmm_files/PhageProteins.hmm.zip -d hmm_files/
rm hmm_files/PhageProteins.hmm.zip
hmmpress hmm_files/PlasmidProteins.hmm
hmmpress hmm_files/PhageProteins.hmm
conda env create -f PPF.yml
conda activate PPF
```
If you dont have a genomad_db directory, you can download it using 
```
genomad download-database .
```
## Execution

### Running the PPF
  The PPF uses 3 scripts, PhagePlasmidFinder.py is the input wrapper script, while PPF1.py and PPF2.py are the scripts that execute genomad and hmmscan, along with the parsing 
  out the data and extracting the scaffolds

### Required Command Line Inputs
```
 -i, --input_fasta           Input fasta file in .fna format
 -g, --genomad-db            Path to the genomad_db
 -o, --output_folder         Name of folder for results
 -j, --phage_proteins        Path to Phage HMMs
 -l, --plasmid_proteins      Path to Plasmid HMMs
 ```
### Example Required Command Line Input
```
  python3 PhagePlasmidFinderDX.py (Input.fasta) (Path/to/genomad_db/) -o (OutputFolder/) -j (Path/to/PhageProteins.hmm) -l (Path/to/PlasmidProteins.hmm) 
```
### Optional Command Line Inputs
  The following inputs are optional commands for controlling the strictness of parameters
```   
  -h, --help                  Opens the help menu
  -s, --splits                Determines number of splits for genomad (default 8)
  -t, --threads               Determines number of threads for genomad (default 10)
  -e, --evalue_cutoff         Set E-value cutoff for hmmscan (default 1e-5)
  -g, --gene_min              Minimum amount of genes for a phage plasmid to be identified (default 10)
  -p, --percent_min           Minimum percent crossover of phages and plasmids for a phage plasmid to be identified (default 0.15)
  -m, --plasmid_threshold     Minimum plasmid_score needed to be have an HMMscan run (default 0.05)
  -c, --phage_threshold       Minimum phage_score sum to have an HMMscan run (default 0.20)
  -x, --extract_toggle        Toggles extraction of scaffolds
```
### Test Run of Phage Plasmid Finder
Run the following command on the provided .fna file
```
python3 PhagePlasmidFinderDX.py Tritonibacter_mobilis_A3R06.fna genomad_db -o Tritonibacter_mobilis_Output -j hmm_files/PhageProteins.hmm -l hmm_files/PlasmidProteins.hmm -s 8 -t 30 -e 1e-5 -p .15 -g 15 -m .05 -c .20 
```
## Output

| Output Directory | Output File | Description |
| --- | --- | --- |
| Main | PlasmidHits.csv | .csv containing the a predicted phage plasmid scaffold, predicted number of genes, percentage of MGE crossover, and fasta description from detected genomad plasmids |
| Main | ProphageHits.csv | .csv containing the a predicted phage plasmid scaffold, predicted number of genes, percentage of MGE crossover, and fasta description from detected genomad prophages |
| Main | error_log.txt | Error log |
| Extras | phage_hmmscan_output.tbl | HMMscan of genomad predicted phages Extras against plasmid HMMs |
| Extras | plasmid_hmmscan_output.tbl | HMMscan of genomad predicted plasmid Extras against phage HMMs |
| Extras | [fasta].phage.csv | Results and overall crossover of genomad predicted phages against plasmid HMMs |
| Extras | [fasta].plasmid.csv | Results and overall crossover of genomad predicted plasmids against phage HMMs |
| Extras | [fasta].Plasmids.fasta | Plasmids ID'd using genomad's nn_classification to be run against phage HMMs |
| Extras/prodigal | PlasmidProdigal.faa | Prodigal output used for finding phage plasmids from genomad predicted plasmids |
| Extras/prodigal | PlasmidProdigal.gff | Prodigal output used for finding phage plasmids from genomad predicted plasmids |
| Extras/prodigal | Phage_Plasmid_hmmscan.tblout | HMMscan of genomad predicted plasmids against phage HMMs |
| Extras/genomad_output | genomad_outputs | Standard genomad outputs |
| PhagePlasmidFasta | [scaffold].fasta | Extracted .fasta files of each scaffold from the PhagePlasmids.csv file |
