import os
import csv
import re
import argparse
import subprocess

def run_prodigal(fasta_file, output_dir):
    """
    Runs Prodigal on the given FASTA file to predict proteins.
    """
    output_gff = os.path.join(output_dir, "PlasmidProdigal.gff")
    output_faa = os.path.join(output_dir, "PlasmidProdigal.faa")
    
    command = ["prodigal", "-i", fasta_file, "-o", output_gff, "-a", output_faa, "-f", "gff"]
    
    try:
        subprocess.run(command, check=True)
        print(f"Prodigal completed successfully. GFF output saved to {output_gff}, FAA output saved to {output_faa}.")
    except subprocess.CalledProcessError as e:
        print(f"Prodigal failed with error: {e}")

def extract_input_name(genomad_output):
    for item in os.listdir(genomad_output):
        if item.endswith("_nn_classification"):
            return item.replace("_nn_classification", "")
    raise ValueError("No valid classification folder found in the Genomad output directory.")

def extract_scaffolds(genomad_output, output_folder, fasta_file, plasmid_threshold=0.80, combined_threshold=0.85, verbose=False):
    try:
        input_name = extract_input_name(genomad_output)
        if verbose:
            print(f"Input name detected: {input_name}")
        
        classification_file = os.path.join(genomad_output, f"{input_name}_nn_classification", f"{input_name}_nn_classification.tsv")
        
        if not os.path.exists(classification_file):
            raise FileNotFoundError(f"Classification file not found: {classification_file}")
        
        os.makedirs(output_folder, exist_ok=True)
        
        combined_output_file = os.path.join(output_folder, f"{input_name}_combined_scaffolds.fasta")
        
        with open(combined_output_file, 'w') as combined_out:
            with open(classification_file, 'r') as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter='\t')
                
                for row in reader:
                    plasmid_score = float(row['plasmid_score'])
                    virus_score = float(row['virus_score'])
                    total_score = plasmid_score + virus_score  # Sum of plasmid and virus scores
                    
                    # Extract if plasmid score > plasmid_threshold
                    if plasmid_score > plasmid_threshold:
                        scaffold_data = extract_scaffold(fasta_file, row['seq_name'])
                        if scaffold_data:
                            combined_out.writelines(scaffold_data)
                            if verbose:
                                print(f"Scaffold {row['seq_name']} added based on plasmid threshold.")
                    
                    # Extract if the sum of plasmid and virus score > combined_threshold
                    if total_score > combined_threshold:
                        scaffold_data = extract_scaffold(fasta_file, row['seq_name'])
                        if scaffold_data:
                            combined_out.writelines(scaffold_data)
                            if verbose:
                                print(f"Scaffold {row['seq_name']} added based on combined plasmid + virus score threshold.")
        
        if verbose:
            print(f"Combined scaffolds saved to {combined_output_file}")
        
        return combined_output_file
    
    except FileNotFoundError as e:
        print(e)
    except ValueError as ve:
        print(ve)
    except Exception as e:
        print(f"An error occurred: {e}")

def extract_scaffold(fasta_file, scaffold_name):
    scaffold_data = []
    found_scaffold = False

    with open(fasta_file, 'r') as file:
        write_scaffold = False

        for line in file:
            if line.startswith('>'):
                if line.startswith(f'>{scaffold_name}'):
                    write_scaffold = True
                    scaffold_data.append(line)
                    found_scaffold = True
                else:
                    write_scaffold = False
            elif write_scaffold:
                scaffold_data.append(line)

    if not found_scaffold:
        print(f"Scaffold {scaffold_name} not found in {fasta_file}")
    
    return scaffold_data if scaffold_data else None

def parse_hmmscan_output(hmmscan_output_file, is_phage):
    """
    Parses the HMMscan output to extract scaffold hits.
    """
    scaffold_hits = {}
    with open(hmmscan_output_file) as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split()
                query_name = fields[2]

                if is_phage:
                    parts = query_name.split('|')
                    if len(parts) > 1:
                        subparts = parts[1].split('_')
                        base_name = f"{parts[0]}|{'_'.join(subparts[:3])}" if len(subparts) >= 4 else query_name
                    else:
                        base_name = query_name
                else:
                    base_name = re.split(r'_\d+$', query_name)[0]

                scaffold_hits.setdefault(base_name, set()).add(query_name)

    return {scaffold: len(genes) for scaffold, genes in scaffold_hits.items()}

def count_genes_per_scaffold(faa_file):
    """
    Counts the number of unique genes (proteins) per scaffold from the Prodigal .faa output.
    """
    gene_counts = {}
    try:
        with open(faa_file, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    scaffold_name = line.split()[0].split('_')[0][1:]
                    protein_name = line.split()[0][1:]  # Get the full protein name

                    # Initialize the set for this scaffold if it doesn't exist
                    if scaffold_name not in gene_counts:
                        gene_counts[scaffold_name] = set()

                    # Add the protein name to the set for this scaffold
                    gene_counts[scaffold_name].add(protein_name)

        # Convert sets to counts
        for scaffold in gene_counts:
            gene_counts[scaffold] = len(gene_counts[scaffold])

    except FileNotFoundError:
        print(f"Error: The file {faa_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return gene_counts

def extract_fasta_description(fasta_file, scaffold_name):
    """
    Extracts the FASTA description for a given scaffold name.
    """
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">") and scaffold_name in line:
                return line.strip()
    return ""

def write_gene_comparison_to_csv(scaffold_hits, gene_counts, fasta_file, output_csv_file, gene_min, percent_min):
    """
    Writes the comparison between total genes and HMMscan hits for each scaffold to a CSV file.
    Also adds the hits/total genes ratio and scaffold description.
    """
    with open(output_csv_file, 'w', newline='') as csvfile:
        fieldnames = ['scaffold', 'total_genes', 'hmmscan_hits', 'hits/total_genes', 'description']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for scaffold in gene_counts:
            total_genes = gene_counts.get(scaffold, 0)
            hmmscan_hits = scaffold_hits.get(scaffold, 0)
            hits_per_gene = hmmscan_hits / total_genes if total_genes > 0 else 0
            if total_genes >= gene_min and hits_per_gene >= percent_min:
                description = extract_fasta_description(fasta_file, scaffold)

                writer.writerow({
                    'scaffold': scaffold,
                    'total_genes': total_genes,
                    'hmmscan_hits': hmmscan_hits,
                    'hits/total_genes': hits_per_gene,
                    'description': description
            })

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Prodigal, HMMscan, and compare gene counts.")
    parser.add_argument('fasta_file', type=str, help='Path to the input FASTA file')
    parser.add_argument('genomad_output', type=str, help='Path to the genomad_output directory')
    parser.add_argument('output_dir', type=str, help='Path to the desired output folder')
    parser.add_argument('hmm_db', type=str, help='Path to the HMM database for HMMscan')
    parser.add_argument('--plasmid_threshold', type=float, default=0.80, help='Plasmid score threshold (default: 0.80)')
    parser.add_argument('--combined_threshold', type=float, default=0.85, help='Plasmid and virus combined threshold (default: 0.85)')
    parser.add_argument('--e_value_threshold', type=float, default=1e-5, help='E-value threshold for HMMscan filtering (default: 1e-5)')
    parser.add_argument('--is_phage', action='store_true', help='Specify if the query is phage-related')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--gene_min', type=int, default=15, help='Minimum gene count per scaffold (default: 15)')
    parser.add_argument('--percent_min', type=float, default=0.15, help='Minimum hits per gene ratio (default: 0.15)')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Extract scaffolds
    combined_fasta = extract_scaffolds(args.genomad_output, args.output_dir, args.fasta_file, args.plasmid_threshold, args.combined_threshold, args.verbose)

    if combined_fasta:
        # Step 2: Run Prodigal
        run_prodigal(combined_fasta, args.output_dir)

        faa_file = os.path.join(args.output_dir, "PlasmidProdigal.faa")
        tblout_file = os.path.join(args.output_dir, "Phage_Plasmid_hmmscan.tblout")

        # Step 3: Run HMMscan with E-value threshold
        subprocess.run(["hmmscan", "--tblout", tblout_file, "-E", str(args.e_value_threshold), args.hmm_db, faa_file], check=True)

        # Step 4: Parse HMMscan results and compare with gene counts
        scaffold_hits = parse_hmmscan_output(tblout_file, args.is_phage)
        gene_counts = count_genes_per_scaffold(faa_file)

        output_csv = os.path.join(args.output_dir, "PhagePlasmids.csv")
        write_gene_comparison_to_csv(scaffold_hits, gene_counts, args.fasta_file, output_csv, args.gene_min, args.percent_min)

        print(f"Gene comparison results saved to {output_csv}")
    else:
        print("No scaffolds were extracted. Prodigal and HMMscan will not run.")
