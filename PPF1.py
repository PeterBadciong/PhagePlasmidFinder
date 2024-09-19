import subprocess
import csv
import argparse
import os
import re

def run_genomad(fasta_file, genomad_db, splits, threads, output_dir):
    os.environ["LANG"] = "en_US.utf-8"
    command = f"genomad end-to-end --splits {splits} --threads {threads} {fasta_file} {output_dir} {genomad_db}"
    subprocess.run(command, shell=True, check=True)

def run_hmmscan(hmm_file, fasta_file, output_file, evalue_cutoff="1e-5"):
    command = f"hmmscan --tblout {output_file} -E {evalue_cutoff} {hmm_file} {fasta_file}"
    subprocess.run(command, shell=True, check=True)

def parse_hmmscan_output(hmmscan_output_file, is_phage):
    scaffold_hits = {}
    with open(hmmscan_output_file) as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split()
                query_name = fields[2]

                # Determine base name for phage or plasmid
                if is_phage:
                    parts = query_name.split('|')
                    if len(parts) > 1:
                        subparts = parts[1].split('_')
                        base_name = f"{parts[0]}|{'_'.join(subparts[:3])}" if len(subparts) >= 4 else query_name
                    else:
                        base_name = query_name
                else:
                    base_name = re.split(r'_\d+$', query_name)[0]

                # Accumulate hits
                scaffold_hits.setdefault(base_name, set()).add(query_name)

    return {k: len(v) for k, v in scaffold_hits.items()}

def read_gene_count_file(gene_count_file):
    scaffold_genes = {}
    with open(gene_count_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            scaffold_name = row['seq_name']
            gene_count = int(row['n_genes'])
            scaffold_genes[scaffold_name] = gene_count
    return scaffold_genes

def extract_fasta_description(fasta_file, scaffold_name):
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">") and scaffold_name in line:
                return line.strip()
    return ""

def combine_results(scaffold_hits, scaffold_genes, fasta_file, gene_min, percent_min):
    combined_results = []
    for base_name, gene_count in scaffold_genes.items():
        hits = scaffold_hits.get(base_name, 0)
        ratio = hits / gene_count if gene_count > 0 else 0
        if gene_count >= 5 and ratio >= .10:
            description = extract_fasta_description(fasta_file, base_name)
            combined_results.append({
                'Scaffold': base_name,
                'TotalGenes': gene_count,
                'TotalHits': hits,
                'HitsPerGene': ratio,
                'Description': description
            })
    return combined_results

def write_output_file(combined_results, output_file):
    fieldnames = ['Scaffold', 'TotalGenes', 'TotalHits', 'HitsPerGene', 'Description']
    with open(output_file, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(combined_results)

def save_hmmscan_output(hmmscan_output_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, os.path.basename(hmmscan_output_file))
    os.rename(hmmscan_output_file, output_path)

def extract_scaffold(fasta_file, scaffold_name, output_file_handle):
    found_scaffold = False
    with open(fasta_file, 'r') as file:
        write_scaffold = False
        scaffold_data = []
        for line in file:
            if line.startswith('>'):
                if line.startswith(f'>{scaffold_name}'):
                    write_scaffold = True
                    scaffold_data.append(line)
                    found_scaffold = True
                elif write_scaffold:
                    output_file_handle.writelines(scaffold_data)
                    scaffold_data = []
                    write_scaffold = False
            elif write_scaffold:
                scaffold_data.append(line)

        if write_scaffold and scaffold_data:
            output_file_handle.writelines(scaffold_data)

    if not found_scaffold:
        print(f"Scaffold {scaffold_name} not found in {fasta_file}")

def run_full_analysis(fasta_file, genomad_db, splits, threads, output_dir, phage_hmm_file, plasmid_hmm_file, evalue_cutoff="1e-5", gene_min=0, percent_min=0, extract_toggle=False):
    os.makedirs(output_dir, exist_ok=True)
    genomad_output_dir = os.path.join(output_dir, 'genomad_output')
    os.makedirs(genomad_output_dir, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    run_genomad(fasta_file, genomad_db, splits, threads, genomad_output_dir)

    phage_fasta_file = os.path.join(genomad_output_dir, f"{base_name}_summary", f"{base_name}_virus.fna")
    plasmid_fasta_file = os.path.join(genomad_output_dir, f"{base_name}_summary", f"{base_name}_plasmid.fna")

    phage_proteins_file = os.path.join(genomad_output_dir, f"{base_name}_summary", f"{base_name}_virus_proteins.faa")
    plasmid_proteins_file = os.path.join(genomad_output_dir, f"{base_name}_summary", f"{base_name}_plasmid_proteins.faa")

    phage_gene_count_file = os.path.join(genomad_output_dir, f"{base_name}_summary", f"{base_name}_virus_summary.tsv")
    plasmid_gene_count_file = os.path.join(genomad_output_dir, f"{base_name}_summary", f"{base_name}_plasmid_summary.tsv")

    phage_hmmscan_output_file = os.path.join(output_dir, "phage_hmmscan_output.tbl")
    plasmid_hmmscan_output_file = os.path.join(output_dir, "plasmid_hmmscan_output.tbl")

    if os.path.exists(phage_proteins_file) and os.path.getsize(phage_proteins_file) > 0:
        run_hmmscan(plasmid_hmm_file, phage_proteins_file, phage_hmmscan_output_file, evalue_cutoff)
        phage_scaffold_hits = parse_hmmscan_output(phage_hmmscan_output_file, is_phage=True)
        phage_scaffold_genes = read_gene_count_file(phage_gene_count_file)
        phage_combined_results = combine_results(phage_scaffold_hits, phage_scaffold_genes, fasta_file, gene_min, percent_min)
        write_output_file(phage_combined_results, os.path.join(output_dir, f"{base_name}_phage.csv"))
    else:
        print(f"Phage proteins file '{phage_proteins_file}' is empty or misformatted. Skipping phage HMMScan analysis.")

    if os.path.exists(plasmid_proteins_file) and os.path.getsize(plasmid_proteins_file) > 0:
        run_hmmscan(phage_hmm_file, plasmid_proteins_file, plasmid_hmmscan_output_file, evalue_cutoff)
        plasmid_scaffold_hits = parse_hmmscan_output(plasmid_hmmscan_output_file, is_phage=False)
        plasmid_scaffold_genes = read_gene_count_file(plasmid_gene_count_file)
        plasmid_combined_results = combine_results(plasmid_scaffold_hits, plasmid_scaffold_genes, fasta_file, gene_min, percent_min)
        write_output_file(plasmid_combined_results, os.path.join(output_dir, f"{base_name}_plasmid.csv"))
    else:
        print(f"Plasmid proteins file '{plasmid_proteins_file}' is empty or misformatted. Skipping plasmid HMMScan analysis.")

    if os.path.exists(phage_hmmscan_output_file):
        save_hmmscan_output(phage_hmmscan_output_file, output_dir)
    if os.path.exists(plasmid_hmmscan_output_file):
        save_hmmscan_output(plasmid_hmmscan_output_file, output_dir)

    if extract_toggle:
        # Collect scaffolds for virus and plasmid
        scaffolds_to_extract = [result['Scaffold'] for result in phage_combined_results + plasmid_combined_results]

        # Extract scaffolds and write to a single FASTA file
        output_fasta_file = os.path.join(output_dir, 'PhagePlasmidScaffolds.fna')
        with open(output_fasta_file, 'w') as output_file:
            for scaffold in set(scaffolds_to_extract):  # Use set to avoid duplicates
                extract_scaffold(phage_fasta_file, scaffold, output_file)  # If phage_proteins_file was processed
                extract_scaffold(plasmid_fasta_file, scaffold, output_file)  # If plasmid_proteins_file was processed

def main():
    parser = argparse.ArgumentParser(description='Run Genomad and HMMScan analysis on FASTA files.')
    parser.add_argument('fasta_file', help='Input FASTA file.')
    parser.add_argument('genomad_db', help='Genomad database path.')
    parser.add_argument('-s', '--splits', type=int, default=1, help='Number of splits for Genomad.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for Genomad.')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for output files.')
    parser.add_argument('-j', '--phage_hmm_file', required=True, help='HMM file for phage analysis.')
    parser.add_argument('-l', '--plasmid_hmm_file', required=True, help='HMM file for plasmid analysis.')
    parser.add_argument('-e', '--evalue_cutoff', default='1e-5', help='E-value cutoff for HMMScan.')
    parser.add_argument('-g', '--gene_min', type=int, default=0, help='Minimum number of genes.')
    parser.add_argument('-p', '--percent_min', type=float, default=0, help='Minimum hits per gene percentage.')
    parser.add_argument('-x', '--extract_toggle', action='store_true', help='Toggle scaffold extraction.')

    args = parser.parse_args()

    run_full_analysis(
        fasta_file=args.fasta_file,
        genomad_db=args.genomad_db,
        splits=args.splits,
        threads=args.threads,
        output_dir=args.output_dir,
        phage_hmm_file=args.phage_hmm_file,
        plasmid_hmm_file=args.plasmid_hmm_file,
        evalue_cutoff=args.evalue_cutoff,
        gene_min=args.gene_min,
        percent_min=args.percent_min,
        extract_toggle=args.extract_toggle
    )

if __name__ == "__main__":
    main()
