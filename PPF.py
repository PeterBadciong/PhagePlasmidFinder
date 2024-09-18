import subprocess
import os
import argparse
import shutil  # Import shutil for moving files

# Setup argparse to take inputs from the command line
parser = argparse.ArgumentParser(description="Run script1 and script2 sequentially with shared inputs.")
parser.add_argument('fasta_file', help='Input FASTA file.')
parser.add_argument('genomad_db', help='Genomad database path.')  # Include genomad_db here
parser.add_argument('-s', '--splits', type=int, default=1, help='Number of splits for Genomad.')
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for Genomad.')
parser.add_argument('-o', '--output_dir', required=True, help='Directory for output files.')
parser.add_argument('-j', '--phage_hmm_file', required=True, help='HMM file for phage analysis.')
parser.add_argument('-l', '--plasmid_hmm_file', required=True, help='HMM file for plasmid analysis.')
parser.add_argument('-e', '--evalue_cutoff', default='1e-5', help='E-value cutoff for HMMScan.')
parser.add_argument('-g', '--gene_min', type=int, default=0, help='Minimum number of genes (default: 15)')
parser.add_argument('-p', '--percent_min', type=float, default=0, help='Minimum hits per gene percentage (default .15)')
parser.add_argument('-x', '--extract_toggle', action='store_true', help='Toggle scaffold extraction.')
parser.add_argument('-m', '--plasmid_threshold', type=float, default=0.75, help='Plasmid score threshold (default: 0.75)')
parser.add_argument('-c', '--combined_threshold', type=float, default=0.85, help='Plasmid and virus combined threshold (default: 0.85)')

args = parser.parse_args()

# Ensure the output directory exists
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Create the Discovery directory inside the output directory
discovery_dir = os.path.join(args.output_dir, 'Discovery')
if not os.path.exists(discovery_dir):
    os.makedirs(discovery_dir)

# Set the path for the error log
error_log_file = os.path.join(args.output_dir, 'error_log.txt')

# Function to log errors and output
def run_with_logging(command):
    try:
        with open(error_log_file, 'a') as log_file:
            log_file.write(f"Running command: {' '.join(command)}\n")
            process = subprocess.run(command, capture_output=True, text=True)
            log_file.write("Standard Output:\n")
            log_file.write(process.stdout + "\n")
            log_file.write("Error Output:\n")
            log_file.write(process.stderr + "\n")
            if process.returncode != 0:
                log_file.write(f"Command failed with return code: {process.returncode}\n")
            log_file.write("\n")
    except Exception as e:
        print(f"Error running command {' '.join(command)}: {e}")
        with open(error_log_file, 'a') as log_file:
            log_file.write(f"Exception occurred: {e}\n")

# Log info to console as well
def log_to_console_and_file(message):
    print(message)
    with open(error_log_file, 'a') as log_file:
        log_file.write(message + "\n")

log_to_console_and_file(f"Writing logs to {error_log_file}")

# Run the first script (PhagePlasmidFinder.py) with arguments from the wrapper script
script1_args = [
    'python3', 'PhagePlasmidFinder.py',
    args.fasta_file,
    args.genomad_db,
    '--plasmid_hmm_file', args.plasmid_hmm_file,
    '--phage_hmm_file', args.phage_hmm_file,
    '-o', args.output_dir,
    '-e', args.evalue_cutoff,
    '-g', str(args.gene_min),
    '-p', str(args.percent_min),
    '-t', str(args.threads)
] + (['-x'] if args.extract_toggle else [])

log_to_console_and_file("Running script1 (PhagePlasmidFinder.py)...")
run_with_logging(script1_args)

# Check if the first script was successful before running the second script
genomad_output = os.path.join(args.output_dir, 'genomad_output')

if os.path.exists(genomad_output):
    # Run the second script (PlasmidBypass.py) with only the required arguments
    script2_args = [
        'python3', 'PlasmidBypass.py',
        args.fasta_file,
        genomad_output,
        args.output_dir,
        args.phage_hmm_file,
        '--plasmid_threshold', str(args.plasmid_threshold),
	'--combined_threshold', str(args.combined_threshold),
	'--gene_min', str(args.gene_min),
	'--percent_min', str(args.percent_min)
    ]

    log_to_console_and_file("Running script2 (PlasmidBypass.py)...")
    run_with_logging(script2_args)

# Move specified output files to the Discovery directory
output_extensions = ['.tbl', '.fasta', 'phage.csv', 'plasmid.csv']
for filename in os.listdir(args.output_dir):
    if any(filename.endswith(ext) for ext in output_extensions):
        source_path = os.path.join(args.output_dir, filename)
        destination_path = os.path.join(discovery_dir, filename)
        shutil.move(source_path, destination_path)
        log_to_console_and_file(f"Moved {filename} to {discovery_dir}")

# Create the prodigal directory inside the output directory
prodigal_dir = os.path.join(args.output_dir, 'prodigal')
if not os.path.exists(prodigal_dir):
    os.makedirs(prodigal_dir)

# Move specified output files to the prodigal directory
output_extensions_prodigal = ['.tblout', '.faa', '.gff']
for filename in os.listdir(args.output_dir):
    if any(filename.endswith(ext) for ext in output_extensions_prodigal):
        source_path = os.path.join(args.output_dir, filename)
        destination_path = os.path.join(prodigal_dir, filename)
        shutil.move(source_path, destination_path)
        log_to_console_and_file(f"Moved {filename} to {prodigal_dir}")
