# Title: Getting genomics co-ordinates for significantly differentially expressed TE's from Telescope outputs
# Author: Dr. Alice M. Godden

import csv

# Function to read TE_rmsk.gtf file and store start and end positions in a dictionary
def read_gtf_file(gtf_file):
    te_positions = {}
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chromosome = parts[0]
            attributes = parts[8].split(';')
            te_id_parts = {}
            for attr in attributes:
                try:
                    key, value = attr.strip().split(' ')
                    te_id_parts[key] = value.strip('"')
                except ValueError:
                    continue
            te_id = te_id_parts.get('gene_id')
            family_id = te_id_parts.get('family_id')
            class_id = te_id_parts.get('class_id')
            if te_id and family_id and class_id:
                start = int(parts[3])
                end = int(parts[4])
                te_key = (te_id, family_id, class_id)
                if te_key in te_positions:
                    te_positions[te_key].append((chromosome, start, end))
                else:
                    te_positions[te_key] = [(chromosome, start, end)]
    return te_positions

# Function to match TE IDs from csv file with TE_rmsk.gtf file and create output
def match_and_write(csv_file, gtf_file, output_file):
    te_positions = read_gtf_file(gtf_file)
    with open(csv_file, 'r') as csvfile, open(output_file, 'w', newline='') as output:
        csv_reader = csv.reader(csvfile)
        csv_writer = csv.writer(output)
        next(csv_reader)  # Skip the header line
        for row in csv_reader:
            te_name = row[0]
            if "#" in te_name:
                gene_id, rest = te_name.split("#", 1)
                if "/" in rest:
                    family_id, class_id = rest.split("/", 1)
                else:
                    family_id, class_id = rest, ""
            else:
                gene_id, family_id, class_id = te_name, "", ""

            te_key = (gene_id, family_id, class_id)
            if te_key in te_positions:
                for chromosome, start, end in te_positions[te_key]:
                    csv_writer.writerow([chromosome, gene_id, family_id, class_id, start, end])
            else:
                print(f"No matching TE found in GTF file for: {te_name}")

# Example usage
if __name__ == "__main__":
    csv_file = "sigDETEs.csv"
    gtf_file = "/gpfs/home/nfv16zpu/bear/ASM1731132v1_rmsk_TE.gtf"
    output_file = "sigDETEs_output.csv"
    match_and_write(csv_file, gtf_file, output_file)
