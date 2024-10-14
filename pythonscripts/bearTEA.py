# Title: Getting genomics co-ordinates for significantly differentially expressed TE's from Telescope outputs
# Author: Dr. Alice M. Godden

import csv

def match_TEs(csv_file, gtf_file, output_file):
    # Open the output file to write results
    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')  # Use tab delimiter for output

        # Read the list of TEs from the CSV file
        with open(csv_file, 'r') as te_file:
            te_reader = csv.reader(te_file)
            for row in te_reader:
                te_name = row[0]  # Assuming TEs are in the first column
                print(f"Searching for TE: {te_name}")  # Optional: Print the TE being searched
                
                # Search for the TE in the GTF file
                with open(gtf_file, 'r') as gtf:
                    for gtf_line in gtf:
                        if te_name in gtf_line:
                            parts = gtf_line.strip().split('\t')
                            # Write the 1st, 4th, and 5th columns to the output
                            writer.writerow([parts[0], parts[3], parts[4]])
    
    print(f"Matching completed. Results saved to {output_file}.")

# Example usage
if __name__ == "__main__":
    csv_file = "sigDETEs.csv"  # Path to the input CSV file
    gtf_file = "ASM1731132v1_rmsk_TE.gtf"  # Path to the input GTF file
    output_file = "matched_TEs_output.txt"  # Path to save the output
    match_TEs(csv_file, gtf_file, output_file)
