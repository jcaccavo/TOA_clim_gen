import sys

def subset_fasta(fasta_file, id_file, output_file):
    # Read IDs from the ID file
    with open(id_file, 'r') as id_f:
        id_set = set(line.strip() for line in id_f)

    # Open the FASTA file and the output file
    with open(fasta_file, 'r') as fasta_f, open(output_file, 'w') as out_f:
        # Initialize variables
        write_sequence = False
        sequence_header = ""
        
        # Iterate through each line of the FASTA file
        for line in fasta_f:
            line = line.strip()
            if line.startswith('>'):  # It's a header line
                sequence_header = line[1:]  # Remove the '>' symbol
                if any(id_part in sequence_header for id_part in id_set):
                    write_sequence = True
                    out_f.write(line + "\n")  # Write the header to the output
                else:
                    write_sequence = False
            elif write_sequence:  # It's part of the sequence we want
                out_f.write(line + "\n")

def main():
    # Check command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python subset_fasta.py <fasta_file> <id_file> <output_file>")
        sys.exit(1)

    # Parse command-line arguments
    fasta_file = sys.argv[1]
    id_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Run the subsetting function
    subset_fasta(fasta_file, id_file, output_file)
    print(f"Subset sequences written to {output_file}")

if __name__ == "__main__":
    main()
