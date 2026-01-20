def update_fasta_headers(fasta_file, header_id_file, output_file):
    """
    Updates the headers of a FASTA file using a list of new header IDs.

    Args:
        fasta_file (str): Path to the input FASTA file.
        header_id_file (str): Path to the file containing the new header IDs (one per line).
        output_file (str): Path to the output FASTA file with updated headers.

    """
    try:
        # Read new header IDs
        with open(header_id_file, 'r') as id_file:
            header_ids = [line.strip() for line in id_file if line.strip()]
        
        # Read the FASTA file and process it
        with open(fasta_file, 'r') as fasta:
            fasta_lines = fasta.readlines()
        
        updated_lines = []
        header_count = 0

        for line in fasta_lines:
            if line.startswith(">"):  # Header line
                if header_count < len(header_ids):
                    updated_lines.append(f">" + header_ids[header_count] + "\n")
                    header_count += 1
                else:
                    raise ValueError("Number of headers in the FASTA file exceeds the number of IDs in the header ID file.")
            else:
                updated_lines.append(line)  # Sequence line
        
        if header_count != len(header_ids):
            raise ValueError("Number of headers in the FASTA file does not match the number of IDs in the header ID file.")
        
        # Write updated lines to the output file
        with open(output_file, 'w') as out_fasta:
            out_fasta.writelines(updated_lines)
        
        print(f"Headers updated successfully! Output written to {output_file}")
    except Exception as e:
        print(f"Error: {e}")

# Example usage
fasta_file = "input.fasta"        # Path to your FASTA file
header_id_file = "header_ids.txt" # Path to the file with new header IDs
output_file = "output.fasta"      # Path to the output file

update_fasta_headers(fasta_file, header_id_file, output_file)
