import os
import csv

def extract_coordinates(pdb_path, output_file, selected_chain='A'):
    """
    Extract alpha carbon coordinates for a specific chain from a PDB file
    """
    try:
        coordinates = []
        
        with open(pdb_path, 'r') as file:
            for line in file:
                if (line.startswith("ATOM") and 
                    line[12:16].strip() == "CA" and 
                    line[21].strip() == selected_chain):
                    
                    try:
                        residue_number = int(line[22:26].strip())
                        
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        coordinates.append([residue_number, x, y, z])
                    
                    except (ValueError, IndexError) as e:
                        print(f"Skipping line: {line.strip()}")
        
        coordinates.sort(key=lambda x: x[0])
        
        # Write to CSV
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['residue_number', 'x', 'y', 'z'])
            csv_writer.writerows(coordinates)
        
        print(f"Extracted {len(coordinates)} coordinates for chain {selected_chain}")
        return True
    
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return False

def process_pdb_files(pdb_directory, output_directory, selected_chain='A'):
    """
    Process all PDB files in a directory, extracting coordinates for a specific chain
    """

    os.makedirs(output_directory, exist_ok=True)
    
    total_files = 0
    processed_files = 0
    
    for pdb_file in os.listdir(pdb_directory):
        if pdb_file.endswith(".pdb"):
            total_files += 1
            
            pdb_path = os.path.join(pdb_directory, pdb_file)
            base_name = os.path.splitext(pdb_file)[0]
            output_file = os.path.join(output_directory, f"{base_name}_chain_{selected_chain}.csv")
            
            if extract_coordinates(pdb_path, output_file, selected_chain):
                processed_files += 1
    
    # Summary
    print("\nProcessing Summary:")
    print(f"Total PDB files found: {total_files}")
    print(f"Successfully processed files: {processed_files}")

pdb_directory = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/2_pdb_files_full"
output_directory = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/3_coordinate_files"

if __name__ == "__main__":
    process_pdb_files(pdb_directory, output_directory, selected_chain='A')