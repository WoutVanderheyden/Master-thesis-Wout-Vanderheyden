import os
import pandas as pd
from Bio import AlignIO
import re

def extract_coordinates_from_alignment(alignment_file, coordinate_dir, output_dir, summary_dir):
    """
    Extract coordinates from CSV files based on alignment data.
    
    Args:
        alignment_file (str): Path to the protein sequence alignment file
        coordinate_dir (str): Directory containing CSV coordinate files
        output_dir (str): Directory to output the extracted coordinates files
        summary_dir (str): Directory to save the summary files
    """
    # Create output and summary directories
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)
        
    # Load alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    
    coordinate_files = {}
    coordinate_data = {}
    seq_id_to_full_id = {}
    
    # Dictionary to track which coordinates are used for each output file
    coordinate_tracking = {}
    
    for record in alignment:
        # Extract first part of sequence ID
        full_id = record.id
        seq_id = full_id.split('|')[0]
 
        normalized_id = seq_id.replace('/', '_').lower()
        seq_id_to_full_id[normalized_id] = full_id
    
    for filename in os.listdir(coordinate_dir):
        if filename.endswith('.csv'):
            ############################# CHANGE THIS BASED ON WHICH STRAIN WE HAVE !!!!!!!! ##################################""""
            match = re.search(r'fold_(a_.*?)_model', filename)
            if match:
                csv_name = match.group(1).lower() 

                for normalized_id, full_id in seq_id_to_full_id.items():
                    if normalized_id in csv_name:
                        print(f"Matched {full_id} to {filename}")
                        #### for now i changed full_id by filename #####
                        coordinate_files[full_id] = os.path.join(coordinate_dir, filename)
                        # Load the coordinate data
                        coordinate_data[full_id] = pd.read_csv(coordinate_files[full_id])
                        # Initialize tracking for this sequence
                        coordinate_tracking[full_id] = {
                            'source_file': filename,
                            'used_coordinates': []
                        }
                        break
    
    # Check if we found all coordinate files
    missing_files = [seq.id for seq in alignment if seq.id not in coordinate_files]
    if missing_files:
        print(f"Warning: Could not find coordinate files for: {', '.join(missing_files)}")

    # Check for CSV files that weren't matched to any sequence
    matched_csv_files = set(info['source_file'] for info in coordinate_tracking.values())
    all_csv_files = [f for f in os.listdir(coordinate_dir) if f.endswith('.csv')]
    unmatched_csv_files = set(all_csv_files) - matched_csv_files
    if unmatched_csv_files:
        print(f"Warning: The following CSV files were not matched to any sequence: {', '.join(unmatched_csv_files)}")
    
    # Initialize dictionaries to store output data and counters for each sequence
    output_data = {seq_id: [] for seq_id in coordinate_files}
    residue_counters = {seq.id: 1 for seq in alignment}
    
    # Get the alignment length
    align_length = alignment.get_alignment_length()
    
    for pos in range(align_length):
        # Check if any sequence has a gap at this position
        has_gap = False
        for record in alignment:
            if record[pos] == '-':
                has_gap = True
                break

        if not has_gap:
            # No gaps at this position, extract coordinates for each sequence
            for record in alignment:
                seq_id = record.id
                if seq_id in coordinate_files:
                    residue_idx = residue_counters[seq_id]

                    # Make sure we're within bounds of the coordinate data
                    if residue_idx <= len(coordinate_data[seq_id]):
                        # Get the coordinate data for this residue
                        coords = coordinate_data[seq_id].iloc[residue_idx - 1]
                        
                        # Add to output data for this sequence
                        output_data[seq_id].append({
                            'residue_number': residue_idx,
                            'x': coords['x'],
                            'y': coords['y'],
                            'z': coords['z']
                        })
                        
                        # Track which coordinates are used
                        coordinate_tracking[seq_id]['used_coordinates'].append(residue_idx)
                    else:
                        print(f"Warning: Residue index {residue_idx} exceeds coordinate data for {seq_id}")
                
                # Increment the residue counter for non-gap positions
                residue_counters[seq_id] += 1
        else:
            # Increment counters only for sequences that don't have a gap at this position
            for record in alignment:
                seq_id = record.id
                if record[pos] != '-':
                    residue_counters[seq_id] += 1

    # Create and save output files for each sequence
    for seq_id, data in output_data.items():
        if data:  # Only create files for sequences with data
            output_df = pd.DataFrame(data)

            # Create output filename using the sequence ID
            # Replace any characters that are not valid in filenames
            safe_seq_id = seq_id.replace('/', '_').replace('|', '_')
            output_file = os.path.join(output_dir, f"{safe_seq_id}.csv")
            
            output_df.to_csv(output_file, index=False)
            print(f"Extracted coordinates for {seq_id} saved to {output_file}")
            
            # Save coordinate tracking info to the summary directory
            tracking_file = os.path.join(summary_dir, f"{safe_seq_id}_coordinate_info.txt")
            with open(tracking_file, 'w') as f:
                f.write(f"Source coordinate file: {coordinate_tracking[seq_id]['source_file']}\n")
                f.write(f"Total coordinates used: {len(coordinate_tracking[seq_id]['used_coordinates'])}\n")
                f.write(f"Residue indices used: {', '.join(map(str, coordinate_tracking[seq_id]['used_coordinates']))}\n")
                f.write(f"First residue: {coordinate_tracking[seq_id]['used_coordinates'][0]}\n")
                f.write(f"Last residue: {coordinate_tracking[seq_id]['used_coordinates'][-1]}\n")
            print(f"Coordinate tracking info saved to {tracking_file}")
    
    # Additionally, save a summary of all coordinate files used to the summary directory
    summary_file = os.path.join(summary_dir, "coordinate_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Coordinate File Usage Summary\n")
        f.write("===========================\n\n")
        for seq_id, tracking in coordinate_tracking.items():
            f.write(f"Sequence ID: {seq_id}\n")
            f.write(f"Source file: {tracking['source_file']}\n")
            f.write(f"Total coordinates used: {len(tracking['used_coordinates'])}\n")
            f.write(f"Residue range: {tracking['used_coordinates'][0]}-{tracking['used_coordinates'][-1]}\n\n")
    print(f"Summary of coordinate usage saved to {summary_file}")
    
    print(f"Processing complete. Files saved to {output_dir}")
    print(f"Summary files saved to {summary_dir}")

if __name__ == "__main__":
    # Example usage
    alignment_file = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/4_align_coordinates/TEST_aligned.fasta"
    coordinate_dir = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/3_coordinate_files"
    output_dir = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/5_aligned_coordinate_files"
    summary_dir = "/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/5_summary"
    
    extract_coordinates_from_alignment(alignment_file, coordinate_dir, output_dir, summary_dir)