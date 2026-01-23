import json
import os

def combine_json_files(input_directory, output_file):
    combined_data = {}

    # Loop through all files in the given directory
    for filename in os.listdir(input_directory):
        if filename.endswith('.json'):
            file_path = os.path.join(input_directory, filename)
            with open(file_path, 'r') as file:
                # Load JSON content
                content = json.load(file)
                # Use the filename without the suffix as the key
                key = os.path.splitext(filename)[0]
                combined_data[key] = content

    # Write the combined data to the output file
    with open(output_file, 'w') as outfile:
        json.dump(combined_data, outfile, indent=4)

# Example usage
input_dir = '../../results/boltz/confidence_scores'
output_json_file = '../../results/boltz/combined_plddt_scores.json'
combine_json_files(input_dir, output_json_file)
