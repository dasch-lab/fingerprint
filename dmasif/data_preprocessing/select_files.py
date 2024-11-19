import random

# Define the file path
file_path = '/disk1/fingerprint/SAbDab_preparation/list_pp/testing_ppi.txt'

# Read the file content
with open(file_path, 'r') as file:
    lines = file.readlines()

# Randomly sample 20 lines
sampled_lines = random.sample(lines, 500)

# Process each line to reformat it
processed_lines = []
for line in sampled_lines:
    parts = line.strip().split('_')
    if len(parts) == 3:
        new_name = f"{parts[0]}_{parts[2]}"  # Keep the first and the last part
        processed_lines.append(new_name + '\n')

# Save the processed lines to a new file or print them
with open('sampled_output.txt', 'w') as output_file:
    output_file.writelines(processed_lines)

# Optionally, print the processed lines
for line in processed_lines:
    print(line.strip())

