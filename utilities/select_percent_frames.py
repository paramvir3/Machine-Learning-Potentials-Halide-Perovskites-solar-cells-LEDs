from ase.io import read, write
import random
from collections import Counter

# Read all frames from data.extxyz
atoms_list = read('data.extxyz', index=':')

# Get all unique chemical elements in the dataset
all_elements = set()
for atoms in atoms_list:
    all_elements.update(atoms.get_chemical_symbols())
all_elements = sorted(list(all_elements))  # Sort for consistency
print(f"All chemical elements found: {all_elements}")

# Calculate target number of frames (2% of total)
total_frames = len(atoms_list)
target_frames = max(1, int(total_frames * 0.02))  # At least 1 frame

# Function to get elements in a subset
def get_elements_in_subset(subset):
    elements = set()
    for atoms in subset:
        elements.update(atoms.get_chemical_symbols())
    return elements

# Initial random selection
random.seed(42)  # For reproducibility
selected_frames = random.sample(atoms_list, target_frames)
selected_elements = get_elements_in_subset(selected_frames)

# Keep adding frames until all elements are represented
remaining_frames = [f for f in atoms_list if f not in selected_frames]
while selected_elements != set(all_elements) and remaining_frames:
    # Find missing elements
    missing_elements = set(all_elements) - selected_elements
    
    # Find a frame with at least one missing element
    for i, frame in enumerate(remaining_frames):
        frame_elements = set(frame.get_chemical_symbols())
        if frame_elements & missing_elements:  # If frame has any missing element
            selected_frames.append(frame)
            remaining_frames.pop(i)
            selected_elements = get_elements_in_subset(selected_frames)
            break

# Write the selected frames to output file
write('test.extxyz', selected_frames)

# Print statistics
final_count = len(selected_frames)
print(f"Original number of frames: {total_frames}")
print(f"Target frames (2%): {target_frames}")
print(f"Actual frames selected: {final_count}")
print(f"Percentage selected: {(final_count/total_frames)*100:.1f}%")
print(f"Elements in selection: {sorted(list(get_elements_in_subset(selected_frames)))}")
