from ase.io import read, write

# Read all frames from data.extxyz
atoms_list = read('data.extxyz', index=':')  # ':' reads all frames

# Select every alternate frame (50% of frames)
selected_frames = atoms_list[::2]

# Write the selected frames to half.extxyz
write('half.extxyz', selected_frames)

print(f"Original number of frames: {len(atoms_list)}")
print(f"Frames written to new file: {len(selected_frames)}")
