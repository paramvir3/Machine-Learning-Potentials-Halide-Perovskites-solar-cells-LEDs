from ase.io import read, write
import numpy as np

def split_extxyz(input_file, test_ratio=0.001, random_seed=None):
    """
    Split an extxyz file into train and test sets.

    Parameters:
        input_file (str): Path to the input extxyz file.
        test_ratio (float): Fraction of frames to use for testing (default: 0.1).
        random_seed (int): Seed for random number generation (optional).

    Returns:
        None
    """
    # Set random seed for reproducibility
    if random_seed is not None:
        np.random.seed(random_seed)

    # Read all frames from the extxyz file
    frames = read(input_file, index=':')
    num_frames = len(frames)

    # Shuffle the frames randomly
    indices = np.arange(num_frames)
    np.random.shuffle(indices)

    # Split indices into test and train sets
    num_test = int(num_frames * test_ratio)
    test_indices = indices[:num_test]
    train_indices = indices[num_test:]

    # Extract test and train frames
    test_frames = [frames[i] for i in test_indices]
    train_frames = [frames[i] for i in train_indices]

    # Write test frames to a new extxyz file
    test_file = input_file.replace('.extxyz', '_test.extxyz')
    write(test_file, test_frames, format='extxyz')
    print(f"Test frames saved to {test_file} ({len(test_frames)} frames)")

    # Write train frames to a new extxyz file
    train_file = input_file.replace('.extxyz', '_train.extxyz')
    write(train_file, train_frames, format='extxyz')
    print(f"Train frames saved to {train_file} ({len(train_frames)} frames)")

# Example usage
input_file = 'data.extxyz'  # Input extxyz file
split_extxyz(input_file, test_ratio=0.001, random_seed=42)
