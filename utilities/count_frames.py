from ase.io import read

def count_extxyz_frames(file_path=None):
    """
    Count the number of frames in an extxyz file.
    
    Parameters:
    file_path (str, optional): Path to the extxyz file. If None, prompts user.
    
    Returns:
    int: Number of frames in the file
    """
    if file_path is None:
        file_path = input("Enter the path to the extxyz file: ")
    
    try:
        # Read all structures from the file
        structures = read(file_path, index=':')
        
        # If a single structure is returned, convert to list
        if not isinstance(structures, list):
            structures = [structures]
        
        num_frames = len(structures)
        print(f"Number of frames in {file_path}: {num_frames}")
        return num_frames
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        return -1
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return -1

if __name__ == "__main__":
    count_extxyz_frames()
