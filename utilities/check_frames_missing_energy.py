def find_frames_without_energy(file_path):
    """
    Analyzes an extxyz file to find frames and line numbers without energy information.
    
    Args:
        file_path (str): Path to the extxyz file
        
    Returns:
        list: List of tuples containing (frame_number, line_number) for frames without energy
    """
    frames_without_energy = []
    current_line = 0
    frame_number = 0
    
    try:
        with open(file_path, 'r') as file:
            while True:
                # Read the number of atoms line (first line of frame)
                num_atoms_line = file.readline()
                current_line += 1
                
                # Check if we've reached end of file
                if not num_atoms_line:
                    break
                    
                # Skip if line is empty
                if not num_atoms_line.strip():
                    continue
                    
                try:
                    num_atoms = int(num_atoms_line.strip())
                except ValueError:
                    print(f"Error: Invalid number of atoms at line {current_line}")
                    break
                    
                # Read the comment/properties line (second line of frame)
                properties_line = file.readline()
                current_line += 1
                
                # Check if energy is missing in properties line
                # Typically energy is indicated by "energy=" or similar keyword
                if 'energy=' not in properties_line.lower():
                    frames_without_energy.append((frame_number, current_line - 1))
                
                # Skip the atom coordinate lines
                for _ in range(num_atoms):
                    file.readline()
                    current_line += 1
                    
                frame_number += 1
                
        return frames_without_energy
        
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        return []
    except Exception as e:
        print(f"Error: An unexpected error occurred: {str(e)}")
        return []

def main():
    # Example usage
    file_path = "./unit_cells.extxyz"  # Replace with your file path
    results = find_frames_without_energy(file_path)
    
    if results:
        print("Frames without energy found at:")
        for frame, line in results:
            print(f"Frame {frame} (Line {line})")
    else:
        print("No frames without energy found or file was empty")

if __name__ == "__main__":
    main()
