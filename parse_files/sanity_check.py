import sys
import glob
import re
import json
import numpy as np


jsonErr=4
valErr=5
matrix_not_exist_error=6
matrix_cond_error=7
atom_position_error=8
duplicate_position_error=9
#This script runs a basic sanity check of data in conf file




# Read JSON data from stdin
try:
    config_json = sys.stdin.read()
    parsed_config = json.loads(config_json)

except json.JSONDecodeError as e:
    print(f"Error parsing JSON input: {e}")
    exit(jsonErr)


#####################################################################################
#check matrices
def check_matrix_condition(matrix, matrix_name="Matrix", det_threshold=1e-12, cond_threshold=1e12):
    """

    :param matrix: 2D list or numpy array representing the matrix
    :param matrix_name: Name of the matrix for error messages
    :param det_threshold: Minimum absolute determinant value (default: 1e-12)
    :param cond_threshold: Maximum condition number (default: 1e12)
    :return: tuple: (is_valid, error_message)
               is_valid: True if matrix passes all checks
               error_message: None if valid, error string if invalid
    """
    try:
        # Convert to numpy array if it's a list
        if isinstance(matrix, list):
            np_matrix = np.array(matrix)
        else:
            np_matrix = matrix
        # Check if it's a square matrix
        if np_matrix.shape[0] != np_matrix.shape[1]:
            return False, f"{matrix_name} is not square: shape {np_matrix.shape}"

        # Check determinant (non-degenerate test)
        det = np.linalg.det(np_matrix)
        if abs(det) < det_threshold:
            return False, f"{matrix_name} is degenerate (determinant ≈ 0): det = {det:.2e}"

        # Check condition number (ill-conditioning test)
        cond_num = np.linalg.cond(np_matrix)
        if cond_num > cond_threshold:
            return False, f"{matrix_name} is ill-conditioned: condition number = {cond_num:.2e}"

        return True, None

    except Exception as e:
        return False, f"Error analyzing {matrix_name}: {str(e)}"



# Check for required matrix fields
if 'lattice_basis' not in parsed_config or not parsed_config['lattice_basis']:
    print("Error: Missing or empty required field 'lattice_basis'")
    exit(matrix_not_exist_error)

if 'space_group_basis' not in parsed_config or not parsed_config['space_group_basis']:
    print("Error: Missing or empty required field 'space_group_basis'")
    exit(matrix_not_exist_error)


# Check matrix conditions
is_valid, error_msg = check_matrix_condition(parsed_config['lattice_basis'], "Lattice basis")
if not is_valid:
    print(f"Error: {error_msg}")
    exit(matrix_cond_error)

is_valid, error_msg = check_matrix_condition(parsed_config['space_group_basis'], "Space group basis")
if not is_valid:
    print(f"Error: {error_msg}")
    exit(matrix_cond_error)

# end checking matrices
#####################################################################################
def check_atom_positions(parsed_config):
    """
    Check that the number of atom positions matches the number of atoms for each type

    Returns:
        tuple: (is_valid, error_message)
    """
    atom_types = parsed_config.get('atom_types', {})
    atom_positions = parsed_config.get('atom_positions', [])

    if not atom_types:
        return False, "No atom types defined"

    if not atom_positions:
        return False, "No atom positions defined"

    # Count positions for each atom type
    position_counts = {}
    for position in atom_positions:
        # Changed from 'type' to 'atom_type'
        atom_type = position.get('atom_type')
        if atom_type:
            position_counts[atom_type] = position_counts.get(atom_type, 0) + 1

    # Check each atom type
    for atom_type, atom_data in atom_types.items():
        # Extract count from the atom data dictionary
        if isinstance(atom_data, dict):
            expected_count = atom_data.get('count', 0)
        else:
            # Fallback for simple integer format
            expected_count = atom_data

        actual_count = position_counts.get(atom_type, 0)

        if actual_count != expected_count:
            return False, f"Atom type '{atom_type}': expected {expected_count} positions, found {actual_count}"

    # Check for positions of undefined atom types
    for atom_type in position_counts:
        if atom_type not in atom_types:
            return False, f"Found positions for undefined atom type '{atom_type}'"

    return True, None
def check_duplicate_positions(parsed_config, tolerance=1e-6):
    """
    Check for duplicate atomic positions after reducing by lattice vectors
    """
    # Get lattice basis and atom positions
    lattice_basis = parsed_config.get('lattice_basis')
    atom_positions = parsed_config.get('atom_positions', [])

    if not lattice_basis:
        return False, "Lattice basis not found"
    if not atom_positions:
        return True, None  # No positions to check

    try:
        # Convert lattice basis to numpy array
        lattice_matrix = np.array(lattice_basis)
        # Store reduced positions for comparison
        reduced_positions = []
        position_info = []  # For error reporting

        for i, position in enumerate(atom_positions):
            # Use position_name if available, otherwise fall back to atom_type, then generic name
            position_name = position.get('position_name')
            atom_type = position.get('atom_type')

            if position_name:
                display_name = position_name
            elif atom_type:
                display_name = atom_type
            else:
                display_name = f'atom_{i}'

            coords = position.get('coordinates')

            if not coords or len(coords) != 3:
                return False, f"Invalid coordinates for {display_name} at index {i}"

            # Convert to numpy array and reduce modulo lattice vectors
            coord_array = np.array(coords, dtype=float)
            # Reduce coordinates to [0, 1) in fractional coordinates
            reduced_coord = coord_array % 1.0
            # Handle numerical precision issues near boundaries
            reduced_coord = np.where(reduced_coord > 1.0 - tolerance, 0.0, reduced_coord)

            reduced_positions.append(reduced_coord)
            position_info.append((display_name, i, coords))

        # Check for duplicates
        for i in range(len(reduced_positions)):
            for j in range(i + 1, len(reduced_positions)):
                pos1 = reduced_positions[i]
                pos2 = reduced_positions[j]
                # Calculate distance considering periodic boundary conditions
                diff = pos1 - pos2
                # Handle wraparound (e.g., 0.999 and 0.001 should be close)
                diff = np.where(diff > 0.5, diff - 1.0, diff)
                diff = np.where(diff < -0.5, diff + 1.0, diff)

                distance = np.linalg.norm(diff)
                if distance < tolerance:
                    name1, idx1, coords1 = position_info[i]
                    name2, idx2, coords2 = position_info[j]
                    return False, (f"Duplicate positions found: {name1} at {coords1} "
                                   f"and {name2} at {coords2} are equivalent after "
                                   f"lattice reduction (distance: {distance:.2e})")

        return True, None
    except Exception as e:
        return False, f"Error checking duplicate positions: {str(e)}"

# Check atom positions match atom counts
is_valid, error_msg = check_atom_positions(parsed_config)
if not is_valid:
    print(f"Error: {error_msg}")
    exit(atom_position_error)

# Debug output before atom position check
# print("DEBUG: atom_types:", parsed_config.get('atom_types'))
# print("DEBUG: atom_positions:", parsed_config.get('atom_positions'))
# if 'atom_positions' in parsed_config:
#     print("DEBUG: Position details:")
#     for i, pos in enumerate(parsed_config['atom_positions']):
#         print(f"  Position {i}: type='{pos.get('type')}', coords={pos.get('coordinates')}")
# Check for duplicate positions
is_valid, error_msg = check_duplicate_positions(parsed_config)
if not is_valid:
    print(f"Error: {error_msg}")
    exit(duplicate_position_error)
#check space group

print("SUCCESS: All sanity checks passed!")