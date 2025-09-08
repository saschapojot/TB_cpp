import numpy as np
import sys
import json
import re
#original file: /home/adada/Documents/pyCode/TB/cd/SymGroup.py
# this script computes space group representations
json_err_code=4
key_err_code=5
val_err_code=6

# Read JSON data from stdin

try:
    config_json = sys.stdin.read()
    parsed_config = json.loads(config_json)

except json.JSONDecodeError as e:
    print(f"Error parsing JSON input: {e}")
    exit(json_err_code)

## assume that everything is under conventional basis

# Fetch all space group related data
try:
    lattice_basis = parsed_config['lattice_basis']
    space_group = parsed_config['space_group']
    space_group_origin = parsed_config['space_group_origin']
    space_group_basis = parsed_config['space_group_basis']

    # Validate data
    if not lattice_basis:
        raise ValueError("lattice_basis is empty")
    if not space_group:
        raise ValueError("space_group is empty")
    if space_group_origin is None:
        raise ValueError("space_group_origin is None")
    if space_group_basis is None:
        raise ValueError("space_group_basis is None")

    # Convert to numpy arrays
    lattice_basis_np = np.array(lattice_basis)
    origin_np = np.array(space_group_origin)
    basis_np = np.array(space_group_basis)

    # print(f"Lattice basis shape: {lattice_basis_np.shape}")
    # print(f"Space group: {space_group}")
    # print(f"Space group origin: {space_group_origin}")
    # print(f"Space group basis shape: {basis_np.shape}")

except KeyError as e:
    print(f"Error: Required key {e} not found in configuration")
    exit(key_err_code)
except ValueError as e:
    print(f"Error with configuration data: {e}")
    exit(val_err_code)

def removeCommentsAndEmptyLines(file):
    """
    Remove comments and empty lines from configuration file

    :param file: conf file path
    :return: list of cleaned lines (comments and empty lines removed)
    """
    with open(file, "r") as fptr:
        lines = fptr.readlines()

    linesToReturn = []
    for oneLine in lines:
        # Remove comments (everything after #) and remove empty lines
        oneLine = re.sub(r'#.*$', '', oneLine).strip()
        if oneLine:  # Only add non-empty lines
            linesToReturn.append(oneLine)

    return linesToReturn

# read space group file
def read_space_group(in_space_group_file, space_group_num):
    """

    :param in_space_group_file: file containing matrices of all space groups
    :param space_group_num: space group number
    :return: space group matrices (affine ) for space_group_num
    """
    contents=removeCommentsAndEmptyLines(in_space_group_file)
    # print(contents)

    # Create regex pattern to match _space_group_num_ followed by number of matrices
    space_group_pattern = r'_(\d+)_\s+(\d+)'
    matrix_elem_pattern = r'([+-]?\d+(?:/\d+)?)'

    for line_num in range(len(contents)):
        match_space_group=re.match(space_group_pattern,contents[line_num])
        if match_space_group:
            #when space group header is matched
            sgn=int(match_space_group.group(1))  # First captured group (space group number)
            if sgn==space_group_num: #if this is the sought space group
                num_matrices = int(match_space_group.group(2))  # Second captured group (number of matrices)

                # Initialize array to store space group matrices
                #first 3 columns is the linear part, last column is the translation part
                space_group_matrices = np.zeros((num_matrices, 3, 4))
                # Read the matrices following the space group header
                for matrix_idx in range(0,num_matrices):
                    matrix_line = contents[line_num + matrix_idx + 1]
                    elements = re.findall(matrix_elem_pattern, matrix_line)
                    if len(elements) != 12:
                        raise ValueError(f"Expected 12 elements, got {len(elements)} in line: {matrix_line}")

                    # Parse 12 elements (3x4 matrix flattened)
                    matrix_elements = []
                    for one_elem in elements:
                        if "/" in one_elem:
                            # Handle fractions like "1/2", "-1/2"
                            numerator, denominator = one_elem.split("/")
                            matrix_elements.append(float(numerator) / float(denominator))
                        else:
                            # Handle integers like "1", "-1", "+1"
                            matrix_elements.append(float(one_elem))

                    # Reshape to 3x4 and store
                    space_group_matrices[matrix_idx] = np.array(matrix_elements).reshape((3, 4))

                return space_group_matrices
    # If space group not found
    raise ValueError(f"Space group {space_group_num} not found in {in_space_group_file}")

def space_group_conv_to_cartesian_basis(space_group_matrices,space_group_basis):
    """

    :param space_group_matrices: space group operators (affine)
    :param space_group_basis: the basis of space_group_matrices, under Cartesian coordinates
    :return: space group operators under Cartesian coordinates
    """


in_space_group_file="./read_only/space_group_matrices_Bilbao.txt"
space_group_matrices=read_space_group(in_space_group_file,225)
# print(space_group_matrices)