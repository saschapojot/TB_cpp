import numpy as np
import sys
import json

# Original file: /home/adada/Documents/pyCode/TB/cd/NbrAtom.py
# This script finds neighbors of atoms in unit cell [0,0,0] to neighboring cells

# Exit codes
json_err_code = 4
key_err_code = 5
val_err_code = 6


# ==============================================================================
# STEP 1: Read and parse JSON input from stdin
# ==============================================================================
try:
    config_json = sys.stdin.read()
    parsed_config = json.loads(config_json)

except json.JSONDecodeError as e:
    print(f"Error parsing JSON input: {e}", file=sys.stderr)
    exit(json_err_code)


# ==============================================================================
# STEP 2: Extract configuration parameters from JSON
# ==============================================================================
try:
    # Primitive cell lattice basis vectors (3x3 matrix)
    # Original variable name: Lv
    lattice_basis_primitive = parsed_config['lattice_basis']
    lattice_basis_primitive = np.array(lattice_basis_primitive)

    # Atom positions in the primitive cell
    # Original variable name: AtLv
    atom_positions_raw = parsed_config["atom_positions"]  # List of dicts

    # Number of neighbor cells to consider in each direction
    # Original variable name: Nbr
    neighbors = parsed_config["neighbors"]

    # Fractional coordinates of Bilbao space group origin
    # Original variable name: originBilbao
    space_group_origin = parsed_config["space_group_origin"]

    # Dimensionality of the system (2D or 3D)
    dim = parsed_config["dim"]

    # Determine neighbor cell range based on dimensionality
    if dim == 3:
        # 3D: study all 3 directions
        neighbor_cell_num_vec = [neighbors, neighbors, neighbors]
    else:
        # 2D: study only first 2 directions
        neighbor_cell_num_vec = [neighbors, neighbors, 0]

    # Extract atom information from list of dictionaries
    atom_position_names = []
    atom_types = []
    atom_coordinates_frac = []

    for atom in atom_positions_raw:
        atom_position_names.append(atom['position_name'])
        atom_types.append(atom['atom_type'])
        atom_coordinates_frac.append(atom['fractional_coordinates'])

    atom_coordinates_frac = np.array(atom_coordinates_frac)  # Convert to numpy array


except KeyError as e:
    print(f"Error: Required key {e} not found in configuration", file=sys.stderr)
    exit(key_err_code)
except ValueError as e:
    print(f"Error with configuration data: {e}", file=sys.stderr)
    exit(val_err_code)


# ==============================================================================
# STEP 3: Generate all atoms in neighboring cells
# ==============================================================================
# Neighbor range: [[-N0, N0], [-N1, N1], [-N2, N2]]
N0, N1, N2 = neighbor_cell_num_vec
total_num_cells = (2*N0+1) * (2*N1+1) * (2*N2+1)
atom_num_in_1_cell = len(atom_position_names)

# Create list of dictionaries for all atoms in the supercell
atoms_in_cells = []

count = 0
for n0 in range(-N0, N0+1):
    for n1 in range(-N1, N1+1):
        for n2 in range(-N2, N2+1):
            for j in range(atom_num_in_1_cell):
                atom_dict = {
                    'cell': [n0, n1, n2],                    # Which unit cell
                    'atom_index': j,                         # Which atom in that cell
                    'position_name': atom_position_names[j],
                    'atom_type': atom_types[j],
                    'frac_coords': atom_coordinates_frac[j].tolist(),  # Convert to list
                    'global_index': count                    # Unique identifier
                }
                atoms_in_cells.append(atom_dict)
                count += 1


# ==============================================================================
# STEP 4: Calculate Cartesian coordinates
# ==============================================================================
# Cartesian coordinates for atoms in the origin cell [0,0,0]
# Formula: cart_coords = frac_coords @ lattice_basis
atom_coordinates_cart_origin = atom_coordinates_frac @ lattice_basis_primitive

# Cartesian coordinates for all atoms in neighboring cells
atoms_cart_coords = np.zeros((len(atoms_in_cells), 3))
for i, atom in enumerate(atoms_in_cells):
    n0, n1, n2 = atom['cell']
    j = atom['atom_index']

    # Cell translation vector: n0*a1 + n1*a2 + n2*a3
    cell_translation = (n0*lattice_basis_primitive[0,:] +
                        n1*lattice_basis_primitive[1,:] +
                        n2*lattice_basis_primitive[2,:])

    # Atom position within the cell
    atom_position = atom_coordinates_frac[j,:] @ lattice_basis_primitive

    # Total position = cell_translation + atom_position_in_cell
    atoms_cart_coords[i,:] = cell_translation + atom_position


# ==============================================================================
# STEP 5: Create atom pairs and calculate distances
# ==============================================================================
# Pair atoms in origin cell [0,0,0] with atoms in all neighboring cells
atom_pairs = []
count = 0
num_digits = 6  # Truncate distances to 6 decimal places

for i in range(atom_num_in_1_cell):  # Atoms in origin cell [0,0,0]
    for j, atom_j in enumerate(atoms_in_cells):  # Atoms in all cells

        # Calculate distance between the two atoms
        displacement = atoms_cart_coords[j] - atom_coordinates_cart_origin[i]
        distance = np.linalg.norm(displacement)
        distance = round(distance, num_digits)  # Truncate to num_digits decimal places

        # Create pair dictionary with complete information
        pair_dict = {
            'pair_index': count,
            'distance': distance,
            'atom_at_center_cell': {
                'cell': [0, 0, 0],
                'atom_index': i,
                'position_name': atom_position_names[i],
                'atom_type': atom_types[i],
                'frac_coords': atom_coordinates_frac[i].tolist(),      # Convert to list
                'cart_coords': atom_coordinates_cart_origin[i].tolist()  # Convert to list
            },
            'atom_at_neighbor_cell': {
                'cell': atom_j['cell'],  # Already a list, no need to copy
                'atom_index': atom_j['atom_index'],
                'position_name': atom_j['position_name'],
                'atom_type': atom_j['atom_type'],
                'frac_coords': atom_j['frac_coords'],  # Already converted to list in STEP 3
                'cart_coords': atoms_cart_coords[j].tolist(),  # Convert to list
                'global_index': atom_j['global_index']
            }
        }
        atom_pairs.append(pair_dict)
        count += 1


# ==============================================================================
# STEP 6: Find unique distances and classify pairs
# ==============================================================================
# Extract all distances from pairs
all_distances = [pair['distance'] for pair in atom_pairs]

# Find unique distances (duplicates automatically removed by set())
# Sort them in ascending order
unique_distances = sorted(set(all_distances))

# Create lookup dictionary: distance -> index
# Example: {1.5: 0, 2.1: 1, 3.0: 2}
distance_to_index = {dist: idx for idx, dist in enumerate(unique_distances)}

# Add distance classification to each pair
for pair in atom_pairs:
    distance = pair['distance']
    pair['distance_index'] = distance_to_index[distance]              # Index of distance bin
    pair['unique_distance_value'] = unique_distances[distance_to_index[distance]]  # Canonical distance value


# ==============================================================================
# STEP 7: Output results as JSON
# ==============================================================================
print(json.dumps(atom_pairs), file=sys.stdout)