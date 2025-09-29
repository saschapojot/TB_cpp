import numpy as np
import sys
import re
import json
import copy

#this script checks if orbitals given by user is complete
# if orbitals are not complete, then orbitals related by symmetry are added

#original file: /home/adada/Documents/pyCode/TB/ck/CheckOrbCpl.py
json_err_code=4
try:
    combined_input_json=sys.stdin.read()
    combined_input=json.loads(combined_input_json)
except json.JSONDecodeError as e:
    print(f"Error parsing JSON input: {e}", file=sys.stderr)
    exit(json_err_code)

orbital_map = {
    # 1s (index 0)
    '1s': 0,

    # 2s, 2p (indices 1-4)
    '2s': 1,
    '2px': 2, '2py': 3, '2pz': 4,

    # 3s, 3p, 3d (indices 5-13)
    '3s': 5,
    '3px': 6, '3py': 7, '3pz': 8,
    '3dxy': 9, '3dyz': 10, '3dxz': 11, '3dx2-y2': 12, '3dz2': 13,

    # 4s, 4p, 4d, 4f (indices 14-29)
    '4s': 14,
    '4px': 15, '4py': 16, '4pz': 17,
    '4dxy': 18, '4dyz': 19, '4dxz': 20, '4dx2-y2': 21, '4dz2': 22,
    '4fxyz': 23, '4fz3': 24, '4fxz2': 25, '4fyz2': 26,
    '4fz(x2-y2)': 27, '4fx(x2-3y2)': 28, '4fy(3x2-y2)': 29,

    # 5s, 5p, 5d, 5f (indices 30-45)
    '5s': 30,
    '5px': 31, '5py': 32, '5pz': 33,
    '5dxy': 34, '5dyz': 35, '5dxz': 36, '5dx2-y2': 37, '5dz2': 38,
    '5fxyz': 39, '5fz3': 40, '5fxz2': 41, '5fyz2': 42,
    '5fz(x2-y2)': 43, '5fx(x2-3y2)': 44, '5fy(3x2-y2)': 45,

    # 6s, 6p, 6d, 6f (indices 46-61)
    '6s': 46,
    '6px': 47, '6py': 48, '6pz': 49,
    '6dxy': 50, '6dyz': 51, '6dxz': 52, '6dx2-y2': 53, '6dz2': 54,
    '6fxyz': 55, '6fz3': 56, '6fxz2': 57, '6fyz2': 58,
    '6fz(x2-y2)': 59, '6fx(x2-3y2)': 60, '6fy(3x2-y2)': 61,

    # 7s, 7p, 7d, 7f (indices 62-77)
    '7s': 62,
    '7px': 63, '7py': 64, '7pz': 65,
    '7dxy': 66, '7dyz': 67, '7dxz': 68, '7dx2-y2': 69, '7dz2': 70,
    '7fxyz': 71, '7fz3': 72, '7fxz2': 73, '7fyz2': 74,
    '7fz(x2-y2)': 75, '7fx(x2-3y2)': 76, '7fy(3x2-y2)': 77,
}

# Extract the two main components
parsed_config = combined_input["parsed_config"]
space_group_representations = combined_input["space_group_representations"]


# Extract space group data
space_group_matrices = np.array(space_group_representations["space_group_matrices"])
space_group_matrices_cartesian = np.array(space_group_representations["space_group_matrices_cartesian"])
space_group_matrices_primitive = np.array(space_group_representations["space_group_matrices_primitive"])

# Extract representation matrices for different orbital types
repr_s, repr_p, repr_d, repr_f = space_group_representations["repr_s_p_d_f"]

repr_s_np = np.array(repr_s)
repr_p_np = np.array(repr_p)
repr_d_np = np.array(repr_d)
repr_f_np = np.array(repr_f)

#num_operations is N in notes
num_operations,_,_=repr_s_np.shape
#combine repr_s_np, repr_p_np, repr_d_np, repr_f_np

#IndSPDF
#1s,2s,2p,3s,...,7f
orbital_nums_spdf = np.array([1, #1s
                              1,3, #2s, 2p
                              1,3,5,#3s, 3p, 3d
                              1,3,5,7,#4s, 4p, 4d, 4f
                              1,3,5,7,#5s, 5p, 5d, 5f
                              1,3,5,7,#6s, 6p, 6d, 6f
                              1,3,5,7,#7s, 7p, 7d, 7f
                              ])
print(f"len(orbital_nums_spdf)={len(orbital_nums_spdf)}",file=sys.stderr)
#78
orbital_max_dim=np.sum(orbital_nums_spdf)
print(f"orbital_max_dim={orbital_max_dim}",file=sys.stderr)

#SymSPDF
spdf_combined=np.zeros((num_operations,orbital_max_dim,orbital_max_dim))





def build_orbital_vectors(parsed_config):
    """
    Build a length 78 orbital vector for each atom in the configuration

    :param parsed_config: Dictionary containing atom types and their orbitals
    :return: Dictionary mapping atom position names to their orbital vectors
    """
    # parsed_config={'name': 'hBN', 'dim': 3,
    # 'spin': 'False', 'neighbors': 1, 'atom_type_num': 2,
    # 'lattice_type': 'primitive',
    # 'lattice_basis': [[0.5, -0.8660254, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 1.0]],
    # 'space_group': 187, 'space_group_origin': [0.0, 0.0, 0.0],
    # 'space_group_basis': [[1.0, 0.0, 0.0], [-0.5, 0.86602540378, 0.0], [0.0, 0.0, 1.0]],
    # 'atom_types': {'B': {'count': 1, 'orbitals': ['2pz', '2s']}, 'N': {'count': 1, 'orbitals': ['2px', '2py']}}, 'atom_positions': [{'position_name': 'B', 'atom_type': 'B', 'coordinates': [0.33333333, 0.66666666, 0.0]}, {'position_name': 'N', 'atom_type': 'N', 'coordinates': [0.66666666, 0.33333333, 0.0]}]}
    # Create orbital vectors for each atom type
    # Build vectors for each atom position
    atom_orbital_vectors = {}

    for atom in parsed_config['atom_positions']:
        position_name = atom['position_name']
        atom_type = atom['atom_type']

        # Get orbitals for this atom type
        orbitals = parsed_config['atom_types'][atom_type]['orbitals']

        # Create orbital vector
        orbital_vector = np.zeros(78)
        for orbital in orbitals:
            if orbital in orbital_map:
                orbital_vector[orbital_map[orbital]] = 1
            else:
                print(f"Warning: Orbital '{orbital}' for atom '{position_name}' not recognized")

        atom_orbital_vectors[position_name] = orbital_vector

    return atom_orbital_vectors


# Fill the diagonal blocks of spdf_combined
for j in range(0,num_operations):
    current_idx = 0
    # Iterate through orbital_nums_spdf to place each block
    for i, block_size in enumerate(orbital_nums_spdf):
        if block_size == 1:  # s orbital
            spdf_combined[j, current_idx:current_idx+1, current_idx:current_idx+1] = repr_s_np[j]
        elif block_size == 3:  # p orbital
            spdf_combined[j, current_idx:current_idx+3, current_idx:current_idx+3] = repr_p_np[j]

        elif block_size == 5:  # d orbital
            spdf_combined[j, current_idx:current_idx+5, current_idx:current_idx+5] = repr_d_np[j]

        elif block_size == 7:  # f orbital
            spdf_combined[j, current_idx:current_idx+7, current_idx:current_idx+7] = repr_f_np[j]
        current_idx += block_size

#IndNonZero
non_zero_spdf_combined=np.sum(np.abs(spdf_combined),axis=0) >1e-6
# print(f"parsed_config={parsed_config}",file=sys.stderr)
atom_orbital_vectors=build_orbital_vectors(parsed_config)
# print(f"atom_orbital_vectors={atom_orbital_vectors}",file=sys.stderr)

# Update atom_orbital_vectors based on symmetry coupling
updated_atom_orbital_vectors = {}
added_orbitals_dict = {}  # Dictionary to store added orbitals for each atom

for atom_name, orbital_vector in atom_orbital_vectors.items():
    # Start with a copy of the original vector
    updated_vector = copy.deepcopy(orbital_vector)
    # Find indices where the orbital vector has 1 (active orbitals)
    active_orbital_indices = np.where(orbital_vector == 1)[0]
    # For each active orbital
    for orbital_idx in active_orbital_indices:
        # Find all orbitals coupled to this one by symmetry
        # Look at column orbital_idx in non_zero_spdf_combined
        coupled_orbital_indices = np.where(non_zero_spdf_combined[:, orbital_idx])[0]
        # Set all coupled positions to 1
        updated_vector[coupled_orbital_indices] = 1
    updated_atom_orbital_vectors[atom_name] = updated_vector

    # Report which orbitals were added
    added_indices = np.where((updated_vector == 1) & (orbital_vector == 0))[0]
    if len(added_indices) > 0:
        added_orbitals = [k for k, v in orbital_map.items() if v in added_indices]
        added_orbitals_dict[atom_name] = added_orbitals  # Store in dictionary
        # print(f"Atom {atom_name}: Added orbitals {added_orbitals} by symmetry", file=sys.stderr)
    else:
        added_orbitals_dict[atom_name] = []  # Empty list if no orbitals added

# Replace the original vectors with updated ones
atom_orbital_vectors = updated_atom_orbital_vectors
# print(f"atom_orbital_vectors={atom_orbital_vectors}", file=sys.stderr)
# print(f"added_orbitals_dict={added_orbitals_dict}", file=sys.stderr)