import numpy as np
import sys

import json

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

num_operations,_,_=repr_s_np.shape
#combine repr_s_np, repr_p_np, repr_d_np, repr_f_np

#IndSPDF
orbital_nums_spdf = np.array([1,1,3,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7])

#94
orbital_max_dim=np.sum(orbital_nums_spdf)
# print(orbital_max_dim,file=sys.stderr)

#SymSPDF
spdf_combined=np.zeros((num_operations,orbital_max_dim,orbital_max_dim))


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


