import re
import subprocess
import sys
import os
import json
import numpy as np


argErrCode=20
if (len(sys.argv)!=2):
    print("wrong number of arguments")
    print("example: python preprocessing.py /path/to/mc.conf")
    exit(argErrCode)

confFileName=str(sys.argv[1])


#########################################################################################
# parse conf file, get parsed_config
confResult = subprocess.run(["python3", "./parse_files/parse_conf.py", confFileName], capture_output=True, text=True)

# Check if the subprocess ran successfully
if confResult.returncode != 0:
    print("Error running parse_conf.py:")
    print(confResult.stderr)
    exit(confResult.returncode)

# Parse the JSON output
try:
    parsed_config = json.loads(confResult.stdout)
    # print(parsed_config)
    print("=" * 60)
    print("COMPLETE PARSED CONFIGURATION")
    print("=" * 60)

    # Print basic configuration
    print(f"Name: {parsed_config['name']}")
    print(f"Dimensions: {parsed_config['dim']}")
    print(f"Spin: {parsed_config['spin']}")
    print(f"Neighbors: {parsed_config['neighbors']}")
    print(f"Atom Type Number: {parsed_config['atom_type_num']}")
    print(f"Lattice Type: {parsed_config['lattice_type']}")
    print(f"Space Group: {parsed_config['space_group']}")

    # Print space group origin
    print(f"Space Group Origin: [{', '.join(map(str, parsed_config['space_group_origin']))}]")

    # Print lattice basis
    print("Lattice Basis:")
    for i, vector in enumerate(parsed_config['lattice_basis']):
        print(f"  Vector {i+1}: [{', '.join(map(str, vector))}]")

    # Print space group basis
    print("Space Group Basis:")
    for i, vector in enumerate(parsed_config['space_group_basis']):
        print(f"  Vector {i+1}: [{', '.join(map(str, vector))}]")

    # Print atom types
    print("\nAtom Types:")
    for atom_type, info in parsed_config['atom_types'].items():
        print(f"  {atom_type}:")
        print(f"    Count: {info['count']}")
        print(f"    Orbitals: {info['orbitals']}")

    # Print atom positions
    print(f"\nAtom Positions (Total: {len(parsed_config['atom_positions'])}):")
    for i, pos in enumerate(parsed_config['atom_positions']):
        print(f"  Position {i+1}:")
        print(f"    Name: {pos['position_name']}")
        print(f"    Atom Type: {pos['atom_type']}")
        print(f"    Coordinates: [{', '.join(map(str, pos['coordinates']))}]")


except json.JSONDecodeError as e:
    print("Error parsing JSON output from parse_conf.py:")
    print(f"JSON Error: {e}")
    print("Raw output was:")
    print(confResult.stdout)
    exit(1)
# get parsed_config
# Convert parsed_config back to JSON string for passing to subprocess
config_json = json.dumps(parsed_config)
##end parsing conf
#########################################################################################


#########################################################################################
# Pass parsed_config to sanity_check.py
print("\n" + "=" * 60)
print("RUNNING SANITY CHECK")
print("=" * 60)

# Run sanity_check.py and pass the JSON data via stdin
sanity_result = subprocess.run(
    ["python3", "./parse_files/sanity_check.py"],
    input=config_json,
    capture_output=True,
    text=True
)
#open the following line for debugging
# print(f"Output: {sanity_result.stdout}")
print(f"Exit code: {sanity_result.returncode}")
# Check sanity check results
if sanity_result.returncode != 0:
    print("Sanity check failed!")
    print(f"return code={sanity_result.returncode}")
    print("Error output:")
    print(sanity_result.stdout)
    exit(sanity_result.returncode)
else:
    print("Sanity check passed!")
    print("Output:")
    print(sanity_result.stdout)
#end sanity check
#########################################################################################


#########################################################################################
# computing space group representations

print("\n" + "=" * 60)
print("COMPUTING SPACE GROUP REPRESENTATIONS")
print("=" * 60)

# Run generate_space_group_representations.py and pass the JSON data via stdin
sgr_result = subprocess.run(
    ["python3", "./symmetry/generate_space_group_representations.py"],
    input=config_json,
    capture_output=True,
    text=True
)

print(f"Exit code: {sgr_result.returncode}")
# Check space group representations results
if sgr_result.returncode != 0:
    print("Space group representations generation failed!")
    print(f"return code={sgr_result.returncode}")
    print("Error output:")
    print(sgr_result.stderr)
    print("Standard output:")
    print(sgr_result.stdout)
    exit(sgr_result.returncode)
else:
    print("Space group representations generated successfully!")

    # Parse the JSON output containing space group representations
    try:
        space_group_representations = json.loads(sgr_result.stdout)
        # print(space_group_representations)
        print("\n" + "=" * 60)
        print("SPACE GROUP REPRESENTATIONS SUMMARY")
        print("=" * 60)

        # Get dimensions
        num_operations = len(space_group_representations["space_group_matrices"])

        print(f"Number of space group operations: {num_operations}")

        # Print dimensions of each representation
        repr_S, repr_P, repr_D, repr_F = space_group_representations["repr_S_P_D_F"]

        print(f"S orbital representations: {len(repr_S)} operations × {len(repr_S[0])}×{len(repr_S[0][0])} matrices")
        print(f"P orbital representations: {len(repr_P)} operations × {len(repr_P[0])}×{len(repr_P[0][0])} matrices")
        print(f"D orbital representations: {len(repr_D)} operations × {len(repr_D[0])}×{len(repr_D[0][0])} matrices")
        print(f"F orbital representations: {len(repr_F)} operations × {len(repr_F[0])}×{len(repr_F[0][0])} matrices")

        # Convert back to NumPy arrays if needed for further processing
        space_group_matrices = np.array(space_group_representations["space_group_matrices"])
        space_group_matrices_cartesian = np.array(space_group_representations["space_group_matrices_cartesian"])
        space_group_matrices_primitive = np.array(space_group_representations["space_group_matrices_primitive"])

        repr_S_np = np.array(repr_S)
        repr_P_np = np.array(repr_P)
        repr_D_np = np.array(repr_D)
        repr_F_np = np.array(repr_F)

        print("\nSpace group representations loaded and converted to NumPy arrays.")
        print(f"Available matrices:")
        print(f"  - space_group_matrices: {space_group_matrices.shape}")
        print(f"  - space_group_matrices_cartesian: {space_group_matrices_cartesian.shape}")
        print(f"  - space_group_matrices_primitive: {space_group_matrices_primitive.shape}")
        print(f"  - S orbital representations: {repr_S_np.shape}")
        print(f"  - P orbital representations: {repr_P_np.shape}")
        print(f"  - D orbital representations: {repr_D_np.shape}")
        print(f"  - F orbital representations: {repr_F_np.shape}")

        # Example: Print first space group operation matrix (identity)
        # print(f"\nFirst space group operation (should be identity):")
        # print("Cartesian representation:")
        # print(space_group_matrices_cartesian[0])

    except json.JSONDecodeError as e:
        print("Error parsing JSON output from space group representations:")
        print(f"JSON Error: {e}")
        print("Raw output was:")
        print(sgr_result.stdout)
        exit(1)

    except KeyError as e:
        print(f"Missing key in space group representations output: {e}")
        print("Available keys:", list(space_group_representations.keys()) if 'space_group_representations' in locals() else "Could not parse JSON")
        exit(1)
# get dict space_group_representations
# end computing space group representations
#########################################################################################