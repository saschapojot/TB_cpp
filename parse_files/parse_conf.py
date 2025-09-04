import re
import sys
import json
import os

# Parse xxx.conf file
fmtErrStr = "format error: "
formatErrCode = 1
valueMissingCode = 2
paramErrCode = 3
fileNotExistErrCode = 4

if len(sys.argv) != 2:
    print("wrong number of arguments.")
    print("usage: python parse_conf.py /path/to/xxx.conf")
    exit(paramErrCode)

conf_file = sys.argv[1]

# Check if file exists
if not os.path.exists(conf_file):
    print(f"file not found: {conf_file}")
    exit(fileNotExistErrCode)

# Pattern definitions
key_value_pattern = r'^([^=\s]+)\s*=\s*([^=]+)\s*$'
float_pattern = r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?"
# Updated to allow optional spaces after semicolon
atom_oribital_pattern = r'^([A-Za-z]+\d*)\s*=\s*(\d+)\s*;\s*((?:s|px|py|pz|dxy|dxz|dyz|dx2-y2|dz2|fxyz|fx3-3xy2|f3x2y-y3|fxz2|fyz2|fx2-y2z|fz3)(?:\s*,\s*(?:s|px|py|pz|dxy|dxz|dyz|dx2-y2|dz2|fxyz|fx3-3xy2|f3x2y-y3|fxz2|fyz2|fx2-y2z|fz3))*)\s*$'
name_pattern = r'^name\s*=\s*([a-zA-Z0-9_-]+)\s*$'
dim_pattern = r"^dim\s*=\s*(\d+)\s*$"
neighbors_pattern = r"^neighbors\s*=\s*(\d+)\s*$"
atom_type_num_pattern = r"^atom_type_num\s*=\s*(\d+)\s*$"
# Updated to handle atom names with numbers (like O1, O2, etc.)
atom_position_pattern_3d = rf'^([a-zA-Z]+\d*)_position_coefs\s*=\s*({float_pattern})\s*,\s*({float_pattern})\s*,\s*({float_pattern})\s*$'
lattice_basis_pattern = rf'^lattice_basis\s*=\s*({float_pattern}\s*,\s*{float_pattern}\s*,\s*{float_pattern})(?:\s*;\s*({float_pattern}\s*,\s*{float_pattern}\s*,\s*{float_pattern})){{2}}\s*$'
space_group_pattern = r"^space_group\s*=\s*(\d+)\s*$"
space_group_origin_pattern = rf"^space_group_origin\s*=\s*({float_pattern})\s*,\s*({float_pattern})\s*,\s*({float_pattern})\s*$"
space_group_basis_pattern = rf"^space_group_basis\s*=\s*({float_pattern}\s*,\s*{float_pattern}\s*,\s*{float_pattern})(?:\s*;\s*({float_pattern}\s*,\s*{float_pattern}\s*,\s*{float_pattern})){{2}}\s*$"
spin_pattern = r'^spin\s*=\s*((?i:true|false))\s*$'
lattice_type_pattern = r'^lattice_type\s*=\s*((?i:primitive|conventional))\s*$'

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

def parseConfContents(file):
    """
    Parse configuration file contents

    :param file: conf file path
    :return: dictionary containing parsed configuration
    """
    linesWithCommentsRemoved = removeCommentsAndEmptyLines(file)

    # Initialize result dictionary
    config = {
        'name': '',
        'dim': '',
        'spin': '',
        'neighbors': '',
        'atom_type_num': '',
        'lattice_type': '',
        'lattice_basis': '',
        'space_group': '',
        'space_group_origin': '',
        'space_group_basis': '',
        'atom_types': {},  # Changed to dictionary to store atom type definitions
        'atom_positions': []  # Store all atom positions with their types
    }

    for oneLine in linesWithCommentsRemoved:
        matchLine = re.match(key_value_pattern, oneLine)

        # If this line has format key=value
        if matchLine:
            # Match name
            name_match = re.match(name_pattern, oneLine)
            if name_match:
                config['name'] = name_match.group(1)
                # print(f"name={config['name']}")
                continue

            # Match dim
            dim_match = re.match(dim_pattern, oneLine)
            if dim_match:
                config['dim'] = int(dim_match.group(1))
                # print(f"dim={config['dim']}")
                continue

            # Match spin
            spin_match = re.match(spin_pattern, oneLine)
            if spin_match:
                config['spin'] = spin_match.group(1)
                # print(f"spin={config['spin']}")
                continue

            # Match neighbors
            match_neighbors = re.match(neighbors_pattern, oneLine)
            if match_neighbors:
                config['neighbors'] = int(match_neighbors.group(1))
                # print(f"neighbors={config['neighbors']}")
                continue

            # Match atom_type_num
            match_atom_type_num = re.match(atom_type_num_pattern, oneLine)
            if match_atom_type_num:
                config['atom_type_num'] = int(match_atom_type_num.group(1))
                # print(f"atom_type_num={config['atom_type_num']}")
                continue

            # Match lattice_type
            match_lattice_type = re.match(lattice_type_pattern, oneLine)
            if match_lattice_type:
                config['lattice_type'] = match_lattice_type.group(1)
                # print(f"lattice_type={config['lattice_type']}")
                continue

            # Match lattice_basis
            match_lattice_basis = re.match(lattice_basis_pattern, oneLine)
            if match_lattice_basis:
                # Extract the full lattice basis value after the = sign
                full_value = oneLine.split('=')[1].strip()
                # Split into 3 vectors
                vectors = []
                for vector in full_value.split(';'):
                    coords = [float(x.strip()) for x in vector.strip().split(',')]
                    vectors.append(coords)
                config['lattice_basis'] = vectors
                # print(f"lattice_basis={vectors}")
                continue

            # Match space group
            match_space_group = re.match(space_group_pattern, oneLine)
            if match_space_group:
                config['space_group'] = int(match_space_group.group(1))
                # print(f"space_group={config['space_group']}")
                continue

            # Match space group origin
            match_space_group_origin = re.match(space_group_origin_pattern, oneLine)
            if match_space_group_origin:
                x_coord = float(match_space_group_origin.group(1))
                y_coord = float(match_space_group_origin.group(2))
                z_coord = float(match_space_group_origin.group(3))
                config['space_group_origin'] = [x_coord, y_coord, z_coord]
                # print(f"space_group_origin={x_coord},{y_coord},{z_coord}")
                continue

            # Match space group basis
            match_space_group_basis = re.match(space_group_basis_pattern, oneLine)
            if match_space_group_basis:
                full_value = oneLine.split('=')[1].strip()
                vectors = []
                for vector in full_value.split(';'):
                    coords = [float(x.strip()) for x in vector.strip().split(',')]
                    vectors.append(coords)
                config['space_group_basis'] = vectors
                # print(f"space_group_basis={vectors}")
                continue

            # Match atom type definitions (A=1;s, B=1;s, O=3;s)
            atom_match = re.match(atom_oribital_pattern, oneLine)
            if atom_match:
                atom_type = atom_match.group(1)  # A, B, O
                atom_count = int(atom_match.group(2))  # 1, 1, 3
                orbitals = atom_match.group(3).strip()  # s

                config['atom_types'][atom_type] = {
                    'count': atom_count,
                    'orbitals': orbitals
                }
                # print(f"atom_type={atom_type}, count={atom_count}, orbitals={orbitals}")
                continue

            # Match position coefficients (A_position_coefs, B_position_coefs, O1_position_coefs, etc.)
            coefs_match = re.match(atom_position_pattern_3d, oneLine)
            if coefs_match:
                position_name = coefs_match.group(1)  # A, B, O1, O2, O3
                x_coord = float(coefs_match.group(2))
                y_coord = float(coefs_match.group(3))
                z_coord = float(coefs_match.group(4))

                # Determine the base atom type (remove numbers: O1 -> O, O2 -> O, etc.)
                base_atom_type = re.sub(r'\d+$', '', position_name)

                position_info = {
                    'position_name': position_name,
                    'atom_type': base_atom_type,
                    'coordinates': [x_coord, y_coord, z_coord]
                }
                config['atom_positions'].append(position_info)
                # print(f"position_name={position_name}, atom_type={base_atom_type}, coordinates={x_coord},{y_coord},{z_coord}")
                continue

            # If no pattern matched, log unrecognized line
            print(f"Unrecognized key-value line: {oneLine}")

        else:
            print("line: " + oneLine + " is discarded.")

    return config

# Parse the configuration and return structured data
parsed_config = parseConfContents(conf_file)

# output the parsed config as JSON
print(json.dumps(parsed_config, indent=2))