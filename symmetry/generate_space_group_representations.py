import numpy as np
import sys
import json
import re
import copy
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

## assume that everything is under primitive cell basis

# Fetch all space group related data
try:
    #primitive cell lattice basis
    lattice_basis_primitive = parsed_config['lattice_basis']
    lattice_basis_primitive=np.array(lattice_basis_primitive)
    # print(f"lattice_basis={lattice_basis}")
    space_group = parsed_config['space_group']
    # print(f"space_group={space_group}, {type(space_group)}")
    space_group_origin = parsed_config['space_group_origin']
    space_group_origin=np.array(space_group_origin)
    # print(f"space_group_origin={space_group_origin}")
    space_group_basis = parsed_config['space_group_basis']
    space_group_basis=np.array(space_group_basis)
    # print(f"space_group_basis={space_group_basis}")

    # Validate data
    # if not lattice_basis:
    #     raise ValueError("lattice_basis is empty")
    # if not space_group:
    #     raise ValueError("space_group is empty")
    # if space_group_origin is None:
    #     raise ValueError("space_group_origin is None")
    # if space_group_basis is None:
    #     raise ValueError("space_group_basis is None")


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

def space_group_to_cartesian_basis(space_group_matrices,space_group_basis):
    """
    GetSymXyz(SymLvSG,LvSG) in cd/SymGroup.py
    :param space_group_matrices: space group operators (affine)
    :param space_group_basis: the basis of space_group_matrices, each row is a basis, under Cartesian coordinates
    :return: space group operators under Cartesian coordinates, dim: num_operators by 3 by 4
    """
    A=space_group_basis
    AT=space_group_basis.T
    AT_inv=np.linalg.inv(AT)

    num_operators=len(space_group_matrices)

    space_group_matrices_cartesian=np.zeros((num_operators,3,4),dtype=float)
    for j in range(0,num_operators):
        space_group_matrices_cartesian[j,:,0:3]=AT @ space_group_matrices[j,:,0:3] @AT_inv
        space_group_matrices_cartesian[j,:,3]=AT@space_group_matrices[j,:,3]

    return space_group_matrices_cartesian



def space_group_to_primitive_cell_basis(space_group_matrices_cartesian,lattice_basis_primitive):
    """

    :param space_group_matrices_cartesian: space group operators (affine) under Cartesian basis
    :param lattice_basis_primitive: primitive cell basis, each row is a basis, under Cartesian coordinates
    :return: space group operators (affine) under primitive cell basis, dim: num_operators by 3 by 4
    """
    B=lattice_basis_primitive

    BT=B.T
    BT_inv=np.linalg.inv(BT)
    num_operators=len(space_group_matrices)
    space_group_matrices_primitive=np.zeros((num_operators,3,4),dtype=float)
    for j in range(0,num_operators):
        space_group_matrices_primitive[j,:,0:3]=BT_inv@space_group_matrices_cartesian[j,:,0:3] @ BT
        space_group_matrices_primitive[j,:,3]=BT_inv @space_group_matrices_cartesian[j,:,3]
    return space_group_matrices_primitive

#TODO: double check this function
def space_group_representation_D_orbitals(R):
    """

    function GetSymD(R) in cd/SymGroup.py
    :param R: linear part of space group, under Cartesian basis
    :return: representation of space group on d orbitals
    """
    [[R_11, R_12, R_13], [R_21, R_22, R_23], [R_31, R_32, R_33]] = R
    RD = np.zeros((5,5))
    sr3 = np.sqrt(3)
    #
    RD[0,0] = R_11*R_22+R_12*R_21
    RD[0,1] = R_21*R_32+R_22*R_31
    RD[0,2] = R_11*R_32+R_12*R_31
    RD[0,3] = 2*R_11*R_12+R_31*R_32
    RD[0,4] = sr3*R_31*R_32
    #
    RD[1,0] = R_12*R_23+R_13*R_22
    RD[1,1] = R_22*R_33+R_23*R_32
    RD[1,2] = R_12*R_33+R_13*R_32
    RD[1,3] = 2*R_12*R_13+R_32*R_33
    RD[1,4] = sr3*R_32*R_33
    #
    RD[2,0] = R_11*R_23+R_13*R_21
    RD[2,1] = R_21*R_33+R_23*R_31
    RD[2,2] = R_11*R_33+R_13*R_31
    RD[2,3] = 2*R_11*R_13+R_31*R_33
    RD[2,4] = sr3*R_31*R_33
    #
    RD[3,0] = R_11*R_21-R_12*R_22
    RD[3,1] = R_21*R_31-R_22*R_32
    RD[3,2] = R_11*R_31-R_12*R_32
    RD[3,3] = (R_11**2-R_12**2 )+1/2*(R_31**2-R_32**2 )
    RD[3,4] = sr3/2*(R_31**2-R_32**2 )
    #
    RD[4,0] = 1/sr3*(2*R_13*R_23-R_11*R_21-R_12*R_22)
    RD[4,1] = 1/sr3*(2*R_23*R_33-R_21*R_31-R_22*R_32)
    RD[4,2] = 1/sr3*(2*R_13*R_33-R_11*R_31-R_12*R_32)
    RD[4,3] = 1/sr3*(2*R_13**2-R_11**2-R_12**2 )+1/sr3/2*(2*R_33**2-R_31**2-R_32**2 )
    RD[4,4] = 1/2*(2*R_33**2-R_31**2-R_32**2 )

    return RD.T

#TODO: double check this function
def space_group_representation_F_orbitals(R):
    """
    function GetSymF(R) in cd/SymGroup.py
    :param R: linear part of space group, under Cartesian basis
    :return: representation of space group on f orbitals
    """
    sr3 = np.sqrt(3); sr5 = np.sqrt(5); sr15 = np.sqrt(15)
    x1x2x3 = np.array([[ 1, 1, 1], # x3
                       [ 2, 2, 2], # y3
                       [ 3, 3, 3], # z3
                       [ 1, 1, 2], # x2y
                       [ 1, 2, 2], # xy2
                       [ 1, 1, 3], # x2z
                       [ 1, 3, 3], # xz2
                       [ 2, 2, 3], # y2z
                       [ 2, 3, 3], # yz2
                       [ 1, 2, 3]  # xyz
                       ],int)
    # 10*10, x3~xyz, a~j
    Rx1x2x3 = np.zeros((10,10))
    for i in range(10):
        n1,n2,n3 = x1x2x3[i]
        Rx1x2x3[i,0] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 1- 1,n3- 1] # a
        Rx1x2x3[i,1] = R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 2- 1,n3- 1] # b
        Rx1x2x3[i,2] = R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 3- 1,n3- 1] # c
        Rx1x2x3[i,3] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 2- 1,n3- 1] \
                       + R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 1- 1,n3- 1] \
                       + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 1- 1,n3- 1] # d
        Rx1x2x3[i,4] = R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 2- 1,n3- 1] \
                       + R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 1- 1,n3- 1] \
                       + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 2- 1,n3- 1] # e
        Rx1x2x3[i,5] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 3- 1,n3- 1] \
                       + R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 1- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 1- 1,n3- 1] # f
        Rx1x2x3[i,6] = R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 3- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 1- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 3- 1,n3- 1] # g
        Rx1x2x3[i,7] = R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 3- 1,n3- 1] \
                       + R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 2- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 2- 1,n3- 1] # h
        Rx1x2x3[i,8] = R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 3- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 2- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 3- 1,n3- 1] # i
        Rx1x2x3[i,9] = R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 3- 1,n3- 1] \
                       + R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 2- 1,n3- 1] \
                       + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 3- 1,n3- 1] \
                       + R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 1- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 2- 1,n3- 1] \
                       + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 1- 1,n3- 1] # j
    # 7*10, fz3~fy(3x2-y2), x3~xyz
    '''           [     "x3",     "y3",     "z3",    "x2y",    "xy2",    "x2z",    "xz2",    "y2z",    "yz2",    "xyz"] '''
    F = np.array([[        0,        0,   1/sr15,        0,        0,-3/2/sr15,        0,-3/2/sr15,        0,        0], # fz3
                  [ -1/2/sr5,        0,        0,        0, -1/2/sr5,        0,    2/sr5,        0,        0,        0], # fxz2
                  [        0, -1/2/sr5,        0, -1/2/sr5,        0,        0,        0,        0,    2/sr5,        0], # fyz2
                  [        0,        0,        0,        0,        0,        0,        0,        0,        0,        1], # fxyz
                  [        0,        0,        0,        0,        0,      1/2,        0,     -1/2,        0,        0], # fz(x2-y2)
                  [  1/2/sr3,        0,        0,        0,   -sr3/2,        0,        0,        0,        0,        0], # fx(x2-3y2)
                  [        0, -1/2/sr3,        0,    sr3/2,        0,        0,        0,        0,        0,        0]  # fy(3x2-y2)
                  ])
    FR = F @ Rx1x2x3 # 7*10, fz3~fy(3x2-y2), a~j
    # 7*10, fz3~fy(3x2-y2), a~j
    '''            [    "a",    "b",    "c",    "d",    "e",    "f",    "g",    "h",    "i",    "j"] '''
    CF = np.array([[      0,      0,   sr15,      0,      0,      0,      0,      0,      0,      0], # fz3
                   [      0,      0,      0,      0,      0,      0,  sr5/2,      0,      0,      0], # fxz2
                   [      0,      0,      0,      0,      0,      0,      0,      0,  sr5/2,      0], # fyz2
                   [      0,      0,      0,      0,      0,      0,      0,      0,      0,      1], # fxyz
                   [      0,      0,      3,      0,      0,      2,      0,      0,      0,      0], # fz(x2-y2)
                   [  2*sr3,      0,      0,      0,      0,      0,  sr3/2,      0,      0,      0], # fx(x2-3y2)
                   [      0, -2*sr3,      0,      0,      0,      0,      0,      0, -sr3/2,      0]  # fy(3x2-y2)
                   ])
    RF = FR @ CF.T
    return RF.T


def space_group_representation_orbitals_all(space_group_matrices_cartesian):
    """
    function GetSymOrb(SymXyz) in cd/SymGroup.py
    :param space_group_matrices_cartesian: space group matrices (affine) under Cartesian basis
    :return: space group representations on atomic orbitals
    """
    num_matrices,_,_=space_group_matrices_cartesian.shape
    print(f"num_matrices={num_matrices}", file=sys.stderr)  # Send to stderr

    # S: s
    repr_S=np.ones((num_matrices,1,1))
    # P: px,py,pz
    repr_P=copy.deepcopy(space_group_matrices_cartesian[:, :3, :3])

    # D: d orbitals (5x5 representation)
    # d_xy,d_yz,d_zx,d_(x^2-y^2 ),d_(3z^2-r^2 )
    repr_D = np.zeros((num_matrices, 5, 5))
    for i in range(num_matrices):
        R = space_group_matrices_cartesian[i, :3, :3]
        repr_D[i] = space_group_representation_D_orbitals(R)
    # F: f orbitals (7x7 representation)
    # fz3,fxz2,fyz2,fxyz,fz(x2-y2),fx(x2-3y2),fy(3x2-y2)
    repr_F = np.zeros((num_matrices, 7, 7))
    for i in range(num_matrices):
        R = space_group_matrices_cartesian[i, :3, :3]
        repr_F[i] = space_group_representation_F_orbitals(R)

    repr_S_P_D_F=[repr_S,repr_P,repr_D,repr_F]

    return repr_S_P_D_F



in_space_group_file="./read_only/space_group_matrices_Bilbao.txt"
space_group_matrices=read_space_group(in_space_group_file,space_group)
space_group_matrices_cartesian=space_group_to_cartesian_basis(space_group_matrices,space_group_basis)
space_group_matrices_primitive=space_group_to_primitive_cell_basis(space_group_matrices_cartesian,lattice_basis_primitive)

repr_S_P_D_F=space_group_representation_orbitals_all(space_group_matrices_cartesian)

space_group_representations={
    "space_group_matrices":space_group_matrices.tolist(),
    "space_group_matrices_cartesian": space_group_matrices_cartesian.tolist(),
    "space_group_matrices_primitive":space_group_matrices_primitive.tolist(),
    "repr_S_P_D_F":  [
        repr_S_P_D_F[0].tolist(),
        repr_S_P_D_F[1].tolist(),
        repr_S_P_D_F[2].tolist(),
        repr_S_P_D_F[3].tolist()
    ]

}

# output as JSON
print(json.dumps(space_group_representations, indent=2))