import re
import sys
import json
import os

#parse xxx.conf file
fmtErrStr = "format error: "
formatErrCode=1
valueMissingCode = 2
paramErrCode = 3
fileNotExistErrCode = 4

if (len(sys.argv) != 2):
    print("wrong number of arguments.")
    print("usage: python parse_conf.py /path/to/xxx.conf")
    exit(paramErrCode)

conf_file = sys.argv[1]

#check if file exists
if not os.path.exists(conf_file):
    print(f"file not found: {conf_file}")
    exit(fileNotExistErrCode)

##################################################
#a line to be read should have pattern:
key_value_pattern =  r'^([^=\s]+)\s*=\s*([^=]+)$'
#matches: key=value
# my-key = value
# my.key= value with spaces
#KEY_NAME =   value123
# complex-key.name = some value here
#does not match:
# = value              # No key
# key =                # No value
# key value            # No equals sign
# my key = value       # Space in key (not allowed)
# key==value          # Multiple equals (first = becomes part of key)
##################################################

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

    :param file: conf file
    :return:
    """
    linesWithCommentsRemoved = removeCommentsAndEmptyLines(file)
    name_str=""
    dim_str=""
    spin_str=""
    neighbors_str=""
    lattice_type_str=""
    lattice_basis_str=""
    space_group_str=""
    atom_type_num_str=""
    atom_orbitals_list_str=""
    atom_positions_list_str=""
    for oneLine in linesWithCommentsRemoved:
        matchLine=re.match(key_value_pattern,oneLine)

        #if this line has format key=value
        if matchLine:
            key = matchLine.group(1).strip()
            value = matchLine.group(2).strip()

        #if this line does not have format key=value, this line is discarded
        else:
            print("line: " + oneLine + " is discarded.")
            continue






