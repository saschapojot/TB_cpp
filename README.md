


###########################
Part I, proprocessing
1. python preprocessing.py ./path/to/xxx.conf
    (i) parse conf file
        ./parse_files/parse_conf.py
    (ii) check if conf file data are valid
        ./parse_files/sanity_check.py
    (iii)   read space group matrices (Bilbao),
            convert space group matrices (affine) from conventional basis to Cartesian basis