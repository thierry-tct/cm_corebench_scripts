import os
import sys
import datetime

try:
    import muteria.common.fs as common_fs
    import muteria.common.matrices as common_matrices
except ImportError:
    sys.path.insert(0, '/home/mfi/mytools/muteria')
    import muteria.common.fs as common_fs
    import muteria.common.matrices as common_matrices
    sys.path.pop(0)

equal_log = common_matrices.OutputLogData.outlogdata_equiv

def load(matrice_file):
    # load matrices and compute mutant killtest mapping
    sm_mat = common_matrices.ExecutionMatrix(matrice_file)
    
    print ("[{}] Getting mutants_to_killingtests ...".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    mutants_to_killingtests = sm_mat.query_active_columns_of_rows()
    
    print ("[{}] Getting tests_to_killed_mutants ...".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    tests_to_killed_mutants = sm_mat.query_active_rows_of_columns()
    
    all_tests = sm_mat.get_nonkey_colname_list()
    
    print ("[{}] Loaded".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    # return data
    return all_tests, mutants_to_killingtests, tests_to_killed_mutants
#~ def load()
