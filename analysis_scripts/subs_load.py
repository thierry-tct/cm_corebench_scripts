import os
import sys

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
    mutants_to_killingtests = sm_mat.query_active_columns_of_rows()
    tests_to_killed_mutants = sm_mat.query_active_rows_of_columns()
    
    all_tests = sm_mat.get_nonkey_colname_list()

    # return data
    return all_tests, mutants_to_killingtests, tests_to_killed_mutants
#~ def load()
