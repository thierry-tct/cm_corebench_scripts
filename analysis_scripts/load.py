
# load data for a particular fault
# Input: fault's res dir tar file
# Output: all tests, tests revealing fault, relevant mutants to tests, mutants to killing test map, tests to killed mutants map

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

def is_relevant(preCommitVersion, preCommitMutant, postCommitVersion, postCommitMutant):
    mutant_relevant = False
    if ( (not equal_log(postCommitMutant, postCommitVersion)) and (not equal_log(postCommitMutant, preCommitMutant))):
        mutant_relevant = True
    return mutant_relevant
#~ def is_relevant()

def get_relevant_mutants_to_relevant_tests(pre_orig_ol, pre_muts_ol, post_orig_ol, post_muts_ol):
    for prog, test_out in pre_orig_ol.get_zip_objective_and_data():
        assert prog == 'program', "invalid program out (pre)"
        orig_pre_test_out = test_out
    for prog, test_out in post_orig_ol.get_zip_objective_and_data():
        assert prog == 'program', "invalid program out (post)"
        orig_post_test_out = test_out
    muts_pre_test_out = {mut: test_out for mut, test_out in pre_muts_ol.get_zip_objective_and_data()}
    muts_post_test_out = {mut: test_out for mut, test_out in post_muts_ol.get_zip_objective_and_data()}

    if set(orig_pre_test_out.keys()) != set(orig_post_test_out.keys()):
        pre_post = set(orig_pre_test_out.keys()) - set(orig_post_test_out.keys())
        post_pre = set(orig_post_test_out.keys()) - set(orig_pre_test_out.keys())
        if len(pre_post) > 0:
            print("WARNING: Some tests (flaky in post) in pre are not in post: {}".format(pre_post))
        if len(post_pre) > 0:
            print("WARNING: Some tests (flaky in pre) in post are not in pre: {}".format(pre_post))

    # Only check mutants from post because a condition to be relevant is killed
    #assert set(muts_pre_test_out.keys()) == set(muts_post_test_out.keys()), "mismatch between pre and post (muts)"
    
    relevant_mutants_to_relevant_tests = {}
    for mut in muts_post_test_out.keys():
        # Check whether relevant
        for pt in muts_post_test_out[mut]:
            if pt not in orig_pre_test_out:
                # Flaky in pre but not in post
                continue
            Out_O2 = orig_post_test_out[pt]
            Out_M2 = muts_post_test_out[mut][pt]
            Out_O1 = orig_pre_test_out[pt]
            if mut in muts_pre_test_out and pt in muts_pre_test_out[mut]:
                Out_M1 = muts_pre_test_out[mut][pt]
            else:
                Out_M1 = Out_O1
            
            if is_relevant(Out_O1, Out_M1, Out_O2, Out_M2):
                if mut not in relevant_mutants_to_relevant_tests:
                    relevant_mutants_to_relevant_tests[mut] = []
                relevant_mutants_to_relevant_tests[mut].append(pt)

    return relevant_mutants_to_relevant_tests
#~ def get_relevant_mutants_to_relevant_tests()

def load(resdir, fault_revealing=True):
    # load fault revealing tests
    fail_tests = None
    if fault_revealing:
        f_file = os.path.join(resdir, "fail_test_checking", "fault_reveling_tests.txt")
        fail_tests = []
        with open(f_file) as f:
            for line in f:
                fail_tests.append(line.strip())

    # load all tests
    post_pf_file = os.path.join(resdir, "post", "RESULTS_DATA", "matrices", "PASSFAIL.csv")
    post_pf_mat = common_matrices.ExecutionMatrix(post_pf_file)
    all_tests_post = post_pf_mat.get_nonkey_colname_list()
    pre_pf_file = os.path.join(resdir, "pre", "RESULTS_DATA", "matrices", "PASSFAIL.csv")
    pre_pf_mat = common_matrices.ExecutionMatrix(pre_pf_file)
    all_tests_pre = pre_pf_mat.get_nonkey_colname_list()
    
    # Conider test thta are not flaky in both pre and post
    all_tests = list(set(all_tests_pre) & set(all_tests_post))

    # load execution outputs
    pre_orig_outlog_file = os.path.join(resdir, "pre", "RESULTS_DATA", "testexecution_outputs", "program_output.json")
    pre_muts_outlog_file = os.path.join(resdir, "pre", "RESULTS_DATA", "testexecution_outputs", "STRONG_MUTATION_output.json")
    post_orig_outlog_file = os.path.join(resdir, "post", "RESULTS_DATA", "testexecution_outputs", "program_output.json")
    post_muts_outlog_file = os.path.join(resdir, "post", "RESULTS_DATA", "testexecution_outputs", "STRONG_MUTATION_output.json")
   
    pre_orig_outlog = common_matrices.OutputLogData(pre_orig_outlog_file)
    pre_muts_outlog = common_matrices.OutputLogData(pre_muts_outlog_file)
    post_orig_outlog = common_matrices.OutputLogData(post_orig_outlog_file)
    post_muts_outlog = common_matrices.OutputLogData(post_muts_outlog_file)
 
    # Compute relevant mutants
    relevant_mutants_to_relevant_tests = get_relevant_mutants_to_relevant_tests(pre_orig_outlog, pre_muts_outlog, post_orig_outlog, post_muts_outlog)

    # load matrices and compute mutant killtest mapping
    post_sm_file = os.path.join(resdir, "post", "RESULTS_DATA", "matrices", "STRONG_MUTATION.csv")
    post_sm_mat = common_matrices.ExecutionMatrix(post_sm_file)
    mutants_to_killingtests = post_sm_mat.query_active_columns_of_rows()
    tests_to_killed_mutants = post_sm_mat.query_active_rows_of_columns()

    # return data
    return all_tests, fail_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants
#~ def load()

