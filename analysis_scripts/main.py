#! /usr/bin/env python

import os
import sys
import multiprocessing
import random
import glob
import shutil

import tqdm
import scipy.stats
import scipy.stats.stats as ss

import load

def error_exit(err):
    print("@Error: "+err)
    exit(1)
#~def error_exit()

def get_samples(all_tests, size_percent=5, count=10000):
    assert size_percent > 0 and size_percent < 100, "invalid size_percent"
    ts_size = len(all_tests) * size_percent * 1.0 / 100
    ts_size = max(1, int(round(ts_size)))

    samples = []
    for c in range(count):
        s = set(random.sample(all_tests, ts_size))
        if s not in samples:
            samples.append(s)
    return samples
#~ def get_samples()

def get_ms(test_suite, mutant_list, tests_to_killed_mutants):
    killed_of_tests = set()
    for t in test_suite:
        killed_of_tests |= set(tests_to_killed_mutants[t])
    killed = set(mutant_list) & killed_of_tests
    ms = len(killed) * 1.0 / len(mutant_list)
    return ms
#~ def get_ms()

def get_corelation(X, Y):
    # Compute kendall and pearson correlation and return
    assert (len(X) == len(Y)), "X and Y must have same length"
    assert len(X) > 1, "Both X and Y must have at least 2 elements"
    
    correlation = {}
    cc, p_value = ss.pearsonr(X, Y) #+[0.0] assume that FD=0 when MS=0
    correlation['pearson'] = {'corr': cc, 'p-value': p_value}
    cc, p_value = ss.kendalltau(X, Y) #+[0.0] assume that FD=0 when MS=0
    correlation['kendall'] = {'corr': cc, 'p-value': p_value}
    return correlation
#~ def get_corelation() 

def load_data(in_top_dir, tmpdir, cache_file):
    # Get the project list
    tars = glob.glob(in_top_dir+'/*.tar.gz')

    def get_pname(in_proj_tar):
        tar_name = os.path.basename(in_proj_tar)
        if tar_name.startswith('ar-') or tar_name.startswith('cr-'):
            pname = tar_name.split('.')[0]
        else:
            error_exit("TODO: implement getting pname for Wei's dataset")
    #~ def get_pname()

    projs = {get_pname(t) for t in tars}

    update_cache = True
    if os.path.isfile(cache_file):
        all_tests, fault_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants = \
                                                                                        load.common_fs.loadJSON(cache_file)
        cache_projs = set(all_tests)
        not_cached = projs - cache_projs
    else:
        not_cached = projs
        all_tests = {}
        fault_tests = {}
        relevant_mutants_to_relevant_tests = {}
        mutants_to_killingtests = {}
        tests_to_killed_mutants = {}

    update_cache = (len(not_cached) > 0)
    
    # Load projects data
    tq_data = tqdm.tqdm(tars)
    for proj_tar in tq_data:
        tq_data.set_description("Loading tar: {}".format(proj_tar))
        pname = get_pname(proj_tar)
        if pname not in not_cached:
            continue
        if os.system('cd {} && tar -xzf {} --exclude {} --exclude {} && test -d res'.format(\
                                                        tmpdir, proj_tar, \
                                                        'post/RESULTS_DATA/other_copied_results/Flakiness', \
                                                        'pre/RESULTS_DATA/other_copied_results/Flakiness')) != 0:
            error_exit("untar failed for {}".format(proj_tar))
        res_folder = os.path.join(tmpdir, 'res')
        ldata = load.load(res_folder, fault_revealing=True)
        shutil.rmtree(res_folder)

        all_tests[pname] = ldata[0]
        fault_tests[pname] = ldata[1]
        relevant_mutants_to_relevant_tests[pname] = ldata[2]
        mutants_to_killingtests[pname] = ldata[3]
        tests_to_killed_mutants[pname] = ldata[4]

    if update_cache:
        load.common_fs.dumpJson([all_tests, fault_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants], cache_file)

    return all_tests, fault_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants
#~ def load_data()

def compute(kwargs):
    #  Get subtest list
    sub_testsuites = get_samples(kwargs['all_tests'], kwargs['size_percent'], kwargs['sample_count'])

    #  Compute mutation scores and fault revelation for each test set
    ms_list = []
    rms_list = []
    find_fault = None
    for sub_ts in sub_testsuites:
        ms = get_ms(sub_ts, list(kwargs['mutants_to_killingtests']), kwargs['tests_to_killed_mutants'])
        rms = get_ms(sub_ts, list(kwargs['relevant_mutants_to_relevant_tests']), kwargs['tests_to_killed_mutants'])
        ms_list.append(ms)
        rms_list.append(rms)
        if kwargs['fault_tests']:
            if find_fault is None:
                find_fault = []
            ft = set(kwargs['fault_tests']) & set(sub_ts)
            find_fault.append(int(len(ft) > 0))

    #  compute corellations and print 
    correlations = {}
    correlations['rMS vs MS'] = get_corelation(rms_list, ms_list)
    if find_fault:
        correlations['Fault vs rMS'] = get_corelation(find_fault, rms_list)
        correlations['Fault vs MS'] = get_corelation(find_fault, ms_list)
    
    #print(correlations)
    return correlations
#~ def compute()

def main():
    if len(sys.argv) != 3:
        error_exit("expected 3 args. got {}". format(len(sys.argv)))
    in_top_dir = os.path.abspath(sys.argv[1])
    out_top_dir = os.path.abspath(sys.argv[2])
    if not os.path.isdir(in_top_dir):
        error_exit("in top dir missing ({})".format(in_top_dir))
    if not os.path.isdir(out_top_dir):
        error_exit("out top dir missing ({})".format(out_top_dir))

    out_folder = os.path.join(out_top_dir, "ANALYSIS_OUTPUT")
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    # CHANGABLE Hardcodec parameters
    test_sizes_percents = list(range(5,91,5))
    sample_count = 10000
    parallel_count = 50

    # load data
    tmpdir = os.path.join(out_folder, "tmp_extracttar_dir.tmp")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)
    print("# LOADING DATA ...")
    cache_file = os.path.join(out_folder, "cache_file.json")
    all_tests, fault_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants = \
								              load_data(in_top_dir, tmpdir, cache_file)
    shutil.rmtree(tmpdir)
    
    # update parallel_count
    parallel_count = min(parallel_count, len(all_tests))
    map_func = map
    if parallel_count > 1:
        parallel_pool = multiprocessing.Pool(parallel_count)
        map_func = parallel_pool.map
    
    # organize the data
    order = sorted(list(all_tests.keys()))
    data = [{
                'sample_count': sample_count,
                'all_tests': all_tests[i],
                'fault_tests': fault_tests[i],
                'relevant_mutants_to_relevant_tests': relevant_mutants_to_relevant_tests[i], 
                'mutants_to_killingtests': mutants_to_killingtests, 
                'tests_to_killed_mutants': tests_to_killed_mutants, 
                'size_percent': None,
            } for i in order
    ]

    # for each test size
    #    compute for each project (in parallel) 
    results = {}
    for size_percent in tqdm.tqdm(test_sizes_percents):
        #print("# EXECUTING FOR TEST SIZE {}% ...".format(size_percent))
        # set size_percent
        for d in data:
            d['size_percent'] = size_percent
        # Compute
        tmp_results = map_func(compute, data)

        # Organize results
        results[size_percent] = {}
        for i in range(len(order)):
            results[size_percent][order[i]] = tmp_results[i]

    # dump results
    raw_res_file = os.path.join(out_folder, "raw_res_file.json")
    load.common_fs.dumpJson(results, raw_res_file, pretty=True)

    print("@DONE!")
#~ def main()

if __name__ == "__main__":
    main()
