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
import numpy as np

import load
import plot

def error_exit(err):
    print("@Error: "+err)
    exit(1)
#~def error_exit()

def allmedian_aggregate (proj2repetlists, percentile=0.5, stopAt=None):
    ''' return a key value, where keys are indexes and values the values
    '''
    projlist = list(proj2repetlists)
    size = len(proj2repetlists[projlist[0]][0])
    if stopAt is not None and size > stopAt:
        size = stopAt
    res = {}
    for i in range(size):
        key = i+1
        for proj, pdat in proj2repetlists.items():
            for rep_dat in pdat:
                res[key].append(rep_dat[i])
        
        res[key] = np.quantile(res[key], percentile)

    return res
#~ def allmedian_aggregate()

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
        return pname
    #~ def get_pname()

    projs = set([get_pname(t) for t in tars])

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
        load.common_fs.dumpJSON([all_tests, fault_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants], cache_file)

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
                'mutants_to_killingtests': mutants_to_killingtests[i], 
                'tests_to_killed_mutants': tests_to_killed_mutants[i], 
                'size_percent': None,
            } for i in order
    ]

    # for each test size
    #    compute for each project (in parallel) 
    raw_results = {}
    print("\n# COMPUTING CORRELATIONS ...\n")
    for size_percent in tqdm.tqdm(test_sizes_percents):
        #print("# EXECUTING FOR TEST SIZE {}% ...".format(size_percent))
        # set size_percent
        for d in data:
            d['size_percent'] = size_percent
        # Compute
        tmp_results = map_func(compute, data)
        print(tmp_results)        

        # Organize raw_results
        raw_results[size_percent] = {}
        for i in range(len(order)):
            raw_results[size_percent][order[i]] = tmp_results[i]

    # dump results
    raw_res_file = os.path.join(out_folder, "raw_res_file.json")
    load.common_fs.dumpJSON(raw_results, raw_res_file, pretty=True)

    results = {}
    results_median = {}
    for size_percent, size_data in raw_results.items():
        results[size_percent] = {}
        for proj, proj_data in sorted(list(size_data.items())):
            for cor_subj, cor_subj_data in proj_data.items():
                if cor_subj not in results[size_percent]:
                    results[size_percent][cor_subj] = {}
                for cor_type, corr in cor_subj_data.items():
                    if cor_type not in results[size_percent][cor_subj]:
                        results[size_percent][cor_subj][cor_type] = []
                    results[size_percent][cor_subj][cor_type].append(corr)
    res_file = os.path.join(out_folder, "res_file.json")
    load.common_fs.dumpJSON(results, res_file, pretty=True)

    exclude_nan = True
    if exclude_nan:
        print("# INFO: Excluding nan correlation for med")
    for size_percent in results:
        results_median[size_percent] = {}
        for cor_subj in results[size_percent]:
            results_median[size_percent][cor_subj] = {}
            for cor_type in results[size_percent][cor_subj]:
                if exclude_nan:
                    c_vals = [x['corr'] for x in results[size_percent][cor_subj][cor_type] if not np.isnan(x['corr'])]
                else:
                    c_vals = [x['corr'] for x in results[size_percent][cor_subj][cor_type]]
                med_vals = { 
                               'min':np.quantile(c_vals, 0),
                               '1st-Qt':np.quantile(c_vals, 0.25),
                               'med':np.quantile(c_vals, 0.5),
                               '3rd-Qt':np.quantile(c_vals, 0.75),
                               'max':np.quantile(c_vals, 1), 
                               'avg':sum(c_vals) * 1.0 / len(c_vals), 
                           }
                results_median[size_percent][cor_subj][cor_type] = med_vals
    meds_res_file = os.path.join(out_folder, "res_file_meds.json")
    load.common_fs.dumpJSON(results_median, meds_res_file, pretty=True)
                
    print("\n# COMPUTING SIMULATIONS ...\n")

    scenario = 'independent_kill'
    nRepeat = 10000

    sim_cache_file = os.path.join(out_folder, "sim_cache_file.{}.json".format(scenario))
    if os.path.isfile (sim_cache_file)
        randomAll_rMS, randomKillable_rMS, randomRelevant_rMS, randomAll_FR, randomKillable_FR, randomRelevant_FR = load.common_fs.loadJSON(sim_cache_file)
    else:
        randomAll_rMS = {}
        randomKillable_rMS = {}
        randomRelevant_rMS = {}
        randomAll_FR = {}
        randomKillable_FR = {}
        randomRelevant_FR = {}
        for ind, proj in enumerate(order):
            print ("processing project {}/{} ...".format(ind+1, len(order)))

            fr_tests = set(fault_tests[proj])
            tests_killing_relevant_muts = set()
            for m, tl in relevant_mutants_to_relevant_tests[proj].items():
                tests_killing_relevant_muts |= set(tl)

            allMuts = list(mutants_to_killingtests[proj])
            killableMutants = [m for m, kt in mutants_to_killingtests[proj].items() if len(kt) > 0]
            relevantMuts = list(relevant_mutants_to_relevant_tests[proj])
            randomAll_rMS[proj] = []
            randomKillable_rMS[proj] = []
            randomRelevant_rMS[proj] = []
            randomAll_FR[proj] = []
            randomKillable_FR[proj] = []
            randomRelevant_FR[proj] = []
            for i in tqdm.tqdm(range(nRepeat)):
                allMuts.shuffle()
                killableMutants.shuffle()
                relevantMuts.shuffle()
                # compute the incremental relevant score and fault detection
                for inList, outTopList_rMS, outTopList_FR in [(allMuts, randomAll_rMS[proj], randomAll_FR[proj]), \
                                                                (killableMutants, randomKillable_rMS[proj], randomKillable_FR[proj]), \
                                                                (relevantMuts, randomRelevant_rMS[proj], randomRelevant_FR[proj])]:
                    tmp_rMS = []
                    tmp_FR = []
                    seen_fr_tests = set()
                    killed_relevant_muts_set = set()
                    for mut in inList:  # TODO: Consider other scenarios
                        kts = mutants_to_killingtests[mut]
                        if len(kts) != 0:
                            # pick a killing test
                            t = random.choice(kts)

                            # set FR and rMS
                            if t in tests_killing_relevant_muts:
                                killed_relevant_muts_set |= set(tests_to_killed_mutants[t])
                            if t in fr_tests:
                                seen_fr_tests.add(t)
                        tmp_FR.append(len(seen_fr_tests))
                        tmp_rMS.append(len(killed_relevant_muts_set) * 1.0 / len(relevantMuts))

                    outTopList_rMS.append(tmp_rMS)
                    outTopList_FR.append(tmp_FR)
        load.common_fs.dumpJSON([randomAll_rMS, randomKillable_rMS, randomRelevant_rMS, randomAll_FR, randomKillable_FR, randomRelevant_FR], sim_cache_file)

    # unifirmization
    minstopat = min([len(rm2tests) for proj, rm2tests in relevant_mutants_to_relevant_tests.items()])
    maxstopat = max([len(rm2tests) for proj, rm2tests in relevant_mutants_to_relevant_tests.items()])

    # XXX Aggregate and Plot the data
    for pc in [0.25, 0.5, 0.75]:
        ## FR
        img_file = os.path.join(out_folder, 'FR_PLOT_{}_{}'.format(scenario, pc))
        allMedToPlot = {'Random': randomAll_FR, 'RandomKillable': randomKillable_FR, 'RandomRelevant': randomRelevant_FR}
        for k,v in allMedToPlot:
            allMedToPlot[k] = allmedian_aggregate (v, percentile=pc, stopAt=minstopat)
        plot.plotTrend(allMedToPlot, img_file, 'Number of Mutants', 'Fault Revelation')
        ## rMS
        img_file = os.path.join(out_folder, 'rMS_PLOT_{}_{}'.format(scenario, pc))
        allMedToPlot = {'Random': randomAll_rMS, 'RandomKillable': randomKillable_rMS, 'RandomRelevant': randomRelevant_rMS}
        for k,v in allMedToPlot:
            allMedToPlot[k] = allmedian_aggregate (v, percentile=pc, stopAt=minstopat)
        plot.plotTrend(allMedToPlot, img_file, 'Number of Mutants', 'Relevant Mutation Score')

    print("@DONE!")
#~ def main()

if __name__ == "__main__":
    main()
