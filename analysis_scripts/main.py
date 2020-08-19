#! /usr/bin/env python

import os
import sys
import multiprocessing
import random
import glob
import shutil
import copy

import tqdm
import scipy.stats
import scipy.stats.stats as ss
import numpy as np
import pandas as pd

import load
import plot

NO_CORRELATION = True
NO_SIMULATION = False
NO_RELEVANT_IN_PLOT = True
STOP_AT_N_MUTANTS = 50 #100 # None

def error_exit(err):
    print("@Error: "+err)
    exit(1)
#~def error_exit()

###### Non Parametic Vargha Delaney A12 ######
# Taken from -- https://gist.github.com/timm/5630491

def a12(lst1,lst2,pairwise=False, rev=True):
    "how often is x in lst1 more than y in lst2?"
    more = same = 0.0
    for i,x in enumerate(lst1):
        second = [lst2[i]] if pairwise else lst2
        for y in second:
            if   x==y : same += 1
            elif rev     and x > y : more += 1
            elif not rev and x < y : more += 1
    return (more + 0.5*same)  / (len(lst1) if pairwise else len(lst1)*len(lst2))

def wilcoxon(list1, list2, isranksum=True):
    if isranksum:
        p_value = scipy.stats.ranksums(list1, list2)
    else:
        p_value = scipy.stats.wilcoxon(list1, list2)
    return p_value
#~ def wilcoxon()

#############################################

def repetavg_and_proj_proportion_aggregate (proj2repetlists, stopAt=None):
    projlist = list(proj2repetlists)
    size = len(proj2repetlists[projlist[0]][0])
    if stopAt is not None and size > stopAt:
        size = stopAt
    res = {}
    for i in range(size):
        key = i+1
        res[key] = 0
        for proj, pdat in proj2repetlists.items():
            plist = []
            for rep_dat in pdat:
                if key not in res:
                    res[key] = []
                plist.append(rep_dat[i])
        
            res[key] += sum(plist) * 1.0 / len(plist)
        res[key] = res[key] * 1.0 / len(proj2repetlists)

    return res
#~ def repetavg_and_proj_sum_aggregate()

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
                if key not in res:
                    res[key] = []
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

        # XXX fix matrices
        #if True:
        #    fix_matrices(os.path.join(res_folder, 'pre', 'RESULTS_DATA'))
        #    fix_matrices(os.path.join(res_folder, 'post', 'RESULTS_DATA'))
        #    if os.system('tar -czf {} {}'.format(proj_tar, res_folder)) != 0:
        #        error_exit ("failed to update tar by fixing matrices")
        #~~~~~~~~~~~~

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

def fix_matrices (results_data_dir):
    import os
    from muteria.drivers import DriversUtils
    sm_mat = os.path.join(results_data_dir, 'matrices', 'STRONG_MUTATION.csv')
    sm_out = os.path.join(results_data_dir, 'testexecution_outputs', 'STRONG_MUTATION_output.json')
    p_mat = os.path.join(results_data_dir, 'matrices', 'PASSFAIL.csv')
    p_out = os.path.join(results_data_dir, 'testexecution_outputs', 'program_output.json')
    DriversUtils.update_matrix_to_cover_when_difference (sm_mat, sm_out, p_mat, p_out)
    print ("Fix success for {}!".format(results_data_dir))
#~ def fix_matrices ()

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
    relmut_pred_file = None
    if len(sys.argv) != 3:
        if len(sys.argv) == 4:
            relmut_pred_file = os.path.abspath(sys.argv[3])
            if not os.path.isfile(relmut_pred_file):
                error_exit("relmut_pred_file not existing")
        else:
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

    # Load predicted if set
    proj_to_pred_mut_to_relscore = None
    if relmut_pred_file is not None:
        pred_raw = load.common_fs.loadJSON(relmut_pred_file)
        proj_to_pred_mut_to_relscore = {}
        for proj, dat in pred_raw.items():
            proj = proj.split("_")[1]
            proj_to_pred_mut_to_relscore[proj] = {}
            if type(dat) == dict:
                dat = dat['score']
            for mut, oracle, pred in dat:
                mut = ":".join(['mart_0', mut])
                if mut not in mutants_to_killingtests[proj]:
                    continue
                proj_to_pred_mut_to_relscore[proj][mut] = pred
                if oracle == 1:
                    if mut not in relevant_mutants_to_relevant_tests[proj]:
                        #print("mutant remote relevant but not local ({}). proj is {}".format(mut, proj))
                        assert False, "mutant remote relevant but not local ({}). proj is {}".format(mut, proj)
                else:
                    if mut in relevant_mutants_to_relevant_tests[proj]:
                        assert False, "reltests is {}. mutant remote non-relevant but not local ({}). proj is {}".format(relevant_mutants_to_relevant_tests[proj][mut], mut, proj)
            
    for pp_1, pp_2 in [('cr-12', 'cr-17'), ('cr-5', 'cr-16')]:
        if pp_1 in proj_to_pred_mut_to_relscore and pp_2 not in proj_to_pred_mut_to_relscore:
            proj_to_pred_mut_to_relscore[pp_2] = copy.deepcopy(proj_to_pred_mut_to_relscore[pp_1])
        if pp_2 in proj_to_pred_mut_to_relscore and pp_1 not in proj_to_pred_mut_to_relscore:
            proj_to_pred_mut_to_relscore[pp_1] = copy.deepcopy(proj_to_pred_mut_to_relscore[pp_2])
    
    # update parallel_count
    parallel_count = min(parallel_count, len(all_tests))
    map_func = map
    if parallel_count > 1:
        parallel_pool = multiprocessing.Pool(parallel_count)
        map_func = parallel_pool.map
    
    if parallel_count > 1:
        parallel_pool = multiprocessing.Pool(parallel_count)
        map_func = parallel_pool.map
    
    # organize the data
    order = sorted(list(all_tests.keys()))
    if proj_to_pred_mut_to_relscore is not None:
        order = list(proj_to_pred_mut_to_relscore)
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

    if not NO_CORRELATION:
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
                

    # Mutants on commit
    mutant_on_commit_file = os.path.join(in_top_dir, 'Corebench_mutantsinpatch.json')
    proj2mutoncommit = None
    if os.path.isfile(mutant_on_commit_file):
        proj2mutoncommit = load.common_fs.loadJSON(mutant_on_commit_file)

    if not NO_SIMULATION:
        print("\n# COMPUTING SIMULATIONS ...\n")

        INDEPENDENT_KILL = 'independent_kill'
        COLLATERALLY_KILL = "collaterally_kill"
        for scenario in [COLLATERALLY_KILL]: #, INDEPENDENT_KILL]:
            nRepeat = 100
            if proj_to_pred_mut_to_relscore is not None:
                nRepeat = 400

            sim_cache_file = os.path.join(out_folder, "sim_cache_file.{}.json".format(scenario))
            if os.path.isfile (sim_cache_file):
                if proj_to_pred_mut_to_relscore is None:
                    randomAll_rMS, randomKillable_rMS, randomOnCommit_rMS, randomRelevant_rMS, randomAll_FR, randomKillable_FR, randomOnCommit_FR, randomRelevant_FR = load.common_fs.loadJSON(sim_cache_file)
                else:
                    randomAll_rMS, randomKillable_rMS, randomOnCommit_rMS, randomRelevant_rMS, predictedRelevant_rMS, randomAll_FR, randomKillable_FR, randomOnCommit_FR, randomRelevant_FR, predictedRelevant_FR = load.common_fs.loadJSON(sim_cache_file)
            else:
                randomAll_rMS = {}
                randomKillable_rMS = {}
                randomOnCommit_rMS = {}
                randomRelevant_rMS = {}
                predictedRelevant_rMS = {}
                randomAll_FR = {}
                randomKillable_FR = {}
                randomOnCommit_FR = {}
                randomRelevant_FR = {}
                predictedRelevant_FR = {}
                for ind, proj in enumerate(order):
                    print ("processing project {}/{} ...".format(ind+1, len(order)))

                    proj_tests_to_killed_relevant_muts = {}
                    for m in relevant_mutants_to_relevant_tests[proj]:
                        kill_tests = mutants_to_killingtests[proj][m]
                        for t in kill_tests:
                            if t not in proj_tests_to_killed_relevant_muts:
                                proj_tests_to_killed_relevant_muts[t] = set()
                            proj_tests_to_killed_relevant_muts[t].add(m)

                    fr_tests = set(fault_tests[proj])
                    tests_killing_relevant_muts = set()
                    for m, tl in relevant_mutants_to_relevant_tests[proj].items():
                        tests_killing_relevant_muts |= set(tl)

                    allMuts = list(mutants_to_killingtests[proj])
                    killableMutants = [m for m, kt in mutants_to_killingtests[proj].items() if len(kt) > 0]
                    onCommitMutants = list(proj2mutoncommit[proj]) if proj2mutoncommit is not None else []
                    relevantMuts = list(relevant_mutants_to_relevant_tests[proj])
                    predictedMuts = list(proj_to_pred_mut_to_relscore[proj])

                    randomAll_rMS[proj] = []
                    randomKillable_rMS[proj] = []
                    randomOnCommit_rMS[proj] = []
                    randomRelevant_rMS[proj] = []
                    predictedRelevant_rMS[proj] = []
                    randomAll_FR[proj] = []
                    randomKillable_FR[proj] = []
                    randomOnCommit_FR[proj] = []
                    randomRelevant_FR[proj] = []
                    predictedRelevant_FR[proj] = []
                    for i in tqdm.tqdm(range(nRepeat)):
                        random.shuffle(allMuts)
                        random.shuffle(killableMutants)
                        random.shuffle(relevantMuts)
                        workList = [(allMuts, randomAll_rMS[proj], randomAll_FR[proj]), \
                                                                        (killableMutants, randomKillable_rMS[proj], randomKillable_FR[proj]), \
                                                                        (onCommitMutants, randomOnCommit_rMS[proj], randomOnCommit_FR[proj]), \
                                                                        (relevantMuts, randomRelevant_rMS[proj], randomRelevant_FR[proj])]
                        if proj_to_pred_mut_to_relscore is not None:
                            predictedMuts.sort(reverse=True, key=lambda m: (proj_to_pred_mut_to_relscore[proj][m], random.random()))
                            workList.append((predictedMuts, predictedRelevant_rMS[proj], predictedRelevant_FR[proj]))
                        # compute the incremental relevant score and fault detection
                        for inList, outTopList_rMS, outTopList_FR in workList:
                            tmp_rMS = []
                            tmp_FR = []
                            seen_fr_tests = set()
                            killed_relevant_muts_set = set()
                            collaterally_killed = set()
                            for mut in inList:  # TODO: Consider other scenarios
                                if scenario == COLLATERALLY_KILL and mut in collaterally_killed:
                                    continue 
                                #print(list(mutants_to_killingtests[proj]))
                                kts = mutants_to_killingtests[proj][mut]
                                if len(kts) != 0:
                                    # pick a killing test
                                    t = random.choice(kts)
                                    if scenario == COLLATERALLY_KILL:
                                        collaterally_killed |= set(tests_to_killed_mutants[proj][t])

                                    # set FR and rMS
                                    if t in proj_tests_to_killed_relevant_muts:
                                        killed_relevant_muts_set |= set(proj_tests_to_killed_relevant_muts[t])
                                    if t in fr_tests:
                                        seen_fr_tests.add(t)
                                tmp_FR.append(int(len(seen_fr_tests) > 0))
                                tmp_rMS.append(len(killed_relevant_muts_set) * 1.0 / len(relevantMuts))

                            outTopList_rMS.append(tmp_rMS)
                            outTopList_FR.append(tmp_FR)
                            
                if proj_to_pred_mut_to_relscore is None:
                    load.common_fs.dumpJSON([randomAll_rMS, randomKillable_rMS, randomOnCommit_rMS, randomRelevant_rMS, randomAll_FR, randomKillable_FR, randomOnCommit_FR, randomRelevant_FR], sim_cache_file)
                else:
                    load.common_fs.dumpJSON([randomAll_rMS, randomKillable_rMS, randomOnCommit_rMS, randomRelevant_rMS, predictedRelevant_rMS, randomAll_FR, randomKillable_FR, randomOnCommit_FR, randomRelevant_FR, predictedRelevant_FR], sim_cache_file)

            with_random_killable = False
            data_lists = (randomAll_rMS, randomAll_FR)
            if not NO_RELEVANT_IN_PLOT:
                data_lists = data_lists + (randomRelevant_rMS, randomRelevant_FR)
            if with_random_killable:
                data_lists = data_lists + (randomKillable_rMS, randomKillable_FR)
            if proj2mutoncommit is not None:
                data_lists = data_lists + (randomOnCommit_rMS, randomOnCommit_FR)
            if proj_to_pred_mut_to_relscore is not None:
                data_lists = data_lists + (predictedRelevant_rMS, predictedRelevant_FR)

            # uniformization
            minstopat = 999999999999
            maxstopat = 0
            for l in data_lists:
                for p,d in l.items():
                    for e_list in d:
                        minstopat = min(minstopat, len(e_list))
                        maxstopat = max(maxstopat, len(e_list))
            print("# minstopat is {}, maxstopat is {}".format(minstopat, maxstopat))

            x_label = 'Number of Mutants'
            if scenario == COLLATERALLY_KILL:
                if STOP_AT_N_MUTANTS is not None:
                    minstopat = STOP_AT_N_MUTANTS
                    assert STOP_AT_N_MUTANTS > 0
                    # Stop at STOP_AT_N_MUTANTS
                    print ("\n# Stop at", STOP_AT_N_MUTANTS)
                    for dl in data_lists:
                        for proj, p_data in dl.items():
                            for ind, rep_dat in enumerate(p_data):
                                to_add = [rep_dat[-1]] * (STOP_AT_N_MUTANTS - len(rep_dat))
                                if to_add:
                                    rep_dat.extend(to_add)
                                else:
                                    del rep_dat[STOP_AT_N_MUTANTS:]
                else:
                    # normalize to 0-100
                    x_label = "Percentage of Mutants"
                    minstopat = 100
                    size_per_proj = normalize_data_x(data_lists, relevant_dat=(None if NO_RELEVANT_IN_PLOT else randomRelevant_FR))
                    avg_proportion = [size_per_proj[p] * 1.0 / len(mutants_to_killingtests[p]) for p in size_per_proj]
                    avg_proportion = sum(avg_proportion) / len(avg_proportion)
                    avg_num = [size_per_proj[p] for p in size_per_proj]
                    avg_num = sum(avg_num) * 1.0 / len(avg_num)
                    print ("\n# ALL MUTANTS >  AVG plot Proportion: {} AVG plot Number {}".format(avg_proportion, avg_num))
                    avg_proportion = [size_per_proj[p] * 1.0 / len([m for m, t in mutants_to_killingtests[p].items() if len(t) > 0]) for p in size_per_proj]
                    avg_proportion = sum(avg_proportion) / len(avg_proportion)
                    avg_num = [size_per_proj[p] for p in size_per_proj]
                    avg_num = sum(avg_num) * 1.0 / len(avg_num)
                    print ("\n# KILLABLE MUTANTS > AVG plot Proportion: {} AVG plot Number {}".format(avg_proportion, avg_num))

            # XXX: Change this if normalize_data_x changes ()
            if not NO_RELEVANT_IN_PLOT:
                stat_dat = stat_test (randomRelevant_FR, randomRelevant_rMS, 'Relevant', randomAll_FR, randomAll_rMS, 'Random')
                stat_file = os.path.join(out_folder, "RelevantVSRandom-stat_test.csv")
                load.common_fs.dumpCSV(pd.DataFrame(stat_dat), stat_file, separator=',')

            if proj2mutoncommit is not None:
                if not NO_RELEVANT_IN_PLOT:
                    stat_dat = stat_test (randomRelevant_FR, randomRelevant_rMS, 'Relevant', randomOnCommit_FR, randomOnCommit_rMS, 'Modification')
                    stat_file = os.path.join(out_folder, "RelevantVSModification-stat_test.csv")
                    load.common_fs.dumpCSV(pd.DataFrame(stat_dat), stat_file, separator=',')

                stat_dat = stat_test (randomOnCommit_FR, randomOnCommit_rMS, 'Modification', randomAll_FR, randomAll_rMS, 'Random')
                stat_file = os.path.join(out_folder, "ModificationVSRandom-stat_test.csv")
                load.common_fs.dumpCSV(pd.DataFrame(stat_dat), stat_file, separator=',')

            if proj_to_pred_mut_to_relscore is not None:
                if not NO_RELEVANT_IN_PLOT:
                    stat_dat = stat_test (randomRelevant_FR, randomRelevant_rMS, 'Relevant', predictedRelevant_FR, predictedRelevant_rMS, 'Prediction')
                    stat_file = os.path.join(out_folder, "RelevantVSPrediction-stat_test.csv")
                    load.common_fs.dumpCSV(pd.DataFrame(stat_dat), stat_file, separator=',')

                stat_dat = stat_test (predictedRelevant_FR, predictedRelevant_rMS, 'Prediction', randomAll_FR, randomAll_rMS, 'Random')
                stat_file = os.path.join(out_folder, "PredictionVSRandom-stat_test.csv")
                load.common_fs.dumpCSV(pd.DataFrame(stat_dat), stat_file, separator=',')

            # XXX Aggregate and Plot the data
            plot_order = []
            if not NO_RELEVANT_IN_PLOT:
                plot_order.append('Relevant')
            if proj_to_pred_mut_to_relscore is not None:
                plot_order.append("Prediction")
            plot_order.append('Random')
            if proj2mutoncommit is not None:
                plot_order.append("Modification")

            ## FR
            img_file = os.path.join(out_folder, 'FR_PLOT_{}'.format(scenario))
            allMedToPlot = {'Random': randomAll_FR}
            if not NO_RELEVANT_IN_PLOT:
                allMedToPlot['Relevant'] = randomRelevant_FR
            if proj2mutoncommit is not None:
                allMedToPlot['Modification'] = randomOnCommit_FR
            if proj_to_pred_mut_to_relscore is not None:
                allMedToPlot['Prediction'] = predictedRelevant_FR
                
            # Save raw data and plot boxes
            raw_FR_plot_data_json_file = os.path.join(out_folder, "raw_FR_plot_data_{}.json".format(scenario))
            load.common_fs.dumpJSON(allMedToPlot, raw_FR_plot_data_json_file)
            apfd_FR_plot_data_img_file = os.path.join(out_folder, "apfd_FR_plot_data_{}".format(scenario))
            apfd_FR_median_file = os.path.join(out_folder, "apfd_FR_medians_{}.json".format(scenario))
            apfd_data = {}
            IS_PERCENTAGE_FAULTS = True
            for tech, t_data in allMedToPlot.items():
                apfd_data[tech] = []
                tmp_acc = {}
                for proj, p_data in t_data.items():
                    for rep_ind, rep_data in enumerate(p_data):
			            assert len(rep_data) <= minstopat
                        if IS_PERCENTAGE_FAULTS:
                            if rep_ind not in tmp_acc:
                                tmp_acc[rep_ind] = list(rep_data)
                            else:
                                # TODO: pairwise sum the lists...
                        else:
                            _apfd_val = np.trapz(rep_data) * 100.0 / (len(rep_data) - 1)
                            apfd_data[tech].append(_apfd_val)
            medians = plot.plotBoxes(apfd_data, plot_order, apfd_FR_plot_data_img_file, plot.colors_bw, ylabel="APFD", yticks_range=range(0,101,20), fontsize=26, title=None)
            load.common_fs.dumpJSON(medians, apfd_FR_median_file)
            
            for k,v in allMedToPlot.items():
                allMedToPlot[k] = repetavg_and_proj_proportion_aggregate (v, stopAt=minstopat)
            plot.plotTrend(allMedToPlot, img_file, x_label, 'Fault Revelation', order=plot_order)
            for pc, pc_name in {0:'min', 0.25: '1stQuantile', 0.5: 'median', 0.75: '3rdQuantile', 1: 'max'}.items():
                ## rMS
                img_file = os.path.join(out_folder, 'rMS_PLOT_{}_{}'.format(scenario, pc_name))
                allMedToPlot = {'Random': randomAll_rMS}
                if not NO_RELEVANT_IN_PLOT:
                    allMedToPlot['Relevant'] = randomRelevant_rMS
                if proj2mutoncommit is not None:
                    allMedToPlot['Modification'] = randomOnCommit_rMS
                if proj_to_pred_mut_to_relscore is not None:
                    allMedToPlot['Prediction'] = predictedRelevant_rMS
                for k,v in allMedToPlot.items():
                    allMedToPlot[k] = allmedian_aggregate (v, percentile=pc, stopAt=minstopat)
                plot.plotTrend(allMedToPlot, img_file, x_label, 'Relevant Mutation Score', order=plot_order)

            # Box
            #groupedData = {i: flattened_data
            #plot.plot_Box_Grouped(groupedData, imagefile, colors_bw, ylabel

    print("@DONE!")
#~ def main()

def normalize_data_x(data_lists, relevant_dat=None, use_med=False):
    # Default use min

    size_per_proj = {}
    for proj in data_lists[0]:
        if relevant_dat is None:
            min_max_len = len(data_lists[0][proj][0])
            for dl in data_lists:
                for rep_dat in dl[proj]:
                    min_max_len = min(min_max_len, len(rep_dat))
        else:
            min_max_len = len(relevant_dat[proj][0])
            for rep_dat in relevant_dat[proj]:
                min_max_len = max(min_max_len, len(rep_dat))
            
        for dl in data_lists:
            for ind, rep_dat in enumerate(dl[proj]):
                to_add = [rep_dat[-1]] * (min_max_len - len(rep_dat))
                if to_add:
                    rep_dat.extend(to_add)
                else:
                    del rep_dat[min_max_len:]
                dl[proj][ind] = normalized_x(rep_dat)
        size_per_proj[proj] = min_max_len
    return size_per_proj
#~ def normalize_data_x()

def normalized_x(arr):
    l = len(arr)
    step = l/100.0

    if step >= 1:
        ret = []
        next_p = 1
        for ind in range(len(arr)):
            if (ind+1) / step >= next_p:
                ret.append(arr[ind])
                next_p += 1
        if len(ret) == 99:
            ret.append(arr[-1])
    else:
        ret = []
        next_p = 1
        ind = 0
        for v in arr:
            while ind < next_p:
                ret.append(v)
                ind += step
            next_p += 1
    if len(ret) < 100:
        print (len(arr), len(ret), 'PB')
        #error_exit('PB')

    return ret[:100]
#~ def normalize()

def stat_test (left_FR, left_rMS, left_name, right_FR, right_rMS, right_name):
    res = {}
    left_fr_list = {} #i: [] for i in range(1, 101)} 
    left_rms_list = {} #i: [] for i in range(1, 101)}
    right_fr_list = {} #i: [] for i in range(1, 101)}
    right_rms_list = {} #i: [] for i in range(1, 101)}
    for in_dat, out_list in [(left_FR, left_fr_list), (left_rMS, left_rms_list), (right_FR, right_fr_list), (right_rMS, right_rms_list)]:
        for proj, dat in in_dat.items():
            for rep_dat in dat:
                for ind, val in enumerate(rep_dat):
                    if ind+1 not in out_list:
                        out_list[ind+1] = []
                    out_list[ind+1].append(val)
    res['proportions'] = range(1,1+len(left_fr_list))
    res['FR-p_value'] = []
    res['rMS-p_value'] = []
    res['FR-A12'] = []
    res['rMS-A12'] = []
    for i in left_fr_list:
        res['FR-p_value'].append(wilcoxon(left_fr_list[i], right_fr_list[i], isranksum=True).pvalue)
        res['rMS-p_value'].append(wilcoxon(left_rms_list[i], right_rms_list[i], isranksum=True).pvalue)
        res['FR-A12'].append(a12(left_fr_list[i],right_fr_list[i], pairwise=False, rev=True))
        res['rMS-A12'].append(a12(left_rms_list[i],right_rms_list[i], pairwise=False, rev=True))

    return res
#~ def stat_test()

if __name__ == "__main__":
    main()
