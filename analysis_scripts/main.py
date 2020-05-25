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
import pandas as pd

import load
import plot

NO_CORRELATION = True
NO_SIMULATION = False

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
                
    if not NO_SIMULATION:
        print("\n# COMPUTING SIMULATIONS ...\n")

        INDEPENDENT_KILL = 'independent_kill'
        COLLATERALLY_KILL = "collaterally_kill"
        for scenario in [COLLATERALLY_KILL]: #, INDEPENDENT_KILL]:
            nRepeat = 100

            sim_cache_file = os.path.join(out_folder, "sim_cache_file.{}.json".format(scenario))
            if os.path.isfile (sim_cache_file):
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
                    relevantMuts = list(relevant_mutants_to_relevant_tests[proj])
                    randomAll_rMS[proj] = []
                    randomKillable_rMS[proj] = []
                    randomRelevant_rMS[proj] = []
                    randomAll_FR[proj] = []
                    randomKillable_FR[proj] = []
                    randomRelevant_FR[proj] = []
                    for i in tqdm.tqdm(range(nRepeat)):
                        random.shuffle(allMuts)
                        random.shuffle(killableMutants)
                        random.shuffle(relevantMuts)
                        # compute the incremental relevant score and fault detection
                        for inList, outTopList_rMS, outTopList_FR in [(allMuts, randomAll_rMS[proj], randomAll_FR[proj]), \
                                                                        (killableMutants, randomKillable_rMS[proj], randomKillable_FR[proj]), \
                                                                        (relevantMuts, randomRelevant_rMS[proj], randomRelevant_FR[proj])]:
                            tmp_rMS = []
                            tmp_FR = []
                            seen_fr_tests = set()
                            killed_relevant_muts_set = set()
                            collaterally_killed = set()
                            for mut in inList:  # TODO: Consider other scenarios
                                if scenario == COLLATERALLY_KILL and mut in collaterally_killed:
                                    continue 
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
                            
                load.common_fs.dumpJSON([randomAll_rMS, randomKillable_rMS, randomRelevant_rMS, randomAll_FR, randomKillable_FR, randomRelevant_FR], sim_cache_file)

            with_random_killable = False
            data_lists = (randomAll_rMS, randomRelevant_rMS, randomAll_FR, randomRelevant_FR)
            if with_random_killable:
                data_lists = data_lists + (randomKillable_rMS, randomKillable_FR)

            # unifirmization
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
                # normalize to 0-100
                x_label = "Percentage of Mutants"
                minstopat = 100
                normalize_data_x(data_lists, relevant_dat=randomRelevant_FR)

            # XXX: Change this if normalize_data_x changes ()
            stat_dat = stat_test (randomRelevant_FR, randomRelevant_rMS, 'Relevant', randomAll_FR, randomAll_rMS, 'Random')
            stat_file = os.path.join(out_folder, "stat_test.csv")
            load.common_fs.dumpCSV(pd.DataFrame(stat_dat), stat_file, separator=',')

            # XXX Aggregate and Plot the data
            order = ['Relevant', 'Random']
            ## FR
            img_file = os.path.join(out_folder, 'FR_PLOT_{}'.format(scenario))
            allMedToPlot = {'Random': randomAll_FR, 'Relevant': randomRelevant_FR}
            for k,v in allMedToPlot.items():
                allMedToPlot[k] = repetavg_and_proj_proportion_aggregate (v, stopAt=minstopat)
            plot.plotTrend(allMedToPlot, img_file, x_label, 'Fault Revelation', order=order)
            for pc in [0.25, 0.5, 0.75]:
                ## rMS
                pc_name = 'median'
                if pc == 0.25:
                    pc_name = '1stQuantile'
                elif pc == 0.75:
                    pc_name = '3rdQuantile'
                elif pc != 0.5:
                    pc_name = pc
                img_file = os.path.join(out_folder, 'rMS_PLOT_{}_{}'.format(scenario, pc_name))
                allMedToPlot = {'Random': randomAll_rMS, 'Relevant': randomRelevant_rMS}
                for k,v in allMedToPlot.items():
                    allMedToPlot[k] = allmedian_aggregate (v, percentile=pc, stopAt=minstopat)
                plot.plotTrend(allMedToPlot, img_file, x_label, 'Relevant Mutation Score', order=order)

            # Box
            #groupedData = {i: flattened_data
            #plot.plot_Box_Grouped(groupedData, imagefile, colors_bw, ylabel

    print("@DONE!")
#~ def main()

def normalize_data_x(data_lists, relevant_dat=None, use_med=False):
    # Default use min

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
    left_fr_list = {i: [] for i in range(1, 101)} 
    left_rms_list = {i: [] for i in range(1, 101)}
    right_fr_list = {i: [] for i in range(1, 101)}
    right_rms_list = {i: [] for i in range(1, 101)}
    for in_dat, out_list in [(left_FR, left_fr_list), (left_rMS, left_rms_list), (right_FR, right_fr_list), (right_rMS, right_rms_list)]:
        for proj, dat in in_dat.items():
            for rep_dat in dat:
                for ind, val in enumerate(rep_dat):
                    out_list[ind+1].append(val)
    res['proportions'] = range(1,101)
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
