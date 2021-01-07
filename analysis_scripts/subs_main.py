#! /usr/bin/env python

import os
import sys
from multiprocessing import Pool
import random
import glob
import shutil
import copy
import time
from datetime import timedelta
import bisect

import tqdm
import scipy.stats
import scipy.stats.stats as ss
import numpy as np
import pandas as pd

import seaborn as sns

import subs_load
import plot
 
#STOP_AT_N_MUTANTS = 100 # None

NUM_REPETITIONS = 1000 #1000

SUB_REPET_NUM = 1

RANDOM = "RANDOM"
PRED_MACHINE_TRANSLATION = "MACHINE-TRANSLATION"
PRED_DECISION_TREES = "DECISION-TREES"

Use_proportion_analysed_mutants = True

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

def inner_stattest(in_data, filename, order=None):
    """ 
        in_data = {'tech1': value_list1, 'tech2': value_list2, ...}
        or
        in_data = {'tech1': {'proj1': value_list1, ...}, 'tech2': {'proj2': value_list2, ...}, ...}
    """
    if order is None:
        order = list(in_data)
        
    statstest_obj = {}
    for pos1,g1 in enumerate(order):
        for pos2, g2 in enumerate(order):
            if pos1 >= pos2:
                continue
            if type(in_data[g1]) == dict:
                tmp_stats = {v:{} for v in list(in_data[g1])}
                for k in tmp_stats:
                    tmp_stats[k]['p_value'] = wilcoxon(in_data[g1][k], in_data[g2][k], isranksum=False)
                    tmp_stats[k]['A12'] = a12(in_data[g1][k], in_data[g2][k], pairwise=True)
            else:
                assert type(in_data[g1]) == list or type(in_data[g1]) == tuple, "invalid data type"
                tmp_stats = {}
                tmp_stats['p_value'] = wilcoxon(in_data[g1], in_data[g2], isranksum=False)
                tmp_stats['A12'] = a12(in_data[g1], in_data[g2], pairwise=True)
            statstest_obj[str((g1, g2))] = tmp_stats
            
    subs_load.common_fs.dumpJSON(statstest_obj, filename, pretty=True)
    
    return statstest_obj
#~ def inner_stattest()

def repetavg_and_proj_proportion_aggregate (proj2repetlists, stopAt=None, agg_median=False):
    projlist = list(proj2repetlists)
    size = len(proj2repetlists[projlist[0]][0])
    if stopAt is not None and size > stopAt:
        size = stopAt
    res = {}
    variance_res = None
    if agg_median:
        variance_res  = {}
    for i in range(size):
        key = i+1
        if agg_median:
            res[key] = []
        else:
            res[key] = 0
        for proj, pdat in proj2repetlists.items():
            plist = []
            for rep_dat in pdat:
                #if key not in res:
                #    res[key] = []
                plist.append(rep_dat[i])
        
            if agg_median:
                res[key].append(sum(plist) * 1.0 / len(plist))
            else:
                res[key] += sum(plist) * 1.0 / len(plist)
        if agg_median:
            variance_res[key] = (min(res[key]), max(res[key]))
            res[key] = np.median(res[key])
        else:
            res[key] = res[key] * 1.0 / len(proj2repetlists)

    return res, variance_res
#~ def repetavg_and_proj_proportion_aggregate()

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

def get_subs_ms(test_suite, tests_to_killed_clusters):
    killed_of_tests = set()
    for t in test_suite:
        killed_of_tests |= set(tests_to_killed_clusters[t])
    clusters = set()
    for t in tests_to_killed_clusters:
        clusters |= set(tests_to_killed_clusters[t])
    subs_ms = len(killed_of_tests) * 1.0 / len(clusters)
    return subs_ms
#~ def get_ms()

def getProjMatricesLabelFiles(in_top_dir, proj, commit=None):
    if commit is not None:
        label_data_folder = os.path.join(in_top_dir, "label_data")
        matrices_folder = os.path.join(in_top_dir, "matrices")
        mat_file = os.path.join(matrices_folder, prog, commit, "STRONG_MUTATION.csv") # TODO, checkz
        label_file = os.path.join(label_data_folder, prog, commit, label.json)
        mut_info_file = os.path.join(matrices_folder, prog, commit, "mutantsInfos.json")
    else:
        mat_file = os.path.join(in_top_dir, 'SEMu-Experiement-data', 'semu_cleaned_data', proj, 'STRONG_MUTATION.csv')
        label_file = os.path.join(in_top_dir, 'SEMu-Experiement-data', 'semu_cleaned_data', proj, 'subsuming-clusters.json')
        mut_info_file = os.path.join(in_top_dir, 'SEMu-Experiement-data', 'semu_cleaned_data', proj, "mutantsInfos.json")
    return mat_file, label_file, mut_info_file
#~ def getProjMatricesLabelFiles()

def load_data(in_top_dir, model_in_dir, cache=True):
    cache_file = os.path.join(model_in_dir, "_cache_file.json")
    
    # Get the project list
    if os.path.isfile(cache_file):
        all_tests, mutants_to_killingtests, tests_to_killed_mutants = subs_load.common_fs.loadJSON(cache_file)
        #cache_projs = set(all_tests)
        print ("# Cache loaded with {} projects.".format(len(all_tests)))
    else:
        #cache_projs = set()
        all_tests = {}
        mutants_to_killingtests = {}
        tests_to_killed_mutants = {}
	
    all_mutants = {}
    machine_translation_mutants = {}
    decision_trees_mutants = {}
    tests_to_killed_subs_cluster = {}
    mutant_to_subs_cluster = {}
    subs_cluster_to_mutants = {}

    machine_translation_muts_json = os.path.join(model_in_dir, "predicted_mutants.json")
    all_muts_json = os.path.join(model_in_dir, "all_mutants.json")
    decision_trees_muts_json = os.path.join(model_in_dir, "projects_probabilities.json")

    update_cache = False
    
    machine_translation_muts_obj = subs_load.common_fs.loadJSON(machine_translation_muts_json)
    decision_trees_muts_obj = subs_load.common_fs.loadJSON(decision_trees_muts_json)
    all_muts_obj = subs_load.common_fs.loadJSON(all_muts_json)
    for obj in (machine_translation_muts_obj, all_muts_obj, decision_trees_muts_obj):
        for pname, pobj in obj.items():
            if type(pobj) == dict:
                new_obj = {}
                for k,v in pobj.items():
                    if 'mart_0:' not in k:
                        new_obj['mart_0:' + k] = v
                    else:
                        new_obj[k] = v
                pobj.clear()
                pobj.update(new_obj)
            else:
                for i in range(len(pobj)):
                    if 'mart_0:' not in pobj[i]:
                        pobj[i] = 'mart_0:' + pobj[i]
    
    assert set(machine_translation_muts_obj) == set(decision_trees_muts_obj), \
                "project mismatch between MT and DT. \n >> MT - DT = {}\n >> DT - MT = {}\n".format(
                        set(machine_translation_muts_obj) - set(decision_trees_muts_obj),
                        set(decision_trees_muts_obj) - set(machine_translation_muts_obj)
                    )
    
    # Load projects data
    tq_data = tqdm.tqdm(list(machine_translation_muts_obj))
    for pname in tq_data:
        tq_data.set_description("Loading {} ...".format(pname))
        if '_' in pname:
            prog, commit = pname.split('_')
        else:
            commit = None
            prog = pname
        #if pname not in not_cached:
            #continue
        
        # get clusters
        sm_mat_file, label_data_file, mut_info_file = getProjMatricesLabelFiles(in_top_dir, prog, commit=commit)
        raw_subs_clust = subs_load.common_fs.loadJSON(label_data_file)['subsume'][1]
        mutant_to_subs_cluster[pname] = {}
        subs_cluster_to_mutants[pname] = {}
        for c_id, c in enumerate(raw_subs_clust):
            subs_cluster_to_mutants[pname][c_id] = list(c)
            for m_id in c:
                assert m_id not in mutant_to_subs_cluster[pname]
                mutant_to_subs_cluster[pname][m_id] = c_id
        
        if pname not in all_tests:
            all_tests[pname], mutants_to_killingtests[pname], tests_to_killed_mutants[pname] = subs_load.load(sm_mat_file)
            update_cache = True
        
        # Add mutants not in matrix but in mutinfo
        minf = subs_load.common_fs.loadJSON(mut_info_file)
        for mut in minf:
            mut = "mart_0:"+mut
            if mut not in mutants_to_killingtests[pname]:
                mutants_to_killingtests[pname][mut] = []

        all_mutants[pname] = all_muts_obj[pname]
        machine_translation_mutants[pname] = machine_translation_muts_obj[pname]
        decision_trees_mutants[pname] = decision_trees_muts_obj[pname]
        
	#assert set(all_mutants[pname]) == set(decision_trees_mutants[pname]), \
        assert len(set(all_mutants[pname]) - set(decision_trees_mutants[pname])) == 0, \
                "mismatch between all_mutants and DT for pname {}.\n DT - all = {}. \n all - DT = {}".format(
                        #pname, set(decision_trees_mutants[pname]) - set(all_mutants[pname]), set(all_mutants[pname]) - set(decision_trees_mutants[pname])) #DBG
                        pname, len(set(decision_trees_mutants[pname]) - set(all_mutants[pname])), len(set(all_mutants[pname]) - set(decision_trees_mutants[pname]))) #DBG
	
        tests_to_killed_subs_cluster[pname] = {}
        for t, kmuts in tests_to_killed_mutants[pname].items():
            tests_to_killed_subs_cluster[pname][t] = set()
            for km in kmuts:
                if km in mutant_to_subs_cluster[pname]:
                    tests_to_killed_subs_cluster[pname][t].add(mutant_to_subs_cluster[pname][km])

    if update_cache:
        subs_load.common_fs.dumpJSON([all_tests, mutants_to_killingtests, tests_to_killed_mutants], cache_file)
        print ("# Cache Written, with {} projects!".format(len(all_tests)))
        
    return all_tests, all_mutants, machine_translation_mutants, decision_trees_mutants, mutants_to_killingtests, \
            tests_to_killed_mutants, tests_to_killed_subs_cluster, mutant_to_subs_cluster, subs_cluster_to_mutants
#~ def load_data()

def main():
    start_time = time.time()
    #relmut_pred_file = None
    if len(sys.argv) != 3:
        #if len(sys.argv) == 4:
        #    relmut_pred_file = os.path.abspath(sys.argv[3])
        #    if not os.path.isfile(relmut_pred_file):
        #        error_exit("relmut_pred_file not existing")
        #else:
        #    error_exit("expected 3 or 2 args. got {}". format(len(sys.argv)))
        error_exit("expected 2 args. got {}". format(len(sys.argv)))
    in_top_dir = os.path.abspath(sys.argv[1])
    out_top_dir = in_top_dir
    model_in_dir = os.path.abspath(sys.argv[2])
    if not os.path.isdir(in_top_dir):
        error_exit("in top dir missing ({})".format(in_top_dir))
    if not os.path.isdir(out_top_dir):
        error_exit("out top dir missing ({})".format(out_top_dir))
    if not os.path.isdir(model_in_dir):
        error_exit("model in dir missing ({})".format(model_in_dir))

    out_folder = os.path.join(out_top_dir, "ANALYSIS_OUTPUT")
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    # load data
    print("# LOADING DATA ...")
    all_tests, all_mutants, machine_translation_mutants, decision_trees_mutants, mutants_to_killingtests, \
        tests_to_killed_mutants, tests_to_killed_subs_cluster, mutant_to_subs_cluster, subs_cluster_to_mutant = \
								              load_data(in_top_dir, model_in_dir, cache=True)

    # Save some data about prediction cluster coverage
    pred_clust_cov = {}
    for proj, mut_list in machine_translation_mutants.items():
        cc_prop = {}
        for c, m_set in subs_cluster_to_mutant[proj].items():
            cc_prop[c] = len(set(mut_list) & set(m_set)) * 100.0 / len(m_set)
        pred_clust_cov[proj] = cc_prop
    pred_clust_cov_file = os.path.join(out_folder, "Pred_machine_translation-subs_cluster_coverage.json")
    subs_load.common_fs.dumpJSON(pred_clust_cov, pred_clust_cov_file, pretty=True)
    
    # Simulation
    print ("# Running Simulations ...")
    for fixed_size in (None,):# 5, 10, 20, 30, "subs_cluster_size"):
        proj2used_size = {}
        sim_res = {}
        other_sim_res = {}
        anal_sim_res = {}
        anal_other_sim_res = {}
        testexec_sim_res = {}
        test_exec_other_sim_res = {}
        mutant_analysis_cost_obj = {}
        test_execution_cost_obj = {}
        tq_data = tqdm.tqdm(list(all_tests))
        for proj in tq_data:
            tq_data.set_description("Simulating for {} ...".format(proj))

            if fixed_size == "subs_cluster_size":
                used_fixed_size = len(subs_cluster_to_mutant[proj])
                proj2used_size[proj] = used_fixed_size
            elif fixed_size is None:
                used_fixed_size = len(machine_translation_mutants[proj])
                proj2used_size[proj] = used_fixed_size
            else:
                assert type(fixed_size) == int
                used_fixed_size = fixed_size
                
            if used_fixed_size > len(machine_translation_mutants[proj]):
                # not enough data to check
                continue
                
            sim_res[proj] = {RANDOM: None, PRED_MACHINE_TRANSLATION: None, PRED_DECISION_TREES: None}
            sim_res[proj][RANDOM], sim_res[proj][PRED_MACHINE_TRANSLATION], \
                                sim_res[proj][PRED_DECISION_TREES], \
                                mutant_analysis_cost_obj[proj], \
                                test_execution_cost_obj[proj] = simulation(NUM_REPETITIONS, all_tests[proj], \
                                                                             all_mutants[proj], \
                                                                             machine_translation_mutants[proj], \
                                                                             decision_trees_mutants[proj], \
                                                                             tests_to_killed_mutants[proj], \
                                                                             tests_to_killed_subs_cluster[proj], \
                                                                             mutants_to_killingtests[proj], \
                                                                             fixed_size=used_fixed_size)
            
            print ("## Doing additional sim ...")
            machine_translation_sMS2size = {}
            mt_sMS2analysed = {}
            mt_sMS2testexec = {}
            mt_analysed2sMS = {}
            mt_testexec2sMS = {}
            for pos,sMS in enumerate(sim_res[proj][PRED_MACHINE_TRANSLATION]):
                anal_here = mutant_analysis_cost_obj[proj][PRED_MACHINE_TRANSLATION][pos]
                testexec_here = test_execution_cost_obj[proj][PRED_MACHINE_TRANSLATION][pos]
                if sMS not in machine_translation_sMS2size:
                    machine_translation_sMS2size[sMS] = []
                    mt_sMS2analysed[sMS] = []
                    mt_sMS2testexec[sMS] = []
                if anal_here not in mt_analysed2sMS:
                    mt_analysed2sMS[anal_here] = []
                if testexec_here not in mt_testexec2sMS:
                    mt_testexec2sMS[testexec_here] = []
                machine_translation_sMS2size[sMS].append(used_fixed_size)
                mt_sMS2analysed[sMS].append(anal_here)
                mt_sMS2testexec[sMS].append(testexec_here)
                mt_analysed2sMS[anal_here].append(sMS)
                mt_testexec2sMS[testexec_here].append(sMS)
            other_sim_res[proj], anal_sim_res[proj], anal_other_sim_res[proj], \
                       testexec_sim_res[proj], test_exec_other_sim_res[proj] = additional_simulation (SUB_REPET_NUM, all_tests[proj], \
                                                                                                        all_mutants[proj], \
                                                                                                        decision_trees_mutants[proj], \
                                                                                                        tests_to_killed_mutants[proj], \
                                                                                                        tests_to_killed_subs_cluster[proj], \
                                                                                                        mutants_to_killingtests[proj], \
                                                                                                        machine_translation_sMS2size, \
                                                                                                        mt_sMS2analysed=mt_sMS2analysed, \
                                                                                                        mt_sMS2testexec=mt_sMS2testexec, \
                                                                                                        mt_analysed2sMS=mt_analysed2sMS, \
                                                                                                        mt_testexec2sMS=mt_testexec2sMS, \
                                                                                                        parallel_count=16)

        # Store sizes
        if len(proj2used_size) > 0:
            saved_size_obj = {
                              'PREDICTED_SIZES': proj2used_size,
                              'TOTAL_SIZES': {p: len(am) for p, am in all_mutants.items()},
                              'SUBSUMING_SIZES': {p: len(sm) for p, sm in mutant_to_subs_cluster.items()}
                             }
            size_file_prefix = os.path.join(out_folder, "used_fixed_size-{}".format("pred_size" if fixed_size is None else fixed_size))
            subs_load.common_fs.dumpJSON(saved_size_obj, size_file_prefix+'.json' , pretty=True)
            size_prop = {'PREDICTED_SIZES': [], 'SUBSUMING_SIZES': []}
            for proj in saved_size_obj['PREDICTED_SIZES']:
                size_prop['PREDICTED_SIZES'].append(saved_size_obj['PREDICTED_SIZES'][proj] * 1.0 / saved_size_obj['TOTAL_SIZES'][proj])
                size_prop['SUBSUMING_SIZES'].append(saved_size_obj['SUBSUMING_SIZES'][proj] * 1.0 / saved_size_obj['TOTAL_SIZES'][proj])
            plot.plotBoxesHorizontal(size_prop, list(size_prop), size_file_prefix, plot.colors_bw, ylabel="Mutants Proportion" , yticks_range=plot.np.arange(0,1.01,0.2))
            
        print("# Plotting ...")
        for fname_prefix, metric, data_obj, is_proportion in [('SELECTION-', 'Subsuming MS', sim_res, True), 
                                 ('SEL-UNUSED-', 'Proportion of Mutant Analysed' if Use_proportion_analysed_mutants else '# Mutant Analysed', \
                                                                                    mutant_analysis_cost_obj, Use_proportion_analysed_mutants), 
                                 ('SEL-UNUSED', '# Tests Executed', test_execution_cost_obj, False), 
                                 ('SELECTION-', 'Selection Size for Same Subsuming MS', other_sim_res, True), 
                                 ('ANALYSIS-', 'Subsuming MS', anal_sim_res, True), 
                                 ('ANALYSIS-', 'Analysed Mutants for Same Subsuming MS', anal_other_sim_res, True), 
                                 ('TESTEXECUTION-', 'Subsuming MS', testexec_sim_res, True), 
                                 ('TESTEXECUTION-', 'Test Execution for same Subsuming MS', test_exec_other_sim_res, False)]:
            # Plot box plot
            print ("@Plot: Plotting {} - {} ...".format(fname_prefix, metric))
            image_file = os.path.join(out_folder, fname_prefix + metric.replace('#', 'num').replace(' ', '_') + '-' + \
                                                                "boxplot_all-{}".format(("pred_size" if fixed_size is None else fixed_size)))
            image_file_agg = os.path.join(out_folder, fname_prefix + metric.replace('#', 'num').replace(' ', '_') + '-' + \
                                                                "merged_boxplot_all-{}".format(("pred_size" if fixed_size is None else fixed_size)))
            median_file_agg = os.path.join(out_folder, fname_prefix + metric.replace('#', 'num').replace(' ', '_') + '-' + \
                                                                "merged_boxplot_all-{}".format(("pred_size" if fixed_size is None else fixed_size)) + "-median.json")
            stattest_json_file_agg = os.path.join(out_folder, fname_prefix + metric.replace('#', 'num').replace(' ', '_') + '-' + \
                                                                "merged_boxplot_all-{}".format(("pred_size" if fixed_size is None else fixed_size)) + "-stat_test.json")

            order = [PRED_MACHINE_TRANSLATION, PRED_DECISION_TREES, RANDOM]
            data_df = []
            merged_dat = {t: [] for t in order}
            max_metric_val = 0
            for proj, p_dat in data_obj.items():
                for tech, t_dat in p_dat.items():
                    for metric_val in t_dat:
                        data_df.append({'Program': proj[:10], metric: metric_val, 'Tech': tech})
                        merged_dat[tech].append(metric_val)
                        if metric_val > max_metric_val:
                            max_metric_val = metric_val
            if len(data_df) > 0:
                if is_proportion:
                    yticks_range = plot.np.arange(0,1.01,0.2)
                else:
                    yticks_range = plot.np.linspace(0, max_metric_val + 1, 10)

                data_df = pd.DataFrame(data_df)
                plot.plt.figure(figsize=(16, 8)) 
                ax = sns.boxplot(x="Program", y=metric, hue="Tech", data=data_df, palette="Set3", medianprops={'linewidth':5}) #, linewidth=2.5)
                plot.plt.yticks(yticks_range) #, fontsize=fontsize)
                plot.plt.savefig(image_file+".pdf", format='pdf') #, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plot.plt.close('all')

                median_list = plot.plotBoxes(merged_dat, order, image_file_agg, plot.colors_bw, ylabel=metric, yticks_range=yticks_range)
                subs_load.common_fs.dumpJSON({order[i]: median_list[i] for i in range(len(median_list))}, median_file_agg, pretty=True)
                
                # Stat_test agg
                stat_test_obj = inner_stattest(merged_dat, stattest_json_file_agg, order=order)
                print ("   :) @@Plot: Done for {} - {}!".format(fname_prefix, metric))
            else:
                print ("   :( @@Plot: No Data for {} - {}!".format(fname_prefix, metric))
                
    print("@DONE (after {}  h:min:sec) !".format(str(timedelta(seconds=(time.time() - start_time)))))
#~ def main()

def simulation(num_repet, test_list, mutant_list, machine_translation_mutant_list,
                  decision_trees_mutant_dict,
                  tests_to_killed_mutants, tests_to_killed_subs_cluster, 
                  mutants_to_killingtests, fixed_size=None):
    ordered_tests_mode = False
    
    if fixed_size is None:
        assert machine_translation_mutant_list is not None, "Must have MT here"
        selection_size = len(machine_translation_mutant_list)
    else:
        selection_size = fixed_size
    random_test_suites = []
    machine_translation_test_suites = []
    decision_trees_test_suites = []
    mutant_analysis_cost = {n: [] for n in (RANDOM, PRED_MACHINE_TRANSLATION, PRED_DECISION_TREES)}
    test_execution_cost = {n: [] for n in (RANDOM, PRED_MACHINE_TRANSLATION, PRED_DECISION_TREES)}
    
    repet_bar = tqdm.tqdm(range(num_repet), desc='Repetitions', leave=False)
    for repet_id in repet_bar:
        # randomly sample
        random_M = set(random.sample(mutant_list, selection_size))
        machine_translation_M = set(random.sample(machine_translation_mutant_list, selection_size)) if machine_translation_mutant_list is not None else None
        decision_trees_M = set(sorted(
                                        random.sample(mutant_list, len(mutant_list)), 
                                        reverse=True, 
                                        key=lambda x: float(decision_trees_mutant_dict[x])
                                    ) [:selection_size])

        random_test_suites.append([])
        if machine_translation_M is not None:
            machine_translation_test_suites.append([])
        decision_trees_test_suites.append([])
        
        if ordered_tests_mode:
            test_order = list(test_list)
            random.shuffle(test_order)
            for t in test_order:
                # get killed mutants
                rand_kill_mut = set(tests_to_killed_mutants[t]) & random_M
                machine_translation_kill_mut = set(tests_to_killed_mutants[t]) & machine_translation_M if machine_translation_M is not None else None
                decision_trees_kill_mut = set(tests_to_killed_mutants[t]) & decision_trees_M
                if len(rand_kill_mut) > 0:
                    random_test_suites[-1].append(t)
                    random_M -= rand_kill_mut
                if machine_translation_M is not None and len(machine_translation_kill_mut) > 0:
                    machine_translation_test_suites[-1].append(t)
                    machine_translation_M -= machine_translation_kill_mut
                if len(decision_trees_kill_mut) > 0:
                    decision_trees_test_suites[-1].append(t)
                    decision_trees_M -= decision_trees_kill_mut
        else:
            tasks = [(RANDOM, random_M, random_test_suites), \
                                                      (PRED_DECISION_TREES, decision_trees_M, decision_trees_test_suites)]
            if machine_translation_M is not None:
                tasks.append((PRED_MACHINE_TRANSLATION, machine_translation_M, machine_translation_test_suites))
            for techname, rem_set, TS_list in tasks:
                analysed_muts_num = 0
                exec_tests_num = 0
                while len(rem_set) > 0:
                    # pick a mutant
                    if techname == PRED_DECISION_TREES: 
                        # decision trees
                        m = max(rem_set, key=lambda x: float(decision_trees_mutant_dict[x]))
                    else:
                        m = random.choice(tuple(rem_set))
                        
                    # generate a test to kill m
                    if m not in mutants_to_killingtests:
                        error_exit("Mutant not in mutants to killingtests. \n mutants_to_killing tests is {}. \nMissing mutants is {}".format(\
                                                                                                                list(mutants_to_killingtests), m))
          
                    analysed_muts_num += 1
          
                    if len(mutants_to_killingtests[m]) > 0:
                        # The tests is executed with all remaining mutant
                        exec_tests_num += len(rem_set)
                        
                        # generate a test
                        t = random.choice(mutants_to_killingtests[m])
                        TS_list[-1].append(t)
                        # remove all collaterally killed mutants
                        rem_set -= set(tests_to_killed_mutants[t]) & rem_set
                    else:
                        rem_set -= {m}
                        
                if Use_proportion_analysed_mutants:
                    mutant_analysis_cost[techname].append(analysed_muts_num * 1.0 / len(mutant_list))
                else:
                    mutant_analysis_cost[techname].append(analysed_muts_num)
                test_execution_cost[techname].append(exec_tests_num)
                        

    # Computer sMS
    rand_sMS = []
    machine_translation_sMS = []
    decision_trees_sMS = []
    for ts in random_test_suites:
        rand_sMS.append(get_subs_ms(ts, tests_to_killed_subs_cluster))
    for ts in machine_translation_test_suites:
        machine_translation_sMS.append(get_subs_ms(ts, tests_to_killed_subs_cluster))
    for ts in decision_trees_test_suites:
        decision_trees_sMS.append(get_subs_ms(ts, tests_to_killed_subs_cluster))

    return rand_sMS, machine_translation_sMS, decision_trees_sMS, mutant_analysis_cost, test_execution_cost
#~ def simulation()
    
def additional_sub_parallel (in_data_valuerange_args_kwargs):
    sMS2selsize = {RANDOM: {}, PRED_DECISION_TREES: {}}
    sMS2analysed = {RANDOM: {}, PRED_DECISION_TREES: {}}
    sMS2testexec = {RANDOM: {}, PRED_DECISION_TREES: {}}
    analysed2sMS = {RANDOM: {}, PRED_DECISION_TREES: {}}
    testexec2sMS = {RANDOM: {}, PRED_DECISION_TREES: {}}
    multi_sizes_bar = tqdm.tqdm(in_data_valuerange_args_kwargs[0], desc='Multiple Sizes', leave=False)
    args = in_data_valuerange_args_kwargs[1]
    kwargs = dict(in_data_valuerange_args_kwargs[2]) # make a new copy
    for fixed_size in multi_sizes_bar:
        kwargs['fixed_size'] = fixed_size
        rand_sMS, machine_translation_sMS, dt_sMS, \
                            mutant_analysis_cost, test_execution_cost = simulation (*args, **kwargs)
        for indat, tech in ((rand_sMS, RANDOM), (dt_sMS, PRED_DECISION_TREES)):
            for pos,sMS in enumerate(indat):
                if sMS not in sMS2selsize[tech]:
                    sMS2selsize[tech][sMS] = []
                    sMS2analysed[tech][sMS] = []
                    sMS2testexec[tech][sMS] = []
                sMS2selsize[tech][sMS].append(fixed_size)
                sMS2analysed[tech][sMS].append(mutant_analysis_cost[tech][pos])
                sMS2testexec[tech][sMS].append(test_execution_cost[tech][pos])
                if mutant_analysis_cost[tech][pos] not in analysed2sMS:
                    analysed2sMS[tech][mutant_analysis_cost[tech][pos]] = []
                analysed2sMS[tech][mutant_analysis_cost[tech][pos]].append(sMS)
                if test_execution_cost[tech][pos] not in testexec2sMS:
                    testexec2sMS[tech][test_execution_cost[tech][pos]] = []
                testexec2sMS[tech][test_execution_cost[tech][pos]].append(sMS)
            
    return sMS2selsize, sMS2analysed, sMS2testexec, analysed2sMS, testexec2sMS
#~ def additional_sub_parallel()

def additional_simulation (num_sub_repet, test_list, mutant_list, 
                              decision_trees_mutant_dict,
                              tests_to_killed_mutants, tests_to_killed_subs_cluster, 
                              mutants_to_killingtests, machine_translation_sMS2size,
                              mt_sMS2analysed=None, 
                              mt_sMS2testexec=None, \
                              mt_analysed2sMS=None, \
                              mt_testexec2sMS=None, \
                              use_raw_number=False, parallel_count=1):
    
    assert parallel_count > 0, "invalid parallel_count"
    
    args = [num_sub_repet, test_list, mutant_list, None, decision_trees_mutant_dict,
                           tests_to_killed_mutants, tests_to_killed_subs_cluster, mutants_to_killingtests]
    kwargs = {}
    
    list_in_data_tqdmrange_args_kwargs = []
    for para_ind in range(parallel_count):
        tqdmrange = range (1 + para_ind, len(mutant_list) + 1, parallel_count)
        list_in_data_tqdmrange_args_kwargs.append((tqdmrange, args, kwargs))
    
    with Pool(parallel_count) as p:
        map_list = p.map(additional_sub_parallel, list_in_data_tqdmrange_args_kwargs)
        
    # Merge maps
    sMS2selsize = {RANDOM: {}, PRED_DECISION_TREES: {}}
    sMS2analysed = {RANDOM: {}, PRED_DECISION_TREES: {}}
    sMS2testexec = {RANDOM: {}, PRED_DECISION_TREES: {}}
    analysed2sMS = {RANDOM: {}, PRED_DECISION_TREES: {}}
    testexec2sMS = {RANDOM: {}, PRED_DECISION_TREES: {}}
    for res_tmp in map_list:
        for partial, aggregated in zip(res_tmp, [sMS2selsize, sMS2analysed, sMS2testexec, analysed2sMS, testexec2sMS]):
            for tech, tdat in partial.items():
                for metric, ss in tdat.items():
                    if metric not in aggregated[tech]:
                        aggregated[tech][metric] = []
                    aggregated[tech][metric] += ss
            
    #for aggregated in [sMS2selsize, sMS2analysed, sMS2testexec, analysed2sMS, testexec2sMS]:
    #    for tech in aggregated:
    #        for sMS in aggregated[tech]:
    #            aggregated[tech][sMS] = list (aggregated[tech][sMS])
            
    sorted_keys_sMS2size = {RANDOM: sorted(list(sMS2selsize[RANDOM])), PRED_DECISION_TREES: sorted(list(sMS2selsize[PRED_DECISION_TREES]))}
    sorted_keys_sMS2analysed = {RANDOM: sorted(list(sMS2analysed[RANDOM])), PRED_DECISION_TREES: sorted(list(sMS2analysed[PRED_DECISION_TREES]))}
    sorted_keys_sMS2testexec = {RANDOM: sorted(list(sMS2testexec[RANDOM])), PRED_DECISION_TREES: sorted(list(sMS2testexec[PRED_DECISION_TREES]))}
    sorted_keys_analysed2sMS = {RANDOM: sorted(list(analysed2sMS[RANDOM])), PRED_DECISION_TREES: sorted(list(analysed2sMS[PRED_DECISION_TREES]))}
    sorted_keys_testexec2sMS = {RANDOM: sorted(list(testexec2sMS[RANDOM])), PRED_DECISION_TREES: sorted(list(testexec2sMS[PRED_DECISION_TREES]))}
    
    def get_other_values (in_sMS, mt_size_list, dict_data, sorted_keys):
        pos_r = max(bisect.bisect_right(sorted_keys[RANDOM], in_sMS) - 1, 0)
        sMS_r = sorted_keys[RANDOM][pos_r]
        pos_d = max(bisect.bisect_right(sorted_keys[PRED_DECISION_TREES], in_sMS) - 1, 0)
        sMS_d = sorted_keys[PRED_DECISION_TREES][pos_d]
        
        min_ss = min (len(dict_data[RANDOM][sMS_r]), len(dict_data[PRED_DECISION_TREES][sMS_d]), len(mt_size_list))
        rand_size = random.sample(dict_data[RANDOM][sMS_r], min_ss)
        dt_size = random.sample(dict_data[PRED_DECISION_TREES][sMS_d], min_ss)
        mt_size = random.sample(mt_size_list, min_ss)
        return rand_size, dt_size, mt_size
    #~def get_other_values ()
    
    sizes = {RANDOM: [], PRED_DECISION_TREES: [], PRED_MACHINE_TRANSLATION: []}
    analysed_sMS = {RANDOM: [], PRED_DECISION_TREES: [], PRED_MACHINE_TRANSLATION: []}
    analysed = {RANDOM: [], PRED_DECISION_TREES: [], PRED_MACHINE_TRANSLATION: []}
    testexec_sMS = {RANDOM: [], PRED_DECISION_TREES: [], PRED_MACHINE_TRANSLATION: []}
    testexec = {RANDOM: [], PRED_DECISION_TREES: [], PRED_MACHINE_TRANSLATION: []}
    
    for out_obj, mt_data, cmp_data, cmp_sorted_keys, is_mutant_proportion in [
                        (sizes, machine_translation_sMS2size, sMS2selsize, sorted_keys_sMS2size, (not use_raw_number)),
                        (analysed, mt_sMS2analysed, sMS2analysed, sorted_keys_sMS2analysed, (not use_raw_number)),
                        (testexec, mt_sMS2testexec, sMS2testexec, sorted_keys_sMS2testexec, False),
                        (analysed_sMS, mt_analysed2sMS, analysed2sMS, sorted_keys_analysed2sMS, False),
                        (testexec_sMS, mt_testexec2sMS, testexec2sMS, sorted_keys_testexec2sMS, False),                        
                                                                    ]:    
        for mt_key, mt_val_list in mt_data.items():
            rand_vals, dt_vals, mt_vals = get_other_values (mt_key, mt_val_list, cmp_data, cmp_sorted_keys)
            if is_mutant_proportion:
                out_obj[PRED_MACHINE_TRANSLATION] += [s * 1.0 / len(mutant_list) for s in mt_vals]
                out_obj[PRED_DECISION_TREES] += [s * 1.0 / len(mutant_list) for s in dt_vals]
                out_obj[RANDOM] += [s * 1.0 / len(mutant_list) for s in rand_vals ]
            else:
                out_obj[PRED_MACHINE_TRANSLATION] += mt_vals
                out_obj[PRED_DECISION_TREES] += dt_vals
                out_obj[RANDOM] += rand_vals
            
    return sizes, analysed_sMS, analysed, testexec_sMS, testexec
#~ def additional_simulation ()
    
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
