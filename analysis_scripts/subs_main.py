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

import seaborn as sns

import subs_load
import plot
 
STOP_AT_N_MUTANTS = 100 # None

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
        mat_file = os.path.join(matrices_folder, prog, commit, "STRONG_MUTATION.csv") # TODO, check
        label_file = os.path.join(label_data_folder, prog, commit, label.json)
    else:
        mat_file = os.path.join(in_top_dir, 'SEMu-Experiement-data', 'semu_cleaned_data', proj, 'STRONG_MUTATION.csv')
        label_file = os.path.join(in_top_dir, 'SEMu-Experiement-data', 'semu_cleaned_data', proj, 'subsuming-clusters.json')
    return mat_file, label_file
#~ def getProjMatricesLabelFiles()

def load_data(in_top_dir, cache_file):
    # Get the project list
    if os.path.isfile(cache_file):
        all_tests, mutants_to_killingtests, tests_to_killed_mutants = subs_load.common_fs.loadJSON(cache_file)
        cache_projs = set(all_tests)
    else:
        cache_projs = set()
        all_tests = {}
        all_mutants = {}
        pred_mutants = {}
        mutants_to_killingtests = {}
        tests_to_killed_mutants = {}
        tests_to_killed_subs_cluster = {}
        mutant_to_subs_cluster = {}
        subs_cluster_to_mutants = {}

    pred_muts_json = os.path.join(in_top_dir, "predicted_mutants.json")
    all_muts_json = os.path.join(in_top_dir, "all_mutants.json")

    #update_cache = (len(not_cached) > 0)
    
    pred_muts_obj = subs_load.common_fs.loadJSON(pred_muts_json)
    all_muts_obj = subs_load.common_fs.loadJSON(all_muts_json)
    for obj in (pred_muts_obj, all_muts_obj):
        for pname, pobj in obj.items():
            for i in range(len(pobj)):
                if 'mart_0:' not in pobj[i]:
                    pobj[i] = 'mart_0:' + pobj[i]
    
    # Load projects data
    tq_data = tqdm.tqdm(list(pred_muts_obj))
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
        sm_mat_file, label_data_file = getProjMatricesLabelFiles(in_top_dir, prog, commit=commit)
        raw_subs_clust = subs_load.common_fs.loadJSON(label_data_file)['subsume'][1]
        mutant_to_subs_cluster[pname] = {}
        subs_cluster_to_mutants[pname] = {}
        for c_id, c in enumerate(raw_subs_clust):
            subs_cluster_to_mutants[pname][c_id] = list(c)
            for m_id in c:
                assert m_id not in mutant_to_subs_cluster[pname]
                mutant_to_subs_cluster[pname][m_id] = c_id
        
        all_tests[pname], mutants_to_killingtests[pname], tests_to_killed_mutants[pname] = subs_load.load(sm_mat_file)

        all_mutants[pname] = all_muts_obj[pname]
        pred_mutants[pname] = pred_muts_obj[pname]
        
        tests_to_killed_subs_cluster[pname] = {}
        for t, kmuts in tests_to_killed_mutants[pname].items():
            tests_to_killed_subs_cluster[pname][t] = set()
            for km in kmuts:
                if km in mutant_to_subs_cluster[pname]:
                    tests_to_killed_subs_cluster[pname][t].add(mutant_to_subs_cluster[pname][km])

    #if update_cache:
        #load.common_fs.dumpJSON([all_tests, fault_tests, relevant_mutants_to_relevant_tests, mutants_to_killingtests, tests_to_killed_mutants], cache_file)
        
    return all_tests, all_mutants, pred_mutants, mutants_to_killingtests, tests_to_killed_mutants, tests_to_killed_subs_cluster, mutant_to_subs_cluster, subs_cluster_to_mutants
#~ def load_data()

def main():
    relmut_pred_file = None
    if len(sys.argv) != 3 and len(sys.argv) != 2:
        if len(sys.argv) == 4:
            relmut_pred_file = os.path.abspath(sys.argv[3])
            if not os.path.isfile(relmut_pred_file):
                error_exit("relmut_pred_file not existing")
        else:
            error_exit("expected 3 or 2 args. got {}". format(len(sys.argv)))
    in_top_dir = os.path.abspath(sys.argv[1])
    if len(sys.argv) == 3:
        out_top_dir = os.path.abspath(sys.argv[2])
    else:
        out_top_dir = in_top_dir
    if not os.path.isdir(in_top_dir):
        error_exit("in top dir missing ({})".format(in_top_dir))
    if not os.path.isdir(out_top_dir):
        error_exit("out top dir missing ({})".format(out_top_dir))

    out_folder = os.path.join(out_top_dir, "ANALYSIS_OUTPUT")
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    # load data
    print("# LOADING DATA ...")
    cache_file = os.path.join(out_folder, "cache_file.json")
    all_tests, all_mutants, pred_mutants, mutants_to_killingtests, tests_to_killed_mutants, tests_to_killed_subs_cluster, mutant_to_subs_cluster, subs_cluster_to_mutant = \
								              load_data(in_top_dir, cache_file)

    num_repet = 1000
    
    # Simulation
    print ("# Running Simulations ...")
    for fixed_size in (None, 10, 20, 30, "subs_cluster_size"):
        sim_res = {}
        tq_data = tqdm.tqdm(list(all_tests))
        for proj in tq_data:
            tq_data.set_description("Loading {} ...".format(proj))

            if fixed_size == "subs_cluster_size":
                used_fixed_size = len(subs_cluster_to_mutant[proj])
            elif fixed_size is None:
                used_fixed_size = len(pred_mutants[proj])
                fixed_size = "pred_size"
            else:
                used_fixed_size = fixed_size
                
            if used_fixed_size > len(pred_mutants[proj]):
                # not enough data to check
                continue
                
            sim_res[proj] = {'RANDOM': None, "PREDICTED": None}
            sim_res[proj]["RANDOM"], sim_res[proj]["PREDICTED"] = simulation(num_repet, all_tests[proj], \
                                                                             all_mutants[proj], \
                                                                             pred_mutants[proj], \
                                                                             tests_to_killed_mutants[proj], \
                                                                             tests_to_killed_subs_cluster[proj], \
                                                                             fixed_size=used_fixed_size)

        print("# Plotting ...")
        # Plot box plot
        image_file = os.path.join(out_folder, "boxplot_all-{}-{}".format(used_fixed_size, fixed_size))
        data_df = []
        for proj, p_dat in sim_res.items():
            for tech, t_dat in p_dat.items():
                for sMS in t_dat:
                    data_df.append({'Program': proj[:6], 'Subsuming MS': sMS, 'Tech': tech})
        if len(data_df) > 0:
            data_df = pd.DataFrame(data_df)
            plt.figure(figsize=(16, 8)) 
            ax = sns.boxplot(x="Program", y="Subsuming MS", hue="Tech", data=data_df, palette="Set3") #, linewidth=2.5)
            plot.plt.savefig(image_file+".pdf", format='pdf') #, bbox_extra_artists=(lgd,), bbox_inches='tight')
            plot.plt.close('all')
    
    print("@DONE!")
#~ def main()

def simulation(num_repet, test_list, mutant_list, pred_mutant_list,
                  tests_to_killed_mutants, tests_to_killed_subs_cluster, fixed_size=None):
    if fixed_size is None:
        selection_size = len(pred_mutant_list)
    else:
        selection_size = fixed_size
    random_test_suites = []
    pred_test_suites = []
    for repet_id in range(num_repet):
        # randomly sample
        random_M = set(random.sample(mutant_list, selection_size))
        pred_M = set(pred_mutant_list)

        test_order = list(test_list)
        random.shuffle(test_order)

        random_test_suites.append([])
        pred_test_suites.append([])
        for t in test_order:
            # get killed mutants
            rand_kill_mut = set(tests_to_killed_mutants[t]) & random_M
            pred_kill_mut = set(tests_to_killed_mutants[t]) & pred_M
            if len(rand_kill_mut) > 0:
                random_test_suites[-1].append(t)
                random_M -= rand_kill_mut
            if len(pred_kill_mut) > 0:
                pred_test_suites[-1].append(t)
                pred_M -= pred_kill_mut

    # Computer sMS
    rand_sMS = []
    pred_sMS = []
    for ts in random_test_suites:
        rand_sMS.append(get_subs_ms(ts, tests_to_killed_subs_cluster))
    for ts in pred_test_suites:
        pred_sMS.append(get_subs_ms(ts, tests_to_killed_subs_cluster))

    return rand_sMS, pred_sMS
#~ def simulation()
    
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
