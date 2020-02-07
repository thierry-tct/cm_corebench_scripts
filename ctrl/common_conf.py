
from __future__ import print_function

import glob
import shutil
import os
import sys
import json
import logging
from collections import defaultdict


from muteria.common.mix import GlobalConstants

from muteria.configmanager.configurations import SessionMode
from muteria.configmanager.configurations import TestcaseToolsConfig
from muteria.configmanager.configurations import CriteriaToolsConfig
from muteria.configmanager.configurations import ToolUserCustom

from muteria.drivers import DriversUtils
from muteria.drivers.testgeneration import TestToolType
from muteria.drivers.criteria import CriteriaToolType
from muteria.drivers.criteria import TestCriteria
from muteria.drivers.testgeneration.testcase_formats.system_devtest import system_test_runner
from muteria.drivers.testgeneration.custom_dev_testcase.system_wrappers.native_code import SystemWrapper

test_alias_to_script = defaultdict(list)
devtestlist = []
bashlist = defaultdict(list)
rootlist = defaultdict(list)
expensivelist = defaultdict(list)
fiemaplist = defaultdict(list)
usingExpensiveTest = False
usingRootTest = False
if 'COREUTILS_TEST_EXPENSIVE' in os.environ:
    ctesetv =  os.getenv('COREUTILS_TEST_EXPENSIVE')
    if ctesetv == '1':
       usingExpensiveTest = True  

if 'COREUTILS_TEST_ROOT' in os.environ:
    rootest =  os.getenv('COREUTILS_TEST_ROOT')
    if rootest == '1':
       usingRootTest = True  

def initdevtestlist(repo_root_dir):
#    print("*********** Init Test List ********** \n")
    global devtestlist
    global bashlist
    global rootlist
    cwd = os.getcwd()
    os.chdir(repo_root_dir)
    bashlist_tmp = read_files_folder("./tests/", suffix=".sh")
    perllist_tmp = read_files_folder("./tests/", suffix=".pl")
    
    for e in bashlist_tmp:
        newe = e #.replace("tests_splitting", 'tests')
        shutil.copyfile(e, newe)
        tfname = os.path.basename(newe)
        ppath = os.path.dirname(newe)
        isroot = False
        isExpensive = False
        count = 0
        def find_test_cases():
            nonlocal count, isroot, isExpensive
            with open(newe, 'r') as f:
                 for line in f.readlines():
                     if "require_root_" in line:
                         if usingRootTest:
                            isroot = True
                         else:
                            return
                     if "very_expensive_" in line:
                        if usingExpensiveTest:
                           isExpensive = True
                        else:
                           return
                     if "fiemap_capable_" in line:
                        return
            with open(newe, 'r') as f:
                 for line in f.readlines():
                     if "$TEST_ID" in line:
                         testcase = os.path.normpath(os.path.join(ppath, str(count)+"_"+tfname))
                         bashlist[testcase] = count
                         rootlist[testcase]= isroot
                         expensivelist[testcase] = isExpensive
                         devtestlist.append(testcase)
                         shutil.copyfile(newe, testcase)
                         count = count + 1
                    
        find_test_cases()
        
    for h in perllist_tmp:
        newh = h #.replace("tests_splitting", "tests")
        shutil.copyfile(h, newh)
        tfname = os.path.basename(newh)
        ppath = os.path.dirname(newh)
        with open(newh, 'r') as f:
            def innerloop():
                for line in f.readlines():
                    if "tests_size_splitting" in line:
                        hl = line.split("=")
                        number = int(hl[1].split(";")[0])
                        for i in range(number):
                            testcase = os.path.normpath(os.path.join(ppath, str(i)+"_"+tfname))
                            bashlist[testcase] = i
                            rootlist[testcase] = False
                            devtestlist.append(testcase)
                            shutil.copyfile(newh, testcase)
                        break

            innerloop()
#    return ['./tests/misc/17_chroot-credentials.sh']
    return devtestlist

def read_files_folder(fpath, suffix='.sh'):
  '''
  read all .path files
  :param fpath:
  :return:
  '''
  files = [f for f in glob.glob(fpath + "**/*" + suffix, recursive=True)]
  return files


# Corebench
def initdevtestlist_bench(test_folder):
    global bashlist
    global roottests
    global test_alias_to_script

    alltests_file = os.path.join(test_folder, "splittests_all.txt")
    roottests_file = os.path.join(test_folder, "splittests_root.txt")
    expensivetests_file = os.path.join(test_folder, "splittests_veryexpensive.txt")
    
    alltests = [v.strip() for v in open(alltests_file)]
    roottests = []
    expensivetests = []
    if os.path.isfile(roottests_file):
        roottests = [v.strip() for v in open(roottests_file)]
    if os.path.isfile(expensivetests_file):
        expensivetests = [v.strip() for v in open(expensivetests_file)]
    if not usingRootTest:
        alltests = list(set(alltests) - set(roottests))
    if not usingExpensiveTest:
        alltests = list(set(alltests) - set(expensivetests))
	
    for testcase in alltests:
        rootlist[testcase] = (testcase in roottests)
        d, f = os.path.split(testcase)
        parts = f.split('_')
        if len(parts) > 1:
            assert parts[0].isdigit(), "split problem, first part not digit"
            bashlist[testcase] = int(parts[0])
            test_alias_to_script[testcase] = os.path.join(d, '_'.join(parts[1:]))
        else:
            bashlist[testcase] = None
            test_alias_to_script[testcase] = testcase

    return alltests
#~ def initdevtestlist_bench

#TCT
# Import coreutils make check test bypasser
import muteria
#py_file_path = os.path.join(os.path.dirname(os.path.dirname(muteria.__file__)), \
#                                                     'miscellaneous', 'coreutils_experiments')
#print(py_file_path, os.listdir(py_file_path))
#sys.path.insert(0, py_file_path)
import bypass_make_check_tests
#sys.path.pop(0)
#~

check_dev_testname_in_run = [True]
run_devtest_folder_level_from_root = [0]

def dev_test_runner(test_name, repo_root_dir, *args, **kwargs):
    #global devtestlist
    global bashlist
    global rootlist
    cwd = os.getcwd()
    os.chdir(repo_root_dir)


#    print("&&&&&& file list   {} \n {} \n {}".format(test_name, args, kwargs))
    if 'env_vars' not in kwargs or kwargs['env_vars'] is None:
       kwargs['env_vars'] = {}

    # TCT
    make_check_test_env_vars = bypass_make_check_tests.get_make_check_tests_env_vars()
    kwargs['env_vars'].update(make_check_test_env_vars)
    
    kwargs['env_vars'].update({"RUN_VERY_EXPENSIVE_TESTS":"no", "RUN_EXPENSIVE_TESTS":"yes"})
    #~

    kwargs['env_vars']["TEST_ID"] = str(bashlist[test_name])
    test_script = test_alias_to_script[test_name] 
    tests_rel_path = 'tests'
    if run_devtest_folder_level_from_root[0] == 1:
        os.chdir('tests')
        test_script = os.path.normpath(os.path.join("..",test_script))
        tests_rel_path = '.'
    elif run_devtest_folder_level_from_root[0] == 2:
        os.chdir(os.path.dirname(test_script))
        test_script = os.path.normpath(os.path.join("..",'..',test_script))
        tests_rel_path = '..'
    else:
        assert run_devtest_folder_level_from_root[0] == 0, "invalid folder level"
#        print('+++++++', os.getcwd(), test_script) #DBG

    kwargs['env_vars']["PATH"] = os.path.join(repo_root_dir, 'src')+":"+os.environ["PATH"]
    old_path = os.environ["PATH"]
    os.environ["PATH"] =  os.path.join(repo_root_dir, 'src')+":"+os.environ["PATH"]

    if test_name.endswith(".pl") or test_name.endswith(".xpl"):
        kwargs['env_vars']["srcdir"] = tests_rel_path
        #os.environ["TEST_ID"] = str(bashlist[test_script])
        #os.environ["srcdir"] ='src'
        #os.environ["PATH"] = os.path.join(os.getcwd(), 'src')+":"+os.environ["PATH"]
        #os.system(" ".join(['perl', "-w", "-I./tests", "-MCoreutils", "-MCuSkip", "-M\"CuTmpdir qw("+test_script+")\"", test_script]))
        retcode = system_test_runner('perl',[ "-w", "-I./"+tests_rel_path, "-MCuSkip", "-MCoreutils", "-MCuTmpdir qw("+test_script+")", test_script], None, repo_root_dir, *args, dbg_log_execution_out=False, **kwargs)
    elif test_name.endswith(".sh") or test_name.endswith(""):
#           print('---- DBG', kwargs['env_vars']["TEST_ID"])
            #TESTS_ENVIRONMENT...
        if rootlist[test_name]:
#           print("SUDO ++++++++++")
            retcode = system_test_runner('sudo', ['-E','bash',test_script ], (test_script if check_dev_testname_in_run[0] else None), repo_root_dir, *args, dbg_log_execution_out=False, **kwargs)
            #retcode = system_test_runner('sudo', ['-E', 'make','check', 'TESTS='+test_script, "SUBDIRS=.","VERBOSE=yes"], test_script, repo_root_dir, *args,
            #                             **kwargs)
            os.system("chmod 777 {}/src".format(repo_root_dir))
        else:
            retcode = system_test_runner('bash', [test_script ], (test_script if check_dev_testname_in_run[0] else None), repo_root_dir, *args, dbg_log_execution_out=False, **kwargs)


    os.environ["PATH"] = old_path


    os.chdir(cwd)

    return retcode


def make_build(repo_root_dir, exe_rel_paths, compiler, flags_list, clean,\
                                                                reconfigure):
    cwd = os.getcwd()
    os.chdir(repo_root_dir)
    if flags_list is None:
         flags_list = [ "-DUSING_MUTERIA" ]
		
    else:	
    	flags_list.append("-DUSING_MUTERIA")
    # try:
    tmp_env = os.environ.copy()
    if compiler is not None:
        tmp_env["CC"] = compiler
    if flags_list is not None:
        tmp_env["CFLAGS"] = " ".join(flags_list)
    args_list=['build_project.sh']
    if reconfigure:
        args_list.append('1')
    else:
        args_list.append('0')

    if clean:
        args_list.append('1')
    else:
        args_list.append('0')

    retcode, out, _ = DriversUtils.execute_and_get_retcode_out_err(prog='bash', args_list=args_list,env=tmp_env, merge_err_to_out=True)
    def print_err(out, msg):
        #out = out.splitlines()
        #print(out, msg)
        logging.error(str(out)+msg)
    if retcode != 0:
        print_err(out, "make")
        os.chdir(cwd)
        return GlobalConstants.COMMAND_FAILURE
    os.chdir(cwd)
    return GlobalConstants.COMMAND_SUCCESS

def build_func(*args, **kwargs):
    return make_build(*args, **kwargs)


this_dir = os.path.dirname(os.path.abspath(__file__))
###

PROGRAMMING_LANGUAGE='C'
##REPOSITORY_ROOT_DIR=os.path.join(os.path.dirname(this_dir), 'coreutils')
##OUTPUT_ROOT_DIR=os.path.join(os.path.dirname(this_dir), 'ctrl', "output")
RUN_MODE=SessionMode.EXECUTE_MODE

##TARGET_SOURCE_INTERMEDIATE_CODE_MAP = { 'src/chroot.c':'src/chroot.o'} #'lib/lib.c':'lib/lib.o',
##REPO_EXECUTABLE_RELATIVE_PATHS = ['src/chroot']
CODE_BUILDER_FUNCTION = build_func

CUSTOM_DEV_TEST_RUNNER_FUNCTION = dev_test_runner
CUSTOM_DEV_TEST_PROGRAM_WRAPPER_CLASS = SystemWrapper

##DEVELOPER_TESTS_LIST = initdevtestlist(REPOSITORY_ROOT_DIR)

#print("2  file list   {}".format(devtestlist))
# custom devtest
dev_test = TestcaseToolsConfig(tooltype=TestToolType.USE_ONLY_CODE, toolname='custom_devtests') #, config_id=0)
dev_test.set_one_test_execution_timeout(5)
dev_test.set_test_oracle_test(False) # To enable wrapper splitting of dev tests

# Consider test execution error as a test failure
#dev_test.set_test_execution_error_as_failure(True) 

# klee tests
klee_test = TestcaseToolsConfig(tooltype=TestToolType.USE_ONLY_CODE, toolname='klee', \
                        tool_user_custom=ToolUserCustom(POST_TARGET_CMD_ORDERED_FLAGS_LIST=[('-sym-args', '2', '2', '2')]))
klee_test.set_one_test_execution_timeout(5)

# semu tests TODO: add sym args from json
semu_sym_args = []
semu_test = TestcaseToolsConfig(tooltype=TestToolType.USE_CODE_AND_TESTS, toolname='semu', \
                        tool_user_custom=ToolUserCustom(
                            PRE_TARGET_CMD_ORDERED_FLAGS_LIST=[
                                #('-semu-forkprocessfor-segv-externalcalls', False),

                                #('-semu-disable-statediff-in-testgen',),
                                #('-semu-continue-mindist-out-heuristic',),
                                #('-semu-use-basicblock-for-distance',),
                                ('-semu-forkprocessfor-segv-externalcalls',),
                                #('-semu-testsgen-only-for-critical-diffs',),
                                #('-semu-consider-outenv-for-diffs',),

                                ('-semu-mutant-max-fork', '0'),
                                ('-semu-checknum-before-testgen-for-discarded', '2'),
                                ('-semu-mutant-state-continue-proba', '0.25'),
                                ('-semu-precondition-length', '0'), # start from top
                                #('-semu-max-total-tests-gen', '1000')
                                ('-semu-max-tests-gen-per-mutant', '5'),
                                ('-solver-backend', 'z3'),
                                ('-max-memory', '150000')
                            ],
                            POST_TARGET_CMD_ORDERED_FLAGS_LIST=semu_sym_args)
                            )
semu_test.set_one_test_execution_timeout(5)


# Shadow TCT
import muteria.drivers.testgeneration.tools_by_languages.c.shadow_se\
                                                .shadow_se as shadow_se_module
def build_func_shadow(*args, **kwargs):
    # TODO: shadow flags 
    if args[3] is None:
        flags = []
    else:
        flags = list(args[3])

    klee_change_macro_file = os.path.dirname(shadow_se_module.__file__)
    flags += ['-I'+klee_change_macro_file]
    args = tuple(args[:3]) + (flags,) + tuple(args[4:])
    #return make_build_func(*args, **kwargs)
    return make_build(*args, **kwargs)
#~ def build_func_shadow()

CODE_BUILDER_FUNCTION = build_func_shadow

shadow_se_test = TestcaseToolsConfig(tooltype=TestToolType.USE_CODE_AND_TESTS, toolname='shadow_se', \
                        tool_user_custom=ToolUserCustom(\
                            PATH_TO_TOOL_BINARY_DIR='/home/shadowvm/shadow/klee-change/Release+Asserts/bin/'
                        ))
shadow_se_test.set_one_test_execution_timeout(5)
#~


# test tool list
TESTCASE_TOOLS_CONFIGS = [
        dev_test, 
#        klee_test,
        semu_test,
        shadow_se_test,
]

ENABLED_CRITERIA = [
        TestCriteria.STATEMENT_COVERAGE,
        TestCriteria.BRANCH_COVERAGE,
        TestCriteria.FUNCTION_COVERAGE,
        TestCriteria.WEAK_MUTATION,
        TestCriteria.MUTANT_COVERAGE,
        TestCriteria.STRONG_MUTATION,
]
CRITERIA_WITH_OUTPUT_SUMMARY = [
        TestCriteria.STRONG_MUTATION,
]

gnucov = CriteriaToolsConfig(tooltype=CriteriaToolType.USE_ONLY_CODE, toolname='gcov', config_id=0)
mart = CriteriaToolsConfig(tooltype=CriteriaToolType.USE_ONLY_CODE, toolname='mart', tool_user_custom=ToolUserCustom(POST_TARGET_CMD_ORDERED_FLAGS_LIST=[('-linking-flags', '-lattr') ])  ,config_id=0)
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA = {}
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA[TestCriteria.STATEMENT_COVERAGE] = [gnucov]
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA[TestCriteria.BRANCH_COVERAGE] = [gnucov]
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA[TestCriteria.FUNCTION_COVERAGE] = [gnucov]
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA[TestCriteria.WEAK_MUTATION] = [mart]
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA[TestCriteria.MUTANT_COVERAGE] = [mart]
CRITERIA_TOOLS_CONFIGS_BY_CRITERIA[TestCriteria.STRONG_MUTATION] = [mart]

CRITERIA_EXECUTION_OPTIMIZERS = {
    "STRONG_MUTATION": "SM_OPTIMIZED_BY_WM",
}

COVER_CRITERIA_ELEMENTS_ONCE = False #True
LOG_DEBUG = True
