from __future__ import print_function

from common_conf import *
import json
import os
from pathlib import Path
workdir = Path(__file__).parent.absolute()
coreutils = os.path.join(workdir.parent.parent.parent.parent, 'shadow-test', 'coreutils')
confprg = os.path.join(workdir, "init.conf")
#print("{}\n".format(confprg))
if not os.path.isfile(confprg):
  print("Lack init.conf file in ctrl folder, please check")
  exit()
f = open(confprg, "r")
config = json.load(f)
REPOSITORY_ROOT_DIR=os.path.join(coreutils, config['id'])
OUTPUT_ROOT_DIR=os.path.join(os.path.dirname(this_dir), 'ctrl', "output")

TARGET_SOURCE_INTERMEDIATE_CODE_MAP = config["TARGET_SOURCE_INTERMEDIATE_CODE_MAP"] #'lib/lib.c':'lib/lib.o',
REPO_EXECUTABLE_RELATIVE_PATHS = config["REPO_EXECUTABLE_RELATIVE_PATHS"]

DEVELOPER_TESTS_LIST = initdevtestlist_bench(os.path.join(REPOSITORY_ROOT_DIR, 'tests'))

# set semu symagrs, read from 'klee_symb_args.json'
with open(os.path.join(workdir, "klee_symb_args.json")) as f:
      coreutils_sym_args = json.load(f)
#      print("===== {}".format(coreutils_sym_args))
      coreutils_tool = os.path.basename(REPO_EXECUTABLE_RELATIVE_PATHS[0])
      if coreutils_tool not in coreutils_sym_args:
            coreutils_tool = ""
      semu_sym_args.extend(coreutils_sym_args[coreutils_tool])


check_dev_testname_in_run = config['check_dev_testname_in_run']

# TCT
#bypass_make_check_tests.compute_make_check_tests_env_vars(REPOSITORY_ROOT_DIR)
diff_make_build_dir = None
if config['make_check_test_rel_path']:
    diff_make_build_dir = os.path.join(REPOSITORY_ROOT_DIR, config['make_check_test_rel_path'])
bypass_make_check_tests.compute_make_check_tests_env_vars(REPOSITORY_ROOT_DIR, diff_make_build_dir)
#~

HASH_OUTLOG=False
