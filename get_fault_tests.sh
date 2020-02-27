#! /bin/bash

# ./get_fault_tests.sh <id> <CPU count>

set -u

error_exit()
{
    echo "ERROR: $1"
    exit 1
}

# id
# cpu
# root
# expensive

[ $# = 2 ] || error_exit "expect 2 argument (id and number of cpu). passed $# agrs"

id=$1
CPU=$2

####

TOPDIR=$(dirname $(readlink -f $0))
mountfold=$(readlink -f $TOPDIR/../..)
workspace_dir=$(readlink -f $TOPDIR/../workspace/$id)
tmpcmd=in_docker_run.sh
in_docker_script=$workspace_dir/$tmpcmd
test -d $workspace_dir || mkdir $workspace_dir || error_exit "failed to create workspace dir"

echo '
#! /bin/bash

set -u

error_exit()
{
    echo "ERROR: $1"
    exit 1
}

pip install -U muteria

#export COREUTILS_TEST_EXPENSIVE=off
#export COREUTILS_TEST_ROOT=on
export COREUTILS_TEST_ROOT=1

TOPDIR=$(dirname $(readlink -f $0))
conf_py=$TOPDIR/ctrl/conf.py
collect=$TOPDIR/res
pass_fail_matrix=$collect/post/RESULTS_DATA/matrices/PASSFAIL.csv

p_id=$(basename $TOPDIR)
exe_dir=/work/executions/cm_corebench_scripts/bug_fixing_exes/$p_id
exe_file=$(ls $exe_dir/old)

if test -f $exe_dir/alias;
then
	alias_=$(cat $exe_dir/alias)
	conf_py=$TOPDIR/../$alias_/ctrl/conf.py
        test -d $collect || cp -r $TOPDIR/../$alias_/res $collect || error_exit "copy re failed"
fi

fail_test_execution=$collect/fail_test_checking
test -d $fail_test_execution || mkdir $fail_test_execution

bug_finding_tests_list=$fail_test_execution/fault_reveling_tests.txt

# temporary
test_list_file=$TOPDIR/test_list.tmp

cat $pass_fail_matrix | head -n1 | tr " " "\n" | sed 1d > $test_list_file || error_exit "failed to get testlist"

custom_exe=$exe_dir/old/$exe_file
echo "tests
$fail_test_execution/old
{\"src/$(basename $custom_exe)\": \"$custom_exe\"}
$test_list_file" | muteria --config $conf_py --lang c customexec || error_exit "run failed"

custom_exe=$exe_dir/new/$exe_file
echo "tests
$fail_test_execution/new
{\"src/$(basename $custom_exe)\": \"$custom_exe\"}
$test_list_file" | muteria --config $conf_py --lang c customexec || error_exit "run failed"

rm -f $test_list_file

echo "import muteria.common.fs as common_matrices" > $test_list_file
echo "_, old_o = common_matrices.OutputLogData(\"$fail_test_execution/old/program_output.json\").get_zip_objective_and_data()" >> $test_list_file
echo "_, new_o = common_matrices.OutputLogData(\"$fail_test_execution/new/program_output.json\").get_zip_objective_and_data()" >> $test_list_file
echo "assert len(old_o) == len(new_o)" >> $test_list_file
echo "diff_tests = []" >> $test_list_file
echo "for tc in old_o:" >> $test_list_file
echo "    eq = common_matrices.OutputLogData.outlogdata_equiv(old_o[tc], new_o[tc]):" >> $test_list_file
echo "    assert eq i not None, \"PB\"" >> $test_list_file
echo "    if not eq:" >> $test_list_file
echo "        diff_tests.append(tc)" >> $test_list_file
echo "with open(\"$bug_finding_tests_list\", \"w\") as f:" >> $test_list_file
echo "    for tc in diff_tests:" >> $test_list_file
echo "        f.write(tc+\"\n\")" >> $test_list_file
echo "print(\"# list printed\")" >> $test_list_file

python $test_list_file || error_exit "python failed"

rm -f $test_list_file
' > $in_docker_script

chmod +x $in_docker_script

# RUN DOCKER
cd $mountfold || error_exit "cd failed to mountfold"
# Check this for ptrace user uid problem: https://seravo.fi/2019/align-user-ids-inside-and-outside-docker-with-subuser-mapping
# Also: https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine
#sudo docker run -it --rm --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -security-opt apparmor=unconfined --mount type=bind,src=$(pwd),dst=/work --user 1000:1000 --privileged \
sudo docker run -it --rm --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --mount type=bind,src=$(pwd),dst=/work --user 1000:1000 --privileged \
							 --cpus=${CPU} maweimarvin/cm bash -c "cd /work/executions/workspace/$id && bash ./${tmpcmd}"

rm $in_docker_script || error_exit "failed to remove in_docker_script"

echo "ALL DONE!"
