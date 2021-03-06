#! /bin/bash

# ./run.sh <id> <CPU count>

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

echo '
#! /bin/bash

set -u

#export COREUTILS_TEST_EXPENSIVE=off
#export COREUTILS_TEST_ROOT=on
export COREUTILS_TEST_ROOT=1

TOPDIR=$(dirname $(readlink -f $0))
conf_py=$TOPDIR/ctrl/conf.py
mut_out=$TOPDIR/ctrl/output
collect=$TOPDIR/res

runner=/home/cmtools/run.sh

# rm PB, rebuild
#cd $TOPDIR/../../../shadow-test/coreutils/$(basename $(pwd))/
#make || CFLAGS="-I/usr/local/lib/python3.6/dist-packages/muteria/drivers/testgeneration/tools_by_languages/c/shadow_se -DUSING_MUTERIA" ./build_project.sh 0 1
#cd -

$runner $conf_py $mut_out $collect

exit $?
' > $in_docker_script

chmod +x $in_docker_script

# RUN DOCKER
cd $mountfold || error_exit "cd failed to mountfold"
# Check this for ptrace user uid problem: https://seravo.fi/2019/align-user-ids-inside-and-outside-docker-with-subuser-mapping
# Also: https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine
#sudo docker run -it --rm --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -security-opt apparmor=unconfined --mount type=bind,src=$(pwd),dst=/work --user 1000:1000 --privileged \
sudo docker run -it --rm --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --mount type=bind,src=$(pwd),dst=/work --user 1000:1000 --privileged \
									 --cpus=${CPU} maweimarvin/cm bash -c "cd /work/executions/workspace/$id && bash ./${tmpcmd}"

# Copy info about changed lines
test -f $workspace_dir/res/klee_changed_src.summary || cp $mountfold/shadow-test/coreutils/$id/klee_changed_src.summary $workspace_dir/res/klee_changed_src.summary \
															|| error_exit "failed to copy change lines summary"

rm $in_docker_script || error_exit "failed to remove in_docker_script"

echo "ALL DONE!"
