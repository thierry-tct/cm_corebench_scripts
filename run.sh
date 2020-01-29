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
export COREUTILS_TEST_ROOT=on

TOPDIR=$(dirname $(readlink -f $0))
conf_py=$TOPDIR/ctrl/conf.py
mut_out=$TOPDIR/ctrl/output
collect=$TOPDIR/res

runner=/home/cmtools/run.sh

$runner $conf_py $mut_out $collect

exit $?
' > $in_docker_script

chmod +x $in_docker_script

# RUN DOCKER
cd $mountfold || error_exit "cd failed to mountfold"
sudo docker run -it --rm --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --mount type=bind,src=$(pwd),dst=/work --user 1000:1000 --privileged \
									 --cpus=${CPU} maweimarvin/cm bash -c "cd /work/executions/workspace/$id && bash ./${tmpcmd}"

rm $in_docker_script || error_exit "failed to remove in_docker_script"

echo "ALL DONE!"
