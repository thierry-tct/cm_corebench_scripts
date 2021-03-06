#! /bin/bash
# Run the following command from within the foder that contains the tar files
# WITH_PREDICTION=cross_validation.json /media/disk2/CONTINUOUS_MUTATION/executions/cm_corebench_scripts/analysis_scripts/run.sh . docker

set -u

error_exit()
{
    echo "@run.sh-ERROR: $1"
    exit 1
}

TOPDIR=$(dirname $(readlink -f $0))

[ $# -eq 1 -o $# -eq 2 ] || error_exit "Invalid number of arguments"

in_top_dir=$(readlink -f $1)
test -d $in_top_dir || error_exit "specified in_top_di missing"
out_top_dir=$in_top_dir/ANALYSIS_OUT
test -d $out_top_dir || mkdir $out_top_dir || error_exit "failed to mae out_top_dir"

in_docker=0

if [ $# -eq 2 ]; then
    [ "$2" = 'docker' ] || error_exit "secon parameter must be docker"
    in_docker=1
fi

pred_file=""
if [ "${WITH_PREDICTION:-}" != "" ]; then
    echo "DBG: With prediction is $WITH_PREDICTION ...."
    pred_file=/work_in/$WITH_PREDICTION
    test -f $in_top_dir/$WITH_PREDICTION || error_exit "pred file $WITH_PREDICTION not existing"
fi

if [ $in_docker -eq 1 ]
then
    docker_image_name="maweimarvin/cm"
    sudo docker run -it --rm --cap-add=SYS_PTRACE --security-opt seccomp=unconfined \
                                         --mount type=bind,src=$in_top_dir,dst=/work_in \
                                         --mount type=bind,src=$out_top_dir,dst=/work_out \
                                         --mount type=bind,src=$TOPDIR,dst=/work_script \
                                         --user 1000:1000 --privileged $docker_image_name bash -c "pip install seaborn; python /work_script/main.py /work_in /work_out $pred_file"    
else
    python3 $TOPDIR/main.py $in_top_dir $out_top_dir $pred_file || error_exit "Execution failed"
fi

echo "@run.sh: DONE!"
