#! /bin/bash

# Run this first (after copying shadow-test from shadow virtual machine (https://srg.doc.ic.ac.uk/projects/shadow/shadow.html))

set -u

error_exit() {
    echo "prepare-error: $1"
    exit 1
}

TOPDIR=$(dirname $(readlink $0))
shadow_data_dir=$TOPDIR/../../shadow-test
repos_topdir=$shadow_data_dir/coreutils

projid=$1

project_repo=$repos_topdir/$projid

#workdir=$TOPDIR/../../

# - Apply patch on analysis
do_analysis_script=$shadow_data_dir/do-analysis.sh
sed -i'' 's/^check_LLVM_support$/#check_LLVM_support/g; s/^check_instrumentation$/#check_instrumentation/g' $do_analysis_script || error_exit "patch do-analysis failed"
grep "^exit 1 # TCT DBG$" || sed -i'' '/^extraHeaders=""$/iexit 1 # TCT DBG/' $do_analysis_script || error_exit "patch do-analysis 2 failed"

# - patch on others 
sed -i'' '2iexit 0' $shadow_data_dir/prepare-native-versions.sh || error_exit "patch failed nat"
sed -i'' '2iexit 0' $shadow_data_dir/prepare-native-and-for-klee.sh || error_exit "patch failed prepare nat and klee"
sed -i'' '2iexit 0' $shadow_data_dir/patch-covering-tcs.sh || error_exit "patch failed patch cov"

# - apply patch on masters
master_configure=$repos_topdir/master.configure.sh
master_make=$repos_topdir/master.make.sh

sed -i'' 's/^  exit 1$/#  exit 1/g' $master_configure
sed -i'' 's/^  exit 1$/#  exit 1/g' $master_make

# - prepare
bash $shadow_data_dir/$projid-patches/build-me.sh || error_exit "build-me failed"

# - get the list of affected files and comment klee_change... headers
affected_summary=$project_repo/klee_changed_src.summary
test -f $affected_summary || error_exit "affected summary missing"
for src in grep "\."
do
    sed -i'' 's|^#include "klee_change_macros.h"$|#include "instrument_wei.h" //#include "klee_change_macros.h"|g' $project_repo/$src || error_exit "failed to sed src; $src"
    sed -i'' 's|^#include "klee_change_functions.h"$|#include "klee_change_functions.h"|g' $project_repo/$src || error_exit "failed to sed 2 src; $src"
    
    # - copy instrument_wei to src folder
    cp $TOPDIR/instrument_wei.h $repos_topdir/$(dirname $src)
done

# - run master.configure.sh and master.make.sh after patching their exit 1
cd $project_repo
$master_configure || error_exit "master configure failed"
$master_make || error_exit "master make failed"

