#! /bin/bash

# Run this first (after copying shadow-test from shadow virtual machine (https://srg.doc.ic.ac.uk/projects/shadow/shadow.html))

set -u

error_exit() {
    echo "prepare-error: $1"
    exit 1
}

TOPDIR=$(dirname $(readlink -f $0))
shadow_data_dir=$TOPDIR/../../shadow-test
repos_topdir=$shadow_data_dir/coreutils

projid=$1

project_repo=$repos_topdir/$projid

#workdir=$TOPDIR/../../

# - Apply patch on analysis
do_analysis_script=$shadow_data_dir/do-analysis.sh
sed -i'' 's/^check_LLVM_support$/#check_LLVM_support/g; s/^check_instrumentation$/#check_instrumentation/g' $do_analysis_script || error_exit "patch do-analysis failed"
grep "^exit 0 # TCT DBG$" $do_analysis_script || sed -i'' '/^extraHeaders=""$/iexit 0 # TCT DBG/' $do_analysis_script || error_exit "patch do-analysis 2 failed"

# - patch on others 
head -n2 $shadow_data_dir/prepare-native-versions.sh | grep '^exit 0$' || sed -i'' '2iexit 0' $shadow_data_dir/prepare-native-versions.sh || error_exit "patch failed nat"
head -n2 $shadow_data_dir/prepare-native-and-for-klee.sh | grep '^exit 0$' || sed -i'' '2iexit 0' $shadow_data_dir/prepare-native-and-for-klee.sh || error_exit "patch failed prepare nat and klee"
head -n2 $shadow_data_dir/patch-covering-tcs.sh | grep '^exit 0$' || sed -i'' '2iexit 0' $shadow_data_dir/patch-covering-tcs.sh || error_exit "patch failed patch cov"

# - apply patch on masters
master_configure=$repos_topdir/master.configure.sh
master_make=$repos_topdir/master.make.sh

sed -i'' 's/^  exit 1$/#  exit 1/g' $master_configure
sed '67q;d' $master_configure | grep "\--skip-po" || sed -i'' '67s/.\/bootstrap/.\/bootstrap || .\/bootstrap --skip-po/g' $master_configure || error_exit "sed skip-po failed"
sed -i'' 's/^  exit 1$/#  exit 1/g' $master_make

# - prepare
cd $shadow_data_dir
echo "@prepare: running build-me..."
bash $shadow_data_dir/$projid-patches/build-me.sh || error_exit "build-me failed"
echo "@prepare: running build-me done!"

# - get the list of affected files and comment klee_change... headers
affected_summary=$project_repo/klee_changed_src.summary
test -f $affected_summary || error_exit "affected summary missing"
for src in `grep "\." $affected_summary`
do
    sed -i'' 's|^#include "klee_change_macros.h"$|#include "instrument_wei.h" //#include "klee_change_macros.h"|g' $project_repo/$src || error_exit "failed to sed src: $src"
    sed -i'' 's|^#include "klee_change_functions.h"$|//#include "klee_change_functions.h"|g' $project_repo/$src || error_exit "failed to sed 2 src; $src"
    
    # - copy instrument_wei to src folder
    cp $TOPDIR/instrument_wei.h $project_repo/$(dirname $src) || error_exit "cp failed on instrument_wei to src: $src"
done

# - run master.configure.sh and master.make.sh after patching their exit 1
cd $project_repo
cp $TOPDIR/build_project.sh . || error_exit "Failed to copy build_project"
$master_configure || error_exit "master configure failed"
$master_make || error_exit "master make failed"


###########
# Make sure that make check runs fine for any test script
tmp_test_script="check_make_check.delete.tmp.sh"
echo "#! /bin/bash" > $tmp_test_script
echo "exit 0" >> $tmp_test_script
make check-TEST TESTS="$tmp_test_script" || error_exit "# PREPARE was fine but make check TESTS failed on test with any name"
rm -rf $tmp_test_script
#~
