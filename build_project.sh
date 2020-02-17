#!/bin/bash -l

#TCT
set -u

error_exit(){
	echo "Build Error: $1"
	exit 1
}

fix_tail() {
	tail_fixing_on=0
	if [ "${CC:-}" = 'wllvm' ]; then
		if [ "$(basename $(readlink -f ${LLVM_COMPILER_PATH:-}/${LLVM_COMPILER}))" = 'llvm-gcc' ]; then
			tail_fixing_on=1
                        export CFLAGS="${CFLAGS:-} -std=c99"
			grep "^#define __MFI_REAL_GNUC__" src/tail.c > /dev/null \
					|| sed -i '1s/^/#define __MFI_REAL_GNUC__  __GNUC__\n#undef __GNUC__\n#define __GNUC__   1\n/'  src/tail.c \
					|| error_exit "Failed to replace src/tail.c's included 'select.h'"
		fi
	fi
}

revert_tail_fix(){
	if [ $tail_fixing_on -eq 1 ]; then
		if grep "^#define __MFI_REAL_GNUC__" src/tail.c > /dev/null
		then
			sed -i '1,3d' src/tail.c || error_exit "Failed to revert src/tail.c's included 'select.h'"
		fi
	fi
}
disable_test_framework_check()
{
    sed -i 's/^fail_[[:blank:]]*([[:blank:]]*).*$/fail_ () { warn_ "$ME_: failed test: $@";}/g' $1
    sed -i 's/^skip_[[:blank:]]*([[:blank:]]*).*$/skip_ () { warn_ "$ME_: skipped test: $@";}/g' $1
    sed -i 's/^fatal_[[:blank:]]*([[:blank:]]*).*$/fatal_ () { warn_ "$ME_: hard error: $@"; }/g' $1
    sed -i 's/^framework_failure_[[:blank:]]*([[:blank:]]*).*$/framework_failure_ () { warn_ "$ME_: set-up failure: $@"; }/g' $1
}


#~


#compile projects
#./bootstrap
if [ $1 -eq 1 ]; then
	./configure  --disable-gcc-warnings  --disable-nls --disable-selinux LDFLAGS="-Wl,--no-as-needed,-ldl"
	if [ "$(basename $(pwd))" = "ar-4-1" ] 
	then
		test -f lib/xsize.h || cp intl/xsize.h lib/xsize.h
		echo "all: ;" > doc/Makefile
	elif [ "$(basename $(pwd))" = "cr-15" -o "$(basename $(pwd))" = "cr-1" -o "$(basename $(pwd))" = "ar-2" ] 
	then
		echo "all: ;" > doc/Makefile
	fi
	#echo "all: ;" > doc/Makefile
	#echo "all: ;" > po/Makefile
#	sed -i '/found non-public submodule commit/,+1d' ./gnulib/top/maint.mk
	sed -i 's/^_GL_WARN_ON_USE (gets, "gets is a security hole - use fgets instead");/\/\/_GL_WARN_ON_USE (gets, "gets is a security hole - use fgets instead");/g' ./lib/stdio.h	
fi

if [ $2 -eq 1 ]; then
	make clean
fi

fix_tail
make -j4 
disable_test_framework_check tests/init.sh
#TCT
revert_tail_fix
#~


