#ifdef USING_MUTERIA
#include "klee_change_macros.h"
#else
#ifndef INSTRUMENT_H_WEI
#define INSTRUMENT_H_WEI
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
unsigned long long  klee_change(unsigned long long    oldV,   unsigned long long   newV){
	char*   version = getenv("CI_VERSION");
	if(NULL == version || strlen(version) != 1 || ! isdigit(version[0])){
		//printf("safty check CI_VERSION=%s\n", version);
		return newV;
	}
	int v = atoi(version); // Using atoi()
        //printf("%s", "klee change\n");
	if(v==1){
                //printf("%s", "return new value\n");
		//printf("new version CI_VERSION=%s\n", version);
		return newV;
	}else if( v==0 ){
                //printf("%s", "return old value\n");
		//printf("old version CI_VERSION=%s\n", version);
		return oldV;
	}else{
		//printf("default CI_VERSION=%s\n", version);
		return newV;
	}
}

int klee_get_true( ){
        //printf("%s", "klee_get_true\n");
	return 1;
}

int klee_get_false( ){
        //printf("%s", "klee_get_false\n");
	return 0;
}

#endif
#endif

