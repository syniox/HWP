// 与内容无关的基础代码

#include "utils.h"
#include <cstdio>

double eps=1e-6;
double inf=1e12;

int rnd(const int l,const int r){
	return rand()%(r-l+1)+l;
}
