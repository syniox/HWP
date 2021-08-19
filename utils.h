// 与内容无关的基础代码

#pragma once

#include <cstdio>
#include <algorithm>
#include <vector>

extern double eps;
extern double inf;

int rnd(const int l,const int r);

template <typename T> void apn(T &x,const T y){
	x=x<=y?x:y;
}
template <typename T> void apx(T &x,const T y){
	x=x>=y?x:y;
}
template <typename T> const T cabs(const T &x){
	return x<0?-x:x;
}

template <typename T1,typename T2>
void replace_with(T1 &str,std::vector<T2> repl,T2 p){
	for(T2 c:repl) std::replace(str.begin(),str.end(),c,p);
}
