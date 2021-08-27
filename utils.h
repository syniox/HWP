// 与内容无关的基础代码 + 一些全局变量

#pragma once

#include <cstdio>
#include <algorithm>
#include <vector>

#include "types.h"

extern const double eps;
extern const double inf;

extern std::vector<std::vector<double>> m_lim; // 模块间的限制，太耗内存？
extern std::vector<std::vector<edg>> e_lim; // 模块与边的限制

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
template <typename T> void inc_swp(T &x,T &y){ //将x和y增续排列
	if(x<y) std::swap(x,y);
}
template <typename T> void flip_vec(std::vector <T> &vt){
	for(T &x:vt) x.flip();
}
template <typename T1,typename T2>
void replace_with(T1 &str,std::vector<T2> repl,T2 p){
	for(T2 c:repl) std::replace(str.begin(),str.end(),c,p);
}

std::string get_line(std::vector<char> repl); //读入一行非空串并替换部分字符为空格
void sanitize_vec(cls_s &cl);
// O(e(cl)*e(cl)) 删除长度为0的边，合并相邻且方向相反的边
double calc_res(std::vector<mdl> m_res,std::vector<mdl> m_org,std::vector<int> ref);
// O(n)  计算一组方案的连线长度和
