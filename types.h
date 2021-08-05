#pragma once

#include <iostream>
// x,y: 平面直角坐标系

struct col_s{
	double r,g,b;
};
struct vec{ // 向量
	double x,y;
	static vec get();
	void flip();
};
struct edg{ // 边，a为起点，b为终点，合法区域在这个边向量的左边
	vec a,b;
	int dr();
	void flip();
	bool ispnt();
};
struct mdl{ // module, 记录该模块长方形的四个顶点，保证连续
	vec v[2];
	vec cntr();//返回中心位置
	void flip();
	void set_inf();
	static mdl build(const vec &ctr,const vec &rct);
};

template <typename T> vec operator * (const vec &v,const T x);
double operator *(const vec &a,const vec &b);
vec operator + (const vec &a,const vec &b);
vec operator - (const vec &a,const vec &b);
bool operator < (const vec &a,const vec &b);
bool operator == (const vec &a,const vec &b);
std::ostream& operator << (std::ostream &out,const vec &v);
std::ostream& operator << (std::ostream &out,const edg &e);
