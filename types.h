#pragma once

#include <iostream>
#include <vector>
// x,y: 平面直角坐标系

struct col_s{
	double r,g,b;
};
struct vec{ // 向量
	double x,y;
	static vec get(); // 输入向量
	void flip(); // 交换x,y
};
struct edg{ // 边，a为起点，b为终点，合法区域在这个边向量的左边
	vec a,b;
	int dr();
	void flip(); // 分别对a,b交换其x,y（get_great_pos的坐标系变换）
	bool ispnt(); // 这条边是不是退化成了一个点(a,b相同)
};
struct mdl{ // module, 记录该模块长方形的四个顶点，保证连续
	vec v[2];
	vec cntr(); // 返回模块中心
	void flip(); // 对每个顶点进行x,y坐标互换
	void set_inf(); // 初始化成无穷大
	static mdl build(const vec &ctr,const vec &rct); // 根据中心坐标ctr和模块长宽rct构建模块
};

using cls_s=std::vector<edg>;

template <typename T> vec operator * (const vec &v,const T x);
double operator *(const vec &a,const vec &b);
vec operator + (const vec &a,const vec &b);
vec operator - (const vec &a,const vec &b);
bool operator < (const vec &a,const vec &b);
bool operator == (const vec &a,const vec &b);
std::ostream& operator << (std::ostream &out,const vec &v);
std::ostream& operator << (std::ostream &out,const edg &e);
