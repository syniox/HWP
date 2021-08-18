#pragma once

#include <cairo/cairo.h>

#include "types.h"

// grey: 矩形对应的位置	red: 合法区域边框
// green: 边界 			blue: 无法放入的矩形
// white: 坐标系		cyan: 被放入的矩形
static const col_s col_grey=col_s{0.4,0.4,0.4},col_red=col_s{0.8,0.1,0.3};
static const col_s col_green=col_s{0.1,0.7,0.1},col_blue=col_s{0.2,0.4,1.0};
static const col_s col_white=col_s{0.8,0.8,0.8},col_cyan=col_s{0.1,0.8,0.8};

struct drawer{
	cairo_surface_t *surface;
	cairo_t *cr;
	std::string oput_str; // 输出文件名
	double d_sf; // 画布边长
	double x_low,x_up,y_low,y_up; // x和y方向上的最大最小值
	drawer(std::string str,double d_sf=2000);
	drawer(const drawer &d); // assert(0),防止意外传递
	~drawer();
	void flush();
	void upd(const vec &a); // 更新最大最小值
	void zoom_out(const double edg=10);
	vec mat2sf(const vec &a)const; // 坐标系坐标转为图上坐标
	void draw_line(vec x,vec y,col_s c=col_red,double width=1,bool mat=1,double rect=0)const; // 画一条x到y的线段 是否是原图坐标
	void draw_grid(int lcnt=10)const; // 画出cr的参考坐标系
	void draw_mdl(mdl m,col_s c=col_grey,int id=-1)const; // 画出模块m
	void draw_cl(const cls_s &cl)const; //画出该闭合回路（debug用）
};

void dbg_cl(const cls_s &cl); // 在dbg.png上画出这个闭合回路cl的形状和位置
