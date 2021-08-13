#include <cairo/cairo.h>

#include "types.h"

// gry: 矩形对应的位置	red: 合法区域边框
// grn: 边界 			blu: 无法放入的矩形
// wht: 坐标系			cyan: 被放入的矩形
static const col_s col_gry=col_s{0.4,0.4,0.4},col_red=col_s{0.8,0.1,0.3};
static const col_s col_grn=col_s{0.1,0.7,0.1},col_blu=col_s{0.2,0.4,1.0};
static const col_s col_wht=col_s{0.8,0.8,0.8},col_cyan=col_s{0.1,0.8,0.8};

struct drawer{
	cairo_surface_t *surface;
	cairo_t *cr;
	std::string oput_str; // 输出文件名
	double d; // 画布边长
	drawer(std::string str,double l=500);
	drawer(const drawer &d); // assert(0),防止意外传递
	~drawer();
	void flush();
	void draw_line(vec x,vec y,col_s c=col_red,double width=1); // 画一条x到y的线段
	void draw_grid(int interval=10); // 画出cr的参考坐标系
	void draw_mdl(mdl m,col_s c=col_gry,int id=-1); // 画出模块m
	void draw_cl(const cls_s cl); //画出该闭合回路（debug用）
};

void dbg_cl(const cls_s &cl); // 在dbg.png上画出这个闭合回路cl的形状和位置