#include <cairo/cairo.h>
#include <cstdio>
#include <cassert>

#include "types.h"
#include "draw.h"

drawer::drawer(std::string str,double l){
	oput_str=str;
	d=l;
	surface=cairo_image_surface_create(CAIRO_FORMAT_ARGB32,d,d);
	cr=cairo_create(surface);
}
drawer::drawer(const drawer &d){assert(0);} // 防止意外传递
drawer::~drawer(){
	cairo_surface_destroy(surface);
}
void drawer::flush(){
	cairo_surface_write_to_png(surface,oput_str.data());
}
void drawer::draw_line(vec x,vec y,col_s c,double width){
	// 画一条x到y的线段
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	cairo_set_line_width(cr,width);
	cairo_move_to(cr,x.x*10,x.y*10);
	cairo_line_to(cr,y.x*10,y.y*10);
	cairo_stroke(cr);
}
void drawer::draw_grid(int interval){
	// 画出cr的参考坐标系
	for(int i=0; i<=d; i+=interval){
		double p=i;
		draw_line((vec){0.0,p},(vec){d,p},col_white,0.8);
		draw_line((vec){p,0.0},(vec){p,d},col_white,0.8);
	}
}
void drawer::draw_mdl(mdl m,col_s c,int id){
	// 画出模块m
	static char ch[10];
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	double x=m.v[0].x,y=m.v[0].y,dx=m.v[1].x-x,dy=m.v[1].y-y;
	if(dx<0) dx=-dx,x-=dx;
	if(dy<0) dy=-dy,y-=dy;
	cairo_rectangle(cr,x*10,y*10,dx*10,dy*10);
	cairo_fill(cr);
	draw_line(vec{x,y},vec{x+dx,y},col_green);
	draw_line(vec{x,y},vec{x,y+dy},col_green);
	draw_line(vec{x+dx,y},vec{x+dx,y+dy},col_green);
	draw_line(vec{x,y+dy},vec{x+dx,y+dy},col_green);
	if(id==-1) return;
	sprintf(ch,"%d",id);
	cairo_set_source_rgba(cr,1,1,1,1);
	cairo_set_font_size(cr,12);
	cairo_move_to(cr,(x+dx/2)*10,(y+dy/2)*10);
	cairo_show_text(cr,ch);
}
void drawer::draw_cl(const cls_s &cl){
	for(edg e:cl) draw_line(e.a,e.b);
}

void dbg_cl(const cls_s &cl){
	// 在dbg.png上画出这个闭合回路cl的形状和位置
	drawer dbg("dbg.png");
	dbg.draw_cl(cl);
	dbg.draw_grid();
	dbg.flush();
}
