#include <cairo/cairo.h>
#include <cstdio>
#include <cassert>

#include "types.h"
#include "draw.h"


//helper

bool in_grid(vec x,double limit){
	return x.x>-limit-1e-6&&x.x<limit+1e-6&&x.y>-limit-1e-6&&x.y<limit+1e-6;
}

//main

drawer::drawer(std::string str,double d_sf,double d_mat){
	oput_str=str;
	this->d_sf=d_sf;
	this->d_mat=d_mat;
	surface=cairo_image_surface_create(CAIRO_FORMAT_ARGB32,d_sf,d_sf);
	cr=cairo_create(surface);
}
drawer::drawer(const drawer &d){assert(0);} // 防止意外传递
drawer::~drawer(){
	cairo_surface_destroy(surface);
}
void drawer::flush(){
	cairo_surface_write_to_png(surface,oput_str.data());
}
vec drawer::mat2sf(const vec &a)const{
	return (vec){a.x/d_mat*d_sf/2+d_sf/2,a.y/d_mat*d_sf/2+d_sf/2};
}
void drawer::draw_line(vec x,vec y,col_s c,double width,bool mat)const{
	// 画一条x到y的线段，mat表示输入是否为电路板上的坐标（而不是画布位置）
	std::cerr<<"draw(org): "<<x<<"->"<<y<<std::endl;
	if(mat) x=mat2sf(x),y=mat2sf(y);
	std::cerr<<"draw: "<<x<<"->"<<y<<"(sf: "<<d_sf<<")"<<std::endl;
	assert(in_grid(x,d_sf)&&in_grid(y,d_sf));
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	cairo_set_line_width(cr,width);
	cairo_move_to(cr,x.x,x.y);
	cairo_line_to(cr,y.x,y.y);
	cairo_stroke(cr);
}
void drawer::draw_grid(int lcnt)const{
	// 画出cr的参考坐标系
	double interval=d_sf/lcnt;
	for(int i=0; i<=d_sf+1e-6; i+=interval){
		double p=i;
		draw_line((vec){0.0,p},(vec){d_sf,p},col_white,0.8,0);
		draw_line((vec){p,0.0},(vec){p,d_sf},col_white,0.8,0);
	}
}
void drawer::draw_mdl(mdl m,col_s c,int id)const{
	// 画出模块m
	static char ch[10];
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	m.v[0]=mat2sf(m.v[0]),m.v[1]=mat2sf(m.v[1]);
	double x=m.v[0].x,y=m.v[0].y,dx=m.v[1].x-x,dy=m.v[1].y-y;
	if(dx<0) dx=-dx,x-=dx;
	if(dy<0) dy=-dy,y-=dy;
	cairo_rectangle(cr,x,y,dx,dy);
	cairo_fill(cr);
	draw_line(vec{x,y},vec{x+dx,y},col_green,1,0);
	draw_line(vec{x,y},vec{x,y+dy},col_green,1,0);
	draw_line(vec{x+dx,y},vec{x+dx,y+dy},col_green,1,0);
	draw_line(vec{x,y+dy},vec{x+dx,y+dy},col_green,1,0);
	if(id==-1) return;
	sprintf(ch,"%d",id);
	cairo_set_source_rgba(cr,1,1,1,1);
	cairo_set_font_size(cr,12);
	cairo_move_to(cr,x+dx/2,y+dy/2);
	cairo_show_text(cr,ch);
}
void drawer::draw_cl(const cls_s &cl)const{
	for(edg e:cl) draw_line(e.a,e.b);
}

void dbg_cl(const cls_s &cl){
	// 在dbg.png上画出这个闭合回路cl的形状和位置
	drawer dbg("dbg.png",600,60);
	dbg.draw_grid();
	dbg.draw_cl(cl);
	dbg.flush();
}
