#include <cairo/cairo.h>
#include <cstdio>
#include <cassert>
#include <cmath>

#include "draw.h"
#include "types.h"
#include "utils.h"


//helper

bool in_grid(vec x,double limit){
	return x.x>-limit-eps&&x.x<limit+eps&&x.y>-limit-eps&&x.y<limit+eps;
}

//main

drawer::drawer(std::string str,double d_sf){
	x_low=y_low=inf,x_up=y_up=-inf;
	oput_str=str;
	this->d_sf=d_sf;
	surface=cairo_image_surface_create(CAIRO_FORMAT_ARGB32,d_sf,d_sf);
	cr=cairo_create(surface);
}
drawer::drawer(const drawer &d){
	std::cerr<<"[warning] copying drawer "<<d.oput_str<<'.'<<std::endl;
	x_low=d.x_low,x_up=d.x_up;
	y_low=d.y_low,y_up=d.y_up;
}
drawer::~drawer(){
	cairo_surface_destroy(surface);
}
void drawer::flush(){
	cairo_surface_write_to_png(surface,oput_str.data());
}
void drawer::upd(const vec &a){
	//std::cerr<<"upd: "<<a<<std::endl;
	apn(x_low,a.x),apx(x_up,a.x);
	apn(y_low,a.y),apx(y_up,a.y);
}
void drawer::zoom_out(const double scale){
	x_low-=scale,x_up+=scale;
	y_low-=scale,y_up+=scale;
}
vec drawer::mat2sf(const vec &a)const{
	double dmax=std::max((x_up-x_low),(y_up-y_low));
	double x_pec=(a.x-x_low)/dmax;
	double y_pec=1-(a.y-y_low)/dmax;
	return (vec){x_pec*d_sf,y_pec*d_sf};
}
vec drawer::sf2mat(const vec &a)const{
	double dmax=std::max((x_up-x_low),(y_up-y_low));
	double x_pec=a.x/d_sf,y_pec=(d_sf-a.y)/d_sf;
	return (vec){x_pec*dmax+x_low,y_pec*dmax+y_low};
}
void drawer::draw_line(vec x,vec y,col_s c,double width,bool mat,double rad)const{
	// 画一条x到y的线段，mat表示输入是否为电路板上的坐标（而不是画布位置）
	// 在向量指向的那一段做一个正方形
	if(mat) x=mat2sf(x),y=mat2sf(y);
	//std::cerr<<"draw: "<<x<<"->"<<y<<"(sf: "<<d_sf<<")"<<std::endl;
	assert(in_grid(x,d_sf)&&in_grid(y,d_sf));
	cairo_set_source_rgba(cr,c.r,c.g,c.b,0.6);
	cairo_set_line_width(cr,width);
	cairo_move_to(cr,x.x,x.y);
	cairo_line_to(cr,y.x,y.y);
	cairo_stroke(cr);
	if(cabs(rad)<eps) return;
	vec vt=(y-x).norm(rad);
	vec vl{-vt.y,vt.x},vr{vt.y,-vt.x};
	vec start=y-vt+vl,end=y-vt+vr;
	cairo_move_to(cr,start.x,start.y);
	cairo_curve_to(cr,start.x,start.y,y.x,y.y,end.x,end.y);
	cairo_stroke(cr);
}
void drawer::draw_grid(int lcnt)const{
	// 画出cr的参考坐标系
	double interval=d_sf/lcnt;
	for(int i=0; i<=d_sf+eps; i+=interval){
		double p=i;
		draw_line((vec){0.0,p},(vec){d_sf,p},col_white,0.8,0);
		draw_line((vec){p,0.0},(vec){p,d_sf},col_white,0.8,0);
	}
}
void drawer::draw_mdl(mdl m,col_s c,std::string id)const{
	// 画出模块m
	cairo_set_source_rgba(cr,c.r,c.g,c.b,0.6);
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
	if(!id.length()) return;
	cairo_set_source_rgba(cr,1,1,1,0.6);
	cairo_set_font_size(cr,12);
	cairo_move_to(cr,x+dx/2,y+dy/2);
	cairo_show_text(cr,id.c_str());
}
void drawer::draw_cl(const cls_s &cl)const{
	for(edg e:cl) draw_line(e.a,e.b,col_red,1,1,6);
}

void dbg_cl(const cls_s &cl,std::initializer_list<mdl> mds){
	// 在dbg.png上画出这个闭合回路cl的形状和位置 包括一些模块
	drawer dbg("dbg.png");
	for(edg e:cl){
		dbg.upd(e.a),dbg.upd(e.b);
	}
	for(mdl m:mds){
		dbg.upd(m.v[0]),dbg.upd(m.v[1]);
	}
	dbg.zoom_out();
	dbg.draw_grid();
	dbg.draw_cl(cl);
	for(mdl m:mds) dbg.draw_mdl(m);
	dbg.flush();
}
