#include <cairo/cairo.h>
#include <cstdio>
#include <cassert>

#include "types.h"
#include "draw.h"


//helper

bool in_grid(vec x,double limit){
	return x.x>-limit&&x.x<limit&&x.y>-limit&&x.y<limit;
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
	return (vec){a.x/d_mat*d_sf/2+d_sf/2,a.y/d_mat*d_sf+d_sf/2};
}
void drawer::draw_line(vec x,vec y,col_s c,double width)const{
	// 画一条x到y的线段
	assert(in_grid(x,d_mat)&&in_grid(y,d_mat));
	std::cerr<<"draw(org): "<<x<<"->"<<y<<std::endl;
	x=mat2sf(x),y=mat2sf(y);
	std::cerr<<"draw: "<<x<<"->"<<y<<std::endl;
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	cairo_set_line_width(cr,width);
	cairo_move_to(cr,x.x*10,x.y*10);
	cairo_line_to(cr,y.x*10,y.y*10);
	cairo_stroke(cr);
}
void drawer::draw_grid(int interval)const{ // 使用了sf坐标而非mat坐标 此处有Bug
	// 画出cr的参考坐标系
	for(int i=0; i<=d_sf; i+=interval){
		double p=i;
		draw_line((vec){0.0,p},(vec){d_sf,p},col_white,0.8);
		draw_line((vec){p,0.0},(vec){p,d_sf},col_white,0.8);
	}
}
void drawer::draw_mdl(mdl m,col_s c,int id)const{
	// 画出模块m
	static char ch[10];
	m.v[0]=mat2sf(m.v[0]),m.v[1]=mat2sf(m.v[1]);
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
void drawer::draw_cl(const cls_s &cl)const{
	for(edg e:cl) draw_line(mat2sf(e.a),mat2sf(e.b));
}

void dbg_cl(const cls_s &cl){
	// 在dbg.png上画出这个闭合回路cl的形状和位置
	drawer dbg("dbg.png");
	dbg.draw_grid();
	dbg.draw_cl(cl);
	dbg.flush();
}
