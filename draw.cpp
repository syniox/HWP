#include "draw.h"

void draw_line(cairo_t *cr,vec a,vec b){
	cairo_set_source_rgba(cr,0.8,0.1,0.3,1.0);
	cairo_set_line_width(cr,1);
	cairo_move_to(cr,a.x*10,a.y*10);
	cairo_line_to(cr,b.x*10,b.y*10);
	cairo_stroke(cr);
}

void draw_mdl(cairo_t *cr,mdl m,double r,double g,double b){
	cairo_set_source_rgba(cr,r,g,b,1.0);
	double x=m.v[0].x,y=m.v[0].y,dx=m.v[2].x-x,dy=m.v[2].y-y;
	if(dx<0) dx=-dx,x-=dx;
	if(dy<0) dy=-dy,y-=dy;
	cairo_rectangle(cr,x*10+0.5,y*10+0.5,dx*10-1,dy*10-1);
	cairo_fill(cr);
}
