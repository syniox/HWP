#include "types.h"
#include "utils.h"

#include <iostream>
#include <cassert>
#include <cmath>
// x,y: 平面直角坐标系

vec vec::get(){
	double x,y;
	std::cin>>x>>y;
	return (vec){x,y};
}
double vec::len2()const{
	return x*x+y*y;
}
void vec::flip(){
	std::swap(x,y);
}
bool vec::ispnt()const{
	return cabs(x)+cabs(y)<eps*2;
}
vec vec::norm(const double l){
	return *this=*this*(l/sqrt(x*x+y*y));
}
double operator * (const vec &a,const vec &b){
	return a.x*b.y-a.y*b.x;
}
vec operator + (const vec &a,const vec &b){
	return (vec){a.x+b.x,a.y+b.y};
}
vec operator - (const vec &a,const vec &b){
	return (vec){a.x-b.x,a.y-b.y};
}
vec vec::operator += (const vec &b){
	return *this=*this+b;
}
vec vec::operator -= (const vec &b){
	return *this=*this-b;
}
bool operator < (const vec &a,const vec &b){
	return a.x==b.x?a.y<b.y:a.x<b.x;
}
bool operator == (const vec &a,const vec &b){
	return cabs(a.x-b.x)<eps&&cabs(a.y-b.y)<eps;
	//return a.x==b.x&&a.y==b.y;
}
std::ostream& operator << (std::ostream &out,const vec &v){
	out<<'('<<v.x<<' '<<v.y<<')';
	return out;
}
std::ostream& operator << (std::ostream &out,const edg &e){
	out<<e.a<<"->"<<e.b;
	return out;
}

int edg::dr()const{ // 返回向量方向，右0上1左2下3
	assert(cabs(a.x-b.x)<eps||cabs(a.y-b.y)<eps);
	assert(!(a==b));
	if(cabs(a.y-b.y)<eps) return (a.x>b.x)<<1;
	else return (a.y>b.y)<<1|1;
}
void edg::flip(){
	a.flip(),b.flip();
}
bool edg::ispnt()const{
	return a==b;
}
edg edg::operator +(const vec &v){
	return (edg){a+v,b+v};
}

vec mdl::cntr()const{
	return (v[0]+v[1])*0.5;
}
void mdl::flip(){
	for(int i=0; i<2; ++i){
		std::swap(v[i].x,v[i].y);
	}
}
void mdl::set_inf(){
	for(int i=0; i<2; ++i){
		v[i].x=v[i].y=-inf;
	}
}
mdl mdl::build(const vec &ctr,const vec &rct){
	double rad_x=rct.x*0.5,rad_y=rct.y*0.5;
	mdl m;
	m.v[0]=(vec){ctr.x-rad_x,ctr.y-rad_y};
	m.v[1]=(vec){ctr.x+rad_x,ctr.y+rad_y};
	return m;
}
