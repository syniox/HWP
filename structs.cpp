#include <cassert>
#include <algorithm>
#include "global.h"
#include "structs.h"

vec vec::get(){
	double x,y;
	std::cin>>x>>y;
	return (vec){x,y};
}
void vec::flip(){
	std::swap(x,y);
}
template <typename T>
vec operator * (const vec &v,const T x){
	return (vec){v.x*x,v.y*x};
}
double operator *(const vec &a,const vec &b){
	return a.x*b.y-a.y*b.x;
}
vec operator + (const vec &a,const vec &b){
	return (vec){a.x+b.x,a.y+b.y};
}
vec operator - (const vec &a,const vec &b){
	return (vec){a.x-b.x,a.y-b.y};
}
bool operator < (const vec &a,const vec &b){
	return a.x==b.x?a.y<b.y:a.x<b.x;
}
bool operator == (const vec &a,const vec &b){
	return a.x==b.x&&a.y==b.y;
}
std::ostream& operator << (std::ostream &out,const vec &v){
	out<<'('<<v.x<<' '<<v.y<<')';
	return out;
}
std::ostream& operator << (std::ostream &out,const edg &e){
	out<<e.a<<"->"<<e.b;
	return out;
}

bool crossed(const edg &a,const edg &b){ // 线段是否相交，重合情况不考虑
	return ((b.a-a.a)*(a.b-a.a))*((b.b-a.a)*(a.b-a.a))<-eps;
}
int edg::dr(){ // 返回向量方向，右0上1左2下3
	assert(a.x==b.x||a.y==b.y);
	assert(!(a==b));
	if(a.y==b.y) return (a.x>b.x)<<1;
	else return (a.x>b.x)<<1|1;
}
void edg::flip(){
	a.flip(),b.flip();
}
bool edg::ispnt(){
	return a==b;
}

vec mdl::cntr(){
	return (v[0]+v[1]+v[2]+v[3])*0.25;
}
void mdl::flip(){
	for(int i=0; i<4; ++i){
		std::swap(v[i].x,v[i].y);
	}
}
void mdl::set_inf(){
	for(int i=0; i<4; ++i){
		v[i].x=v[i].y=-1e12;
	}
}
mdl mdl::build(const vec &ctr,const vec &rct){
	double rad_x=rct.x*0.5,rad_y=rct.y*0.5;
	mdl m;
	m.v[0]=(vec){ctr.x-rad_x,ctr.y-rad_y};
	m.v[1]=(vec){ctr.x-rad_x,ctr.y+rad_y};
	m.v[2]=(vec){ctr.x+rad_x,ctr.y+rad_y};
	m.v[3]=(vec){ctr.x+rad_x,ctr.y-rad_y};
	return m;
}
