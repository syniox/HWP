#include <iostream>
#include <vector>
#include <map>
static const int eps=1e-6;
using std::cin;
using std::vector;
using std::map;
//x,y: 平面直角坐标系

struct vec{//向量
	double x,y;
	static vec get();
};
struct edg{//边，a为起点，b为终点，合法区域在这个边向量的左边
	vec a,b;
};
struct rect{
	vec v[4];
};

using cls_s=vector<edg>;

std::vector<cls_s> clss;

vec vec::get(){
	double x,y;
	cin>>x>>y;
	return (vec){x,y};
}
double operator *(const vec &a,const vec &b){
	return a.x*b.y-a.y*b.x;
}
vec operator - (const vec &a,const vec &b){
	return (vec){a.x-b.x,a.y-b.y};
}
bool operator < (const vec &a,const vec &b){
	return a.x==b.x?a.y<b.y:a.x<b.x;
}

bool crossed(const edg &a,const edg &b){//重合情况不考虑
	return ((b.a-a.a)*(a.b-a.a)>0)==((b.b-a.a)*(a.b-a.a)>0);
}

bool isvalid(rect a,cls_s &cl){
	for(edg v:cl){
		for(int i=0; i<4; ++i){
			if(crossed((edg){a.v[i],a.v[(i+1)&3]},v))
				return 0;
		}
	}
	return 1;
}

void addegbuk(vector<edg> &eg,map<vec,std::vector<int>> &pntidx,vec a,vec b){
	eg.push_back((edg){a,b});
	pntidx[a].push_back(eg.size());
	pntidx[b].push_back(eg.size());
}

void get_cls(const int edgcnt){
	vector<edg> eg;
	map<vec,std::vector<int>> pntidx;
	vector<bool> vis;
	for(int i=1,dir; i<=edgcnt; ++i){
		vec a=vec::get(),b=vec::get();
		cin>>dir;
		if(dir) std::swap(a,b);
		if(a.x==b.x||a.y==b.y){
			addegbuk(eg,pntidx,a,b);
		}else{
			vec c;
			if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
			if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
			if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
			if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
			addegbuk(eg,pntidx,a,c);
			addegbuk(eg,pntidx,c,b);
		}
	}
}

vec getgreatpos(rect a,std::vector<edg> &cl){

}

int main(){
	int edgcnt,modcnt;
	cin>>edgcnt>>modcnt;
	get_cls(edgcnt);
	while(modcnt--){

	}
	return 0;
}
