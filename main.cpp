#include <iostream>
#include <vector>
#include <map>
static const int eps=1e-6;
using std::cin;
using std::vector;
using std::map;

struct vec{//向量
	double x,y;
	static vec get();
};
struct edg{//边，a为起点，b为终点，合法区域在这个边向量的左边（叉积大于0）
	vec a,b;
	bool contain(const vec &x);
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
bool edg::contain(const vec &x){
	return (b-a)*(x-a)>eps;
}

bool isvalid(rect a,cls_s &cl){
	for(edg v:cl){
		for(int i=0; i<4; ++i){
			if(!v.contain(a.v[i])) return 0;
		}
	}
	return 1;
}

void get_cls(const int edgcnt){
	vector<edg> eg;
	map<vec,std::vector<int>> pntidx;
	eg.reserve(edgcnt);
	for(int i=1,dir; i<=edgcnt; ++i){
		vec a=vec::get(),b=vec::get();
		cin>>dir;
		if(dir) std::swap(a,b);
		if(a.x==b.x||a.y==b.y){
			eg.push_back((edg){a,b});
			map[a].push_back(i);
			map[b].push_back(i);
		}else{

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
