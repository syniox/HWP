#include <iostream>
#include <vector>
#include <cassert>
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

vector<cls_s> clss;//闭包集合

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
bool operator != (const vec &a,const vec &b){
	return a.x!=b.x||a.y!=b.y;
}

bool crossed(const edg &a,const edg &b){//线段是否相交，重合情况不考虑
	return ((b.a-a.a)*(a.b-a.a))*((b.b-a.a)*(a.b-a.a))<0;
}

int getdr(const edg &e){//右0上1左2下3
	assert(e.a.x==e.b.x||e.a.y==e.b.y);
	if(e.a.y==e.b.y) return (e.a.x>e.b.x)<<1;
	else return (e.a.x>e.b.x)<<1|1;
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
}

void get_cls(const int edgcnt){//根据题目给出的边构建闭包
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
			vec c;//斜边预处理
			if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
			if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
			if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
			if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
			addegbuk(eg,pntidx,a,c);
			addegbuk(eg,pntidx,c,b);
		}
	}

	vis.resize(eg.size());
	for(int i=0; i<(int)eg.size(); ++i){
		if(vis[i]) continue;
		clss.push_back(cls_s());
		clss.rbegin()->push_back(eg[i]);
		for(int id=i,cnt=0; cnt<1e5; ++cnt){//防止数据不合法？
			int curdr=getdr(eg[id]),res=id;
			for(int j:pntidx[eg[id].b]){//WIP
				if(eg[j].a!=eg[id].b) continue;
				if(res==id||((curdr-getdr(eg[res])+4)&3)>((curdr-getdr(eg[j])+4)&3))
					res=j;
			}
			assert(res!=id);
			if(res==i) break;
			id=res;
			clss.rbegin()->push_back(eg[id]);
		}
	}
}

template <typename T>
const T cabs(const T &x){
	return x<0?-x:x;
}
double getdis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}

vec getgreatpos(vec rct,cls_s &cl,bool &rot){
	//rct描述长宽，cl表示搜寻的闭包
	//rot记录是否旋转90度，函数返回左下角位置

}

int main(){
	int edgcnt,modcnt;
	cin>>edgcnt>>modcnt;
	get_cls(edgcnt);
	while(modcnt--){
		vec v=vec::get(),tgt=vec::get();
		for(cls_s cl:clss){

		}
	}
	return 0;
}
