#include <vector>
#include <map>
#include <cassert>
#include <cairo/cairo.h>
#include "structs.h"
#include "global.h"
#include "draw.h"

using std::cin; using std::cout; using std::cerr; using std::endl;
using std::min; using std::max;
using std::vector;
using std::map;
// x,y: 平面直角坐标系

using cls_s=vector<edg>;

void add_edg(cls_s &eg,map<vec,std::vector<int>> &vec_idx,vec a,vec b){
	eg.push_back((edg){a,b});
	vec_idx[a].push_back(eg.size()-1);
}

void get_cls(const int edgcnt,vector<cls_s> clss,cairo_t *cr){ // 根据题目给出的边构建rectilinear block
	vector<edg> eg;
	map<vec,vector<int>> vec_idx;
	vector<bool> vis;
	for(int i=1,dir; i<=edgcnt; ++i){
		vec a=vec::get(),b=vec::get();
		cin>>dir;
		if(dir) std::swap(a,b);
		draw_line(cr,a,b);
		if(a.x==b.x||a.y==b.y){
			add_edg(eg,vec_idx,a,b);
		}else{
			vec c; // 斜边拆分
			if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
			if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
			if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
			if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
			add_edg(eg,vec_idx,a,c);
			add_edg(eg,vec_idx,c,b);
		}
	}
	vis.resize(eg.size());
	for(int i=0; i<(int)eg.size(); ++i){
		if(vis[i]) continue;
		clss.push_back(cls_s());
		clss.rbegin()->push_back(eg[i]);
		for(int id=i,cnt=0; cnt<1e5; ++cnt){ // 防止数据不合法？
			int curdr=eg[id].dr(),res=id;
			vis[id]=1;
			for(int j:vec_idx[eg[id].b]){ // WIP
				assert(eg[j].a==eg[id].b);
				if(vis[j]&&j!=i) continue;
				if(res==id||((curdr-eg[res].dr()+4)&3)<((curdr-eg[j].dr()+4)&3))
					res=j;
			}
			assert(res!=id);
			if(res==i) break;
			id=res;
			clss.rbegin()->push_back(eg[id]);
		}
	}
}
