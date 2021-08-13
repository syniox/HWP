#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cairo/cairo.h>
#include <map>

#include "types.h"
#include "draw.h"

static const double eps=1e-6;
using std::cin; using std::cout; using std::cerr; using std::endl;
using std::min; using std::max;
using std::vector;
using std::map;
// x,y: 平面直角坐标系
// n: 模块数 e: 边数 e(cl): 某个闭合回路的边数
// 时间： O(e*log(e)) + O(n*e*e)

cls_s org_edg;

template <typename T> inline void apn(T &x,const T y){
	x=x<=y?x:y;
}
template <typename T> inline void apx(T &x,const T y){
	x=x>=y?x:y;
}
template <typename T> const T cabs(const T &x){
	return x<0?-x:x;
}

void add_edg(vector<edg> &eg,map<vec,std::vector<int>> &vec_idx,vec a,vec b){
	eg.push_back((edg){a,b});
	vec_idx[a].push_back(eg.size()-1);
}

void get_cls(vector<cls_s> &clss,const int edgcnt){
	// O(e*log(e)) 根据题目给出的边构建出若干个闭合回路
	vector<edg> eg;
	map<vec,vector<int>> vec_idx;
	vector<bool> vis;
	// 输入边并进行存储
	for(int i=1,dir; i<=edgcnt; ++i){
		vec a=vec::get(),b=vec::get();
		cin>>dir;
		if(dir) std::swap(a,b);
		org_edg.push_back(edg{a,b});
		if(a.x==b.x||a.y==b.y){
			add_edg(eg,vec_idx,a,b);
		}else{
			vec c; // 把斜边拆分成横边和竖边
			if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
			if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
			if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
			if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
			add_edg(eg,vec_idx,a,c);
			add_edg(eg,vec_idx,c,b);
		}
	}
	// 构建闭合回路
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

double get_dis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}
template <typename T> void flip_vec(vector <T> &vt){
	for(T &x:vt) x.flip();
}
bool have_crs(double l1,double r1,double l2,double r2,bool inc){
	// (l1,r1)和(l2,r2)是否有相交的部分，inc(inclusive): 是否包含边界
	if(l1>r1) std::swap(l1,r1);
	if(l2>r2) std::swap(l2,r2);
	return max(l1,l2)-inc*eps<min(r1,r2);
}
bool on_line(const edg &e,const vec &p){
	if(e.a.x==e.b.x) return cabs(p.x-e.a.x)<eps;
	return cabs(p.y-e.a.y)<eps;
}
bool on_edge(const edg &e,const vec &p,bool inc){
	// 判断一个向量是否在这条边上，inc(inclusive): 在边界上算不算
	int wgt=inc?1:-1;
	if(e.a.x==e.b.x){
		double l=std::min(e.a.y,e.b.y),r=std::max(e.a.y,e.b.y);
		return cabs(p.x-e.a.x)<eps&&p.y>=l-wgt*eps&&p.y<=r+wgt*eps;
	}
	double l=std::min(e.a.x,e.b.x),r=std::max(e.a.x,e.b.x);
	return cabs(p.y-e.a.y)<eps&&p.x>=l-wgt*eps&&p.x<=r+wgt*eps;
}
bool on_edge(const edg &e,const mdl &m){
	// 模块m的边是否和边e有相交的部分
	int p=0;
	for(; p<2&&!on_line(e,m.v[p]); ++p);
	if(p==2) return 0;
	if(e.a.x==e.b.x) return have_crs(e.a.y,e.b.y,m.v[0].y,m.v[1].y,0);
	return have_crs(e.a.x,e.b.x,m.v[0].x,m.v[1].y,0);
}

void update_gpos(mdl &gpos,double l,double r,const vec rct,const vec tgt,const double line_y){
	// 对某一个合法段[l,r] 找到他的最优解并看看能不能更新当前最优答案gpos
	if(r-l<rct.x) return;
	double x,rad_x=rct.x*0.5;
	if(r<tgt.x+rad_x) x=r-rad_x;
	else if(l>tgt.x-rad_x) x=l+rad_x;
	else x=tgt.x;
	mdl neo=mdl::build((vec){x,line_y},rct);
	if(get_dis(neo.cntr(),tgt)<get_dis(gpos.cntr(),tgt)){
		gpos=neo;
	}
}

mdl get_great_pos_cl(const cls_s &cl,const vector<edg> ebuk[4],const vec rct,const vec tgt,const bool fliped){
	// O(e(cl)*e) 闭合回路cl的最优解
	mdl gpos;
	gpos.set_inf();
	// 以某条边为基准，看看模块至少有一个角在这贴边上时代价最少能做到多少
	for(edg cur_e:cl){
		if(cur_e.dr()&1) continue;
		bool rvld=(cur_e.dr()==2)^fliped;
		double pa=min(cur_e.a.x,cur_e.b.x),a=pa-rct.x;
		double pb=max(cur_e.a.x,cur_e.b.x),b=pb+rct.x;
		double line_y=cur_e.a.y+rct.y*(0.5-rvld);
		double other_y=cur_e.a.y+(rvld?-rct.y:rct.y);
		// 寻找两边可以最多向外延伸多少
		for(int i=0; i<2; ++i){
			for(edg e:ebuk[i<<1|1]){
				double cx=e.a.x;
				if(have_crs(cur_e.a.y,other_y,e.a.y,e.b.y,0)){
					if(cx<pa+eps) apx(a,cx);
					if(cx>pb-eps) apn(b,cx);
				}
			}
		}
		double ed=a;
		// 遍历与这条边方向相反的边，找出合法的段，计算最小代价并更新答案
		for(edg e:ebuk[cur_e.dr()^2]){
			if(cabs(e.a.y-line_y)*2>rct.y-eps&&cabs(e.a.y-cur_e.a.y)>eps) continue;
			if(!cur_e.dr()) std::swap(e.a,e.b);
			if(e.a.x>ed){
				update_gpos(gpos,ed,min(e.a.x,b),rct,tgt,line_y);
			}
			apx(ed,e.b.x);
		}
		update_gpos(gpos,ed,b,rct,tgt,line_y);
	}
	return gpos;
}

mdl get_great_pos_basic(const vector<cls_s> &clss,int &best_cl,const vec rct,const vec tgt,const bool fliped){
	// O(e*e) 一个坐标系和模块方向的最优解（横向，坐标系是否经过变换）
	vector<edg> ebuk[4];
	for(cls_s cl:clss){
		for(edg e:cl) ebuk[e.dr()].push_back(e);
	}
	for(int i=0; i<4; ++i){
		std::sort(ebuk[i].begin(),ebuk[i].end(),
				[](const edg &a,const edg &b){ return min(a.a.x,a.b.x)<min(b.a.x,b.b.x); });
	}
	mdl gpos; // great pos
	gpos.set_inf();
	for(int i=0; i<(int)clss.size(); ++i){
		mdl pos=get_great_pos_cl(clss[i],ebuk,rct,tgt,fliped);
		if(get_dis(pos.cntr(),tgt)<get_dis(gpos.cntr(),tgt)){
			gpos=pos,best_cl=i;
		}
	}
	return gpos;
}

mdl get_great_pos(vector<cls_s> &clss,int &best_cl,vec rct,vec tgt){
	// O(e*e) 寻找模块摆放的最优位置
	// 函数返回模块最后占用的位置
	// TODO: 把排序函数从basic中提出来
	// TODO: flip_vec &cl safe?
	mdl mpos[4];
	int bcl[4]={0};
	mpos[0]=get_great_pos_basic(clss,bcl[0],rct,tgt,0); // 横着的原矩阵 横边
	rct.flip();
	mpos[1]=get_great_pos_basic(clss,bcl[1],rct,tgt,0); // 竖着的原矩阵 横边
	for(cls_s &cl:clss) flip_vec(cl);
	tgt.flip();
	mpos[2]=get_great_pos_basic(clss,bcl[2],rct,tgt,1); // 横着的原矩阵 竖边 坐标系颠倒
	rct.flip();
	mpos[3]=get_great_pos_basic(clss,bcl[3],rct,tgt,1); // 竖着的原矩阵 竖边 坐标系颠倒
	for(cls_s &cl:clss) flip_vec(cl);
	tgt.flip();
	mpos[2].flip(),mpos[3].flip();

	mdl pos=mpos[0];
	best_cl=bcl[0];
	double res=get_dis(mpos[0].cntr(),tgt);
	for(int i=1; i<4; ++i){
		double tmp=get_dis(mpos[i].cntr(),tgt);
		if(tmp<res) res=tmp,pos=mpos[i],best_cl=bcl[i];
	}
	return pos;
}

void sanitize_vec(cls_s &cl){
	//O(e(cl)*e(cl)) 删除长度为0的边，合并相邻且方向相反的边
	// TODO: optimize
	for(int cnt=-1; cnt; ){
		cnt=0;
		for(int i=cl.size()-1; i>=0; --i){
			if(cl[i].ispnt()){
				++cnt;
				cl.erase(cl.begin()+i);
			}
		}
		for(int i=0; i<(int)cl.size(); ++i){
			int a=i,b=(i+1)%cl.size();
			if(((cl[a].dr()^cl[b].dr())&1)==0){
				++cnt;
				assert(cl[a].b==cl[b].a);
				cl[a].b=cl[b].b;
				cl.erase(cl.begin()+b);
				break;
			}
		}
	}
}

void insert_mdl(cls_s &cl,mdl md){
	// O(e**2) 将模块摆放的区域设为不可用区域（假设该模块紧贴边缘）
	// 找到一条和模块相邻的边，在他们的公共位置上找一个断点，把它作为绘画的起点和终点，
	// 画出模块的轮廓，然后对无效的边界进行整理
	using cls_i=cls_s::iterator;
	for(cls_i it=cl.begin(); it!=cl.end(); ++it){
		// 找一条与这个模块相邻的边
		edg e=*it;
		if(!on_edge(e,md)) continue;
		bool fliped=0;
		if(e.dr()&1){
			fliped=1;
			flip_vec(cl);
			e.flip();
			md.flip();
		}
		// 寻找切分点和模块的另一边的位置
		double bk_x=1e9,st_x=-1,ed_x=-1,other_y=-1; // 断点，另一个y
		for(int i=0; i<2; ++i){
			if(cabs(md.v[i].y-e.a.y)>eps){
				other_y=md.v[i].y;
			}
			if((md.v[i].x>md.v[i^1].x)==(e.a.x>e.b.x)){
				st_x=md.v[i].x;
				ed_x=md.v[i^1].x;
			}
			apn(bk_x,md.v[i].x);
		}
		assert(st_x>=0&&ed_x>=0&&other_y>=0);
		apx(bk_x,min(e.a.x,e.b.x));
		vec pb=vec{bk_x,e.a.y},p1=vec{st_x,e.a.y};
		vec p2=vec{st_x,other_y},p3=vec{ed_x,other_y},p4=vec{ed_x,e.a.y};
		it->b=pb;
		it=cl.insert(++it,(edg){pb,p1}),it=cl.insert(++it,(edg){p1,p2});
		it=cl.insert(++it,(edg){p2,p3}),it=cl.insert(++it,(edg){p3,p4});
		it=cl.insert(++it,(edg){p4,pb}),it=cl.insert(++it,(edg){pb,e.b});
		if(fliped){
			flip_vec(cl);
		}
		sanitize_vec(cl);
		return;
	}
	assert(0);
}

double calc_res(vector<mdl> m1,vector<mdl> m2){
	// O(n)  计算一组方案的连线长度和
	assert(m1.size()==m2.size());
	int sz=m1.size();
	double res=0;
	for(int i=0; i<sz; ++i){
		vec v=m1[i].cntr()-m2[i].cntr();
		res+=cabs(v.x)+cabs(v.y);
	}
	return res;
}

int main(){
	// 输入格式1：
	// 第一行输入空白区域数量n和模块数量m
	// 接下来n行，每行第一个数是e，表示该空白区域的边界点数；接下来2e个数，依次表示边界上逆时针顺序的点的坐标；
	// 接下来m行，每行4个数w, h, x, y，表示矩形的最佳位置是(x,y)到(x+w, y+h)画出的矩形。
	// 输入2：
	// 第一行输入有多少条边界e和多少个模块m
	// 接下来e行输入每条边（向量）的起止坐标x1,y1,x2,y2和空白区域的位置在向量的左边还是右边，左边为0右边为1
	// 接下来m行每行输入4个数，代表该模块的长和宽，最优位置中心的x坐标和y坐标
	// 输出：
	// 共m行，每行输出该模块摆放位置的对角端点

	drawer dw_ans("oput.png");
	dw_ans.draw_grid();

	int edgcnt,mdlcnt;
	cin>>edgcnt>>mdlcnt;
	vector<cls_s> org_cls;
	get_cls(org_cls,edgcnt);
	for(int i=1; i<=mdlcnt; ++i){
		vec tgt=vec::get(),v=vec::get(); // 返回最优位置的中心？
		tgt=tgt+v*0.5;
		int best_cl=0;
		mdl mpos=get_great_pos(org_cls,best_cl,v,tgt);
		if(get_dis(mpos.cntr(),tgt)>1e12){
			dw_ans.draw_mdl(mdl::build(tgt,v),col_blu,i);
			cout<<i<<": "<<"cannot be put."<<endl;
		}else{
			insert_mdl(org_cls[best_cl],mpos);
			cout<<i<<": "<<mpos.v[0]<<' '<<mpos.v[1]<<endl;
			dw_ans.draw_mdl(mdl::build(tgt,v),col_cyan,i);
			dw_ans.draw_mdl(mpos,col_gry,i);
		}
	}
	for(edg e:org_edg){
		dw_ans.draw_line(e.a,e.b);
	}
	dw_ans.flush();
	return 0;
}
