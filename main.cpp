#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <sstream>
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

void san_str(string &str){
	std::replace(str.begin(),str.end(),'(',' ');
	std::replace(str.begin(),str.end(),')',' ');
	std::replace(str.begin(),str.end(),',',' ');
}

void get_cls(vector<cls_s> &clss){
	// O(e*log(e)) 根据题目给出的边构建出若干个闭合回路
	for(int i=0; i<(int)clss.size(); ++i){
		std::string str;
		getline(cin,str);
		san_str(str);
		std::istringstream is(str);
		double x0,y0,x1,y1,x2,y2;
		is>>x0>>y0;
		for(x1=x0,y1=y0; ss>>x2>>y2; x1=x2,y1=y2){
			clss[i].push_back((edg){(vec){x1,y1},(vec){x2,y2}});
		}
		clss[i].push_back((vec){x2,y2},(vec){x0,y0});
	}
}

double get_dis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}
template <typename T> void flip_vec(vector <T> &vt){
	for(T &x:vt) x.flip();
}
double get_overlap(double l1,double r1,double l2,double r2){
	// (l1,r1)和(l2,r2)是否有相交的部分，inc(inclusive): 是否包含边界
	if(l1>r1) std::swap(l1,r1);
	if(l2>r2) std::swap(l2,r2);
	return min(r1,r2)-max(l1,l2);
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
	if(e.a.x==e.b.x) return get_overlap(e.a.y,e.b.y,m.v[0].y,m.v[1].y)>eps;
	return get_overlap(e.a.x,e.b.x,m.v[0].x,m.v[1].x)>eps;
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

void toposort(vector<mdl> &idx,vector<vector<int>> &mdl_nxt){
	int n=idx.size();
	vector<bool> vt(n);
	for(int i=0; i<n; ++i){

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
		double pa=min(cur_e.a.x,cur_e.b.x),a=pa-rct.x+eps*2;
		double pb=max(cur_e.a.x,cur_e.b.x),b=pb+rct.x-eps*2;
		double line_y=cur_e.a.y+rct.y*(0.5-rvld);
		double other_y=cur_e.a.y+(rvld?-rct.y:rct.y);
		// 寻找两边可以最多向外延伸多少
		for(int i=0; i<2; ++i){
			for(edg e:ebuk[i<<1|1]){
				double cx=e.a.x;
				if(get_overlap(cur_e.a.y,other_y,e.a.y,e.b.y)>eps){
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
	if(m1.size()!=m2.size()) return 1e18;
	int sz=m1.size();
	double res=0;
	for(int i=0; i<sz; ++i){
		vec v=m1[i].cntr()-m2[i].cntr();
		res+=cabs(v.x)+cabs(v.y);
	}
	return res;
}

vector<mdl> solve_seq(vector<cls_s> clss,const vector<int> &seq,const vector<mdl> &mds){
	// O(n*e*e) 对一个给定的模块摆放优先顺序做贪心最优解
	vector<mdl> res(seq.size());
	for(int id:seq){
		vec tgt=mds[id].cntr();
		vec v=(vec){cabs(mds[id].v[0].x-mds[id].v[1].x),cabs(mds[id].v[0].y-mds[id].v[1].y)};
		int best_cl=0;
		mdl mpos=get_great_pos(clss,best_cl,v,tgt);
		if(get_dis(mpos.cntr(),tgt)>1e12){
			res[id].set_inf();
		}else{
			insert_mdl(clss[best_cl],mpos);
			res[id]=mpos;
		}
	}
	return res;
}

int main(){
	// 输入格式1：（华为）
	// 第一行输入闭合回路个数n，模块个数m
	// 接下来n行，每一行按顺序输入该闭合回路上的所有点(x,y)
	// 再下来m行，每行输入模块代号，长，宽，x和y坐标
	// 输入格式2：
	// 第一行输入空白区域数量n和模块数量m
	// 接下来n行，每行第一个数是e，表示该空白区域的边界点数；接下来2e个数，依次表示边界上逆时针顺序的点的坐标；
	// 接下来m行，每行4个数w, h, x, y，表示矩形的最佳位置是(x,y)到(x+w, y+h)画出的矩形。

	drawer dw_ans("oput.png");
	dw_ans.draw_grid();

	int clcnt,mdlcnt;
	cin>>clcnt>>mdlcnt;
	vector<cls_s> org_cls(clcnt);
	vector<mdl> org_mdl(mdlcnt); // 模组最佳位置
	vector<string> mdl_name(mdlcnt); // 模组代号
	vector<string,int> mdl_idx(mdlcnt); // 模组编号
	vector<vector<int>> mdl_nxt(mdlcnt); // 连接点在它上面的模组
	get_cls(org_cls);
	for(int i=0; i<mdlcnt; ++i){
		cin>>mdl_id[i];
		vec v=vec::get(),tgt;
		string str;
		if(str[0]=='('){
			san_str(str);
			std::istringstream is(str);
			double x,y;
			is>>x>>y;
			tgt=(vec){x,y};
		}else{
			int fa=mdl_idx[str];
			mdl_nxt[fa].push_back(i);
		}
		cin>>str;
		tgt=tgt+v*0.5;
		org_mdl[i]=mdl::build(tgt,v);
	}
	vector<int> idx(mdlcnt);
	for(int i=0; i<mdlcnt; ++i) idx[i]=i;
	std::random_shuffle(idx.begin(),idx.end());
	vector<mdl> res_mdl=solve_seq(org_cls,idx,org_mdl);
	double res_len=calc_res(res_mdl,org_mdl);
	for(int times=5; times--; ){
		random_shuffle(idx.begin(),idx.end());
		vector<mdl> cur_mdl=solve_seq(org_cls,idx,org_mdl);
		double cur_len=calc_res(cur_mdl,org_mdl);
		if(res_len>cur_len){
			res_len=cur_len;
			res_mdl.swap(cur_mdl);
		}
	}
	for(int i=0; i<mdlcnt; ++i){
		if(res_mdl[i].v[0].x<0){
			dw_ans.draw_mdl(org_mdl[i],col_blue,i+1);
			cerr<<"Cannot put "<<i+1<<'.'<<endl;
		}
		dw_ans.draw_mdl(org_mdl[i],col_cyan,i+1);
		cerr<<i+1<<": "<<res_mdl[i].v[0]<<' '<<res_mdl[i].v[1]<<endl;
		dw_ans.draw_mdl(res_mdl[i],col_grey,i+1);
	}
	cerr<<"total length: "<<res_len<<endl;
	for(cls_s cl:org_cls){
		dw_ans.draw_cl(cl);
	}
	dw_ans.flush();
	return 0;
}
