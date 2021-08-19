#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cairo/cairo.h>
#include <sstream>
#include <cstdlib>
#include <map>

#include "utils.h"
#include "types.h"
#include "draw.h"

// x,y: 平面直角坐标系
// n: 模块数 e: 边数 e(cl): 某个闭合回路的边数
// 时间： O(e*log(e)) + O(n*e*e)

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

void get_cls(std::vector<cls_s> &clss,drawer &dw_ans){
	// O(e) 输入边并进行存储
	for(cls_s &cl:clss){
		std::string str;
		for(; !str.length(); getline(std::cin,str));
		replace_with(str,std::vector<char>{'[',']','(',')',','},' ');
		std::istringstream is(str);
		for(double x1,y1,x2,y2; is>>x1>>y1>>x2>>y2; ){
			vec a=(vec){x1,y1},b=(vec){x2,y2};
			dw_ans.upd(a),dw_ans.upd(b);
			if(x1==x2||y1==y2){
				cl.push_back((edg){a,b});
			}else{
				vec c; // 把斜边拆分成横边和竖边
				if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
				if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
				if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
				if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
				cl.push_back((edg){a,c});
				cl.push_back((edg){c,b});
			}
		}
		sanitize_vec(cl);
	}
}

double get_dis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}
template <typename T> void flip_vec(std::vector <T> &vt){
	for(T &x:vt) x.flip();
}
double get_overlap(double l1,double r1,double l2,double r2){
	// (l1,r1)和(l2,r2)交的部分的长度
	if(l1>r1) std::swap(l1,r1);
	if(l2>r2) std::swap(l2,r2);
	return std::min(r1,r2)-std::max(l1,l2);
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

mdl get_great_pos_cl(const cls_s &cl,const std::vector<edg> ebuk[4],const vec rct,const vec tgt,const bool fliped){
	// O(e(cl)*e) 闭合回路cl的最优解
	mdl gpos;
	gpos.set_inf();
	// 以某条边为基准，看看模块至少有一个角在这贴边上时代价最少能做到多少
	for(edg cur_e:cl){
		if(cur_e.dr()&1) continue;
		bool rvld=(cur_e.dr()==2)^fliped;
		double pa=std::min(cur_e.a.x,cur_e.b.x),a=pa-rct.x+eps*2;
		double pb=std::max(cur_e.a.x,cur_e.b.x),b=pb+rct.x-eps*2;
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
				update_gpos(gpos,ed,std::min(e.a.x,b),rct,tgt,line_y);
			}
			apx(ed,e.b.x);
		}
		update_gpos(gpos,ed,b,rct,tgt,line_y);
	}
	return gpos;
}

mdl get_great_pos_basic(const std::vector<cls_s> &clss,int &best_cl,const vec rct,const vec tgt,const bool fliped){
	// O(e*e) 一个坐标系和模块方向的最优解（横向，坐标系是否经过变换）
	std::vector<edg> ebuk[4];
	for(cls_s cl:clss){
		for(edg e:cl) ebuk[e.dr()].push_back(e);
	}
	for(int i=0; i<4; ++i){
		std::sort(ebuk[i].begin(),ebuk[i].end(),
				[](const edg &a,const edg &b){ return std::min(a.a.x,a.b.x)<std::min(b.a.x,b.b.x); });
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

mdl get_great_pos(std::vector<cls_s> &clss,int &best_cl,vec rct,vec tgt){
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
		double bk_x=inf,st_x=inf,ed_x=inf,other_y=inf; // 断点，另一个y
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
		assert(st_x!=inf&&ed_x!=inf&&other_y!=inf);
		apx(bk_x,std::min(e.a.x,e.b.x));
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
	//error occured
	for(edg e:cl){
		std::cerr<<e.a<<e.b<<std::endl;
	}
	std::cerr<<"md:"<<md.v[0]<<md.v[1]<<std::endl;
	dbg_cl(cl,{md});
	assert(0);
}

double calc_res(std::vector<mdl> m1,std::vector<mdl> m2){
	// O(n)  计算一组方案的连线长度和
	if(m1.size()!=m2.size()) return inf;
	int sz=m1.size();
	double res=0;
	for(int i=0; i<sz; ++i){
		vec v=m1[i].cntr()-m2[i].cntr();
		res+=cabs(v.x)+cabs(v.y);
	}
	return res;
}

void topo_rand(std::vector<int> &idx,std::vector<int> &mdl_ref){
	// O(n*n) 引入随机化的拓扑排序 TODO 优化时间复杂度
	int n=idx.size();
	assert(n==(int)mdl_ref.size());
	std::vector<std::vector<int>> mdl_nxt(n);
	std::vector<int> que;
	for(int i=0; i<n; ++i){
		if(mdl_ref[i]==-1) que.push_back(i);
		else mdl_nxt[mdl_ref[i]].push_back(i);
	}
	for(int i=0; i<n; ++i){
		assert(!que.empty()); // 还未实现破环成链流程
		int que_idx=rnd(0,que.size()-1),x=que[que_idx]; // 目前使用确定性随机化算法便于调试
		que.erase(que.begin()+que_idx);
		for(int t:mdl_nxt[x]){
			que.push_back(t);
		}
		idx[i]=x;
	}
}

std::vector<mdl> solve_seq(std::vector<cls_s> clss,const std::vector<int> &seq,const std::vector<mdl> &mds,const std::vector<int> &mdl_ref){
	// O(n*e*e) 对一个给定的模块摆放优先顺序做贪心最优解
	std::vector<mdl> res(seq.size());
	for(int id:seq){
		vec v,tgt;
		if(mdl_ref[id]==-1){
			tgt=mds[id].cntr();
			v=(vec){cabs(mds[id].v[0].x-mds[id].v[1].x),cabs(mds[id].v[0].y-mds[id].v[1].y)};
		}else{
			tgt=res[mdl_ref[id]].cntr();
			v=mds[id].v[0];
		}
		int best_cl=0;
		mdl mpos=get_great_pos(clss,best_cl,v,tgt);
		if(get_dis(mpos.cntr(),tgt)>inf){
			res[id].set_inf();
		}else{
			insert_mdl(clss[best_cl],mpos);
			res[id]=mpos;
		}
	}
	return res;
}

int main(){
	// 输入格式1：(当前使用格式)
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

	int clcnt,mdlcnt;
	std::cin>>clcnt>>mdlcnt;
	std::vector<cls_s> org_cls(clcnt);
	std::vector<mdl> org_mdl(mdlcnt);
	std::vector<std::string> mdl_name(mdlcnt);
	std::vector<int> mdl_ref(mdlcnt,-1);
	std::map<std::string,int> mdl_idx;
	get_cls(org_cls,dw_ans);
	//std::cerr<<"---get_md---"<<std::endl;
	for(int i=0; i<mdlcnt; ++i){
		std::cin>>mdl_name[i];
		vec v=vec::get(),tgt;
		std::string str;
		std::cin>>str;
		if(str[0]=='('){
			std::string str2;
			getline(std::cin,str2);
			str+=str2;
			replace_with(str,std::vector<char>{'(',')',','},' ');
			std::istringstream is(str);
			double x,y;
			is>>x>>y;
			tgt=(vec){x,y};
			org_mdl[i]=mdl::build(tgt,v);
			dw_ans.upd(org_mdl[i].v[0]);
			dw_ans.upd(org_mdl[i].v[1]);
		}else{
			mdl_ref[i]=mdl_idx[str];
			org_mdl[i]=(mdl){{v,vec()}};
		}
	}
	std::vector<int> idx(mdlcnt);
	topo_rand(idx,mdl_ref);
	std::vector<mdl> res_mdl=solve_seq(org_cls,idx,org_mdl,mdl_ref);
	double res_len=calc_res(res_mdl,org_mdl);
	for(int times=5; times--; ){
		topo_rand(idx,mdl_ref);
		std::vector<mdl> cur_mdl=solve_seq(org_cls,idx,org_mdl,mdl_ref);
		double cur_len=calc_res(cur_mdl,org_mdl);
		if(res_len>cur_len){
			res_len=cur_len;
			res_mdl.swap(cur_mdl);
		}
	}
	dw_ans.zoom_out();
	dw_ans.draw_grid();
	for(int i=0; i<mdlcnt; ++i){
		if(res_mdl[i].v[0].x<=-inf){
			if(mdl_ref[i]==-1) dw_ans.draw_mdl(org_mdl[i],col_blue,mdl_name[i]);
			std::cerr<<"Cannot put "<<i+1<<'.'<<std::endl;
		}
		if(mdl_ref[i]==-1) dw_ans.draw_mdl(org_mdl[i],col_cyan,mdl_name[i]);
		std::cerr<<i+1<<": "<<mdl_name[i]<<','<<res_mdl[i].v[0]<<' '<<res_mdl[i].v[1]<<std::endl;
		dw_ans.draw_mdl(res_mdl[i],col_grey,mdl_name[i]);
	}
	//std::cerr<<"edge:"<<dw_ans.sf2mat((vec){0,0})<<dw_ans.sf2mat((vec){dw_ans.d_sf,dw_ans.d_sf})<<std::endl;
	std::cerr<<"total length: "<<res_len<<std::endl;
	for(cls_s cl:org_cls){
		dw_ans.draw_cl(cl);
	}
	dw_ans.flush();
	return 0;
}
