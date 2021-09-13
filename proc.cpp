// 给定一个模块接入顺序，寻找最优方案的主流程

#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

#include "types.h"
#include "utils.h"
#include "draw.h"


static std::vector<edg> e_vec;
static std::vector<std::array<int,4>> e_mdl;

double get_dis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}
double get_overlap(double l1,double r1,double l2,double r2){
	// (l1,r1)和(l2,r2)交的部分的长度或是相差的距离（负数时）
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
		// pa, pb: 边的极大极小值
		// a,b: 模块依附在这条边上时的最左最右值
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
			if(cabs(e.a.y-line_y)*2>rct.y-eps&&cabs(e.a.y-cur_e.a.y)>eps) continue; // 非法段连续
			// 非法段断开（则中间的为合法段）
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
	// 预处理 把边分类按坐标排序
	for(cls_s cl:clss){
		for(edg e:cl) ebuk[e.dr()].push_back(e);
	}
	for(int i=0; i<4; ++i){
		std::sort(ebuk[i].begin(),ebuk[i].end(),
				[](const edg &a,const edg &b){ return std::min(a.a.x,a.b.x)<std::min(b.a.x,b.b.x); });
	}
	mdl gpos; // great pos
	gpos.set_inf();
	// 遍历所有合法闭环，寻找最优解 TODO: 直接把边拆散时是否边数会过多?
	for(int i=0; i<(int)clss.size(); ++i){
		mdl pos=get_great_pos_cl(clss[i],ebuk,rct,tgt,fliped);
		if(get_dis(pos.cntr(),tgt)<get_dis(gpos.cntr(),tgt)){
			gpos=pos,best_cl=i;
		}
	}
	return gpos;
}

mdl get_great_pos(std::vector<cls_s> &clss,int &best_cl,vec rct,vec tgt,std::vector<edg> e_lim){
	// O(e*e) 寻找模块摆放的最优位置
	// 遍历矩阵和依附边的朝向（交换矩阵xy或坐标系的xy），函数返回模块最后占用的位置
	// TODO: 把排序函数从basic中提出来
	// TODO: flip_vec &cl safe?
	mdl mpos[4];
	int bcl[4]={0};
	mpos[0]=get_great_pos_basic(clss,bcl[0],rct,tgt,0); // 横着的原矩阵 横边
	rct.flip();
	mpos[1]=get_great_pos_basic(clss,bcl[1],rct,tgt,0); // 竖着的原矩阵 横边
	flip_vec(e_lim);
	for(cls_s &cl:clss) flip_vec(cl);
	tgt.flip();
	mpos[2]=get_great_pos_basic(clss,bcl[2],rct,tgt,1); // 横着的原矩阵 竖边 坐标系颠倒
	rct.flip();
	mpos[3]=get_great_pos_basic(clss,bcl[3],rct,tgt,1); // 竖着的原矩阵 竖边 坐标系颠倒
	flip_vec(e_lim);
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

bool insert_mdl(cls_s &cl,mdl md,const int id){
	// O(e**2) 将模块摆放的区域设为不可用区域（假设该模块紧贴边缘）
	// 找到一条和模块相邻的边，在他们的公共位置上找一个断点，把它作为绘画的起点和终点，
	// 画出模块的轮廓，然后对无效的边界进行整理
	// 返回是否能与某条边相连
	// TODO 相邻的边不一定在这上面？
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
		// 边(st_x,ed_x)在bk_x处断开，模块放上去后另一条边的y坐标是other_y
		double bk_x=inf,st_x=inf,ed_x=inf,other_y=inf;
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
		// 把那条边从切分点破开，绕模块一圈后从断点继续出发
		vec pb=vec{bk_x,e.a.y},p1=vec{st_x,e.a.y};
		vec p2=vec{st_x,other_y},p3=vec{ed_x,other_y},p4=vec{ed_x,e.a.y};
		it->b=pb;
		it=cl.insert(++it,{pb,p1});
		e_vec.push_back({p1,p2}),it=cl.insert(++it,{p1,p2});
		e_vec.push_back({p2,p3}),it=cl.insert(++it,{p2,p3});
		e_vec.push_back({p3,p4}),it=cl.insert(++it,{p3,p4});
		it=cl.insert(++it,{p4,pb});
		it=cl.insert(++it,{pb,e.b});
		int ecnt=e_vec.size()-1;
		if(fliped){
			e_mdl[id]={ecnt-0,ecnt-1,ecnt-2,0};
			flip_vec(cl);
		}else{
			e_mdl[id]={ecnt-1,ecnt-0,0,ecnt-2};
		}
		sanitize_vec(cl);
		return 1;
	}
	return 0;
}

std::vector<mdl> solve_seq(std::vector<cls_s> clss,const std::vector<int> &seq,
		const std::vector<mdl> &mds,const std::vector<int> &mdl_ref,
		std::vector<edg> e_vec_org,const int mdlcnt){
	// O(n*e*e) 对一个给定的模块摆放优先顺序做贪心最优解
	e_vec.swap(e_vec_org);
	e_mdl.resize(mdlcnt);
	std::vector<mdl> res(seq.size());
	for(mdl &m: res) m.set_inf();
	for(int id:seq){
		vec v,tgt;
		//寻找模块连接位置
		if(mdl_ref[id]==-1){
			tgt=mds[id].cntr();
			v=(vec){cabs(mds[id].v[0].x-mds[id].v[1].x),cabs(mds[id].v[0].y-mds[id].v[1].y)};
		}else{
			tgt=res[mdl_ref[id]].cntr();
			v=mds[id].v[0];
		}
		//寻找离最优位置最近的点
		int best_cl=0;
		mdl mpos=get_great_pos(clss,best_cl,v,tgt,{}); // TODO 使用lim_e
		if(get_dis(mpos.cntr(),tgt)>inf){
			res[id].set_inf();
		}else{
			if(!insert_mdl(clss[best_cl],mpos,id)){ // best_cl?
				clss.push_back(cls_s());
				cls_s &t=clss[clss.size()-1];
				add_bevel(t,mpos.v[0],mpos.v[1]);
				add_bevel(t,mpos.v[1],mpos.v[0]);
			}
			res[id]=mpos;
		}
	}
	return res;
}
