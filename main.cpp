#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cairo/cairo.h>
#include "types.h"
#include <map>
static const double eps=1e-6;
using std::cin; using std::cout; using std::cerr; using std::endl;
using std::min; using std::max;
using std::vector;
using std::map;
// x,y: 平面直角坐标系

using cls_s=vector<edg>;
cls_s org_edg;

// gry: 矩形对应的位置	red: 合法区域边框
// grn: 边界 			blu: 无法放入的矩形
// cyan: 被放入的矩形
static const col_s col_gry=col_s{0.4,0.4,0.4},col_red=col_s{0.8,0.1,0.3};
static const col_s col_grn=col_s{0.1,0.7,0.1},col_blu=col_s{0.2,0.4,1.0};
static const col_s col_cyan=col_s{0.1,0.8,0.8};

template <typename T> inline void apn(T &x,const T y){
	x=x<y?x:y;
}
template <typename T> inline void apx(T &x,const T y){
	x=x>y?x:y;
}
template <typename T> const T cabs(const T &x){
	return x<0?-x:x;
}

void add_edg(vector<edg> &eg,map<vec,std::vector<int>> &vec_idx,vec a,vec b){
	eg.push_back((edg){a,b});
	vec_idx[a].push_back(eg.size()-1);
}

void draw_line(cairo_t *cr,vec x,vec y,col_s c=col_red){
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	cairo_set_line_width(cr,1);
	cairo_move_to(cr,x.x*10,x.y*10);
	cairo_line_to(cr,y.x*10,y.y*10);
	cairo_stroke(cr);
}
void draw_mdl(cairo_t *cr,mdl m,col_s c=col_gry,int id=-1){
	static char ch[10];
	cairo_set_source_rgba(cr,c.r,c.g,c.b,1.0);
	double x=m.v[0].x,y=m.v[0].y,dx=m.v[1].x-x,dy=m.v[1].y-y;
	if(dx<0) dx=-dx,x-=dx;
	if(dy<0) dy=-dy,y-=dy;
	cairo_rectangle(cr,x*10,y*10,dx*10,dy*10);
	cairo_fill(cr);
	draw_line(cr,vec{x,y},vec{x+dx,y},col_grn);
	draw_line(cr,vec{x,y},vec{x,y+dy},col_grn);
	draw_line(cr,vec{x+dx,y},vec{x+dx,y+dy},col_grn);
	draw_line(cr,vec{x,y+dy},vec{x+dx,y+dy},col_grn);
	if(id==-1) return;
	sprintf(ch,"%d",id);
	cairo_set_source_rgba(cr,1,1,1,1);
	cairo_set_font_size(cr,12);
	cairo_move_to(cr,(x+dx/2)*10,(y+dy/2)*10);
	cairo_show_text(cr,ch);
}
void dbg_cl(const cls_s &cl){
	const char* oput_png="dbg.png";
	cairo_surface_t *surface;
	surface=cairo_image_surface_create(CAIRO_FORMAT_ARGB32,500,500);
	cairo_t *cr=cairo_create(surface);
	for(edg e:cl){
		draw_line(cr,e.a,e.b);
	}
	cairo_surface_write_to_png(surface,oput_png);
	cairo_surface_destroy(surface);
}

void get_cls(vector<cls_s> &clss,const int edgcnt,cairo_t *cr){ // 根据题目给出的边构建rectilinear block
	vector<edg> eg;
	map<vec,vector<int>> vec_idx;
	vector<bool> vis;
	for(int i=1,dir; i<=edgcnt; ++i){
		vec a=vec::get(),b=vec::get();
		cin>>dir;
		if(dir) std::swap(a,b);
		org_edg.push_back(edg{a,b});
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

double get_dis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}
template <typename T> void flip_vec(vector <T> &vt){
	for(T &x:vt) x.flip();
}
bool on_edge(const edg &e,const vec &p,bool inc=1){ // inclusive: 在边界上算不算
	int wgt=inc?1:-1;
	if(e.a.x==e.b.x){
		double l=std::min(e.a.y,e.b.y),r=std::max(e.a.y,e.b.y);
		return cabs(p.x-e.a.x)<eps&&p.y>=l-wgt*eps&&p.y<=r+wgt*eps;
	}
	double l=std::min(e.a.x,e.b.x),r=std::max(e.a.x,e.b.x);
	return cabs(p.y-e.a.y)<eps&&p.x>=l-wgt*eps&&p.x<=r+wgt*eps;
}
bool on_edge(const edg &e,const mdl &m){
	int cnt=0;
	for(int i=0; i<2; ++i){
		cnt+=on_edge(e,m.v[i],0);
	}
	for(int i=0; i<2; ++i){
		cnt+=on_edge(e,(vec){m.v[i].x,m.v[i^1].y},0);
	}
	return cnt;
}

void update_gpos(mdl &gpos,double l,double r,const vec rct,const vec tgt,const double line_y){
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
	// O(e(cl)*e) 一个闭合回路的最优解
	mdl gpos;
	gpos.set_inf();
	for(edg cur_e:cl){
		if(cur_e.dr()&1) continue;
		bool rvld=(cur_e.dr()==2)^fliped;
		double pa=min(cur_e.a.x,cur_e.b.x),a=pa-rct.x;
		double pb=max(cur_e.a.x,cur_e.b.x),b=pb+rct.x;
		double line_y=cur_e.a.y+rct.y*(0.5-rvld);
		double other_y=cur_e.a.y+(rvld?-rct.y:rct.y);
		for(int i=0; i<2; ++i){
			for(edg e:ebuk[i<<1|1]){
				double cx=e.a.x;
				vec p1=(vec){cx,cur_e.a.y},p2=(vec){cx,other_y};
				if(on_edge((edg){p1,p2},e.a)||on_edge((edg){p1,p2},e.b)){
					if(cx<pa+eps) apx(a,cx);
					if(cx>pb-eps) apn(b,cx);
				}
			}
		}
		double ed=pa;
		for(edg e:ebuk[cur_e.dr()^2]){
			if(cabs(e.a.y-line_y)*2>rct.y) continue;
			if(!cur_e.dr()) std::swap(e.a,e.b);
			if(e.a.x<ed+eps){
				apx(ed,e.b.x);
				continue;
			}
			update_gpos(gpos,ed,e.a.x,rct,tgt,line_y);
		}
		update_gpos(gpos,ed,pb,rct,tgt,line_y);
	}
	return gpos;
}

mdl get_great_pos_basic(const vector<cls_s> &clss,int &best_cl,const vec rct,const vec tgt,const bool fliped){
	// O(e^2) 一个坐标系的最优解（横向，坐标系是否经过变换）
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

mdl get_great_pos(vector<cls_s> &clss,int &best_cl,vec rct,vec tgt){// 寻找某个闭包的最优位置
	// cl表示搜寻的闭包
	// 函数返回模块最后占用的位置
	mdl mpos[4];
	int bcl[4]={0};
	mpos[0]=get_great_pos_basic(clss,bcl[0],rct,tgt,0); // 横着的原矩阵 横向rb
	rct.flip();
	mpos[1]=get_great_pos_basic(clss,bcl[1],rct,tgt,0); // 竖着的原矩阵 横向rb
	for(cls_s cl:clss) flip_vec(cl);
	mpos[2]=get_great_pos_basic(clss,bcl[2],rct,tgt,1); // 横着的原矩阵 竖向rb（坐标系颠倒）
	rct.flip();
	mpos[3]=get_great_pos_basic(clss,bcl[3],rct,tgt,1); // 竖着的原矩阵 竖向rb（坐标系颠倒）
	for(cls_s cl:clss) flip_vec(cl);
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

void insert_mdl(cls_s &cl,mdl md){ // 将该区域设为不可用区域（假设该模块紧贴边缘）
	using cls_i=cls_s::iterator;
	for(cls_i it=cl.begin(); it!=cl.end(); ++it){
		edg e=*it;
		if(!on_edge(e,md)) continue;
		bool fliped=0;
		if(e.dr()&1){
			fliped=1;
			flip_vec(cl);
			e.flip();
			md.flip();
		}
		double bk_x=1e9,st_x,ed_x,other_y; // 断点，另一个y
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

	cairo_surface_t *surface;
	surface=cairo_image_surface_create(CAIRO_FORMAT_ARGB32,500,500);
	cairo_t *cr=cairo_create(surface);

	int edgcnt,mdlcnt;
	cin>>edgcnt>>mdlcnt;
	vector<cls_s> org_cls;
	get_cls(org_cls,edgcnt,cr);
	for(int i=1; i<=mdlcnt; ++i){
		vec tgt=vec::get(),v=vec::get(); // 返回最优位置的中心？
		tgt=tgt+v*0.5;
		int best_cl=0;
		mdl mpos=get_great_pos(org_cls,best_cl,v,tgt);
		if(get_dis(mpos.cntr(),tgt)>1e12){
			draw_mdl(cr,mdl::build(tgt,v),col_blu,i);
			cout<<i<<": "<<"cannot be put."<<endl;
		}else{
			insert_mdl(org_cls[best_cl],mpos);
			cout<<i<<": "<<mpos.v[0]<<' '<<mpos.v[1]<<endl;
			draw_mdl(cr,mdl::build(tgt,v),col_cyan,i);
			draw_mdl(cr,mpos,col_gry,i);
		}
	}
	for(edg e:org_edg){
		draw_line(cr,e.a,e.b);
	}
	cairo_surface_write_to_png(surface,"test.png");
	cairo_surface_destroy(surface);
	return 0;
}
