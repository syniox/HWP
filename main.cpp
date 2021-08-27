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
#include "proc.h"

// x,y: 平面直角坐标系
// n: 模块数 e: 边数 e(cl): 某个闭合回路的边数
// 时间： O(e*log(e)) + O(n*e*e)
// assumption1: 模块的大小大于间隔大小

//TODO 用list代替vector

static int ecnt; // TODO 重新初始化
static std::vector<edg> e_vec; // 存储边

void add_bevel(cls_s &cl,const vec a,const vec b){
	vec c; // 把斜边拆分成横边和竖边
	if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
	else if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
	else if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
	else if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
	else assert(0);
	cl.push_back({a,c,-1}); // 取消斜边计数?
	cl.push_back({c,b,-1});
}

void get_cls(std::vector<cls_s> &clss,std::vector<cls_s> &input,drawer &dw_ans,double thr=4){
	// O(e) 输入边并进行存储
	int n=clss.size();
	for(int i=0; i<n; ++i){
		std::string str;
		for(; !str.length(); getline(std::cin,str));
		replace_with(str,{'[',']','(',')',','},' ');
		std::istringstream is(str);
		for(double x1,y1,x2,y2; is>>x1>>y1>>x2>>y2; ){
			vec a=(vec){x1,y1},b=(vec){x2,y2};
			dw_ans.upd(a),dw_ans.upd(b);
			if(x1==x2||y1==y2){
				input[i].push_back({a,b,++ecnt});
				clss[i].push_back({a,b,ecnt});
			}else{
				input[i].push_back({a,b,-1});
				for(vec vt=(b-a).norm(thr); (b-a).len2()>vt.len2(); a+=vt){
					add_bevel(clss[i],a,a+vt);
				}
				add_bevel(clss[i],a,b);
			}
		}
		sanitize_vec(clss[i]);
	}
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
	// debug
	std::cerr<<"rand result: ";
	for(int i:idx){
		std::cerr<<i<<' ';
	}
	std::cerr<<std::endl;
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
	std::vector<cls_s> org_cls(clcnt),input_cls(clcnt);
	std::vector<mdl> org_mdl(mdlcnt);
	std::vector<std::string> mdl_name(mdlcnt),ref_name(mdlcnt);
	std::vector<int> mdl_ref(mdlcnt,-1);
	std::map<std::string,int> mdl_idx;
	//std::cerr<<"---get_cls---"<<std::endl;
	get_cls(org_cls,input_cls,dw_ans);
	//std::cerr<<"---get_md---"<<std::endl;
	for(int i=0; i<mdlcnt; ++i){
		std::cin>>mdl_name[i];
		mdl_idx[mdl_name[i]]=i;
		vec v=vec::get(),tgt;
		std::string str;
		std::cin>>str;
		if(str[0]=='('){
			std::string str2;
			getline(std::cin,str2);
			str+=str2;
			replace_with(str,{'(',')',','},' ');
			std::istringstream is(str);
			double x,y;
			is>>x>>y;
			tgt=(vec){x,y};
			org_mdl[i]=mdl::build(tgt,v);
			dw_ans.upd(org_mdl[i].v[0]);
			dw_ans.upd(org_mdl[i].v[1]);
		}else{
			ref_name[i]=str;
			org_mdl[i]={v,vec()};
		}
	}
	for(int i=0; i<mdlcnt; ++i){
		if(ref_name[i]!=""){
			mdl_ref[i]=mdl_idx[ref_name[i]];
		}
	}
	//std::cerr<<"---get_lim---"<<std::endl;

	std::vector<int> idx(mdlcnt);
	topo_rand(idx,mdl_ref);
	std::vector<mdl> res_mdl=solve_seq(org_cls,idx,org_mdl,mdl_ref,e_vec,mdlcnt);
	double res_len=calc_res(res_mdl,org_mdl,mdl_ref);
	for(int times=20; times--; ){
		topo_rand(idx,mdl_ref);
		std::vector<mdl> cur_mdl=solve_seq(org_cls,idx,org_mdl,mdl_ref,e_vec,mdlcnt);
		double cur_len=calc_res(cur_mdl,org_mdl,mdl_ref);
		if(res_len>cur_len){
			std::cerr<<"better."<<std::endl; // debug
			res_len=cur_len;
			res_mdl.swap(cur_mdl);
		}
	}
	dw_ans.zoom_out();
	dw_ans.draw_grid();
	for(int i=0; i<mdlcnt; ++i){
		if(res_mdl[i].v[0].x<=-inf){
			std::cerr<<"Cannot put "<<i+1<<'.'<<std::endl;
			continue;
		}
		std::cerr<<i<<": "<<mdl_name[i]<<','<<res_mdl[i].v[0]<<' '<<res_mdl[i].v[1]<<std::endl;
		dw_ans.draw_mdl(res_mdl[i],col_grey,mdl_name[i]);
	}
	for(int i=0; i<mdlcnt; ++i){
		dw_ans.draw_pnt(res_mdl[i].cntr());
		if(mdl_ref[i]==-1){
			dw_ans.draw_line(res_mdl[i].cntr(),org_mdl[i].cntr(),col_blue,0.9,1,4);
			dw_ans.draw_pnt(org_mdl[i].cntr());
		}else{
			dw_ans.draw_line(res_mdl[i].cntr(),res_mdl[mdl_ref[i]].cntr(),col_cyan,0.9,1,4);
			dw_ans.draw_pnt(res_mdl[mdl_ref[i]].cntr());
		}
	}
	//std::cerr<<"edge:"<<dw_ans.sf2mat({0,0})<<dw_ans.sf2mat({dw_ans.d_sf,dw_ans.d_sf})<<std::endl;
	std::cerr<<"total length: "<<res_len<<std::endl;
#ifdef DEBUG
	for(cls_s cl:org_cls){
		dw_ans.draw_cl(cl,5);
	}
#else
	for(cls_s cl:input_cls){
		dw_ans.draw_cl(cl,0);
	}
#endif
	dw_ans.flush();
	return 0;
}
