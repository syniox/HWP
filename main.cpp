#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include <cctype>
#include <map>

#include "utils.h"
#include "types.h"
#include "draw.h"
#include "proc.h"

// x,y: 平面直角坐标系
// n: 模块数 e: 边数 e(cl): 某个闭合回路的边数
// 时间： O(e*log(e)) + O(n*e*e)
// assumption1: 模块的大小大于间隔大小
// assumption2: 间隔要求只适用于模块在另一个方向上的投影与它有交的情况

//TODO 用list代替vector

static std::vector<edg> e_vec; // 存储边

void get_cls(std::vector<cls_s> &clss,std::vector<cls_s> &input,drawer &dw_ans,double thr=4){
	// O(e) 输入边并进行存储
	int n=clss.size();
	for(int i=0; i<n; ++i){
		std::string str=get_line({'[',']','(',')',','});
		std::istringstream is(str);
		for(double x1,y1,x2,y2; is>>x1>>y1>>x2>>y2; ){
			vec a=(vec){x1,y1},b=(vec){x2,y2};
			dw_ans.upd(a),dw_ans.upd(b);
			input[i].push_back({a,b});
			if(x1==x2||y1==y2){
				clss[i].push_back({a,b});
				e_vec.push_back({a,b});
			}else{
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
	// 第一行输入空白区域数量n 模块数量m 限制条数l
	// 接下来n行，每行第一个数是e，表示该空白区域的边界点数；接下来2e个数，依次表示边界上逆时针顺序的点的坐标；
	// 接下来m行，每行4个数w, h, x, y，表示矩形的最佳位置是(x,y)到(x+w, y+h)画出的矩形。
	// 接下来l行，每输入两个字符串，字母开头代表模块，数字开头代表边
	// 输出：
	// 前m行每行输出该模块摆放位置的对角端点
	// 最后一行输出连线的总长度

	drawer dw_ans("oput.png");

	int clcnt,mdlcnt,limcnt;
	std::cin>>clcnt>>mdlcnt>>limcnt;
	std::vector<cls_s> org_cls(clcnt),input_cls(clcnt);
	std::vector<mdl> org_mdl(mdlcnt);
	std::vector<std::string> mdl_name(mdlcnt),ref_name(mdlcnt);
	std::vector<int> mdl_ref(mdlcnt,-1);
	std::map<std::string,int> mdl_idx;
	get_cls(org_cls,input_cls,dw_ans);
	e_lim.resize(mdlcnt);
	m_lim.resize(mdlcnt);
	for(auto &v:m_lim){
		v.resize(mdlcnt);
	}
	//---get-mdl---
	for(int i=0; i<mdlcnt; ++i){
		std::cin>>mdl_name[i];
		mdl_idx[mdl_name[i]]=i;
		vec v=vec::get(),tgt;
		std::string str;
		std::cin>>str;
		if(str[0]=='('){ // 连接点为固定坐标
			str+=get_line({});
			replace_with(str,{'(',')',','},' ');
			std::istringstream is(str);
			double x,y;
			is>>x>>y;
			tgt=(vec){x,y};
			org_mdl[i]=mdl::build(tgt,v);
			//更新坐标范围
			dw_ans.upd(org_mdl[i].v[0]);
			dw_ans.upd(org_mdl[i].v[1]);
		}else{ // 连接点在模块上
			ref_name[i]=str;
			org_mdl[i]={v,vec()};
		}
	}
	for(int i=0; i<mdlcnt; ++i){
		if(ref_name[i]!=""){
			mdl_ref[i]=mdl_idx[ref_name[i]];
		}
	}
	//---get-lim---
	while(limcnt--){
		// 输入限制并进行类型判断
		static vec dv[]={{1,0},{0,1},{-1,0},{0,-1}};
		std::istringstream is(get_line({'(',')',',','V','H'}));
		std::string str1,str2;
		double dis;
		is>>str1>>str2>>dis;
		if(isdigit(str1[0])^isdigit(str2[0])){ // 边对模块的间隔要求
			if(isdigit(str1[0])) std::swap(str1,str2);
			int m_id=mdl_idx[str1],e_id;
			sscanf(str2.c_str(),"%d",&e_id); // TODO 与标准接轨
			int dr=e_vec[e_id].dr();
			e_lim[m_id].push_back(e_vec[e_id]+(dv[(dr+1)&3])*dis); // 对边
			// 邻边 有没有办法简化？
			if(dr&1){ // 竖向
				// 强制定义边的方向, 下else同理
				double low_y=e_vec[e_id].a.y,high_y=e_vec[e_id].b.y;
				double low_x=e_vec[e_id].a.x,high_x=(e_vec[e_id].a+dv[(dr+1)&3]).x;
				inc_swp(low_y,high_y);
				inc_swp(low_x,high_x);
				e_lim[m_id].push_back({{high_x,low_y},{low_x,low_y}});
				e_lim[m_id].push_back({{low_x,high_y},{high_x,high_y}});
			}else{ // 横向
				double low_x=e_vec[e_id].a.x,high_x=e_vec[e_id].b.x;
				double low_y=e_vec[e_id].a.y,high_y=(e_vec[e_id].a+dv[(dr+1)&3]).y;
				inc_swp(low_y,high_y);
				inc_swp(low_x,high_x);
				e_lim[m_id].push_back({{low_x,low_y},{low_x,high_y}});
				e_lim[m_id].push_back({{high_x,high_y},{high_x,low_y}});
			}
		}else if(isdigit(str1[0])){ // 模块对模块的间隔要求
			int id1=mdl_idx[str1],id2=mdl_idx[str2];
			m_lim[id1][id2]=m_lim[id2][id1]=dis;
		}else{
			std::cerr<<"[Error] limits between edge!"<<std::endl;
			return 1;
		}
	}
	//---get-ans---
	std::vector<int> idx(mdlcnt);
	topo_rand(idx,mdl_ref);
	std::vector<mdl> res_mdl=solve_seq(org_cls,idx,org_mdl,mdl_ref,e_vec,mdlcnt);
	double res_len=calc_res(res_mdl,org_mdl,mdl_ref);
	for(int times=20; times--; ){
		topo_rand(idx,mdl_ref);
		std::vector<mdl> cur_mdl=solve_seq(org_cls,idx,org_mdl,mdl_ref,e_vec,mdlcnt);
		double cur_len=calc_res(cur_mdl,org_mdl,mdl_ref);
		if(res_len>cur_len){ // 当前结果优于最优结果
			std::cerr<<"better."<<std::endl; // debug
			res_len=cur_len;
			res_mdl.swap(cur_mdl);
		}
	}
	//---oput---
	dw_ans.zoom_out();
	dw_ans.draw_grid();
	// 画模块
	for(int i=0; i<mdlcnt; ++i){
		if(res_mdl[i].v[0].x<=-inf){
			std::cerr<<"Cannot put "<<i+1<<'.'<<std::endl;
			continue;
		}
		std::cerr<<i<<": "<<mdl_name[i]<<','<<res_mdl[i].v[0]<<' '<<res_mdl[i].v[1]<<std::endl;
		dw_ans.draw_mdl(res_mdl[i],col_grey,mdl_name[i]);
	}
	// 画模块的连接线
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
