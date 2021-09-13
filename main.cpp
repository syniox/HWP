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

void get_cls(std::string str,std::vector<cls_s> &clss,std::vector<cls_s> &input,drawer &dw_ans,double thr=4){
	// O(e) 给一个闭包串，进行存储
	if(!validstr(str)) return;
	replace_with(str,{'[',']','(',')',','},' ');
	clss.push_back(cls_s());
	input.push_back(cls_s());
	int i=clss.size()-1;
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
	// 输入格式见readme
	// 未能处理距离限制
	drawer dw_ans("oput.png");

	int mdlcnt=0;
	std::vector<cls_s> org_cls,input_cls;
	std::vector<mdl> org_mdl;
	std::vector<std::string> mdl_name,ref_name;
	std::vector<int> mdl_ref;
	std::map<std::string,int> mdl_idx;
	std::string inbuf;
	// 输入闭合回路 TODO 判断合法性
	while(std::getline(std::cin,inbuf)&&inbuf.find("器件")==std::string::npos){
		get_cls(inbuf,org_cls,input_cls,dw_ans);
	}
	//输入模块
	while(std::getline(std::cin,inbuf)){
		if(!validstr(inbuf)) continue;
		int i=mdlcnt++;
		e_lim.push_back(std::vector<edg>());
		m_lim.push_back(std::vector<double>(mdlcnt));
		mdl_name.push_back(std::string());
		ref_name.push_back(std::string());
		org_mdl.push_back({vec(),vec()});
		for(auto &v:m_lim){
			v.push_back(0);
		}
		std::istringstream is(inbuf);
		is>>mdl_name[i];
		vec v,tgt;
		is>>v.x>>v.y;
		std::string str;
		std::getline(is,str);
		for(; str[0]==' '; str.erase(0,0)){
			assert(!str.empty());
		}
		if(str[0]=='('){ // 固定坐标
			replace_with(str,{'(',')',','},' ');
			std::istringstream istr(str);
			double x,y;
			istr>>x>>y;
			tgt=(vec){x,y};
			org_mdl[i]=mdl::build(tgt,v);
			//更新坐标范围
			dw_ans.upd(org_mdl[i].v[0]);
			dw_ans.upd(org_mdl[i].v[1]);
		}else{ // 模块上
			ref_name[i]=str;
			org_mdl[i]={v,vec()};
		}
	}
	for(int i=0; i<mdlcnt; ++i){
		if(ref_name[i]!=""){
			mdl_ref[i]=mdl_idx[ref_name[i]];
		}
	}
	// ---输入模块的距离限制---
	/*
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
	*/
	// ---寻找排列并计算方案---
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
	// ---输出方案---
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
