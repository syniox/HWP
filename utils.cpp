// 与内容无关的基础代码 + 全局变量

#include <iostream>
#include <cstdio>
#include <cassert>

#include "utils.h"

const double eps=1e-6;
const double inf=1e12;

std::vector<std::vector<double>> m_lim; // 模块间的限制
std::vector<std::vector<edg>> e_lim; // 模块与边的限制

int rnd(const int l,const int r){
	return rand()%(r-l+1)+l;
}

std::string get_line(std::vector<char> repl){
	// 读入一行非汉字字符串并替换部分字符为空格
	std::string str;
	for(;;){
		std::getline(std::cin,str);
		if(!str.length()) continue;
		bool f=0;
		for(char c:str){
			f|=c<0;
		}
		if(!f) break;
	}
	replace_with(str,repl,' ');
	return str;
}

void sanitize_vec(cls_s &cl){
	// O(e(cl)*e(cl)) 删除长度为0的边，合并相邻且方向相反的边
	// TODO: optimize(SPFA?)
	for(int cnt=-1; cnt; ){
		cnt=0;
		// 删除长度为0的边
		for(int i=cl.size()-1; i>=0; --i){
			if(cl[i].ispnt()){
				++cnt;
				cl.erase(cl.begin()+i);
			}
		}
		// 合并一次方向相反的边
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

double calc_res(std::vector<mdl> m_res,std::vector<mdl> m_org,std::vector<int> ref){
	// O(n)  计算一组方案的连线长度和
	if(m_res.size()!=m_org.size()) return inf;
	int sz=m_res.size();
	double res=0;
	for(int i=0; i<sz; ++i){
		vec v;
		if(ref[i]==-1) v=m_res[i].cntr()-m_org[i].cntr();
		else v=m_res[i].cntr()-m_res[ref[i]].cntr();
		res+=cabs(v.x)+cabs(v.y);
	}
	return res;
}

void add_bevel(cls_s &cl,const vec a,const vec b){
	vec c; // 把斜边拆分成横边和竖边
	if(a.x<b.x&&a.y<b.y) c=(vec){a.x,b.y};
	else if(a.x>b.x&&a.y<b.y) c=(vec){b.x,a.y};
	else if(a.x>b.x&&a.y>b.y) c=(vec){a.x,b.y};
	else if(a.x<b.x&&a.y>b.y) c=(vec){b.x,a.y};
	else assert(0);
	cl.push_back({a,c}); // 取消斜边计数?
	cl.push_back({c,b});
}
