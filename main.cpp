#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <map>
static const int eps=1e-6;
using std::cin; using std::cout; using std::cerr; using std::endl;
using std::min; using std::max;
using std::vector;
using std::map;
// x,y: 平面直角坐标系

struct vec{ // 向量
	double x,y;
	static vec get();
	inline void flip();
};
struct edg{ // 边，a为起点，b为终点，合法区域在这个边向量的左边
	vec a,b;
	inline int dr();
	inline void flip();
};
struct mdl{ // module, 记录该模块长方形的四个顶点，保证连续
	vec v[4];
	inline vec cntr();//返回中心位置
	inline void flip();
	inline void set_inf();
	static mdl build(const vec &ctr,const vec &rct);
};

using cls_s=vector<edg>;

vector<cls_s> clss; // 闭包集合

vec vec::get(){
	double x,y;
	cin>>x>>y;
	return (vec){x,y};
}
inline void vec::flip(){
	std::swap(x,y);
}
template <typename T>
vec operator * (const vec &v,const T x){
	return (vec){v.x*x,v.y*x};
}
double operator *(const vec &a,const vec &b){
	return a.x*b.y-a.y*b.x;
}
vec operator + (const vec &a,const vec &b){
	return (vec){a.x+b.x,a.y+b.y};
}
vec operator - (const vec &a,const vec &b){
	return (vec){a.x-b.x,a.y-b.y};
}
bool operator < (const vec &a,const vec &b){
	return a.x==b.x?a.y<b.y:a.x<b.x;
}
bool operator == (const vec &a,const vec &b){
	return a.x==b.x&&a.y==b.y;
}
std::ostream& operator << (std::ostream &out,const vec &v){
	out<<'('<<v.x<<' '<<v.y<<')';
	return out;
}
std::ostream& operator << (std::ostream &out,const edg &e){
	out<<e.a<<"->"<<e.b;
	return out;
}

bool crossed(const edg &a,const edg &b){ // 线段是否相交，重合情况不考虑
	return ((b.a-a.a)*(a.b-a.a))*((b.b-a.a)*(a.b-a.a))<-eps;
}
int edg::dr(){ // 返回向量方向，右0上1左2下3
	assert(a.x==b.x||a.y==b.y);
	if(a.y==b.y) return (a.x>b.x)<<1;
	else return (a.x>b.x)<<1|1;
}
inline void edg::flip(){
	a.flip(),b.flip();
}

inline vec mdl::cntr(){
	return (v[0]+v[1]+v[2]+v[3])*0.25;
}
inline void mdl::flip(){
	for(int i=0; i<4; ++i){
		std::swap(v[i].x,v[i].y);
	}
}
inline void mdl::set_inf(){
	for(int i=0; i<4; ++i){
		v[i].x=v[i].y=-1e12;
	}
}
inline mdl mdl::build(const vec &ctr,const vec &rct){
	mdl m;
	m.v[0]=(vec){ctr.x-rct.x,ctr.y-rct.y};
	m.v[1]=(vec){ctr.x-rct.x,ctr.y+rct.y};
	m.v[2]=(vec){ctr.x+rct.x,ctr.y+rct.y};
	m.v[3]=(vec){ctr.x+rct.x,ctr.y-rct.y};
	return m;
}


void add_edg(vector<edg> &eg,map<vec,std::vector<int>> &vec_idx,vec a,vec b){
	eg.push_back((edg){a,b});
	vec_idx[a].push_back(eg.size()-1);
}

void get_cls(const int edgcnt){ // 根据题目给出的边构建rectilinear block
	vector<edg> eg;
	map<vec,vector<int>> vec_idx;
	vector<bool> vis;
	for(int i=1,dir; i<=edgcnt; ++i){
		vec a=vec::get(),b=vec::get();
		cin>>dir;
		if(dir) std::swap(a,b);
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
			//cerr<<eg[id]<<endl;
			vis[id]=1;
			for(int j:vec_idx[eg[id].b]){ // WIP
				//cerr<<eg[j]<<endl;
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

template <typename T>
const T cabs(const T &x){
	return x<0?-x:x;
}
double get_dis(const vec &a,const vec &b){
	return cabs(a.x-b.x)+cabs(a.y-b.y);
}
template <typename T>
void flip_vec(vector <T> vt){
	for(T &x:vt) x.flip();
}

void update_gpos(mdl &gpos,double l,double r,const vec rct,const vec tgt,const double line_y){
	if(r-l<rct.x) return;
	double x,rad_x=rct.x*0.5;
	if(r<tgt.x-rad_x) x=r-rad_x;
	else if(l>tgt.x+rad_x) x=l+rad_x;
	else x=tgt.x;
	mdl neo=mdl::build((vec){x,line_y},rct);
	if(get_dis(neo.cntr(),tgt)<get_dis(gpos.cntr(),tgt)){
		gpos=neo;
	}
}

mdl get_great_pos_basic(const vec rct,cls_s &cl,const vec tgt,const bool fliped){//横向，坐标系是否经过变换
	using pdd=std::pair<double,double>;
	vector<pdd> invalid_seg,neo_seg;
	mdl gpos; // great pos
	for(edg cur_e:cl){
		if(cur_e.dr()&1) continue; // 保证横向
		invalid_seg.clear();
		neo_seg.clear();
		for(edg e:cl){
			if((e.dr()&1)||cur_e.dr()==e.dr()) continue;
			if(cabs(cur_e.a.y-e.a.y)<rct.y){ // 忽略刚好重合的情况
				int a=e.a.x,b=e.b.x;
				if(a>b) std::swap(a,b);
				invalid_seg.push_back(std::make_pair(a,b));
			}
		}
		{
			double end=-1e18,start;
			for(pdd pr:invalid_seg){
				if(end<pr.first){
					if(end>-1e18) neo_seg.push_back(std::make_pair(start,end));
					start=pr.first,end=pr.second;
				}else{
					end=std::max(end,pr.second);
				}
			}
		}
		vector<pdd>::iterator it1=neo_seg.begin(),it2=it1;
		++it2;
		double line_y=cur_e.a.y+rct.y*(0.5-((cur_e.dr()==2)^fliped));
		double a=min(cur_e.a.x,cur_e.b.x);
		double b=max(cur_e.a.x,cur_e.b.x);
		update_gpos(gpos,a,min(b,it1->first),rct,tgt,line_y);
		for(; it2!=neo_seg.end(); ++it1,++it2){
			double l=max(a,it1->second);
			double r=min(b,it2->first);
			update_gpos(gpos,l,r,rct,tgt,line_y);
		}
	}
	return gpos;
}

mdl get_great_pos(vec rct,cls_s &cl,vec tgt){// 寻找某个闭包的最优位置
	// cl表示搜寻的闭包
	// 函数返回模块最后占用的位置
	mdl mpos[4];
	mpos[0]=get_great_pos_basic(rct,cl,tgt,0); // 横着的原矩阵 横向rb
	rct.flip();
	mpos[1]=get_great_pos_basic(rct,cl,tgt,0); // 竖着的原矩阵 横向rb
	flip_vec(cl);
	mpos[2]=get_great_pos_basic(rct,cl,tgt,1); // 横着的原矩阵 竖向rb（坐标系颠倒）
	rct.flip();
	mpos[3]=get_great_pos_basic(rct,cl,tgt,1); // 竖着的原矩阵 竖向rb（坐标系颠倒）
	flip_vec(cl);
	mpos[2].flip(),mpos[3].flip();
	mdl pos=mpos[0];
	double res=get_dis(mpos[0].cntr(),tgt);
	for(int i=1; i<4; ++i){
		double tmp=get_dis(mpos[i].cntr(),tgt);
		if(tmp<res) res=tmp,pos=mpos[i];
	}
	return pos;
}

bool on_edge(const edg &e,const vec &p){
	if(e.a.x==e.b.x) return cabs(p.x-e.a.x)<eps;
	return cabs(p.y-e.a.y)<eps;
}
bool on_edge(const edg &e,const mdl &m){
	int cnt=0;
	for(int i=0; i<4; ++i){
		cnt+=on_edge(e,m.v[i]);
	}
	assert(cnt==0||cnt==2);
	return cnt==2;
}

void insert_mdl(cls_s &cl,mdl md){// 将该区域设为不可用区域（假设该模块紧贴边缘）
	using cls_i=cls_s::iterator;
	for(cls_i it=cl.begin(); it!=cl.end(); ++it){
		edg e=*it;
		if(!on_edge(e,md)) continue;
		bool fliped=0;
		if(e.dr()&1){
			fliped=1;
			flip_vec(cl);
			md.flip();
			e.flip();
		}
		double st_x,ed_x,other_y;
		for(int i=0; i<4; ++i){
			if(cabs(md.v[i].y-e.a.y)>eps){
				other_y=md.v[i].y;
			}
			if(md.v[i].x-e.a.x<md.v[i^3].x-e.a.x){
				st_x=md.v[i].x;
				ed_x=md.v[i^3].x;
			}
		}
		vec p0=(vec){st_x,e.a.y},p1=(vec){st_x,other_y};
		vec p2=(vec){ed_x,other_y},p3=(vec){ed_x,e.a.y};
		it->b=p0;
		++it;
		it=cl.insert(it,(edg){p0,p1});
		it=cl.insert(it,(edg){p1,p2});
		it=cl.insert(it,(edg){p2,p3});
		it=cl.insert(it,(edg){p3,e.b});
		if(fliped){
			flip_vec(cl);
		}
		return;
	}
	assert(0);
	/*
	vector <edg*> bkg; // TODO safe enough?
	for(edg &e:cl){
		if(on_edge(e,md)) bkg.push_back(&e);
	}
	assert(bkg.size()>0&&bkg.size()<=4);
	if(bkg.size()==1){
		edg &e=**(bkg.begin());
	}else if(bkg.size()==2){
		edg &e1=**(bkg.begin()),&e2=**(bkg.begin());
	}else{
		cerr<<"unimplemented"<<endl;
		exit(1);
	}
	*/
	/*
	int sz=cl.size(),l=-1,r=-1;
	for(int i=0; i<sz; ++i){
		if(on_edge(cl[i],md)){
			assert(l==-1);
			l=i;
			for(r=l+1; r<sz&&on_edge(cl[r],md); ++r);
			i=r--;
		}
	}
	assert(l!=-1);
	if(l==r){
		cls_i it=cl.begin()+l;
	}else if(r==l+1){

	}else{
		cerr<<"unimplemented"<<endl;
		exit(1);
	}
	*/
}

int main(){
	// 输入格式：
	// 第一行输入空白区域数量n和模块数量m
	// 接下来n行，每行第一个数是e，表示该空白区域的边界点数；接下来2e个数，依次表示边界上逆时针顺序的点的坐标；
	// 接下来m行，每行4个数w, h, x, y，表示矩形的最佳位置是(x,y)到(x+w, y+h)画出的矩形。
	// 输入：
	// 第一行输入有多少条边界e和多少个模块m
	// 接下来e行输入每条边（向量）的起止坐标x1,y1,x2,y2和空白区域的位置在向量的左边还是右边，左边为0右边为1
	// 接下来m行每行输入4个数，代表该模块的长和宽，最优位置中心的x坐标和y坐标
	// 输出：
	// 共m行，每行输出该模块摆放位置的对角端点
	int edgcnt,mdlcnt;
	cin>>edgcnt>>mdlcnt;
	get_cls(edgcnt);
	/*for(cls_s cl:clss){
		for(edg e:cl){
			cout<<e.a<<"->"<<e.b<<endl;
		}
		cout<<"---"<<endl;
	}*/
	for(int i=1; i<=mdlcnt; ++i){
		vec v=vec::get(),tgt=vec::get();//返回最优位置的中心？
		double res=1e18;
		cls_s *best_cl;
		mdl mpos;
		for(cls_s &cl:clss){
			mdl cur=get_great_pos(v,cl,tgt);
			double tmp_res=get_dis(cur.cntr(),tgt);
			if(tmp_res<res){
				best_cl=&cl;
				mpos=cur;
			}
		}
		assert(res<1e18);
		insert_mdl(*best_cl,mpos);
		cout<<i<<": "<<mpos.v[0]<<' '<<mpos.v[2]<<endl;
	}
	return 0;
}
