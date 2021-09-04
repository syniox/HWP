// 给定一个模块接入顺序，寻找最优方案的主流程

#pragma once

#include <vector>

#include "types.h"

std::vector<mdl> solve_seq(std::vector<cls_s> clss,const std::vector<int> &seq,
		const std::vector<mdl> &mds,const std::vector<int> &mdl_ref,
		std::vector<edg> e_mdl_org,const int mdlcnt);
// O(n*e*e) 对一个给定的模块摆放优先顺序做贪心最优解
