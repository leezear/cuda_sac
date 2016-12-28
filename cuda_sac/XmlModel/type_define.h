/*
* type_define.h
*
*  Created on: 2016年6月21日
*      Author: leezear
*/

#ifndef TYPE_DEFINE_H_
#define TYPE_DEFINE_H_

#include <vector>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <string>
#include <string.h>
#include <math.h>
#include <iostream>

typedef unsigned int u32;

using namespace std;

namespace cp {
namespace model {

const int MAX_VALUE = INT_MAX - 4096;
const int MAX_OPT = INT_MIN + 4096;

/**
*  @brief  启发式类型
*/
enum Heuristics
{
	HEU_LIFO,		///<无启发式
	HEU_DOM_WDEG	///<dom/wdeg
};

/*!
*\brief 整型元组类型
*/
typedef vector<int> IntTuple;
//typedef
//typedef tuple<int, int> IntTuple;
/**
* \brief 表约束语义
*/
enum Semantices {
	SEM_SUPPORT,		///<支持
	SEM_CONFLICT		///<冲突
};

/**
* \brief 表约束语义
*/
enum PropagateProblem {
	MAIN_PROBLEM,		///<main problem
	SUB_PROBLEM		///<sub problem
};

/**
* \brief SAC consistent test
*/
enum SacConsistentType {
	SAC_CONSISTENT = 1, SAC_INCONSISTENT = -1, SAC_UNTEST = 0
};

/**
* \brief 关系类型
*/
//enum RelationType {
//	REL_TABULAR,			///<表约束关系
//	REL_PREDICATION		///<表达式关系
//};
enum RelationType
{
	REL_TABULAR, REL_PREDICATION
};
/**
* \brief 关系维度
*/
enum RelationArity {
	REL_BINARY,				///<二维
	REL_NONBINARY			///<多维
};

/**
* \brief some elements in prediction expression
*/
enum PredictionElement {
	PE_OP, PE_VAL, PE_VAR
};

//const int SAC_CONSISTENT = -1024;
/**
* \brief Features of Constraints Network
*/
struct NetworkFeatures {
	u32 max_dom_size;
	u32 vars_size;
	u32 cons_size;
	u32 arity = 2;
	RelationType rel_type;
};

enum PredicateOperator {
	PO_NON = INT_MIN, 		///<空		-2147483648
	PO_L_PAR = INT_MIN + 1,	///<左括号	-2147483647
	PO_R_PAR = INT_MIN + 2,	///<右括号	-2147483646
	PO_COMMA = INT_MIN + 3,	///<逗号		-2147483645
	PO_ADD = INT_MIN + 10,	///<加法		-2147483638
	PO_SUB = INT_MIN + 11,	///<减法		-2147483637
	PO_MUL = INT_MIN + 12,	///<乘法		-2147483636
	PO_DIV = INT_MIN + 13,	///<除法		-2147483635
	PO_MOD = INT_MIN + 14,	///<取余		-2147483634
	PO_EQ = INT_MIN + 30,	///<"="		-2147483618
	PO_NE = INT_MIN + 31,	///<"!="	-2147483617
	PO_LT = INT_MIN + 32,	///<"<"		-2147483616
	PO_LE = INT_MIN + 33,	///<"<="	-2147483615
	PO_GT = INT_MIN + 34,	///<">"		-2147483614
	PO_GE = INT_MIN + 35,	///<">="	-2147483613
	PO_ABS = INT_MIN + 60,	///<取模		-2147483588
	PO_AND = INT_MIN + 101,	///<逻辑与	-2147483547
	PO_OR = INT_MIN + 102	///<逻辑或	-2147483546
};

enum OpType {
	OT_PAR = 0, OT__SIG_OP = 1, OT_BI_OP = 2,
};

/**
* \brief Max number
*/
#define MAX(x,y) (x)>(y)?(x):(y)

/**
* \brief Get the max number of blocks
*/

/// u32,u32 tuple
struct bituple {
	u32 x;
	u32 y;

	bool operator==(const bituple &rhs) {
		return (x == rhs.x) && (y == rhs.y);
	}
};

//tuple<u32, u32> BiU32Tuple;

/**
* \brief propagate stack
*/
struct PStack {
	int ps[200];
	u32 num_prs = 0;
	int top = -1;
	int bottom = 0;
};

inline int GetTopNum(int num_elements, int num_threads) {
	return (num_elements + (num_threads - 1)) / num_threads;
}

inline int func_or(int *x) {
	return x[0] || x[1];
}

inline int func_and(int *x) {
	return x[0] && x[1];
}

inline int func_eq(int *x) {
	return x[0] == x[1];
}

inline int func_abs(int *x) {
	return (int)fabsf(x[0]);
}

inline int func_sub(int *x) {
	return x[0] - x[1];
}

inline int func_div(int *x) {
	return x[0] / x[1];
}

inline int func_mod(int *x) {
	return x[0] % x[1];
}

inline int func_ne(int *x) {
	return x[0] != x[1];
}

inline int func_add(int *x) {
	return x[0] + x[1];
}

inline int func_mul(int *x) {
	return x[0] + x[1];
}

inline int func_lt(int *x) {
	return x[0] < x[1];
}

inline int func_le(int *x) {
	return x[0] <= x[1];
}

inline int func_gt(int *x) {
	return x[0] > x[1];
}

inline int func_ge(int *x) {
	return x[0] >= x[1];
}

}
}

#endif /* TYPE_DEFINE_H_ */