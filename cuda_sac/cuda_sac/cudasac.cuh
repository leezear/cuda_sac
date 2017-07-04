/*
* cudasac.cuh
*
*  Created on: 2017年5月4日
*      Author: leezear
*/

#ifndef CUDASAC_CUH_
#define CUDASAC_CUH_
//CUDA runtimes
//#define __CUDA_INTERNAL_COMPILATION__
//#include "math_functions.h"
//#undef __CUDA_INTERNAL_COMPILATION__
// 定义无符号整型掩码
#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_launch_parameters.h>
// Utilities and system includes
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <iostream>
#include <math.h>
#include "XModel.h"
#include "math_functions.hpp"
#include "cuda_runtime_api.h"

using namespace cudacp;

__constant__ __device__ u32 D_U32_MASK1[32] = { 0x80000000, 0x40000000, 0x20000000,
0x10000000, 0x08000000, 0x04000000, 0x02000000, 0x01000000, 0x00800000,
0x00400000, 0x00200000, 0x00100000, 0x00080000, 0x00040000, 0x00020000,
0x00010000, 0x00008000, 0x00004000, 0x00002000, 0x00001000, 0x00000800,
0x00000400, 0x00000200, 0x00000100, 0x00000080, 0x00000040, 0x00000020,
0x00000010, 0x00000008, 0x00000004, 0x00000002, 0x00000001 };

__constant__ __device__ u32 D_U32_MASK0[32] = { 0x7FFFFFFF, 0xBFFFFFFF, 0xDFFFFFFF,
0xEFFFFFFF, 0xF7FFFFFF, 0xFBFFFFFF, 0xFDFFFFFF, 0xFEFFFFFF, 0xFF7FFFFF,
0xFFBFFFFF, 0xFFDFFFFF, 0xFFEFFFFF, 0xFFF7FFFF, 0xFFFBFFFF, 0xFFFDFFFF,
0xFFFEFFFF, 0xFFFF7FFF, 0xFFFFBFFF, 0xFFFFDFFF, 0xFFFFEFFF, 0xFFFFF7FF,
0xFFFFFBFF, 0xFFFFFDFF, 0xFFFFFEFF, 0xFFFFFF7F, 0xFFFFFFBF, 0xFFFFFFDF,
0xFFFFFFEF, 0xFFFFFFF7, 0xFFFFFFFB, 0xFFFFFFFD, 0xFFFFFFFE };

const u32 U32_MASK1[32] = { 0x80000000, 0x40000000, 0x20000000, 0x10000000,
0x08000000, 0x04000000, 0x02000000, 0x01000000, 0x00800000, 0x00400000,
0x00200000, 0x00100000, 0x00080000, 0x00040000, 0x00020000, 0x00010000,
0x00008000, 0x00004000, 0x00002000, 0x00001000, 0x00000800, 0x00000400,
0x00000200, 0x00000100, 0x00000080, 0x00000040, 0x00000020, 0x00000010,
0x00000008, 0x00000004, 0x00000002, 0x00000001 };
const u32 U32_MASK0[32] = { 0x7FFFFFFF, 0xBFFFFFFF, 0xDFFFFFFF, 0xEFFFFFFF,
0xF7FFFFFF, 0xFBFFFFFF, 0xFDFFFFFF, 0xFEFFFFFF, 0xFF7FFFFF, 0xFFBFFFFF,
0xFFDFFFFF, 0xFFEFFFFF, 0xFFF7FFFF, 0xFFFBFFFF, 0xFFFDFFFF, 0xFFFEFFFF,
0xFFFF7FFF, 0xFFFFBFFF, 0xFFFFDFFF, 0xFFFFEFFF, 0xFFFFF7FF, 0xFFFFFBFF,
0xFFFFFDFF, 0xFFFFFEFF, 0xFFFFFF7F, 0xFFFFFFBF, 0xFFFFFFDF, 0xFFFFFFEF,
0xFFFFFFF7, 0xFFFFFFFB, 0xFFFFFFFD, 0xFFFFFFFE };
//////定义一些常量//////
static const int U32_SIZE = sizeof(u32); 	// =4
static const int U32_BIT = U32_SIZE * 8; 	// =32
static const int U32_POS = 5;				//位数
static const int U32_MOD_MASK = 31;
static const int WORKSIZE = 32;
static const int WORKSIZE_LARGE = 256;
static const int WARPSIZE = 32;
//// 一个变量论域占几个32位
//int H_BITDOM_INTSIZE;
//__constant__ int D_BITDOM_INTSIZE;
//论域个数
int H_DS_SIZE;
__constant__ int D_DS_SIZE;
//变量个数
int H_VS_SIZE;
__constant__ int D_VS_SIZE;
//约束个数
int H_CS_SIZE;
__constant__ int D_CS_SIZE;
//最大论域变量的模
int H_MDS;
__constant__ int D_MDS;
//论域的int大小
int H_BD_INTSIZE;
__constant__ int D_BD_INTSIZE;
// 一个约束的bitsup的int长度
int H_BITSUP_INTSIZE;
__constant__ int D_BITSUP_INTSIZE;

///////////定义一些全局变量数组////////////

// 除法向上取整
__inline__ __device__ __host__ int GetTopNum(int num_elements,
	int num_threads) {
	return (int)ceil((float)num_elements / num_threads);
}

#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

#ifndef GetOffSet
#define GetOffSet(x)(U32_BIT - (x & U32_MOD_MASK))
#endif

#ifndef PROPFAILED
#define PROPFAILED -1
#endif

unsigned int nextPow2(unsigned int x) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

inline int GetBitSubDomStartIndexHost(int x, int a) {
	return (x * H_MDS + a) * H_VS_SIZE;
}

inline int GetBitSubDomIndexHost(int x, int a, int y) {
	return (x * H_MDS + a) * H_VS_SIZE + y;
}

inline int2 GetBitSupIndexByTupleHost(int cid, int2 t) {
	return make_int2(cid * H_BITSUP_INTSIZE + t.x, cid * H_BITSUP_INTSIZE + t.y);
}

__device__ int2 GetBitSupIndexByTupleDevice(int cid, int2 t) {
	return make_int2(cid * D_BITSUP_INTSIZE + t.x, cid * D_BITSUP_INTSIZE + t.y);
}

static __inline__ __device__ int GetBitSupIdxDevice(int cid, int a) {
	return cid * D_BITSUP_INTSIZE + a;
}

static __inline__ __device__ int GetBitSubDomIdxDevice(int x, int a, int y) {
	return (x * D_MDS + a) * D_VS_SIZE + y;
}

static __inline__ __device__ int GetSVarPreIdxDevice(int x, int a, int y) {
	return (x * D_MDS + a) * D_VS_SIZE + y;
}

inline int GetSVarPreIdxHost(int x, int a, int y) {
	return (x * H_MDS + a) * H_VS_SIZE + y;
}

static __inline__ __device__ void DelValDevice(u32* bitDom, int* mVarPre, int x, int a) {
	u32 les = bitDom[x] & D_U32_MASK0[a];
	if (bitDom[x] != les) {
		atomicAnd(&bitDom[x], D_U32_MASK0[a]);
		mVarPre[x] = 1;

		if (bitDom[x] == 0)
			mVarPre[x] = INT_MIN;
	}
}

// 主问题论域
u32* h_bitDom;
u32* d_bitDom;
// 子问题论域
u32* h_bitSubDom;
u32* d_bitSubDom;
// bit支持
uint2* h_bitSup;
uint2* d_bitSup;
// 标记主问题变量域是否删减，初始化全部为1
int* h_MVarPre;
int* d_MVarPre;

// 记录每个变量的论域长度
int* h_var_size;
int* d_var_size;

//	存储约束的scope，
// 	类型int3
// 	scope.x: x.id
//	scope.y: y.id
//	scope.z: c.id
int3* h_scope;
int3* d_scope;

//标记主问题约束是否需检查，初始化全部为1
int* h_MConPre;
int* d_MConPre;
//主问题约束传播队列(压缩版)
int* h_MConEvt;
int* d_MConEvt;
//主问题约束传播队列
int* h_MCon;
int* d_MCon;

//标记子问题变量域是否删减，初始化全部为1
int* h_SVarPre;
int* d_SVarPre;
//记录子问题变量域
int3* h_SVar;
int3* d_SVar;
//标记子问题约束是否需检查，初始化全部为1
int* h_SConPre;
int* d_SConPre;
//子问题约束传播队列(压缩版)
int3* h_SConEvt;
int3* d_SConEvt;
//子问题约束传播队列
int3* h_SCon;
int3* d_SCon;

//主问题约束在GPU上执行的规格
int H_MCCount;
int D_MCCount;

//主问题变量在GPU上执行的规格
int H_MVCount;
int D_MVCount;

__device__ __inline__ int pow2i(int e) {
	return 1 << e;
}

__global__ void showConstraints(int3* scope, int len) {
	const int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < len) {
		printf("i: %d, %d, %d\n", idx, scope[idx].x, scope[idx].y,
			scope[idx].z);
	}
}

__global__ void showVariables(u32* bitDom, int* MVarPre, int* var_size,
	int len) {
	const int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < len) {
		printf("i: %d, %x, %d, %d\n", idx, bitDom[idx], MVarPre[idx],
			__popc(bitDom[idx]));
	}
}

__global__ void ShowSubConEvt(int3* ConEvt)
{
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	__shared__ int3 s_cevt;

	//if (tid == 0)
	//{
	//	if()
	//	s_cevt = ConEvt[bid];
	//	printf("(%d, %d, %d)\n", s_cevt.x, s_cevt.y, s_cevt.z);
	//}
}

__global__ void showSubVar(u32* bitSubDom, int* var_size, int vs_size) {
	const int val = blockIdx.x;
	const int var = blockIdx.y;
	const int subVar = threadIdx.x;

	if (val < var_size[var]) {
		if (subVar < vs_size) {
			const int idx = GetBitSubDomIdxDevice(var, val, subVar);
			//if (var == 7 && val == 11)
			//	printf("var_size[%d] = %d, (%d, %d, %d): %x, idx = %d\n", var,
			//		var_size[var], var, val, subVar, bitSubDom[idx], idx);
		}
	}
}

//通过已改变变量生成约束队列
__global__ void GenConPre(int *VarPre, int* BlocksCount, int3* scp, int len) {
	//拿到约束id
	const int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < len) {
		int3 sp = scp[idx];
		__shared__ int failed;
		int pred;
		//		printf("c_id = %d: (%d = %d, %d = %d)\n", sp.z, sp.x, VarPre[sp.x],
		//				sp.y, VarPre[sp.y]);
		if (threadIdx.x == 0)
			failed = 0;
		__syncthreads();

		const int2 v_pre = make_int2(VarPre[sp.x], VarPre[sp.y]);
		if (v_pre.x == 0 && v_pre.y == 0)
			pred = 0;
		else if (v_pre.x == INT32_MIN || v_pre.y == INT32_MIN)
			failed = INT32_MIN;
		else
			pred = 1;
		__syncthreads();

		int BC = (failed == INT32_MIN) ? INT32_MIN : __syncthreads_count(pred);

		if (threadIdx.x == 0)
			BlocksCount[blockIdx.x] = BC;
	}
}

__global__ void CompactQ(int *VarPre, int* ConEvt, int* BOffset, int3* scp,
	int len) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	extern __shared__ int warpTotals[];
	//一个线程块内有128个线程
	//一个块内有4个线程束
	if (idx < len) {
		int3 sp = scp[idx];
		//获得判定
		const int2 v_pre = make_int2(VarPre[sp.x], VarPre[sp.y]);
		int pred = v_pre.x || v_pre.y;

		//warp index
		//线程束索引
		int w_i = threadIdx.x / WARPSIZE;
		//thread index within a warp
		//线程束内线程索引
		int w_l = idx % WARPSIZE;
		//thread mask (ERROR IN THE PAPERminus one is required)
		//线程掩码
		//INT_MAX = 1111 1111 1111 1111 1111 1111 1111 1111
		//若线程内id=0，右移32-0-1 = 31位 右侧剩下1位
		//若线程内id=5，右移32-5-1 = 26位 右侧剩下6位
		//若线程内id=31，右移32-31-1 = 0位 右侧剩下32位
		//线程束内threid|  31~~~~~~0
		//ballot对应位置|   1......1
		int t_m = INT_MAX >> (WARPSIZE - w_l - 1);
		//balres = number whose ith bit is one if the ith's thread pred is true masked up to the current index in warp
		//线程内局部变量pred = 1，与掩码按位与但过滤掉超过该线程id的记录，只保留变量前面的判定
		int b = __ballot(pred) & t_m;
		//popc count the number of bit one. simply count the number predicated true BEFORE MY INDEX
		//计算只计算当前线程索引对应的前N个的位数之和
		//即为线程束内排他扫描
		int t_u = __popc(b);

		//由每个线程束最后一个线程写入共享内存，对应id为线程束ID，将本线程ID加回，
		//将包含求和的最终值写入共享内存，不包含求和的值没有被覆盖
		//warpTotals长度为4
		if (w_l == WARPSIZE - 1)
			warpTotals[w_i] = t_u + pred;

		__syncthreads();

		//线程束id为0，线程束内线程id，若blockDim.x = 128，则w_l < 128/32 = 4
		//线程块内第一个线程束的前（4）个线程束工作，w_l < 活动线程束数（4），即每个线程束被一个线程运行
		if (w_i == 0 && w_l < blockDim.x / WARPSIZE) {
			int w_i_u = 0;
			for (int j = 0; j <= 5; ++j) {
				//# of the ones in the j'th digit of the warp offsets
				//0->5 6个位置：
				//000 001
				//000 010
				//000 100
				//001 000
				//010 000
				//100 000
				int b_j = __ballot(warpTotals[w_l] & pow2i(j));
				w_i_u += (__popc(b_j & t_m)) << j;
				//printf("indice %i t_m=%i,j=%i,b_j=%i,w_i_u=%i\n",w_l,t_m,j,b_j,w_i_u);
			}
			warpTotals[w_l] = w_i_u;
		}
		__syncthreads();

		if (pred) {
			const int evt_idx = t_u + warpTotals[w_i] + BOffset[blockIdx.x];
			ConEvt[evt_idx] = scp[idx].z;
			//			printf("ConEvt[%d] = %d\n", evt_idx, ConEvt[evt_idx]);
		}
	}
}

__global__ void MVarPreSetZero(int* mVarPre) {
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < D_VS_SIZE)
		mVarPre[tid] = 0;
}

__global__ void SVarPreSetZero(int* sVarPre) {
	const int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < D_VS_SIZE*D_VS_SIZE*D_MDS)
		sVarPre[tid] = 0;
}

__global__ void CsCheckMain(int* mConEvt, int* mVarPre, int3* scope,
	u32* bitDom, uint2* bitSup) {
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	__shared__ int3 s_scp;
	__shared__ uint2 s_bitDom;

	uint2 l_sum = make_uint2(0, 0);
	uint2 l_res = make_uint2(0, 0);

	if (tid == 0) {
		const int c_id = mConEvt[bid];
		s_scp = scope[c_id];
		s_bitDom.x = bitDom[s_scp.x];
		s_bitDom.y = bitDom[s_scp.y];
	}
	__syncthreads();

	if (tid < D_MDS)
		l_sum = bitSup[GetBitSupIdxDevice(s_scp.z, tid)];
	__syncthreads();

	//归约
	l_sum.x &= s_bitDom.y;
	l_sum.y &= s_bitDom.x;

	//投票并反转
	l_res.x = __brev(__ballot(l_sum.x));
	l_res.y = __brev(__ballot(l_sum.y));

	__syncthreads();
	if (tid == 0) {
		//存入全局内存,并记录改变
		l_res.x &= s_bitDom.x;
		if (s_bitDom.x != l_res.x) {
			atomicAnd(&bitDom[s_scp.x], l_res.x);
			mVarPre[s_scp.x] = 1;

			if (bitDom[s_scp.x] == 0)
				mVarPre[s_scp.x] = INT_MIN;
		}

		l_res.y &= s_bitDom.y;
		if (s_bitDom.y != l_res.y) {
			atomicAnd(&bitDom[s_scp.y], l_res.y);
			mVarPre[s_scp.y] = 1;

			if (bitDom[s_scp.y] == 0)
				mVarPre[s_scp.y] = INT_MIN;
		}
	}
}

__global__ void CsCheckSub(int3* sConEvt, int* sVarPre, int* mVarPre, int3* scope, u32* bitSubDom, u32* bitDom, uint2* bitSup) {
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	__shared__ int3 s_scp;
	__shared__ int3 s_cevt;
	__shared__ uint2 s_bitSubDom;
	__shared__ int2 s_bitSubDomIdx;
	__shared__ int2 s_sVarPreIdx;
	uint2 l_sum = make_uint2(0, 0);
	uint2 l_res = make_uint2(0, 0);
	if (tid == 0) {
		s_cevt = sConEvt[bid];
		s_scp = scope[s_cevt.z];
		s_bitSubDomIdx.x = GetBitSubDomIdxDevice(s_cevt.x, s_cevt.y, s_scp.x);
		s_bitSubDomIdx.y = GetBitSubDomIdxDevice(s_cevt.x, s_cevt.y, s_scp.y);
		s_sVarPreIdx.x = GetSVarPreIdxDevice(s_cevt.x, s_cevt.y, s_scp.x);
		s_sVarPreIdx.y = GetSVarPreIdxDevice(s_cevt.x, s_cevt.y, s_scp.y);
		s_bitSubDom.x = bitSubDom[s_bitSubDomIdx.x];
		s_bitSubDom.y = bitSubDom[s_bitSubDomIdx.y];
		//printf("c_evt[%d] = (%d, %d, %d) = (%d, %d), bitSubDom.x = %8x, bitSubDom.y = %8x\n", bid, s_cevt.x, s_cevt.y, s_cevt.z, s_scp.x, s_scp.y, s_bitSubDom.x, s_bitSubDom.y);
	}
	__syncthreads();

	if (tid < D_MDS)
		l_sum = bitSup[GetBitSupIdxDevice(s_scp.z, tid)];
	__syncthreads();

	//归约
	l_sum.x &= s_bitSubDom.y;
	l_sum.y &= s_bitSubDom.x;

	//投票并反转
	l_res.x = __brev(__ballot(l_sum.x));
	l_res.y = __brev(__ballot(l_sum.y));

	__syncthreads();
	if (tid == 0) {
		//存入全局内存,并记录改变
		l_res.x &= s_bitSubDom.x;
		if (s_bitSubDom.x != l_res.x) {
			atomicAnd(&bitSubDom[s_bitSubDomIdx.x], l_res.x);
			sVarPre[s_sVarPreIdx.x] = 1;

			if (bitSubDom[s_bitSubDomIdx.x] == 0) {
				sVarPre[s_sVarPreIdx.x] = INT_MIN;
				//printf("bid = %d, (%d, %d), c_id = %d should be delete!\n", bid, s_cevt.x, s_cevt.y, s_cevt.z);
				DelValDevice(bitDom, mVarPre, s_cevt.x, s_cevt.y);
			}
		}
	}
	__syncthreads();

	if (tid == 0) {
		l_res.y &= s_bitSubDom.y;
		if (s_bitSubDom.y != l_res.y) {
			atomicAnd(&bitSubDom[s_bitSubDomIdx.y], l_res.y);
			sVarPre[s_sVarPreIdx.y] = 1;

			if (bitSubDom[s_bitSubDomIdx.y] == 0) {
				sVarPre[s_sVarPreIdx.y] = INT_MIN;
				//printf("(%d, %d), c_id = %d should be delete!\n", s_cevt.x, s_cevt.y, s_cevt.z);
				DelValDevice(bitDom, mVarPre, s_cevt.x, s_cevt.y);
			}
		}
	}
	__syncthreads();
}

__global__ void CsCheckMainLocal(int* mConEvt, int* mVarPre, int3* scope,
	u32* bitDom, uint2* bitSup) {
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int c_id = mConEvt[bid];
	const int3 l_scp = scope[c_id];
	const uint2 l_bitDom = make_uint2(bitDom[l_scp.x], bitDom[l_scp.y]);

	uint2 l_sum = make_uint2(0, 0);
	uint2 l_res = make_uint2(0, 0);
	__syncthreads();

	if (tid < D_MDS)
		l_sum = bitSup[GetBitSupIdxDevice(l_scp.z, tid)];
	__syncthreads();

	//归约
	l_sum.x &= l_bitDom.y;
	l_sum.y &= l_bitDom.x;
	__syncthreads();
	//投票并反转
	l_res.x = __brev(__ballot(l_sum.x));
	l_res.y = __brev(__ballot(l_sum.y));
	__syncthreads();

	if (tid == 0) {
		//		printf("c = %d, x = %x, y = %x\n", s_scp.z, l_res.x, l_res.y);
		//存入全局内存,并记录改变
		l_res.x &= l_bitDom.x;
		if (l_bitDom.x != l_res.x) {
			atomicAnd(&bitDom[l_scp.x], l_res.x);
			mVarPre[l_scp.x] = 1;

			if (bitDom[l_scp.x] == 0)
				mVarPre[l_scp.x] = INT_MIN;
		}

		l_res.y &= l_bitDom.y;
		if (l_bitDom.y != l_res.y) {
			atomicAnd(&bitDom[l_scp.y], l_res.y);
			mVarPre[l_scp.y] = 1;

			if (bitDom[l_scp.x] == 0)
				mVarPre[l_scp.x] = INT_MIN;
		}
	}
}

__global__ void UpdateSubDom(u32* bitDom, u32* bitSubDom, int* SVarPre,
	int* var_size, int vs_size) {
	const int val = blockIdx.x;
	const int var = blockIdx.y;
	const int subVar = threadIdx.x;
	__shared__ int pred;

	//记录子问题标题是否被删除
	if (subVar == 0)
		pred = bitDom[var] & D_U32_MASK1[val];
	__syncthreads();

	if (val < var_size[var]) {
		if (subVar < vs_size) {
			const int sub_idx = GetBitSubDomIdxDevice(var, val, subVar);
			const int subpre_idx = GetSVarPreIdxDevice(var, val, subVar);

			//子问题被删除
			if (!pred) {
				bitSubDom[sub_idx] = 0;
			}
			else {
				if (bitSubDom[sub_idx] != bitDom[subVar]) {
					//若有变量值发生变化
					SVarPre[subpre_idx] = 1;
					bitSubDom[sub_idx] &= bitDom[subVar];
					//printf("(%d, %d, %d) = %8x, %8x, %d\n", var, val, subVar, bitSubDom[sub_idx], bitDom[subVar], SVarPre[subpre_idx]);
				}
				else {
					//若有变量值没发生变化
					SVarPre[subpre_idx] = 0;
				}
			}
			//		printf("(%d, %d): %d = %8x, %d\n", var, val, subVar, bitSubDom[sub_idx],
			//				SVarPre[subpre_idx]);
			//		printf("(%d, %d, %d): %x, idx = %d\n", var, val, subVar, bitSubDom[idx],
			//				idx);
			//		printf("(%d, %d, %d): %x, idx = %d\n", i, j, k, h_bitSubDom[idx], idx);
		}
	}
}

//int CompactConsQueMain(int mcc_blocks, int work_size,
//	thrust::device_vector<int> MCCBCount,
//	thrust::device_vector<int> MCCBOffset);
float SACGPU();

float BuidBitModel32bit(XModel *xm) {
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	//拷贝一些常量
	H_DS_SIZE = xm->feature.ds_size;
	H_VS_SIZE = xm->feature.vs_size;
	H_CS_SIZE = xm->feature.cs_size;
	H_MDS = xm->feature.max_dom_size;
	H_BD_INTSIZE = GetTopNum(H_MDS, U32_BIT);	// = 1
	H_BITSUP_INTSIZE = H_MDS;
	cudaMemcpyToSymbol(D_DS_SIZE, &H_DS_SIZE, U32_SIZE);
	cudaMemcpyToSymbol(D_VS_SIZE, &H_VS_SIZE, U32_SIZE);
	cudaMemcpyToSymbol(D_CS_SIZE, &H_CS_SIZE, U32_SIZE);
	cudaMemcpyToSymbol(D_BD_INTSIZE, &H_BD_INTSIZE, U32_SIZE);
	cudaMemcpyToSymbol(D_MDS, &H_MDS, U32_SIZE);
	cudaMemcpyToSymbol(D_BITSUP_INTSIZE, &H_BITSUP_INTSIZE, U32_SIZE);

	H_MCCount = GetTopNum(H_CS_SIZE, WORKSIZE);
	H_MVCount = GetTopNum(H_VS_SIZE, WORKSIZE);
	cudaDeviceSynchronize();
	//初始化scope
	h_scope = (int3*)malloc(H_CS_SIZE * sizeof(int3));
	cudaMalloc(&d_scope, H_CS_SIZE * sizeof(int3));

	for (int i = 0; i < H_CS_SIZE; ++i) {
		XCon *c = xm->cons[i];
		h_scope[i].x = c->scope[0];
		h_scope[i].y = c->scope[1];
		h_scope[i].z = c->id;
	}

	cudaMemcpy(d_scope, h_scope, H_CS_SIZE * sizeof(int3),
		cudaMemcpyHostToDevice);
	//	showConstraints<<<H_MCCount, WORKSIZE>>>(d_scope, H_CS_SIZE);
	cudaDeviceSynchronize();
	//变量论域（0, 32]
	h_bitDom = (u32*)malloc(H_VS_SIZE * sizeof(u32));
	h_MVarPre = (int*)malloc(H_VS_SIZE * sizeof(int));
	h_var_size = (int*)malloc(H_VS_SIZE * sizeof(u32));
	cudaMalloc(&d_bitDom, H_VS_SIZE * sizeof(u32));
	cudaMalloc(&d_MVarPre, H_VS_SIZE * sizeof(u32));
	cudaMalloc(&d_var_size, H_VS_SIZE * sizeof(u32));

	for (int i = 0; i < H_VS_SIZE; ++i) {
		XVar* v = xm->vars[i];
		XDom* d = xm->doms[v->dom_id];
		const int dom_size = d->size;
		h_bitDom[i] = UINT32_MAX << GetOffSet(dom_size);
		h_MVarPre[i] = 1;
		h_var_size[i] = d->size;
	}

	cudaMemcpy(d_bitDom, h_bitDom, H_VS_SIZE * sizeof(u32),
		cudaMemcpyHostToDevice);
	cudaMemcpy(d_MVarPre, h_MVarPre, H_VS_SIZE * sizeof(int),
		cudaMemcpyHostToDevice);
	cudaMemcpy(d_var_size, h_var_size, H_VS_SIZE * sizeof(u32),
		cudaMemcpyHostToDevice);
	//	showVariables<<<H_MVCount, WORKSIZE>>>(d_bitDom, d_MVarPre, d_var_size,
	//			H_VS_SIZE);
	cudaDeviceSynchronize();
	//创建bitSubDom
	h_bitSubDom = (u32*)malloc(H_VS_SIZE * H_MDS * H_VS_SIZE * sizeof(u32));
	cudaMalloc(&d_bitSubDom, H_VS_SIZE * H_MDS * H_VS_SIZE * sizeof(u32));

	for (int i = 0; i < H_VS_SIZE; ++i) {
		for (int j = 0; j < H_MDS; ++j) {
			//若(i, j)子问题
			for (int k = 0; k < H_VS_SIZE; ++k) {
				const int idx = GetBitSubDomIndexHost(i, j, k);
				if (j < h_var_size[i]) {
					//bitDom[k]有问题
					h_bitSubDom[idx] = h_bitDom[k];
				}
				else {
					h_bitSubDom[idx] = 0;
				}
			}

			//最后将(i, j)的bitDom 改掉
			//获取i,j, i的起始地址，并修改
			const int iji_idx = GetBitSubDomIndexHost(i, j, i);
			h_bitSubDom[iji_idx] = U32_MASK1[j];
		}
	}
	//for (int i = 0; i < H_VS_SIZE; ++i) {
	//	for (int j = 0; j < H_MDS; ++j) {
	//		//若(i, j)子问题
	//		for (int k = 0; k < H_VS_SIZE; ++k) {
	//			const int idx = GetBitSubDomIndexHost(i, j, k);
	//			printf("(%d, %d, %d): %x, idx = %d\n", i, j, k, h_bitSubDom[idx], idx);
	//		}
	//	}
	//}

	cudaMemcpy(d_bitSubDom, h_bitSubDom,
		H_VS_SIZE * H_MDS * H_VS_SIZE * sizeof(u32),
		cudaMemcpyHostToDevice);
	// dim.x dim.y dim.z
	dim3 SubProDim(H_MDS, H_VS_SIZE);
	int topnum_var_threads = GetTopNum(H_VS_SIZE, WORKSIZE) * WORKSIZE;
	showSubVar << <SubProDim, topnum_var_threads >> > (d_bitSubDom, d_var_size,
		H_VS_SIZE);
	cudaDeviceSynchronize();
	//创建bitSup
	h_bitSup = (uint2*)malloc(H_CS_SIZE * H_MDS * sizeof(uint2));
	cudaMalloc(&d_bitSup, H_CS_SIZE * H_MDS * sizeof(uint2));

	for (int i = 0; i < H_CS_SIZE; ++i) {
		XCon* c = xm->cons[i];
		XRel* r = xm->rels[c->rel_id];
		XVar* v[2] = { xm->vars[c->scope[0]], xm->vars[c->scope[1]] };
		//		XDom* d[2] = { xm->doms[v[0]->dom_id], xm->doms[v[1]->dom_id] };
		//初始化位矩阵
		for (int j = 0; j < H_MDS; ++j) {
			const int idx = c->id * H_BITSUP_INTSIZE + j;
			//支持取0x0000..., 冲突取0xFFF...
			h_bitSup[idx].x = (r->sem == SEM_CONFLICT) ? h_bitDom[v[1]->id] : 0;
			h_bitSup[idx].y = (r->sem == SEM_CONFLICT) ? h_bitDom[v[0]->id] : 0;
		}

		//向位矩阵中填充值
		for (int j = 0; j < r->size; ++j) {
			const int2 t = make_int2(r->tuples[j][0], r->tuples[j][1]);
			const int2 idx = GetBitSupIndexByTupleHost(c->id, t);
			if (r->sem == SEM_SUPPORT) {
				h_bitSup[idx.x].x |= U32_MASK1[t.y & U32_MOD_MASK];
				h_bitSup[idx.y].y |= U32_MASK1[t.x & U32_MOD_MASK];
			}
			else {
				h_bitSup[idx.x].x &= U32_MASK0[t.y & U32_MOD_MASK];
				h_bitSup[idx.y].y &= U32_MASK0[t.x & U32_MOD_MASK];
			}
		}
		//		for (int j = 0; j < H_MDS; ++j) {
		//			const int idx = c->id * H_BITSUP_INTSIZE + j;
		//			printf("bitSup[%d, x, %d] = %x, bitSup[%d, y, %d] = %x\n", i, j, h_bitSup[idx].x, i, j,
		//					h_bitSup[idx].y);
		//		}
	}

	cudaMemcpy(d_bitSup, h_bitSup, H_CS_SIZE * H_MDS * sizeof(uint2),
		cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	//生成主问题约束队列
	h_MCon = (int*)malloc(H_CS_SIZE * sizeof(int));
	h_MConPre = (int*)malloc(H_CS_SIZE * sizeof(int));
	h_MConEvt = (int*)malloc(H_CS_SIZE * sizeof(int));
	cudaMalloc(&d_MCon, H_CS_SIZE * sizeof(int));
	cudaMalloc(&d_MConPre, H_CS_SIZE * sizeof(int));
	cudaMalloc(&d_MConEvt, H_CS_SIZE * sizeof(int));

	for (int i = 0; i < H_CS_SIZE; ++i) {
		h_MCon[i] = i;
		h_MConPre[i] = 1;
		h_MConEvt[i] = i;
	}
	cudaMemcpy(d_MCon, h_MCon, H_CS_SIZE * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_MConPre, h_MConPre, H_CS_SIZE * sizeof(int),
		cudaMemcpyHostToDevice);
	cudaMemcpy(d_MConEvt, h_MConEvt, H_CS_SIZE * sizeof(int),
		cudaMemcpyHostToDevice);
	//生成子问题约束队列
	h_SConPre = (int*)malloc(H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int));
	h_SConEvt = (int3*)malloc(H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3));
	h_SCon = (int3*)malloc(H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3));
	h_SVar = (int3*)malloc(H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3));
	h_SVarPre = (int*)malloc(H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int));

	cudaMalloc(&d_SConPre, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int));
	cudaMalloc(&d_SConEvt, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3));
	cudaMalloc(&d_SCon, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3));
	cudaMalloc(&d_SVar, H_VS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3));
	cudaMalloc(&d_SVarPre, H_VS_SIZE * H_VS_SIZE * H_MDS * sizeof(int));

	for (int i = 0; i < H_VS_SIZE; ++i) {
		for (int j = 0; j < H_MDS; ++j) {
			for (int k = 0; k < H_CS_SIZE; ++k) {
				const int idx = (i * H_MDS + j) * H_CS_SIZE + k;
				h_SCon[idx].x = i;
				h_SCon[idx].y = j;
				h_SCon[idx].z = k;

				h_SConEvt[idx].x = i;
				h_SConEvt[idx].y = j;
				h_SConEvt[idx].z = k;

				h_SConPre[idx] = 0;
			}

			for (int k = 0; k < H_VS_SIZE; ++k) {
				//子问题(i, j) k为变量id
				//const int idx = (i * H_MDS + j) * H_VS_SIZE + k;
				const int idx = GetSVarPreIdxHost(i, j, k);
				h_SVar[idx].x = i;
				h_SVar[idx].y = j;
				h_SVar[idx].z = k;

				h_SVarPre[idx] = 0;
			}
		}
	}

	cudaMemcpy(d_SConPre, h_SConPre, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SConEvt, h_SConEvt, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SCon, h_SCon, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SVar, h_SVar, H_VS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SVarPre, h_SVarPre, H_VS_SIZE * H_VS_SIZE * H_MDS * sizeof(int), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	//获得主问题压缩的BLOCK数
	//	H_MCCBLOCK = GetTopNum(CS_SIZE, num_threads);
	//	MCC_BCount.resize(MCC_BLOCK, 0);
	//	MCC_BOffset.resize(MCC_BLOCK, 0);
	//	MCC_BlocksCount = thrust::raw_pointer_cast(MCC_BCount.data());
	//	MCC_BlocksOffset = thrust::raw_pointer_cast(MCC_BOffset.data());

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return elapsedTime;
}

__global__ void GenSubConPre(int* SVarPre, int* SCCBCount, int3* scp, int len) {
	//获取子问题，及
	const int val = blockIdx.x;
	const int var = blockIdx.y;
	const int sub_con = threadIdx.x;
	int3 sp = scp[sub_con];
	__shared__ int failed;
	int pred = 0;

	if (sub_con == 0)
		failed = 0;
	__syncthreads();

	if (sub_con < len) {
		const int2 v_pre = make_int2(SVarPre[GetSVarPreIdxDevice(var, val, sp.x)],
			SVarPre[GetSVarPreIdxDevice(var, val, sp.y)]);
		if (v_pre.x == 0 && v_pre.y == 0)
			pred = 0;
		else if (v_pre.x == INT32_MIN || v_pre.y == INT32_MIN)
			failed = INT32_MIN;
		else
			pred = 1;
	}
	__syncthreads();

	int BC = (failed == INT32_MIN) ? INT32_MIN : __syncthreads_count(pred);

	if (threadIdx.x == 0) {
		SCCBCount[blockIdx.y * gridDim.x + blockIdx.x] = BC;
		//printf("(%d, %d) idx = %d, BC = %d\n", var, val, (blockIdx.y * gridDim.x + blockIdx.x), BC);
	}
}

//__global__ void CompactSubQ(int* SVarPre, int3* ConEvt, int* BOffset, int3* scp, int len) {
//	//获取子问题，及
//	const int val = blockIdx.x;
//	const int var = blockIdx.y;
//	const int sub_con = threadIdx.x;
//	const int blockIdx1D = blockIdx.y * gridDim.x + blockIdx.x;
//	const int g_idx = blockIdx1D * blockDim.x + threadIdx.x;
//	extern __shared__ int warpTotals[];
//	int pred = 0;
//	int3 sp;
//
//	if (sub_con < len) {
//		sp = scp[sub_con];
//		//获得判定
//		const int2 v_pre = make_int2(SVarPre[GetSVarPreIdxDevice(var, val, sp.x)],
//			SVarPre[GetSVarPreIdxDevice(var, val, sp.y)]);
//		pred = v_pre.x || v_pre.y;
//
//
//		//	printf("scp[%d, %d, %d], gidx = %d, pred = %d, blockIdx1D = %d\n", var, val,
//		//			sub_con, g_idx, pred, blockIdx1D);
//
//		//warp index
//		//线程束索引
//		int w_i = sub_con / WARPSIZE;
//		//thread index within a warp
//		//线程束内线程索引
//		int w_l = g_idx % WARPSIZE;
//		//thread mask (ERROR IN THE PAPERminus one is required)
//		//线程掩码
//		//INT_MAX = 1111 1111 1111 1111 1111 1111 1111 1111
//		//若线程内id=0，右移32-0-1 = 31位 右侧剩下1位
//		//若线程内id=5，右移32-5-1 = 26位 右侧剩下6位
//		//若线程内id=31，右移32-31-1 = 0位 右侧剩下32位
//		//线程束内threid|  31~~~~~~0
//		//ballot对应位置|   1......1
//		int t_m = INT_MAX >> (WARPSIZE - w_l - 1);
//		//balres = number whose ith bit is one if the ith's thread pred is true masked up to the current index in warp
//		//线程内局部变量pred = 1，与掩码按位与但过滤掉超过该线程id的记录，只保留变量前面的判定
//		int b = __ballot(pred) & t_m;
//		//popc count the number of bit one. simply count the number predicated true BEFORE MY INDEX
//		//计算只计算当前线程索引对应的前N个的位数之和
//		//即为线程束内排他扫描
//		int t_u = __popc(b);
//
//		//由每个线程束最后一个线程写入共享内存，对应id为线程束ID，将本线程ID加回，
//		//将包含求和的最终值写入共享内存，不包含求和的值没有被覆盖
//		//warpTotals长度为4
//		if (w_l == WARPSIZE - 1)
//			//	{
//			warpTotals[w_i] = t_u + pred;
//		//		printf("warpTotals[%d] = %d\n", w_i, warpTotals[w_i]);
//		//	}
//		__syncthreads();
//
//		//线程束id为0，线程束内线程id，若blockDim.x = 128，则w_l < 128/32 = 4
//		//线程块内第一个线程束的前（4）个线程束工作，w_l < 活动线程束数（4），即每个线程束被一个线程运行
//		if (w_i == 0 && w_l < blockIdx1D / WARPSIZE) {
//			int w_i_u = 0;
//			for (int j = 0; j <= 5; ++j) {
//				//# of the ones in the j'th digit of the warp offsets
//				//0->5 6个位置：
//				//000 001
//				//000 010
//				//000 100
//				//001 000
//				//010 000
//				//100 000
//				int b_j = __ballot(warpTotals[w_l] & pow2i(j));
//				w_i_u += (__popc(b_j & t_m)) << j;
//				//printf("indice %i t_m=%i,j=%i,b_j=%i,w_i_u=%i\n",w_l,t_m,j,b_j,w_i_u);
//			}
//			warpTotals[w_l] = w_i_u;
//		}
//		__syncthreads();
//		if (pred) {
//			const int evt_idx = t_u + warpTotals[w_i] + BOffset[blockIdx1D];
//			//		printf("warpTotals[%d] = %d\n", w_i,  warpTotals[0]);
//			//		printf("ConEvt[%d] = %d, %d, %d\n", evt_idx, 0, 0, 0);
//			ConEvt[evt_idx].x = sp.x;
//			ConEvt[evt_idx].y = sp.y;
//			ConEvt[evt_idx].z = sp.z;
//			//printf("ConEvt[%d] = %d, %d, %d\n", evt_idx, ConEvt[evt_idx].x,
//			//	ConEvt[evt_idx].y, ConEvt[evt_idx].z);
//		}
//	}
//}

// 一维
__global__ void CompactSubQ(int* SVarPre, int3* SConEvt, int3* SCon, int* BOffset, int3* scp, int len) {
	//获取子问题，及
	const int bid = blockIdx.x;
	const int val = bid % D_MDS;
	const int var = bid / D_MDS;
	const int sub_con = threadIdx.x;
	const int idx = threadIdx.x + blockIdx.x*blockDim.x;
	const int g_idx = threadIdx.x + blockIdx.x*len;
	extern __shared__ int warpTotals[];
	int pred = 0;

	if (sub_con < len) {
		int3 sp = scp[sub_con];
		//获得判定
		const int2 v_pre = make_int2(SVarPre[GetSVarPreIdxDevice(var, val, sp.x)],
			SVarPre[GetSVarPreIdxDevice(var, val, sp.y)]);
		pred = v_pre.x || v_pre.y;

		//if (var == 0 && val < 2)
		//	printf("scp[%d, %d, %d], idx = %d, pred = %d, bid = %d\n", var, val,
		//		sub_con, idx, pred, bid);

		//warp index
		//线程束索引
		int w_i = threadIdx.x / WARPSIZE;
		//thread index within a warp
		//线程束内线程索引
		int w_l = idx % WARPSIZE;
		//thread mask (ERROR IN THE PAPERminus one is required)
		//线程掩码
		//INT_MAX = 1111 1111 1111 1111 1111 1111 1111 1111
		//若线程内id=0，右移32-0-1 = 31位 右侧剩下1位
		//若线程内id=5，右移32-5-1 = 26位 右侧剩下6位
		//若线程内id=31，右移32-31-1 = 0位 右侧剩下32位
		//线程束内threid|  31~~~~~~0
		//ballot对应位置|   1......1
		int t_m = INT_MAX >> (WARPSIZE - w_l - 1);
		//balres = number whose ith bit is one if the ith's thread pred is true masked up to the current index in warp
		//线程内局部变量pred = 1，与掩码按位与但过滤掉超过该线程id的记录，只保留变量前面的判定
		int b = __ballot(pred) & t_m;
		//popc count the number of bit one. simply count the number predicated true BEFORE MY INDEX
		//计算只计算当前线程索引对应的前N个的位数之和
		//即为线程束内排他扫描
		int t_u = __popc(b);

		//由每个线程束最后一个线程写入共享内存，对应id为线程束ID，将本线程ID加回，
		//将包含求和的最终值写入共享内存，不包含求和的值没有被覆盖
		//warpTotals长度为4
		if (w_l == WARPSIZE - 1)
			warpTotals[w_i] = t_u + pred;
		__syncthreads();

		//线程束id为0，线程束内线程id，若blockDim.x = 128，则w_l < 128/32 = 4
		//线程块内第一个线程束的前（4）个线程束工作，w_l < 活动线程束数（4），即每个线程束被一个线程运行
		if (w_i == 0 && w_l < blockDim.x / WARPSIZE) {
			int w_i_u = 0;
			for (int j = 0; j <= 5; ++j) {
				//# of the ones in the j'th digit of the warp offsets
				//0->5 6个位置：
				//000 001
				//000 010
				//000 100
				//001 000
				//010 000
				//100 000
				int b_j = __ballot(warpTotals[w_l] & pow2i(j));
				w_i_u += (__popc(b_j & t_m)) << j;
				//printf("indice %i t_m=%i,j=%i,b_j=%i,w_i_u=%i\n",w_l,t_m,j,b_j,w_i_u);
			}
			warpTotals[w_l] = w_i_u;
		}
		__syncthreads();
		if (pred) {
			const int evt_idx = t_u + warpTotals[w_i] + BOffset[blockIdx.x];
			//		printf("warpTotals[%d] = %d\n", w_i,  warpTotals[0]);
			//		printf("ConEvt[%d] = %d, %d, %d\n", evt_idx, 0, 0, 0);
			SConEvt[evt_idx].x = SCon[g_idx].x;
			SConEvt[evt_idx].y = SCon[g_idx].y;
			SConEvt[evt_idx].z = SCon[g_idx].z;
			//printf("ConEvt[%d] = %d, %d, %d\n", evt_idx, SConEvt[evt_idx].x,
			//	SConEvt[evt_idx].y, SConEvt[evt_idx].z);
		}
	}
}


int CompactConsQueMain(int mcc_blocks, int mvc_blocks, int work_size, int work_size_large, thrust::device_vector<int> MCCBCount, thrust::device_vector<int> MCCBOffset) {
	int* d_MCCBCount_ptr;
	int* d_MCCBOffset_ptr;
	d_MCCBCount_ptr = thrust::raw_pointer_cast(MCCBCount.data());
	d_MCCBOffset_ptr = thrust::raw_pointer_cast(MCCBOffset.data());
	//以约束数量启动
	//P1
	GenConPre << <mcc_blocks, work_size >> > (d_MVarPre, d_MCCBCount_ptr, d_scope,
		H_CS_SIZE);
	cudaDeviceSynchronize();
	//P2
	thrust::exclusive_scan(MCCBCount.begin(), MCCBCount.end(),
		MCCBOffset.begin());
	cudaDeviceSynchronize();
	int total = MCCBOffset[MCCBCount.size() - 1]
		+ MCCBCount[MCCBCount.size() - 1];
	//	std::cout << "total = " << total << std::endl;
	if (total < 0)
		return PROPFAILED;

	////P3
	//每个约束一个线程进行归约,共享内存大小 = 一个块内线程束的个数,用来装载线程束计算结果
	CompactQ << <mcc_blocks, work_size, sizeof(int) * (work_size / WARPSIZE) >> > (
		d_MVarPre, d_MConEvt, d_MCCBOffset_ptr, d_scope, H_CS_SIZE);
	cudaDeviceSynchronize();

	MVarPreSetZero << <mvc_blocks, work_size_large >> > (d_MVarPre);
	cudaDeviceSynchronize();
	return total;
}

int CompactConsQueSub(dim3 scc_blocks, int svc_blocks, const int vs_size, int cs_size, int con_thrds, int work_size_large, thrust::device_vector<int> SCCBCount, thrust::device_vector<int> SCCBOffset) {
	int* d_SCCBCount_ptr;
	int* d_SCCBOffset_ptr;
	d_SCCBCount_ptr = thrust::raw_pointer_cast(SCCBCount.data());
	d_SCCBOffset_ptr = thrust::raw_pointer_cast(SCCBOffset.data());
	//以约束数量启动
	//P1
	GenSubConPre << <scc_blocks, con_thrds >> > (d_SVarPre, d_SCCBCount_ptr, d_scope, cs_size);
	cudaDeviceSynchronize();

	//P2
	thrust::exclusive_scan(SCCBCount.begin(), SCCBCount.end(), SCCBOffset.begin());
	cudaDeviceSynchronize();

	int total = SCCBOffset[SCCBCount.size() - 1] + SCCBCount[SCCBCount.size() - 1];
	//std::cout << "total = " << total << std::endl;

	//for (size_t i = 0; i < SCCBCount.size(); i++)
	//	std::cout << i << ":" << SCCBCount[i] << std::endl;

	if (total < 0)
		return PROPFAILED;

	//////P3
	////每个约束一个线程进行归约,共享内存大小 = 一个块内线程束的个数,用来装载线程束计算结果
	//CompactSubQ << <scc_blocks, con_thrds, sizeof(int) * (con_thrds / WARPSIZE) >> > (
	//	d_SVarPre, d_SConEvt, d_SCCBOffset_ptr, d_scope, cs_size);

	CompactSubQ << < H_VS_SIZE * H_MDS, con_thrds, sizeof(int) * (con_thrds / WARPSIZE) >> > (d_SVarPre, d_SConEvt, d_SCon, d_SCCBOffset_ptr, d_scope, cs_size);
	cudaDeviceSynchronize();

	//cudaMemcpy(h_SConEvt, d_SConEvt, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3), cudaMemcpyDeviceToHost);
	//for (size_t i = 0; i < total; i++)
	//{
	//	printf("h_SConEvt[%d] = (%2d, %2d, %3d)\n", i, h_SConEvt[i].x, h_SConEvt[i].y, h_SConEvt[i].z);
	//}

	SVarPreSetZero << <svc_blocks, work_size_large >> > (d_SVarPre);
	cudaDeviceSynchronize();
	//std::cout << "------------------------------" << std::endl;
	//cudaMemcpy(h_SVarPre, d_SVarPre, H_VS_SIZE * H_VS_SIZE * H_MDS * sizeof(int), cudaMemcpyDeviceToHost);
	//for (size_t i = 0; i < H_VS_SIZE; i++)
	//{
	//	for (size_t j = 0; j < H_MDS; j++)
	//	{
	//		for (size_t k = 0; k < H_VS_SIZE; k++)
	//		{
	//			const int idx = GetSVarPreIdxHost(i, j, k);
	//			if (h_SVarPre[idx] == 1)
	//				printf("h_SVarPre[%2d, %2d, %2d], %d = %d\n", i, j, k, idx, h_SVarPre[idx]);
	//		}

	//	}
	//}
	//std::cout << "------------------------------" << std::endl;
	return total;
}
void ConstraintsCheckMain(int c_total) {
	CsCheckMain << <c_total, WORKSIZE >> > (d_MConEvt, d_MVarPre, d_scope, d_bitDom, d_bitSup);
	//	CsCheckMainLocal<<<c_total, WORKSIZE>>>(d_MConEvt, d_MVarPre, d_scope,
	//			d_bitDom, d_bitSup);
	//	CsCheckMain<<<1, WORKSIZE, (sizeof(int3) + sizeof(uint2))>>>(d_MConEvt,
	//			d_MVarPre, d_scope, d_bitDom, d_bitSup);
	cudaDeviceSynchronize();
}

void ConstraintsCheckSub(int c_total) {
	CsCheckSub << <c_total, WORKSIZE >> > (d_SConEvt, d_SVarPre, d_MVarPre, d_scope, d_bitSubDom, d_bitDom, d_bitSup);
	//ShowSubConEvt << <c_total, WORKSIZE >> > (d_SConEvt);
	cudaDeviceSynchronize();
}

float SACGPU() {
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	float elapsedTime;
	//1. 在主问题上执行AC
	//1.1. 主问题运行规格
	//以约束个数运行的程序块
	const int H_MCCBLOCK = GetTopNum(H_CS_SIZE, WORKSIZE);
	thrust::device_vector<int> d_MCCBCount(H_MCCBLOCK, 0);
	thrust::device_vector<int> d_MCCBOffset(H_MCCBLOCK, 0);
	const int H_MVCBLOCK = GetTopNum(H_VS_SIZE, WORKSIZE_LARGE);

	//1.2. 流压缩取得约束队列
	int mc_total = CompactConsQueMain(H_MCCBLOCK, H_MVCBLOCK, WORKSIZE, WORKSIZE_LARGE, d_MCCBCount, d_MCCBOffset);
	//std::cout << "mc_total = " << mc_total << std::endl;
	do {
		////1.3. 约束检查
		ConstraintsCheckMain(mc_total);
		////1.4. 流压缩取得约束队列
		mc_total = CompactConsQueMain(H_MCCBLOCK, H_MVCBLOCK, WORKSIZE, WORKSIZE_LARGE, d_MCCBCount, d_MCCBOffset);
		//std::cout << "mc_total = " << mc_total << std::endl;
		//		showVariables<<<H_MVCount, WORKSIZE>>>(d_bitDom, d_MVarPre, d_var_size, H_VS_SIZE);
	} while (mc_total > 0);
	//showVariables << <H_MVCount, WORKSIZE >> > (d_bitDom, d_MVarPre, d_var_size, H_VS_SIZE);
	//1.5. 失败返回
	if (mc_total == PROPFAILED) {
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);

		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		printf("Failed\n");
		return elapsedTime;
	}

	//2. 在子问题上执行AC
	//2.1. 更新子问题论域
	dim3 subblock_dim(H_MDS, H_VS_SIZE);
	const int var_threads = GetTopNum(H_VS_SIZE, WORKSIZE)*WORKSIZE;
	UpdateSubDom << <subblock_dim, var_threads >> > (d_bitDom, d_bitSubDom, d_SVarPre, d_var_size, H_VS_SIZE);

	//确定子问题运行规格 = 子问题为线程块，子问题约束为线程
	//规格
	const int H_SCCBLOCK = H_VS_SIZE * H_MDS;
	thrust::device_vector<int> d_SCCBCount(H_SCCBLOCK, 0);
	thrust::device_vector<int> d_SCCBOffset(H_SCCBLOCK, 0);
	const int H_SVCBLOCK = GetTopNum(H_VS_SIZE * H_VS_SIZE * H_MDS, WORKSIZE_LARGE);
	const int con_threads = GetTopNum(H_CS_SIZE, WORKSIZE)*WORKSIZE;
	//2.2. 压缩子问题队列
	int sc_total = CompactConsQueSub(subblock_dim, H_SVCBLOCK, H_VS_SIZE, H_CS_SIZE, con_threads, WORKSIZE_LARGE, d_SCCBCount, d_SCCBOffset);
	//std::cout << "sc_total = " << sc_total << std::endl;
	//2.3. 子问题上约束检查
	while (sc_total > 0)
	{
		ConstraintsCheckSub(sc_total);
		//2.4. 主问题上流压缩取得约束队列
		int mc_total = CompactConsQueMain(H_MCCBLOCK, H_MVCBLOCK, WORKSIZE, WORKSIZE_LARGE, d_MCCBCount, d_MCCBOffset);
		//std::cout << "mc_total = " << mc_total << std::endl;
		while (mc_total > 0) {
			////2.5. 约束检查
			ConstraintsCheckMain(mc_total);
			////2.6. 流压缩取得约束队列
			mc_total = CompactConsQueMain(H_MCCBLOCK, H_MVCBLOCK, WORKSIZE, WORKSIZE_LARGE, d_MCCBCount, d_MCCBOffset);
			//std::cout << "mc_total = " << mc_total << std::endl;
		}

		//1.5. 失败返回
		if (mc_total == PROPFAILED) {
			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);

			cudaEventElapsedTime(&elapsedTime, start, stop);
			cudaEventDestroy(start);
			cudaEventDestroy(stop);
			printf("Failed\n");
			return elapsedTime;
		}

		sc_total = CompactConsQueSub(subblock_dim, H_SVCBLOCK, H_VS_SIZE, H_CS_SIZE, con_threads, WORKSIZE_LARGE, d_SCCBCount, d_SCCBOffset);
		//std::cout << "sc_total = " << sc_total << std::endl;
	}

	////检查一下D_SVarPre
	//cudaMemcpy(h_SVarPre, d_SVarPre, H_VS_SIZE * H_VS_SIZE * H_MDS * sizeof(int), cudaMemcpyDeviceToHost);

	//for (size_t i = 0; i < H_VS_SIZE; i++)
	//{
	//	for (size_t j = 0; j < H_MDS; j++)
	//	{
	//		for (size_t k = 0; k < H_VS_SIZE; k++)
	//		{
	//			const int idx = GetSVarPreIdxHost(i, j, k);
	//			if (h_SVarPre[idx] == 1)
	//				printf("h_SVarPre[%2d, %2d, %2d], %d = %d\n", i, j, k, idx, h_SVarPre[idx]);
	//		}

	//	}
	//}
	//cudaMemcpy(h_SConEvt, d_SConEvt, H_CS_SIZE * H_VS_SIZE * H_MDS * sizeof(int3), cudaMemcpyDeviceToHost);
	//for (size_t i = 0; i < sc_total; i++)
	//{
	//	printf("h_SConEvt[%d] = (%2d, %2d, %3d)\n", i, h_SConEvt[i].x, h_SConEvt[i].y, h_SConEvt[i].z);
	//}
		//2.4. 主问题上流压缩取得约束队列

	////
	//////1.5. 失败返回
	////	if (mc_total == PROPFAILED) {
	////		cudaEventRecord(stop, 0);
	////		cudaEventSynchronize(stop);
	////
	////		cudaEventElapsedTime(&elapsedTime, start, stop);
	////		cudaEventDestroy(start);
	////		cudaEventDestroy(stop);
	////		printf("Failed\n");
	////		return elapsedTime;
	////	}
	//showVariables << <H_MVCount, WORKSIZE >> > (d_bitDom, d_MVarPre, d_var_size, H_VS_SIZE);
	//更新子问题论域
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return elapsedTime;
}

void DelGPUModel() {
	free(h_scope);
	free(h_var_size);
	free(h_bitDom);
	free(h_MVarPre);
	free(h_bitSubDom);
	free(h_bitSup);
	free(h_MCon);
	free(h_MConEvt);
	free(h_MConPre);
	free(h_SConPre);
	free(h_SConEvt);
	free(h_SCon);
	free(h_SVar);
	free(h_SVarPre);

	cudaFree(d_scope);
	cudaFree(d_var_size);
	cudaFree(d_bitDom);
	cudaFree(d_MVarPre);
	cudaFree(d_bitSubDom);
	cudaFree(d_bitSup);
	cudaFree(d_MCon);
	cudaFree(d_MConEvt);
	cudaFree(d_MConPre);
	cudaFree(d_SConPre);
	cudaFree(d_SConEvt);
	cudaFree(d_SCon);
	cudaFree(d_SVar);
	cudaFree(d_SVarPre);
}
#endif /* CUDASAC_CUH_ */
