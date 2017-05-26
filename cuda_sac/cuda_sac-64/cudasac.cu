//#pragma once
//// CUDA Runtime
//#include <cuda_runtime.h>
//#include <device_functions.h>
//#include <device_launch_parameters.h>
//// Utilities and system includes
//#include <helper_cuda.h>
//#include <helper_functions.h>
//#include <thrust/device_vector.h>
//#include <thrust/scan.h>
//#include <iostream>
//#include <math.h>
//#include "XModel.h"
//#include "math_functions.hpp"
//#include "cuda_runtime_api.h"
//
//using namespace cudacp;
//
//#ifndef MIN
//#define MIN(x,y) ((x < y) ? x : y)
//#endif
//// num_threads
//static const int num_threads = 128;
//static const int U32_SIZE = sizeof(u32); ///<4
//static const int U32_BIT = U32_SIZE * 8;	///<32
//static const int U32_POS = 5;
//static const int U32_MOD_MASK = 31;
//
//struct int_predicate
//{
//	__host__ __device__	bool operator()(const int x)
//	{
//		return x > 0;
//	}
//};
//
//// 一个bitDom[x]的长度
//__constant__ int D_BITDOM_INTSIZE;
//// 整个bitDom的长度
//
//__constant__ int D_BITDOMS_INTSIZE;
//
//// bit支持， uint2 bitSup[c][a][idx].x = bitSup[c,x,a,idx]
//// bit支持， uint2 bitSup[c][a][idx].y = bitSup[c,y,a,idx]
////__device__ uint2*** d_bitSup;
////__host__ uint2*** h_bitSup;
////__device__ u32** d_bitDom;
////__host__ u32** h_bitDom;
////__device__ u32**
//__constant__ int D_NUM_BD_BLOCK;
//__constant__ int D_NUM_CS_SIZE_BLOCKS;
////////////////////////////////////////////////////////////////////////////
////	一些GPU常量
////////////////////////////////////////////////////////////////////////////
//__device__ __managed__ int NUM_BD_BLOCK;
//__device__ __managed__ int NUM_CS_SIZE_BLOCKS;
//// 一个变量论域的int长度
//__device__ __managed__ int BITDOM_INTSIZE;
//// 整个变量集合的论域的int长度 
//__device__ __managed__ int BITDOMS_INTSIZE;
//// 一个约束的bitsup的int长度 
//__device__ __managed__ int BITSUP_INTSIZE;
//// 整个约束集合的bitsup的int长度
//__device__ __managed__ int BITSUPS_INTSIZE;
////子问题论域总长度
//__device__ __managed__ int BITSUBDOMS_INTSIZE;
////变量个数
//__device__ __managed__ int VS_SIZE;
////约束个数
//__device__ __managed__ int CS_SIZE;
////主问题约束压缩BLOCK数
//__device__ __managed__ int MCC_BLOCK;
////////////////////////////////////////////////////////////////////////////
////	一些GPU变量
////////////////////////////////////////////////////////////////////////////
//__device__ __managed__ int M_Qsize;
//
////////////////////////////////////////////////////////////////////////////
////  GPU约束记录信息，不可更改
////////////////////////////////////////////////////////////////////////////
////	每个变量的int大小
//__device__ __managed__ int *vars_size;
//// 存储约束的scope，类型int3，scope.x: x.id; scope.y: y.id; scope.z: c.id
//__device__ __managed__ int3* scope;
//// 最大dom
//__device__ __managed__ int MAX_DOM_SIZE;
//// subCon长度
//__device__ __managed__ int SUBCON_SIZE;
//
//
////__device__ __managed__ int BITDOM_SIZE;
////__device__ __managed__ int
//// 主问题数据结构，使用UM
//// 表示约束网络论域
//__device__ __managed__ u32* bitDom;
//// 表示约束，不可修改
//__device__ __managed__ uint2* bitSup;
//////类似队列，存储约束id
////__device__ __managed__ int *mainCon;
//////子问题数据结构
////表示子问题的约束网络论域。
//__device__ __managed__ u32* bitSubDom;
//////类似队列，存储子问题约束id subCon.x: variable，subCon.y: value，subCon.z: c.id
////__device__ __managed__ ushort3* subCon;
////标记主问题变量域是否删减，初始化全部为1
//__device__ __managed__ int* M_VarPre;
////标记主问题约束是否需检查，初始化全部为1
//__device__ __managed__ int* M_ConPre;
////主问题约束传播队列(压缩版)
//__device__ __managed__ int* M_ConEvt;
////主问题约束传播队列
//__device__ __managed__ int* M_Con;
//
////标记子问题变量域是否删减，初始化全部为1
//__device__ __managed__ int* S_VarPre;
////记录子问题变量域
//__device__ __managed__ uint3* S_Var;
////标记子问题约束是否需检查，初始化全部为1
//__device__ __managed__ int* S_ConPre;
////子问题约束传播队列(压缩版)
//__device__ __managed__ int3* S_ConEvt;
////子问题约束传播队列
//__device__ __managed__ int3* S_Con;
////
//int* MCC_BlocksCount;
//int* MCC_BlocksOffset;
//
//thrust::device_vector<int> MCC_BCount;
//thrust::device_vector<int> MCC_BOffset;
////
////__device__ __managed__ ushort4* subVar;
//////标记子问题发生改动的变量id，初始化全部为0
////__device__ __managed__ unsigned short* subEvtVar;
//////标记子问题发生改动的约束id，初始化全部为1
////__device__ __managed__ int* subEvtCon;
//
//static const u32 U32_MASK1[32] = {
//	0x80000000, 0x40000000, 0x20000000, 0x10000000,
//	0x08000000, 0x04000000, 0x02000000, 0x01000000,
//	0x00800000, 0x00400000, 0x00200000, 0x00100000,
//	0x00080000, 0x00040000, 0x00020000, 0x00010000,
//	0x00008000, 0x00004000, 0x00002000, 0x00001000,
//	0x00000800, 0x00000400, 0x00000200, 0x00000100,
//	0x00000080, 0x00000040, 0x00000020, 0x00000010,
//	0x00000008, 0x00000004, 0x00000002, 0x00000001,
//};
//
//static const u32 U32_MASK0[32] = {
//	0x7FFFFFFF, 0xBFFFFFFF, 0xDFFFFFFF, 0xEFFFFFFF,
//	0xF7FFFFFF, 0xFBFFFFFF, 0xFDFFFFFF, 0xFEFFFFFF,
//	0xFF7FFFFF, 0xFFBFFFFF, 0xFFDFFFFF, 0xFFEFFFFF,
//	0xFFF7FFFF, 0xFFFBFFFF, 0xFFFDFFFF, 0xFFFEFFFF,
//	0xFFFF7FFF, 0xFFFFBFFF, 0xFFFFDFFF, 0xFFFFEFFF,
//	0xFFFFF7FF, 0xFFFFFBFF, 0xFFFFFDFF, 0xFFFFFEFF,
//	0xFFFFFF7F, 0xFFFFFFBF, 0xFFFFFFDF, 0xFFFFFFEF,
//	0xFFFFFFF7, 0xFFFFFFFB, 0xFFFFFFFD, 0xFFFFFFFE,
//};
//
////__forceinline__ int  GetBitDomIndex(int var_id)
////{
////	return var_id * BITDOM_INTSIZE;
////}
//
//// 根据x和index获得bitDom位置
//#define GetBitDomIndex(x, i) (x * BITDOM_INTSIZE + i)
//// 根据落在最后文字的值的个数获取bit表示的偏移量
//#define GetOffSet(x)(U32_BIT - (x & U32_MOD_MASK))
//
//#define GetBitSubDomStartIndex(x,a)((x * MAX_DOM_SIZE + a) * BITDOMS_INTSIZE)
//#define GetBitSubDomIndex(x, a, y, i)(GetBitSubDomStartIndex(x,a) + GetBitDomIndex(y, i))
//
//__device__ bool IsGtZero(int x)
//{
//	return x > 0;
//}
//
//__inline__ __device__ __host__ int GetTopNum(int num_elements, int num_threads)
//{
//	return (num_elements + (num_threads - 1)) / num_threads;
//}
////************************************
//// Method:    intsizeof
//// FullName:  intsizeof
//// Access:    public 
//// Returns:   int
//// Qualifier: 获取用bit表示的int长度
//// Parameter: const int x
////************************************
//inline int intsizeof(const int x)
//{
//	return (int)ceil((float)x / U32_BIT);
//}
//
//__device__ __inline__ int pow2i(int e)
//{
//	return 1 << e;
//}
////__global__ void enforceACMain(u32* bitDom, u32* bitSup, u32* M_Con, u32* M_ConEvt, u32* M_ConPre, u32* M_VarPre)
////{
////
////}
//
////通过已改变变量生成约束队列
//__global__ void GenConPre(int *VarPre, int* BlocksCount, int3* scp, int len)
//{
//	const int idx = blockDim.x*blockIdx.x + threadIdx.x;
//	if (idx < len)
//	{
//		int3 sp = scp[idx];
//		int pred;
//		if (VarPre[sp.x] == 1 || VarPre[sp.y] == 1)
//			pred = 1;
//		else
//			pred = 0;
//
//		int BC = __syncthreads_count(pred);
//
//		if (threadIdx.x == 0)
//		{
//			BlocksCount[threadIdx.x] = BC;
//		}
//	}
//}
//
//__global__ void CompactQ(int *VarPre, int* ConEvt, int* BOffset, int3* scp, int len)
//{
//	int idx = threadIdx.x + blockIdx.x*blockDim.x;
//	extern __shared__ int warpTotals[];
//	//一个线程块内有128个线程
//	//一个块内有4个线程束
//	if (idx < len)
//	{
//		int3 sp = scp[idx];
//		int pred;
//		//获得判定
//		if (VarPre[sp.x] == 1 || VarPre[sp.y] == 1)
//			pred = 1;
//		else
//			pred = 0;
//
//		//warp index
//		//线程束索引
//		int w_i = threadIdx.x / warpSize;
//		//thread index within a warp
//		//线程束内线程索引
//		int w_l = idx % warpSize;
//		//thread mask (ERROR IN THE PAPERminus one is required)
//		//线程掩码
//		//INT_MAX = 1111 1111 1111 1111 1111 1111 1111 1111 
//		//若线程内id=0，右移32-0-1 = 31位 右侧剩下1位
//		//若线程内id=5，右移32-5-1 = 26位 右侧剩下6位
//		//若线程内id=31，右移32-31-1 = 0位 右侧剩下32位
//		//线程束内threid|  31~~~~~~0
//		//ballot对应位置|   1......1
//		int t_m = INT_MAX >> (warpSize - w_l - 1);
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
//		if (w_l == warpSize - 1)
//			warpTotals[w_i] = t_u + pred;
//
//		__syncthreads();
//
//		//线程束id为0，线程束内线程id，若blockDim.x = 128，则w_l < 128/32 = 4
//		//线程块内第一个线程束的前（4）个线程束工作，w_l < 活动线程束数（4），即每个线程束被一个线程运行
//		if (w_i == 0 && w_l < blockDim.x / warpSize)
//		{
//			int w_i_u = 0;
//			for (int j = 0; j <= 5; ++j)
//			{
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
//
//		if (pred)
//			ConEvt[t_u + warpTotals[w_i] + BOffset[blockIdx.x]] = scp[idx].z;
//
//	}
//}
//
//void CompactQueueMain()
//{
//	//以约束数量启动
//	//P1
//	GenConPre << <MCC_BLOCK, num_threads >> >(M_VarPre, MCC_BlocksCount, scope, CS_SIZE);
//	cudaDeviceSynchronize();
//	//P2
//	thrust::exclusive_scan(MCC_BCount.begin(), MCC_BCount.end(), MCC_BOffset.begin());
//	cudaDeviceSynchronize();
//	//P3
//	//每个约束一个线程进行归约,共享内存大小 = 一个块内线程束的个数,用来装载线程束计算结果
//	CompactQ << <MCC_BLOCK, num_threads, sizeof(int)*(num_threads / warpSize) >> >(M_VarPre, M_Con, MCC_BlocksOffset, scope, CS_SIZE);
//}
//
//#define GetBitSupIndexByINTPrstn(cid,x_val,y_val) (cid * BITSUP_INTSIZE + x_val * BITDOM_INTSIZE + y_val)
//
//__inline__ __device__ __host__ int2 GetBitSupIndexByTuple(int cid, int2 t)
//{
//	return make_int2(cid * BITSUP_INTSIZE + t.x * BITDOM_INTSIZE + (t.y >> U32_POS), cid * BITSUP_INTSIZE + t.y * BITDOM_INTSIZE + (t.x >> U32_POS));
//}
//
//__inline__ __device__ __host__ int GetBitSupIndexById(int cid)
//{
//	return cid * BITSUP_INTSIZE;
//}
//
//void DelGPUModel();
//
//void BuildBitModel(XModel *xm)
//{
//#pragma region 计算常量
//	BITDOM_INTSIZE = intsizeof(xm->feature.max_dom_size);
//	MAX_DOM_SIZE = xm->feature.max_dom_size;
//	VS_SIZE = xm->feature.vs_size;
//	CS_SIZE = xm->feature.cs_size;
//	BITDOMS_INTSIZE = BITDOM_INTSIZE * VS_SIZE;
//	BITSUP_INTSIZE = MAX_DOM_SIZE * BITDOM_INTSIZE;
//	BITSUPS_INTSIZE = BITSUP_INTSIZE * CS_SIZE;
//	BITSUBDOMS_INTSIZE = VS_SIZE * MAX_DOM_SIZE * BITDOMS_INTSIZE;
//	SUBCON_SIZE = VS_SIZE * MAX_DOM_SIZE * CS_SIZE;
//#pragma endregion 计算常量
//#pragma region 约束网络信息
//	//cudaMallocManaged(&vars_size, sizeof(int) * VS_SIZE);
//	//// 初始化变量域大小
//	//for (int i = 0; i < xm->feature.vs_size; ++i)
//	//{
//	//	XVar* v = xm->vars[i];
//	//	XDom* d = xm->doms[v->dom_id];
//	//	vars_size[i] = d->size;
//	//}
//
//	// 初始化scope
//	cudaMallocManaged(&scope, sizeof(int3) * CS_SIZE);
//	for (int i = 0; i < CS_SIZE; ++i)
//	{
//		XCon *c = xm->cons[i];
//		scope[i].x = c->scope[0];
//		scope[i].y = c->scope[1];
//		scope[i].z = c->id;
//	}
//
//	////显示
//	//for (int i = 0; i < CS_SIZE; ++i)
//	//{
//	//	printf("scope[%d] = {%d, %d}\n", scope[i].z, scope[i].x, scope[i].y);
//	//}
//#pragma endregion 约束网络信息
//#pragma region 拷贝bitDom
//	cudaMallocManaged(&bitDom, sizeof(u32) * BITDOMS_INTSIZE);
//	cudaMallocManaged(&M_VarPre, sizeof(int) * VS_SIZE);
//
//	for (int i = 0; i < VS_SIZE; ++i)
//	{
//		XVar* v = xm->vars[i];
//		XDom* d = xm->doms[v->dom_id];
//		const int dom_size = d->size;
//		// 当前变量的实际INT长度
//		const int dom_int_size = intsizeof(dom_size);
//
//		for (int j = 0; j < BITDOM_INTSIZE; ++j)
//		{
//			const int idx = GetBitDomIndex(i, j);
//			//printf("idx = %d\n", idx);
//			// 三种情况
//			if (j < dom_int_size - 1)
//				bitDom[idx] = UINT32_MAX;
//			else if (j == dom_int_size - 1)
//				bitDom[idx] = UINT32_MAX << GetOffSet(dom_size);
//			else
//				bitDom[idx] = 0;
//		}
//
//		M_VarPre[i] = 1;
//	}
//
//	//for (int i = 0; i < VS_SIZE; ++i)
//	//{
//	//	for (int j = 0; j < BITDOM_INTSIZE; ++j)
//	//	{
//	//		int idx = GetBitDomIndex(i, j);
//	//		printf("var = %d, j = %d, idx = %d, bitDom = %x, pre= %x\n", i, j, idx, bitDom[idx], M_VarPre[i]);
//	//	}
//	//}
//#pragma endregion 拷贝bitDom
//#pragma region 创建bitSubDom
//	cudaMallocManaged(&bitSubDom, sizeof(u32)*BITDOMS_INTSIZE*VS_SIZE*MAX_DOM_SIZE);
//	for (int i = 0; i < VS_SIZE; ++i)
//	{
//		for (int j = 0; j < MAX_DOM_SIZE; ++j)
//		{
//			const int start_idx = GetBitSubDomStartIndex(i, j);
//			for (int k = 0; k < BITSUBDOMS_INTSIZE; ++k)
//				bitSubDom[start_idx + k] = bitDom[k];
//			//最后将(i,j)的bitDom 改掉
//			//获取i,j,i的起始地址，
//			const int ijistart = start_idx + i*BITDOM_INTSIZE;
//			for (int k = 0; k < BITDOM_INTSIZE; ++k)
//				// j在索引K的范围内:j/32,将第j%32位置为1
//				if (k == j >> U32_POS)
//					bitSubDom[ijistart + k] = U32_MASK1[j&U32_MOD_MASK];
//			//其它位置为0
//				else
//					bitSubDom[ijistart + k] = 0;
//		}
//	}
//
//	//for (int i = 0; i < VS_SIZE; ++i)
//	//{
//	//	for (int j = 0; j < MAX_DOM_SIZE; ++j)
//	//	{
//	//		printf("sub problem:(%d, %d): ", i, j);
//	//		const int start_idx = GetBitSubDomStartIndex(i, j);
//	//		for (int k = 0; k < BITDOMS_INTSIZE; ++k)
//	//		{
//	//			printf("%x ", bitSubDom[start_idx + k]);
//	//		}
//	//		printf("\n");
//	//	}
//	//}
//#pragma endregion 创建bitSubDom
//#pragma region 拷贝bitSup
//	cudaMallocManaged(&bitSup, sizeof(uint2) * BITSUPS_INTSIZE);
//	for (int i = 0; i < CS_SIZE; ++i)
//	{
//		XCon* c = xm->cons[i];
//		XRel* r = xm->rels[c->rel_id];
//		XVar* v[2] = { xm->vars[c->scope[0]], xm->vars[c->scope[1]] };
//		XDom* d[2] = { xm->doms[v[0]->dom_id], xm->doms[v[1]->dom_id] };
//
//		//初始化位矩阵
//		for (int j = 0; j < MAX_DOM_SIZE; ++j)
//		{
//			for (int k = 0; k < BITDOM_INTSIZE; ++k)
//			{
//				const int idx = GetBitSupIndexByINTPrstn(c->id, j, k);
//				if (j < d[0]->size && (k < (d[1]->size >> U32_POS)))
//				{
//					//支持取0x0000..., 冲突取0xFFF...
//					bitSup[idx].x = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
//					bitSup[idx].y = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
//				}
//				else if (k == (d[1]->size >> U32_POS))
//				{
//					bitSup[idx].x = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
//					bitSup[idx].y = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
//					bitSup[idx].x <<= U32_BIT - (d[1]->size & U32_MOD_MASK);
//					bitSup[idx].y <<= U32_BIT - (d[1]->size & U32_MOD_MASK);
//				}
//				else
//				{
//					bitSup[idx].x = 0;
//					bitSup[idx].y = 0;
//				}
//			}
//		}
//		//向位矩阵中填充值
//		for (int j = 0; j < r->size; ++j)
//		{
//			const int2 t = make_int2(r->tuples[j][0], r->tuples[j][1]);
//			//printf("c_id= %d, %d, %d\n", c->id, t.x, t.y);
//			const int2 idx = GetBitSupIndexByTuple(c->id, t);
//			//printf("idx = %d, %d\n", idx.x, idx.y);
//			if (r->sem == SEM_SUPPORT)
//			{
//				bitSup[idx.x].x |= U32_MASK1[t.y & U32_MOD_MASK];
//				bitSup[idx.y].y |= U32_MASK1[t.x & U32_MOD_MASK];
//			}
//			else
//			{
//				bitSup[idx.x].x &= U32_MASK0[t.y & U32_MOD_MASK];
//				bitSup[idx.y].y &= U32_MASK0[t.x & U32_MOD_MASK];
//			}
//		}
//		//// 初始化位矩阵
//		//for (int j = 0; j < MAX_DOM_SIZE; ++j)
//		//{
//		//	printf("c_id = %d, j = %d: ", i, j);
//		//	for (int k = 0; k < BITDOM_INTSIZE; ++k)
//		//	{
//		//		const int idx = GetBitSupIndexByINTPrstn(c->id, j, k);
//		//		printf("%x, %x", bitSup[idx].x, bitSup[idx].y);
//		//	}
//		//	printf("\n");
//		//}
//	}
//#pragma endregion 拷贝bitSup
//#pragma region 生成约束
//	cudaMallocManaged(&M_Con, sizeof(int)*CS_SIZE);
//	cudaMallocManaged(&M_ConEvt, sizeof(int) * CS_SIZE);
//	cudaMallocManaged(&M_ConPre, sizeof(int)*CS_SIZE);
//
//	for (int i = 0; i < CS_SIZE; ++i)
//	{
//		M_Con[i] = i;
//		M_ConEvt[i] = i;
//		M_ConPre[i] = 1;
//		//printf("i = %d , M_Con = %d, M_ConEvt = %d, M_ConPre = %d\n", i, M_Con[i], M_ConEvt[i], M_ConPre[i]);
//	}
//#pragma endregion 生成约束
//#pragma region 子问题约束队列
//	cudaMallocManaged(&S_ConPre, sizeof(int)*SUBCON_SIZE);
//	cudaMallocManaged(&S_ConEvt, sizeof(int3)*SUBCON_SIZE);
//	cudaMallocManaged(&S_Con, sizeof(int3)*SUBCON_SIZE);
//	cudaMallocManaged(&S_VarPre, sizeof(int)*VS_SIZE*MAX_DOM_SIZE*VS_SIZE);
//	cudaMallocManaged(&S_Var, sizeof(int3)*VS_SIZE*MAX_DOM_SIZE*VS_SIZE);
//
//	for (int i = 0; i < VS_SIZE; ++i)
//	{
//		const XVar* v = xm->vars[i];
//		for (int j = 0; j < MAX_DOM_SIZE; ++j)
//		{
//			for (int k = 0; k < CS_SIZE; ++k)
//			{
//				// 子问题(i, j) k为约束id
//				const int idx = (i*MAX_DOM_SIZE + j)*CS_SIZE + k;
//				//i*xm->feature.max_dom_size*xm->feature.cs_size + j*xm->feature.cs_size + k;
//
//				S_Con[idx].x = i;
//				S_Con[idx].y = j;
//				S_Con[idx].z = k;
//
//				S_ConEvt[idx].x = i;
//				S_ConEvt[idx].y = j;
//				S_ConEvt[idx].z = k;
//
//				S_ConPre[idx] = 1;
//				//printf("S_Con = (%d, %d, %d), S_ConEvt = (%d, %d, %d), pre = %d\n", S_Con[idx].x, S_Con[idx].y, S_Con[idx].z, S_ConEvt[idx].x, S_ConEvt[idx].y, S_ConEvt[idx].z, S_ConPre[idx]);
//			}
//
//			for (int k = 0; k < VS_SIZE; ++k)
//			{
//				//子问题(i, j) k为变量id
//				const int idx = (i*MAX_DOM_SIZE + j)*VS_SIZE + k;
//				S_Var[idx].x = i;
//				S_Var[idx].y = j;
//				S_Var[idx].z = k;
//
//				S_VarPre[idx] = 1;
//
//				//printf("S_Var = (%d, %d, %d), S_VarPre = %d\n", S_Var[idx].x, S_Var[idx].y, S_Var[idx].z, S_VarPre[idx]);
//			}
//		}
//	}
//#pragma endregion 子问题约束队列
//
//#pragma region 程序运行规格
//	//获得主问题压缩的BLOCK数
//	MCC_BLOCK = GetTopNum(CS_SIZE, num_threads);
//	MCC_BCount.resize(MCC_BLOCK, 0);
//	MCC_BOffset.resize(MCC_BLOCK, 0);
//	MCC_BlocksCount = thrust::raw_pointer_cast(MCC_BCount.data());
//	MCC_BlocksOffset = thrust::raw_pointer_cast(MCC_BOffset.data());
//#pragma endregion
//
//}
//
//__global__ void ConCheckMain(int* ConEvt, int* btSp, int2* scp)
//{
//	const int c_id = blockIdx.x;
//	//获取约束在bitSup的开始索引
//	const int start_idx = GetBitSupIndexById(c_id);
//	const int2 sp = scp[c_id];
//	extern __shared__ int2[];
//}
//
//void ConstraintCheckMain()
//{
//	//num_threads最好可变
//	ConCheckMain << <CS_SIZE, num_threads >> >();
//}
//
//float SACGPU()
//{
//	//1. 在主问题上执行AC
//	//1.1. 流压缩
//	CompactQueueMain();
//	//1.2. 约束检查
//	ConstraintCheckMain();
//}
//
//void DelGPUModel()
//{
//	cudaFree(scope);
//	cudaFree(bitDom);
//	cudaFree(M_VarPre);
//	cudaFree(bitSubDom);
//	cudaFree(bitSup);
//	cudaFree(M_Con);
//	cudaFree(M_ConEvt);
//	cudaFree(M_ConPre);
//	cudaFree(S_ConPre);
//	cudaFree(S_ConEvt);
//	cudaFree(S_Con);
//	cudaFree(S_Var);
//	cudaFree(S_VarPre);
//}
//
//
