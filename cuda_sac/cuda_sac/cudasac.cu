////#pragma once
//// CUDA Runtime
////#include <cuda_runtime.h>
////#include <device_functions.h>
////#include <device_launch_parameters.h>
//// Utilities and system includes
////#include <helper_cuda.h>
////#include <helper_functions.h>
////#include <thrust/device_vector.h>
////#include <thrust/scan.h>
////#include <iostream>
////#include <math.h>
////#include "XModel.h"
////#include "math_functions.hpp"
////#include "cuda_runtime_api.h"
////
////using namespace cudacp;
////
////#ifndef MIN
////#define MIN(x,y) ((x < y) ? x : y)
////#endif
////// num_threads
////static const int num_threads = 128;
////static const int U32_SIZE = sizeof(u32); ///<4
////static const int U32_BIT = U32_SIZE * 8;	///<32
////static const int U32_POS = 5;
////static const int U32_MOD_MASK = 31;
////
////struct int_predicate
////{
////	__host__ __device__	bool operator()(const int x)
////	{
////		return x > 0;
////	}
////};
////
////// һ��bitDom[x]�ĳ���
////__constant__ int D_BITDOM_INTSIZE;
////// ����bitDom�ĳ���
////
////__constant__ int D_BITDOMS_INTSIZE;
////
////// bit֧�֣� uint2 bitSup[c][a][idx].x = bitSup[c,x,a,idx]
////// bit֧�֣� uint2 bitSup[c][a][idx].y = bitSup[c,y,a,idx]
//////__device__ uint2*** d_bitSup;
//////__host__ uint2*** h_bitSup;
//////__device__ u32** d_bitDom;
//////__host__ u32** h_bitDom;
//////__device__ u32**
////__constant__ int D_NUM_BD_BLOCK;
////__constant__ int D_NUM_CS_SIZE_BLOCKS;
//////////////////////////////////////////////////////////////////////////////
//////	һЩGPU����
//////////////////////////////////////////////////////////////////////////////
////__device__ __managed__ int NUM_BD_BLOCK;
////__device__ __managed__ int NUM_CS_SIZE_BLOCKS;
////// һ�����������int����
////__device__ __managed__ int BITDOM_INTSIZE;
////// �����������ϵ������int���� 
////__device__ __managed__ int BITDOMS_INTSIZE;
////// һ��Լ����bitsup��int���� 
////__device__ __managed__ int BITSUP_INTSIZE;
////// ����Լ�����ϵ�bitsup��int����
////__device__ __managed__ int BITSUPS_INTSIZE;
//////�����������ܳ���
////__device__ __managed__ int BITSUBDOMS_INTSIZE;
//////��������
////__device__ __managed__ int VS_SIZE;
//////Լ������
////__device__ __managed__ int CS_SIZE;
//////������Լ��ѹ��BLOCK��
////__device__ __managed__ int MCC_BLOCK;
//////////////////////////////////////////////////////////////////////////////
//////	һЩGPU����
//////////////////////////////////////////////////////////////////////////////
////__device__ __managed__ int M_Qsize;
////
//////////////////////////////////////////////////////////////////////////////
//////  GPUԼ����¼��Ϣ�����ɸ���
//////////////////////////////////////////////////////////////////////////////
//////	ÿ��������int��С
////__device__ __managed__ int *vars_size;
////// �洢Լ����scope������int3��scope.x: x.id; scope.y: y.id; scope.z: c.id
////__device__ __managed__ int3* scope;
////// ���dom
////__device__ __managed__ int MAX_DOM_SIZE;
////// subCon����
////__device__ __managed__ int SUBCON_SIZE;
////
////
//////__device__ __managed__ int BITDOM_SIZE;
//////__device__ __managed__ int
////// ���������ݽṹ��ʹ��UM
////// ��ʾԼ����������
////__device__ __managed__ u32* bitDom;
////// ��ʾԼ���������޸�
////__device__ __managed__ uint2* bitSup;
////////���ƶ��У��洢Լ��id
//////__device__ __managed__ int *mainCon;
////////���������ݽṹ
//////��ʾ�������Լ����������
////__device__ __managed__ u32* bitSubDom;
////////���ƶ��У��洢������Լ��id subCon.x: variable��subCon.y: value��subCon.z: c.id
//////__device__ __managed__ ushort3* subCon;
//////���������������Ƿ�ɾ������ʼ��ȫ��Ϊ1
////__device__ __managed__ int* M_VarPre;
//////���������Լ���Ƿ����飬��ʼ��ȫ��Ϊ1
////__device__ __managed__ int* M_ConPre;
//////������Լ����������(ѹ����)
////__device__ __managed__ int* M_ConEvt;
//////������Լ����������
////__device__ __managed__ int* M_Con;
////
//////���������������Ƿ�ɾ������ʼ��ȫ��Ϊ1
////__device__ __managed__ int* S_VarPre;
//////��¼�����������
////__device__ __managed__ uint3* S_Var;
//////���������Լ���Ƿ����飬��ʼ��ȫ��Ϊ1
////__device__ __managed__ int* S_ConPre;
//////������Լ����������(ѹ����)
////__device__ __managed__ int3* S_ConEvt;
//////������Լ����������
////__device__ __managed__ int3* S_Con;
//////
////int* MCC_BlocksCount;
////int* MCC_BlocksOffset;
////
////thrust::device_vector<int> MCC_BCount;
////thrust::device_vector<int> MCC_BOffset;
//////
//////__device__ __managed__ ushort4* subVar;
////////��������ⷢ���Ķ��ı���id����ʼ��ȫ��Ϊ0
//////__device__ __managed__ unsigned short* subEvtVar;
////////��������ⷢ���Ķ���Լ��id����ʼ��ȫ��Ϊ1
//////__device__ __managed__ int* subEvtCon;
////
////static const u32 U32_MASK1[32] = {
////	0x80000000, 0x40000000, 0x20000000, 0x10000000,
////	0x08000000, 0x04000000, 0x02000000, 0x01000000,
////	0x00800000, 0x00400000, 0x00200000, 0x00100000,
////	0x00080000, 0x00040000, 0x00020000, 0x00010000,
////	0x00008000, 0x00004000, 0x00002000, 0x00001000,
////	0x00000800, 0x00000400, 0x00000200, 0x00000100,
////	0x00000080, 0x00000040, 0x00000020, 0x00000010,
////	0x00000008, 0x00000004, 0x00000002, 0x00000001,
////};
////
////static const u32 U32_MASK0[32] = {
////	0x7FFFFFFF, 0xBFFFFFFF, 0xDFFFFFFF, 0xEFFFFFFF,
////	0xF7FFFFFF, 0xFBFFFFFF, 0xFDFFFFFF, 0xFEFFFFFF,
////	0xFF7FFFFF, 0xFFBFFFFF, 0xFFDFFFFF, 0xFFEFFFFF,
////	0xFFF7FFFF, 0xFFFBFFFF, 0xFFFDFFFF, 0xFFFEFFFF,
////	0xFFFF7FFF, 0xFFFFBFFF, 0xFFFFDFFF, 0xFFFFEFFF,
////	0xFFFFF7FF, 0xFFFFFBFF, 0xFFFFFDFF, 0xFFFFFEFF,
////	0xFFFFFF7F, 0xFFFFFFBF, 0xFFFFFFDF, 0xFFFFFFEF,
////	0xFFFFFFF7, 0xFFFFFFFB, 0xFFFFFFFD, 0xFFFFFFFE,
////};
////
//////__forceinline__ int  GetBitDomIndex(int var_id)
//////{
//////	return var_id * BITDOM_INTSIZE;
//////}
////
////// ����x��index���bitDomλ��
////#define GetBitDomIndex(x, i) (x * BITDOM_INTSIZE + i)
////// ��������������ֵ�ֵ�ĸ�����ȡbit��ʾ��ƫ����
////#define GetOffSet(x)(U32_BIT - (x & U32_MOD_MASK))
////
////#define GetBitSubDomStartIndex(x,a)((x * MAX_DOM_SIZE + a) * BITDOMS_INTSIZE)
////#define GetBitSubDomIndex(x, a, y, i)(GetBitSubDomStartIndex(x,a) + GetBitDomIndex(y, i))
////
////__device__ bool IsGtZero(int x)
////{
////	return x > 0;
////}
////
////__inline__ __device__ __host__ int GetTopNum(int num_elements, int num_threads)
////{
////	return (num_elements + (num_threads - 1)) / num_threads;
////}
//////************************************
////// Method:    intsizeof
////// FullName:  intsizeof
////// Access:    public 
////// Returns:   int
////// Qualifier: ��ȡ��bit��ʾ��int����
////// Parameter: const int x
//////************************************
////inline int intsizeof(const int x)
////{
////	return (int)ceil((float)x / U32_BIT);
////}
////
////__device__ __inline__ int pow2i(int e)
////{
////	return 1 << e;
////}
//////__global__ void enforceACMain(u32* bitDom, u32* bitSup, u32* M_Con, u32* M_ConEvt, u32* M_ConPre, u32* M_VarPre)
//////{
//////
//////}
////
//////ͨ���Ѹı��������Լ������
////__global__ void GenConPre(int *VarPre, int* BlocksCount, int3* scp, int len)
////{
////	const int idx = blockDim.x*blockIdx.x + threadIdx.x;
////	if (idx < len)
////	{
////		int3 sp = scp[idx];
////		int pred;
////		if (VarPre[sp.x] == 1 || VarPre[sp.y] == 1)
////			pred = 1;
////		else
////			pred = 0;
////
////		int BC = __syncthreads_count(pred);
////
////		if (threadIdx.x == 0)
////		{
////			BlocksCount[threadIdx.x] = BC;
////		}
////	}
////}
////
////__global__ void CompactQ(int *VarPre, int* ConEvt, int* BOffset, int3* scp, int len)
////{
////	int idx = threadIdx.x + blockIdx.x*blockDim.x;
////	extern __shared__ int warpTotals[];
////	һ���߳̿�����128���߳�
////	һ��������4���߳���
////	if (idx < len)
////	{
////		int3 sp = scp[idx];
////		int pred;
////		����ж�
////		if (VarPre[sp.x] == 1 || VarPre[sp.y] == 1)
////			pred = 1;
////		else
////			pred = 0;
////
////		warp index
////		�߳�������
////		int w_i = threadIdx.x / warpSize;
////		thread index within a warp
////		�߳������߳�����
////		int w_l = idx % warpSize;
////		thread mask (ERROR IN THE PAPERminus one is required)
////		�߳�����
////		INT_MAX = 1111 1111 1111 1111 1111 1111 1111 1111 
////		���߳���id=0������32-0-1 = 31λ �Ҳ�ʣ��1λ
////		���߳���id=5������32-5-1 = 26λ �Ҳ�ʣ��6λ
////		���߳���id=31������32-31-1 = 0λ �Ҳ�ʣ��32λ
////		�߳�����threid|  31~~~~~~0
////		ballot��Ӧλ��|   1......1
////		int t_m = INT_MAX >> (warpSize - w_l - 1);
////		balres = number whose ith bit is one if the ith's thread pred is true masked up to the current index in warp
////		�߳��ھֲ�����pred = 1�������밴λ�뵫���˵��������߳�id�ļ�¼��ֻ��������ǰ����ж�
////		int b = __ballot(pred) & t_m;
////		popc count the number of bit one. simply count the number predicated true BEFORE MY INDEX
////		����ֻ���㵱ǰ�߳�������Ӧ��ǰN����λ��֮��
////		��Ϊ�߳���������ɨ��
////		int t_u = __popc(b);
////
////		��ÿ���߳������һ���߳�д�빲���ڴ棬��ӦidΪ�߳���ID�������߳�ID�ӻأ�
////		��������͵�����ֵд�빲���ڴ棬��������͵�ֵû�б�����
////		warpTotals����Ϊ4
////		if (w_l == warpSize - 1)
////			warpTotals[w_i] = t_u + pred;
////
////		__syncthreads();
////
////		�߳���idΪ0���߳������߳�id����blockDim.x = 128����w_l < 128/32 = 4
////		�߳̿��ڵ�һ���߳�����ǰ��4�����߳���������w_l < ��߳�������4������ÿ���߳�����һ���߳�����
////		if (w_i == 0 && w_l < blockDim.x / warpSize)
////		{
////			int w_i_u = 0;
////			for (int j = 0; j <= 5; ++j)
////			{
////				# of the ones in the j'th digit of the warp offsets
////				0->5 6��λ�ã�
////				000 001
////				000 010
////				000 100
////				001 000
////				010 000
////				100 000
////				int b_j = __ballot(warpTotals[w_l] & pow2i(j));
////				w_i_u += (__popc(b_j & t_m)) << j;
////				printf("indice %i t_m=%i,j=%i,b_j=%i,w_i_u=%i\n",w_l,t_m,j,b_j,w_i_u);
////			}
////			warpTotals[w_l] = w_i_u;
////		}
////		__syncthreads();
////
////		if (pred)
////			ConEvt[t_u + warpTotals[w_i] + BOffset[blockIdx.x]] = scp[idx].z;
////
////	}
////}
////
////void CompactQueueMain()
////{
////	//��Լ����������
////	//P1
////	GenConPre << <MCC_BLOCK, num_threads >> >(M_VarPre, MCC_BlocksCount, scope, CS_SIZE);
////	cudaDeviceSynchronize();
////	//P2
////	thrust::exclusive_scan(MCC_BCount.begin(), MCC_BCount.end(), MCC_BOffset.begin());
////	cudaDeviceSynchronize();
////	//P3
////	//ÿ��Լ��һ���߳̽��й�Լ,�����ڴ��С = һ�������߳����ĸ���,����װ���߳���������
////	CompactQ << <MCC_BLOCK, num_threads, sizeof(int)*(num_threads / warpSize) >> >(M_VarPre, M_Con, MCC_BlocksOffset, scope, CS_SIZE);
////}
////
////#define GetBitSupIndexByINTPrstn(cid,x_val,y_val) (cid * BITSUP_INTSIZE + x_val * BITDOM_INTSIZE + y_val)
////
////__inline__ __device__ __host__ int2 GetBitSupIndexByTuple(int cid, int2 t)
////{
////	return make_int2(cid * BITSUP_INTSIZE + t.x * BITDOM_INTSIZE + (t.y >> U32_POS), cid * BITSUP_INTSIZE + t.y * BITDOM_INTSIZE + (t.x >> U32_POS));
////}
////
////__inline__ __device__ __host__ int GetBitSupIndexById(int cid)
////{
////	return cid * BITSUP_INTSIZE;
////}
////
////void DelGPUModel();
////
////void BuildBitModel(XModel *xm)
////{
////#pragma region ���㳣��
////	BITDOM_INTSIZE = intsizeof(xm->feature.max_dom_size);
////	MAX_DOM_SIZE = xm->feature.max_dom_size;
////	VS_SIZE = xm->feature.vs_size;
////	CS_SIZE = xm->feature.cs_size;
////	BITDOMS_INTSIZE = BITDOM_INTSIZE * VS_SIZE;
////	BITSUP_INTSIZE = MAX_DOM_SIZE * BITDOM_INTSIZE;
////	BITSUPS_INTSIZE = BITSUP_INTSIZE * CS_SIZE;
////	BITSUBDOMS_INTSIZE = VS_SIZE * MAX_DOM_SIZE * BITDOMS_INTSIZE;
////	SUBCON_SIZE = VS_SIZE * MAX_DOM_SIZE * CS_SIZE;
////#pragma endregion ���㳣��
////#pragma region Լ��������Ϣ
////	//cudaMallocManaged(&vars_size, sizeof(int) * VS_SIZE);
////	//// ��ʼ���������С
////	//for (int i = 0; i < xm->feature.vs_size; ++i)
////	//{
////	//	XVar* v = xm->vars[i];
////	//	XDom* d = xm->doms[v->dom_id];
////	//	vars_size[i] = d->size;
////	//}
////
////	// ��ʼ��scope
////	cudaMallocManaged(&scope, sizeof(int3) * CS_SIZE);
////	for (int i = 0; i < CS_SIZE; ++i)
////	{
////		XCon *c = xm->cons[i];
////		scope[i].x = c->scope[0];
////		scope[i].y = c->scope[1];
////		scope[i].z = c->id;
////	}
////
////	////��ʾ
////	//for (int i = 0; i < CS_SIZE; ++i)
////	//{
////	//	printf("scope[%d] = {%d, %d}\n", scope[i].z, scope[i].x, scope[i].y);
////	//}
////#pragma endregion Լ��������Ϣ
////#pragma region ����bitDom
////	cudaMallocManaged(&bitDom, sizeof(u32) * BITDOMS_INTSIZE);
////	cudaMallocManaged(&M_VarPre, sizeof(int) * VS_SIZE);
////
////	for (int i = 0; i < VS_SIZE; ++i)
////	{
////		XVar* v = xm->vars[i];
////		XDom* d = xm->doms[v->dom_id];
////		const int dom_size = d->size;
////		// ��ǰ������ʵ��INT����
////		const int dom_int_size = intsizeof(dom_size);
////
////		for (int j = 0; j < BITDOM_INTSIZE; ++j)
////		{
////			const int idx = GetBitDomIndex(i, j);
////			//printf("idx = %d\n", idx);
////			// �������
////			if (j < dom_int_size - 1)
////				bitDom[idx] = UINT32_MAX;
////			else if (j == dom_int_size - 1)
////				bitDom[idx] = UINT32_MAX << GetOffSet(dom_size);
////			else
////				bitDom[idx] = 0;
////		}
////
////		M_VarPre[i] = 1;
////	}
////
////	//for (int i = 0; i < VS_SIZE; ++i)
////	//{
////	//	for (int j = 0; j < BITDOM_INTSIZE; ++j)
////	//	{
////	//		int idx = GetBitDomIndex(i, j);
////	//		printf("var = %d, j = %d, idx = %d, bitDom = %x, pre= %x\n", i, j, idx, bitDom[idx], M_VarPre[i]);
////	//	}
////	//}
////#pragma endregion ����bitDom
////#pragma region ����bitSubDom
////	cudaMallocManaged(&bitSubDom, sizeof(u32)*BITDOMS_INTSIZE*VS_SIZE*MAX_DOM_SIZE);
////	for (int i = 0; i < VS_SIZE; ++i)
////	{
////		for (int j = 0; j < MAX_DOM_SIZE; ++j)
////		{
////			const int start_idx = GetBitSubDomStartIndex(i, j);
////			for (int k = 0; k < BITSUBDOMS_INTSIZE; ++k)
////				bitSubDom[start_idx + k] = bitDom[k];
////			//���(i,j)��bitDom �ĵ�
////			//��ȡi,j,i����ʼ��ַ��
////			const int ijistart = start_idx + i*BITDOM_INTSIZE;
////			for (int k = 0; k < BITDOM_INTSIZE; ++k)
////				// j������K�ķ�Χ��:j/32,����j%32λ��Ϊ1
////				if (k == j >> U32_POS)
////					bitSubDom[ijistart + k] = U32_MASK1[j&U32_MOD_MASK];
////			//����λ��Ϊ0
////				else
////					bitSubDom[ijistart + k] = 0;
////		}
////	}
////
////	//for (int i = 0; i < VS_SIZE; ++i)
////	//{
////	//	for (int j = 0; j < MAX_DOM_SIZE; ++j)
////	//	{
////	//		printf("sub problem:(%d, %d): ", i, j);
////	//		const int start_idx = GetBitSubDomStartIndex(i, j);
////	//		for (int k = 0; k < BITDOMS_INTSIZE; ++k)
////	//		{
////	//			printf("%x ", bitSubDom[start_idx + k]);
////	//		}
////	//		printf("\n");
////	//	}
////	//}
////#pragma endregion ����bitSubDom
////#pragma region ����bitSup
////	cudaMallocManaged(&bitSup, sizeof(uint2) * BITSUPS_INTSIZE);
////	for (int i = 0; i < CS_SIZE; ++i)
////	{
////		XCon* c = xm->cons[i];
////		XRel* r = xm->rels[c->rel_id];
////		XVar* v[2] = { xm->vars[c->scope[0]], xm->vars[c->scope[1]] };
////		XDom* d[2] = { xm->doms[v[0]->dom_id], xm->doms[v[1]->dom_id] };
////
////		//��ʼ��λ����
////		for (int j = 0; j < MAX_DOM_SIZE; ++j)
////		{
////			for (int k = 0; k < BITDOM_INTSIZE; ++k)
////			{
////				const int idx = GetBitSupIndexByINTPrstn(c->id, j, k);
////				if (j < d[0]->size && (k < (d[1]->size >> U32_POS)))
////				{
////					//֧��ȡ0x0000..., ��ͻȡ0xFFF...
////					bitSup[idx].x = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
////					bitSup[idx].y = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
////				}
////				else if (k == (d[1]->size >> U32_POS))
////				{
////					bitSup[idx].x = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
////					bitSup[idx].y = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
////					bitSup[idx].x <<= U32_BIT - (d[1]->size & U32_MOD_MASK);
////					bitSup[idx].y <<= U32_BIT - (d[1]->size & U32_MOD_MASK);
////				}
////				else
////				{
////					bitSup[idx].x = 0;
////					bitSup[idx].y = 0;
////				}
////			}
////		}
////		//��λ���������ֵ
////		for (int j = 0; j < r->size; ++j)
////		{
////			const int2 t = make_int2(r->tuples[j][0], r->tuples[j][1]);
////			//printf("c_id= %d, %d, %d\n", c->id, t.x, t.y);
////			const int2 idx = GetBitSupIndexByTuple(c->id, t);
////			//printf("idx = %d, %d\n", idx.x, idx.y);
////			if (r->sem == SEM_SUPPORT)
////			{
////				bitSup[idx.x].x |= U32_MASK1[t.y & U32_MOD_MASK];
////				bitSup[idx.y].y |= U32_MASK1[t.x & U32_MOD_MASK];
////			}
////			else
////			{
////				bitSup[idx.x].x &= U32_MASK0[t.y & U32_MOD_MASK];
////				bitSup[idx.y].y &= U32_MASK0[t.x & U32_MOD_MASK];
////			}
////		}
////		//// ��ʼ��λ����
////		//for (int j = 0; j < MAX_DOM_SIZE; ++j)
////		//{
////		//	printf("c_id = %d, j = %d: ", i, j);
////		//	for (int k = 0; k < BITDOM_INTSIZE; ++k)
////		//	{
////		//		const int idx = GetBitSupIndexByINTPrstn(c->id, j, k);
////		//		printf("%x, %x", bitSup[idx].x, bitSup[idx].y);
////		//	}
////		//	printf("\n");
////		//}
////	}
////#pragma endregion ����bitSup
////#pragma region ����Լ��
////	cudaMallocManaged(&M_Con, sizeof(int)*CS_SIZE);
////	cudaMallocManaged(&M_ConEvt, sizeof(int) * CS_SIZE);
////	cudaMallocManaged(&M_ConPre, sizeof(int)*CS_SIZE);
////
////	for (int i = 0; i < CS_SIZE; ++i)
////	{
////		M_Con[i] = i;
////		M_ConEvt[i] = i;
////		M_ConPre[i] = 1;
////		//printf("i = %d , M_Con = %d, M_ConEvt = %d, M_ConPre = %d\n", i, M_Con[i], M_ConEvt[i], M_ConPre[i]);
////	}
////#pragma endregion ����Լ��
////#pragma region ������Լ������
////	cudaMallocManaged(&S_ConPre, sizeof(int)*SUBCON_SIZE);
////	cudaMallocManaged(&S_ConEvt, sizeof(int3)*SUBCON_SIZE);
////	cudaMallocManaged(&S_Con, sizeof(int3)*SUBCON_SIZE);
////	cudaMallocManaged(&S_VarPre, sizeof(int)*VS_SIZE*MAX_DOM_SIZE*VS_SIZE);
////	cudaMallocManaged(&S_Var, sizeof(int3)*VS_SIZE*MAX_DOM_SIZE*VS_SIZE);
////
////	for (int i = 0; i < VS_SIZE; ++i)
////	{
////		const XVar* v = xm->vars[i];
////		for (int j = 0; j < MAX_DOM_SIZE; ++j)
////		{
////			for (int k = 0; k < CS_SIZE; ++k)
////			{
////				// ������(i, j) kΪԼ��id
////				const int idx = (i*MAX_DOM_SIZE + j)*CS_SIZE + k;
////				//i*xm->feature.max_dom_size*xm->feature.cs_size + j*xm->feature.cs_size + k;
////
////				S_Con[idx].x = i;
////				S_Con[idx].y = j;
////				S_Con[idx].z = k;
////
////				S_ConEvt[idx].x = i;
////				S_ConEvt[idx].y = j;
////				S_ConEvt[idx].z = k;
////
////				S_ConPre[idx] = 1;
////				//printf("S_Con = (%d, %d, %d), S_ConEvt = (%d, %d, %d), pre = %d\n", S_Con[idx].x, S_Con[idx].y, S_Con[idx].z, S_ConEvt[idx].x, S_ConEvt[idx].y, S_ConEvt[idx].z, S_ConPre[idx]);
////			}
////
////			for (int k = 0; k < VS_SIZE; ++k)
////			{
////				//������(i, j) kΪ����id
////				const int idx = (i*MAX_DOM_SIZE + j)*VS_SIZE + k;
////				S_Var[idx].x = i;
////				S_Var[idx].y = j;
////				S_Var[idx].z = k;
////
////				S_VarPre[idx] = 1;
////
////				//printf("S_Var = (%d, %d, %d), S_VarPre = %d\n", S_Var[idx].x, S_Var[idx].y, S_Var[idx].z, S_VarPre[idx]);
////			}
////		}
////	}
////#pragma endregion ������Լ������
////
////#pragma region �������й��
////	//���������ѹ����BLOCK��
////	MCC_BLOCK = GetTopNum(CS_SIZE, num_threads);
////	MCC_BCount.resize(MCC_BLOCK, 0);
////	MCC_BOffset.resize(MCC_BLOCK, 0);
////	MCC_BlocksCount = thrust::raw_pointer_cast(MCC_BCount.data());
////	MCC_BlocksOffset = thrust::raw_pointer_cast(MCC_BOffset.data());
////#pragma endregion
////
////}
////
////__global__ void ConCheckMain(int* ConEvt, int* btSp, int2* scp)
////{
////	const int c_id = blockIdx.x;
////	//��ȡԼ����bitSup�Ŀ�ʼ����
////	const int start_idx = GetBitSupIndexById(c_id);
////	const int2 sp = scp[c_id];
////	extern __shared__ int2[];
////}
////
////void ConstraintCheckMain()
////{
////	//num_threads��ÿɱ�
////	ConCheckMain << <CS_SIZE, num_threads >> >();
////}
////
////float SACGPU()
////{
////	//1. ����������ִ��AC
////	//1.1. ��ѹ��
////	CompactQueueMain();
////	//1.2. Լ�����
////	ConstraintCheckMain();
////}
////
////void DelGPUModel()
////{
////	cudaFree(scope);
////	cudaFree(bitDom);
////	cudaFree(M_VarPre);
////	cudaFree(bitSubDom);
////	cudaFree(bitSup);
////	cudaFree(M_Con);
////	cudaFree(M_ConEvt);
////	cudaFree(M_ConPre);
////	cudaFree(S_ConPre);
////	cudaFree(S_ConEvt);
////	cudaFree(S_Con);
////	cudaFree(S_Var);
////	cudaFree(S_VarPre);
////}
////
////
