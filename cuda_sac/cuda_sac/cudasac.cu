#pragma once
//// CUDA Runtime
//#include <cuda_runtime.h>
//#include <device_functions.h>
//#include <device_launch_parameters.h>
//// Utilities and system includes
//#include <helper_cuda.h>
//#include <helper_functions.h>
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
//
//static const int WORK_SIZE = 64;
//static const int U32_SIZE = sizeof(u32); ///<4
//static const int U32_BIT = U32_SIZE * 8;	///<32
//static const int U32_POS = 5;
//static const int U32_MOD_MASK = 31;
//
///// A single variable domain needs how many int, it can also be used in locate the index by variable's id;
//__constant__ int D_BITDOM_INTSIZE;
//__constant__ int D_NUM_BD_BLOCK;
//__constant__ int D_NUM_CS_SIZE_BLOCKS;
////////////////////////////////////////////////////////////////////////////
////	һЩGPU����
////////////////////////////////////////////////////////////////////////////
//__device__ __managed__ int NUM_BD_BLOCK;
//__device__ __managed__ int NUM_CS_SIZE_BLOCKS;
//// һ�����������int����
//__device__ __managed__ int BITDOM_INTSIZE;
//// �����������ϵ������int���� 
//__device__ __managed__ int BITDOMS_INTSIZE;
//// һ��Լ����bitsup��int���� 
//__device__ __managed__ int BITSUP_INTSIZE;
//// ����Լ�����ϵ�bitsup��int����
//__device__ __managed__ int BITSUPS_INTSIZE;
////�����������ܳ���
//__device__ __managed__ int BITSUBDOMS_INTSIZE;
//
////////////////////////////////////////////////////////////////////////////
// //GPUԼ����¼��Ϣ�����ɸ���
////////////////////////////////////////////////////////////////////////////
////	ÿ��������int��С
//__device__ __managed__ int *vars_size;
//// �洢Լ����scope������int3��scope.x: x.id; scope.y: y.id; scope.z: c.id
//__device__ __managed__ int3* scope;
//// ���dom
//__device__ __managed__ int MAX_DOM_SIZE;
//// subCon����
//__device__ __managed__ int SUBCON_SIZE;
//
//
////__device__ __managed__ int BITDOM_SIZE;
////__device__ __managed__ int
//// ���������ݽṹ��ʹ��UM
//// ��ʾԼ����������
//__device__ __managed__ u32* bitDom;
//// ��ʾԼ���������޸�
//__device__ __managed__ u32* bitSup;
////���ƶ��У��洢Լ��id
//__device__ __managed__ int *mainCon;
////���������ݽṹ
////��ʾ�������Լ����������
//__device__ __managed__ u32* bitSubDom;
////���ƶ��У��洢������Լ��id subCon.x: variable��subCon.y: value��subCon.z: c.id
//__device__ __managed__ ushort3* subCon;
////�����Ҫ��Լ��Լ��id����ʼȫ��Ϊ1
//__device__ __managed__ int* mainEvtCon;
////��Ƿ����Ķ��ı���id����ʼ��ȫ��Ϊ0
//__device__ __managed__ int* mainEvtVar;
//
//__device__ __managed__ ushort4* subVar;
////��������ⷢ���Ķ��ı���id����ʼ��ȫ��Ϊ0
//__device__ __managed__ unsigned short* subEvtVar;
////��������ⷢ���Ķ���Լ��id����ʼ��ȫ��Ϊ1
//__device__ __managed__ int* subEvtCon;
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
//#define GetBitDomIndex(i,j) (i * BITDOM_INTSIZE + j)
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
//// Qualifier: ��ȡ��bit��ʾ��int����
//// Parameter: const int x
////************************************
//inline int intsizeof(const int x)
//{
//	return (int)ceil((float)x / U32_BIT);
//}
//
//#define GetBitSupIndexByINTPrstn(cid,x_val,y_val) (cid * BITSUP_INTSIZE + x_val * MAX_DOM_SIZE + y_val)
//#define GetBitSupIndexByTuple(cid,t) (cid * BITSUP_INTSIZE + t.x  + t.y / U32_BIT * MAX_DOM_SIZE)
//void DelGPUModel();
//
//void BuildBitModel(XModel *xm)
//{
//#pragma region ���㳣��
//	BITDOM_INTSIZE = intsizeof(xm->feature.max_dom_size);
//	MAX_DOM_SIZE = xm->feature.max_dom_size;
//	BITDOMS_INTSIZE = BITDOM_INTSIZE * xm->feature.vs_size;
//	BITSUP_INTSIZE = xm->feature.max_dom_size * BITDOM_INTSIZE;
//	BITSUPS_INTSIZE = BITSUP_INTSIZE * xm->feature.cs_size;
//	BITSUBDOMS_INTSIZE = xm->feature.vs_size * xm->feature.max_dom_size * BITDOMS_INTSIZE;
//	SUBCON_SIZE = xm->feature.vs_size * xm->feature.ds_size * xm->feature.cs_size;
//#pragma endregion ���㳣��
//#pragma region Լ��������Ϣ
//	cudaMallocManaged(&vars_size, sizeof(int) * xm->feature.vs_size);
//	// ��ʼ���������С
//	for (int i = 0; i < xm->feature.vs_size; ++i)
//	{
//		XVar* v = xm->vars[i];
//		XDom* d = xm->doms[v->dom_id];
//		vars_size[i] = d->size;
//	}
//
//	// ��ʼ��scope
//	cudaMallocManaged(&scope, sizeof(int3) * xm->feature.cs_size);
//	for (int i = 0; i < xm->feature.cs_size; ++i)
//	{
//		XCon *c = xm->cons[i];
//		scope[i].x = c->scope[0];
//		scope[i].y = c->scope[1];
//		scope[i].z = c->id;
//	}
//#pragma endregion Լ��������Ϣ
//#pragma region ����bitDom
//	cudaMallocManaged(&bitDom, sizeof(int) * BITDOMS_INTSIZE);
//	cudaMallocManaged(&mainEvtVar, sizeof(int) * xm->feature.vs_size);
//
//	for (int i = 0; i < xm->feature.vs_size; ++i)
//	{
//		XVar* v = xm->vars[i];
//		XDom* d = xm->doms[v->dom_id];
//		int dom_size = d->size;
//		int dom_int_size = intsizeof(dom_size);
//
//		for (int j = 0; j < BITDOM_INTSIZE; ++j)
//		{
//			int idx = GetBitDomIndex(i, j);
//
//			if (j < dom_int_size - 1)
//			{
//				bitDom[idx] = UINT32_MAX;
//			}
//			else if (j == dom_int_size - 1)
//			{
//				int offset = U32_BIT - (dom_size & U32_MOD_MASK);
//				bitDom[idx] = UINT32_MAX << offset;
//			}
//			else
//			{
//				bitDom[idx] = 0;
//			}
//		}
//
//		mainEvtVar[i] = 1;
//	}
//
//	//for (int i = 0; i < xm->feature.vs_size; ++i)
//	//{
//	//	for (int j = 0; j < BITDOM_INTSIZE; ++j)
//	//	{
//	//		int idx = GetBitDomIndex(i, j);
//	//		printf("var = %d, j = %d, idx = %d, pre= %x\n", i, j, idx, bitDom[idx]);
//	//	}
//	//}
//#pragma endregion ����bitDom
//#pragma region ����bitSubDom
//
//#pragma endregion ����bitSubDom
//#pragma region ����bitCon
//	cudaMallocManaged(&bitSup, sizeof(int) * BITSUPS_INTSIZE);
//	for (int i = 0; i < xm->feature.cs_size; ++i)
//	{
//		XCon* c = xm->cons[i];
//		XRel* r = xm->rels[c->rel_id];
//		XVar* v[2] = { xm->vars[c->scope[0]],xm->vars[c->scope[1]] };
//		XDom* d[2] = { xm->doms[v[0]->dom_id], xm->doms[v[1]->dom_id] };
//
//		// ��ʼ��λ����
//		for (int j = 0; j < BITDOM_INTSIZE; ++j)
//		{
//			for (int k = 0; k < xm->feature.max_dom_size; ++k)
//			{
//				const int idx = GetBitSupIndexByINTPrstn(c->id, j, k);
//				if (k < d[0]->size && (j < (d[1]->size / U32_BIT)))
//				{
//					// ֧��ȡ0x0000...,��ͻȡ0xFFF...
//					bitSup[idx] = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
//				}
//				else if (j == (d[1]->size / U32_BIT))
//				{
//					bitSup[idx] = (r->sem == SEM_CONFLICT) ? UINT32_MAX : 0;
//					bitSup[idx] <<= U32_BIT - (d[1]->size & U32_MOD_MASK);
//				}
//				else
//				{
//					bitSup[idx] = 0;
//				}
//			}
//		}
//
//		// ��λ���������ֵ
//		for (int j = 0; j < r->size; ++j)
//		{
//			int2 t;
//			t.x = r->tuples[j][0];
//			t.y = r->tuples[j][1];
//
//			const int idx = GetBitSupIndexByTuple(c->id, t);
//
//			if (r->sem == SEM_SUPPORT)
//			{
//				bitSup[idx] |= U32_MASK1[t.y & U32_MOD_MASK];
//			}
//			else
//			{
//				bitSup[idx] &= U32_MASK0[t.y & U32_MOD_MASK];
//			}
//		}
//		//// ��ʼ��λ����
//		//for (int j = 0; j < BITDOM_INTSIZE; ++j)
//		//{
//		//	for (int k = 0; k < xm->feature.max_dom_size; ++k)
//		//	{
//		//		const int idx = GetBitSupIndexByINTPrstn(c->id, j, k);
//		//		printf("c_id = %d, j = %d, k = %d, idx = %d, pre= %x\n", i, j, k, idx, bitSup[idx]);
//		//	}
//		//}
//	}
//#pragma endregion ����bitCon
//#pragma region mainCon
//
//	cudaMallocManaged(&mainCon, sizeof(int)*xm->feature.cs_size);
//	cudaMallocManaged(&mainEvtCon, sizeof(int) * xm->feature.cs_size);
//
//	for (int i = 0; i < xm->feature.cs_size; ++i)
//	{
//		mainCon[i] = i;
//		mainEvtCon[i] = 0;
//	}
//
//#pragma endregion mainCon
//#pragma region subCon
//	cudaMallocManaged(&subCon, sizeof(int3)*SUBCON_SIZE);
//	for (int i = 0; i < xm->feature.vs_size; ++i)
//	{
//		XVar* v = xm->vars[i];
//		for (int j = 0; j < xm->feature.max_dom_size; ++j)
//		{
//			for (int k = 0; k < xm->feature.cs_size; ++k)
//			{
//				const int idx = i*xm->feature.max_dom_size*xm->feature.cs_size + j*xm->feature.cs_size + k;
//
//				subCon[idx].x = i;
//				subCon[idx].y = j;
//				subCon[idx].z = k;
//			}
//		}
//	}
//#pragma endregion subCon
//#pragma region ������������������
//#pragma endregion ������������������
//}
//
//void DelGPUModel()
//{
//	cudaFree(vars_size);
//	cudaFree(scope);
//	cudaFree(bitDom);
//	cudaFree(bitSup);
//}



