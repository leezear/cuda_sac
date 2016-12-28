/*
 * XMLModel.h
 *
 *  Created on: 2016年6月21日
 *      Author: leezear
 */

#ifndef XMLMODEL_H_
#define XMLMODEL_H_
#include "type_define.h"

namespace cp {
namespace model {

struct XMLDomain {
	int id;
	u32 size;
	int *values;
};

struct XMLVariable {
	int id;
	int dm_id;
};

struct XMLRelation {
	int id;
	int size = 0;
	int arity;
	Semantices semantices;
	int num_parameters = 0;
	int *parameters;
	RelationType type;
	PStack *propagators_stack;
	int **tuples;
};

struct XMLConstraint {
	int id;
	int arity;
	int *scope;
	int re_id = -1;
	int pre_id = -1;
	int pars_length;
	int* pars;
};

struct XMLPredicate {
	int id;
	int *parameters;
	int num_parameters;
	PStack propagators_stack;
};

struct XMLModel {
	int ds_size;
	int vs_size;
	int rs_size;
	int ps_size;
	int cs_size;
	int max_d_size;
	int fs_size = 0;
	NetworkFeatures features;
	XMLDomain *domains;
	XMLVariable *variables;
	PStack *f_stacks;
	XMLRelation *relations;
	XMLConstraint *constraints;
	XMLPredicate *predicates;
};

XMLDomain CreateDomain(int id, int size, char* values_str);
XMLVariable CreateVariable(int id, int dm_id);
XMLRelation CreateTabular(int id, int size, int arity, Semantices semantices, char *tuple_str);
XMLConstraint CreateConstraint(const int id, const int arity, int *scope, const int re_id);
XMLConstraint CreateConstraint(const int id, const int arity, int *scope, const int pre_id, char* pars_chr, const int pars_len);
int *CreateScope(char * scope_str, const int arity);
XMLPredicate CreatePredicate(int id, char* pas_chr, char* prs_chr);
bool CreateDomains(XMLModel *model, int size);
bool CreateVariables(XMLModel *model, int size);
bool CreateRelations(XMLModel *model, int size);
bool CreateConstraints(XMLModel *model, int size);
bool CreatePredicates(XMLModel* model, int size);
XMLModel *CreateModel();
XMLModel *CreateModelFromINT();
bool DestroyEXTModel(XMLModel *model);
bool DestroyINTModel(XMLModel *model);
bool DestroyDomains(XMLModel *model);
bool DestroyVariables(XMLModel *model);
bool DestroyRelations(XMLModel *model);
bool DestroyConstraints(XMLModel *model);
bool DestroyPredicates(XMLModel *model);
bool GenerateValues(char* values_str, int values[]);
void GenerateTuples(char *tuple_str, int **tuples, const int size, const int arity);
void GeneratePostfixExpression(PStack *ps);
int GetOperator(std::string s);
PStack GeneratePStack(char* fun_chr);
int* GenerateParameters(char* pars_chr, int pars_len);
bool GetPar(int par_in, int *par);

void push(PStack *ps, int e);
int pop(PStack *ps);

} /* namespace model */
} /* namespace cudacsp */

#endif /* XMLMODEL_H_ */
