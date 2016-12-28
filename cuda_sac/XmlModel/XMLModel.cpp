/*
 * XMLModel.cpp
 *
 *  Created on: 2016年4月7日
 *      Author: leezear
 */

#include "XMLModel.h"

namespace cp {
namespace model {

XMLDomain CreateDomain(int id, int size, char* values_str) {
	XMLDomain domain;
	domain.id = id;
	domain.size = size;
	domain.values = new int[size];
	GenerateValues(values_str, domain.values);

	return domain;
}

XMLVariable CreateVariable(int id, int dm_id) {

	XMLVariable variable;
	variable.id = id;
	variable.dm_id = dm_id;

	return variable;
}

XMLRelation CreateTabular(int id, int size, int arity, Semantices semantices, char *tuple_str) {
	XMLRelation relation;
	relation.id = id;
	relation.size = size;
	relation.arity = arity;
	relation.semantices = semantices;
	relation.tuples = new int*[size];

	for(int i = 0; i < size; ++i) {
		relation.tuples[i] = new int[arity];
	}

	GenerateTuples(tuple_str, relation.tuples, size, arity);

	return relation;
}

XMLPredicate CreatePredicate(int id, char* pas_chr, char* prs_chr) {
	XMLPredicate predicate;
	int max_par_index;
	predicate.id = id;
	std::string par_str;
	par_str = pas_chr;
	par_str = par_str.substr(par_str.rfind("X") + 1);
	max_par_index = atoi(par_str.c_str());
	++max_par_index;
	predicate.num_parameters = max_par_index;
	predicate.parameters = new int[max_par_index];
	for(int i = 0; i < max_par_index; ++i) {
		predicate.parameters[i] = i;
	}
	predicate.propagators_stack = GeneratePStack(prs_chr);
	GeneratePostfixExpression(&predicate.propagators_stack);
	return predicate;
}

XMLConstraint CreateConstraint(const int id, const int arity, int *scope, const int re_id) {
	XMLConstraint constraint;
	constraint.id = id;
	constraint.arity = arity;
	constraint.scope = scope;
	constraint.re_id = re_id;

	return constraint;
}

XMLConstraint CreateConstraint(const int id, const int arity, int *scope, const int pre_id, char* pars_chr,
		const int pars_len) {
	XMLConstraint constraint;
	constraint.id = id;
	constraint.arity = arity;
	constraint.scope = scope;
	constraint.pre_id = pre_id;
	constraint.pars_length = pars_len;
	constraint.pars = GenerateParameters(pars_chr, pars_len);
	return constraint;
}

int *CreateScope(char * scope_str, const int arity) {
	char* ptr;
	int *scope = new int[arity];
	char* context = scope_str;
	const char seps[] = " V";

	for(int i = 0; i < arity; ++i) {
		ptr = strtok_s(context, seps, &context);
		scope[i] = atoi(ptr);
	}
	return scope;
}

bool DestroyEXTModel(XMLModel *model) {
	DestroyDomains(model);
	DestroyVariables(model);
	DestroyRelations(model);
	DestroyConstraints(model);
	delete model;
	model = NULL;
	return true;
}

bool DestroyINTModel(XMLModel *model) {
	DestroyDomains(model);
	DestroyVariables(model);
	DestroyPredicates(model);
	DestroyConstraints(model);
	delete model;
	model = NULL;
	return true;
}

bool CreateDomains(XMLModel *model, int size) {
	model->ds_size = size;
	model->domains = new XMLDomain[size];
	return true;
}

bool CreateVariables(XMLModel *model, int size) {
	model->vs_size = size;
	model->variables = new XMLVariable[size];
	return true;
}

bool CreateRelations(XMLModel *model, int size) {
	model->rs_size = size;
	model->relations = new XMLRelation[size];
	return true;
}

bool CreatePredicates(XMLModel* model, int size) {
	model->ps_size = size;
	model->predicates = new XMLPredicate[size];
	return true;
}

bool CreateConstraints(XMLModel *model, int size) {
	model->cs_size = size;
	model->constraints = new XMLConstraint[size];
	return true;
}

bool DestroyConstraints(XMLModel *model) {
	if(model->features.rel_type == REL_TABULAR) {
		for(int i = 0; i < model->cs_size; ++i) {
			delete[] model->constraints[i].scope;
			model->constraints[i].scope = NULL;
		}
	}
	else if(model->features.rel_type == REL_PREDICATION) {
		if(model->constraints[0].pre_id != -1) {
			for(int i = 0; i < model->cs_size; ++i) {
				delete[] model->constraints[i].pars;
				model->constraints[i].pars = NULL;
			}
		}
	}
	delete[] model->constraints;
	model->constraints = NULL;
	return true;
}

bool DestroyDomains(XMLModel *model) {
	for(int i = 0; i < model->ds_size; ++i) {
		delete[] model->domains[i].values;
		model->domains[i].values = NULL;
	}

	delete[] model->domains;
	model->domains = NULL;

	return true;
}

bool DestroyVariables(XMLModel *model) {
	delete[] model->variables;
	model->variables = NULL;
	return true;
}

bool DestroyRelations(XMLModel *model) {
	for(int i = 0; i < model->rs_size; ++i) {
		const int num_tuples = model->relations[i].size;

		for(int j = 0; j < num_tuples; ++j) {
			delete[] model->relations[i].tuples[j];
			model->relations[i].tuples[j] = NULL;
		}

		delete[] model->relations[i].tuples;
		model->relations[i].tuples = NULL;
	}

	delete[] model->relations;
	model->relations = NULL;

	return true;
}

bool DestroyPredicates(XMLModel *model) {
	for(int i = 0; i < model->ps_size; ++i) {
		delete[] model->predicates[i].parameters;
		model->predicates[i].parameters = NULL;
	}

	delete[] model->predicates;
	model->predicates = NULL;

	return true;
}

void GenerateTuples(char *tuple_str, int **tuples, const int num_tuples, const int arity) {
	char* ptr;
	char* context = tuple_str;
	const char seps[] = " |";
	int i = 0, j = 0;

	for(i = 0; i < num_tuples; ++i) {
		for(j = 0; j < arity; ++j) {
			ptr = strtok_s(context, seps, &context);
			tuples[i][j] = atoi(ptr);
		}
	}
}

bool GenerateValues(char* values_str, int *values) {
	std::string s = values_str;
	s += " ";
	std::string tmp = values_str;
	int start = 0;
	int end = 0;
	int pre_blankpos = 0;
	std::string::size_type blankpos = s.find(" ", 0);
	std::string::size_type pointpos = 0;
	int g_count = 0;
	int value_int;

	if(blankpos == std::string::npos) {
		sscanf_s(values_str, "%d..%d", &start, &end);
		for(int i = start; i <= end; ++i) {
			values[g_count] = i;
			++g_count;
		}
	}
	else {
		while(blankpos != std::string::npos) {
			tmp = s.substr(0, blankpos - pre_blankpos);
			s = s.substr(blankpos + 1);
			pointpos = tmp.find(".", 0);
			blankpos = s.find(" ", 0);

			if(pointpos != std::string::npos) {
				sscanf_s(tmp.c_str(), "%d..%d", &start, &end);

				for(int i = start; i <= end; ++i) {
					values[g_count] = i;
					++g_count;
				}
			}
			else {
				sscanf_s(tmp.c_str(), "%d", &value_int);
				values[g_count] = value_int;
				++g_count;
			}
		}
	}

	return true;
}

PStack GeneratePStack(char* fun_chr) {
	PStack ps;
	int startpos = 0;
	std::string s = fun_chr;
	std::string tmp;
	int op;
	u32 i = 0;
	int j = -1;
	for(i = 0; i < s.length(); ++i) {
		switch(s[i]) {
		case '(':
			tmp = s.substr(startpos, i - startpos);
			op = GetOperator(tmp);
			if(op != PO_NON) {
				ps.ps[++j] = op;
			}
			ps.ps[++j] = PO_L_PAR;
			startpos = i + 1;
			break;
		case ')':
			tmp = s.substr(startpos, i - startpos);
			op = GetOperator(tmp);
			if(op != PO_NON) {
				ps.ps[++j] = op;
			}
			ps.ps[++j] = PO_R_PAR;
			startpos = i + 1;
			break;
		case ',':

			tmp = s.substr(startpos, i - startpos);
			op = GetOperator(tmp);
			if(op != PO_NON) {
				ps.ps[++j] = op;
			}
			startpos = i + 1;
			break;
		case ' ':
			startpos = i + 1;
			break;
		default:
			break;
		}
	}
	ps.num_prs = j + 1;
	ps.top = j;
	ps.bottom = 0;

	return ps;
}
int GetOperator(std::string s) {
	int value;

	if(s == "(") {
		return PO_L_PAR;
	}
	else if(s == ")") {
		return PO_R_PAR;
	}
	else if(s == "div") {
		return PO_DIV;
	}
	else if(s == "mod") {
		return PO_MOD;
	}
	else if(s == "sub") {
		return PO_SUB;
	}
	else if(s == "mul") {
		return PO_MUL;
	}
	else if(s == "abs") {
		return PO_ABS;
	}
	else if(s == "eq") {
		return PO_EQ;
	}
	else if(s == "ne") {
		return PO_NE;
	}
	else if(s == "lt") {
		return PO_LT;
	}
	else if(s == "le") {
		return PO_LE;
	}
	else if(s == "gt") {
		return PO_GT;
	}
	else if(s == "ge") {
		return PO_GE;
	}
	else if(s == "and") {
		return PO_AND;
	}
	else if(s == "or") {
		return PO_OR;
	}
	else if(s[0] == 'X') {
		sscanf_s(s.c_str(), "X%d", &value);
		return value;
	}

	return PO_NON;
}

int GetParameters(const char* par_chr) {
	int value;
	std::string s = par_chr;

	if(s[0] == 'V') {
		sscanf_s(s.c_str(), "V%d", &value);
		value += MAX_VALUE;
	}
	else {
		value = atoi(par_chr);
	}
	return value;
}

void GeneratePostfixExpression(PStack *ps) {
	PStack res_stk;
	int op;
	u32 i = 0, j = 0;
	int last_lpar_idx = 0;
	u32 idx = 0;

	while(i < ps->num_prs) {
		op = ps->ps[i];

		if(op == PO_L_PAR) {
			last_lpar_idx = i;
		}
		else if(op == PO_R_PAR) {
			for(j = last_lpar_idx; j < i; ++j) {
				if(ps->ps[j] > MAX_OPT) {
					push(&res_stk, ps->ps[j]);
					ps->ps[j] = PO_NON;
				}
			}

			ps->ps[last_lpar_idx] = PO_NON;
			idx = last_lpar_idx - 1;
			op = ps->ps[idx];
			push(&res_stk, op);

			while(op != PO_L_PAR && last_lpar_idx > 0) {
				--last_lpar_idx;
				op = ps->ps[last_lpar_idx];
			}
		}
		++i;
	}

	ps->num_prs = res_stk.num_prs;
	ps->top = res_stk.top;
	ps->bottom = res_stk.bottom;

	for(i = 0; i < ps->num_prs; ++i)
		ps->ps[i] = res_stk.ps[i];
}

int* GenerateParameters(char* pars_chr, int pars_len) {
	int* pas = new int[pars_len];
	std::string s = pars_chr;
	s += " ";
	int startpos = 0;
	int j = -1;
	std::string tmp;
//	int op;
	int par;

	for(u32 i = 0; i < s.length(); ++i) {
		switch(s[i]) {
		case ' ':
			tmp = s.substr(startpos, i - startpos);
			par = GetParameters(tmp.c_str());
			pas[++j] = par;
			startpos = i + 1;
			break;
		default:
			break;
		}
	}

	return pas;
}

void push(PStack *ps, int e) {
	if(ps->top <= 200) {
		++ps->num_prs;
		++ps->top;
		ps->ps[ps->top] = e;
	}
}

int pop(PStack *ps) {
	int res = 0;
	if(ps->num_prs >= 0) {
		--ps->num_prs;
		res = ps->ps[ps->top];
		--ps->top;
	}
	return res;
}

} /* namespace model */
} /* namespace cudacsp */
