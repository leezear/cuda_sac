#include "XModel.h"
#include <string>

namespace cudacp
{

XDom::XDom(int id_, int size_, char* val_str_) :
	id(id_), size(size_)
{
	values = new int[size];
	GenerateValues(val_str_);
}

XDom::~XDom()
{
	for (int i = 0; i < size; ++i)
	{
		delete[] values;
		values = NULL;
	}
}

void XDom::GenerateValues(char* values_str)
{
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

	if (blankpos == std::string::npos) {
		sscanf_s(values_str, "%d..%d", &start, &end);
		for (int i = start; i <= end; ++i) {
			values[g_count] = i;
			++g_count;
		}
	}
	else {
		while (blankpos != std::string::npos) {
			tmp = s.substr(0, blankpos - pre_blankpos);
			s = s.substr(blankpos + 1);
			pointpos = tmp.find(".", 0);
			blankpos = s.find(" ", 0);

			if (pointpos != std::string::npos) {
				sscanf_s(tmp.c_str(), "%d..%d", &start, &end);

				for (int i = start; i <= end; ++i) {
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
}

XRel::XRel(const int id_, const int arity_, const int size_, const Semantices sem_, char* tuples_str_)
	:id(id_), arity(arity_), size(size_), sem(sem_)
{
	tuples = new int *[size];

	for (int i = 0; i < size; ++i)
		tuples[i] = new int[arity];

	GenerateTuples(tuples_str_);
}

XRel::~XRel()
{
	for (int i = 0; i < size; ++i)
	{
		delete[] tuples[i];
		tuples[i] = NULL;
	}

	delete[] tuples;
	tuples = NULL;
}

void XRel::GenerateTuples(char *tuples_str_)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < arity; ++j)
			tuples[i][j] = atoi(strtok_s(tuples_str_, " |", &tuples_str_));
}

XCon::XCon(const int id_, const int rel_id_, const int arity_, char * scope_str_)
	:id(id_), rel_id(rel_id_), arity(arity_)
{
	int *scope = new int[arity];

	for (int i = 0; i < arity; ++i)
		scope[i] = atoi(strtok_s(scope_str_, " V", &scope_str_));
}

XCon::~XCon()
{
	delete[] scope;
	scope = NULL;
}

};