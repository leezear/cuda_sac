#include "XModel.h"
#include <string>

namespace cudacp
{

XDom::XDom(int id_, int size_, char* val_str_) :
	id(id_), size(size_)
{
	values = new int[size];
	GenerateValues(val_str_, values);
}

void XDom::GenerateValues(char* values_str, int* values)
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

};