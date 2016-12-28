#pragma once
namespace cudacp
{

typedef unsigned int u32;

enum Semantices {
	//³åÍ»
	SEM_CONFLICT = 0,
	//Ö§³Ö
	SEM_SUPPORT = 1
};

class XDom
{
public:
	int id;
	int size;
	int* values;
	XDom(int id, int size, char* values_str);
private:
	void GenerateValues(char* values_str, int *values);
};

class XVar
{
public:
	int id;
	int dom_id;
};

class XRel
{
public:
	int id;
	int arity;
	int size;
	int **tuple;
	Semantices sem;
};

class XCon
{
public:
	int id;
	int rel_id;
	int arity;
	int* scope;
};


class XModel
{
public:
	XModel()
	{};
	virtual ~XModel() {};
};
}


