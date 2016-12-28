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
	XDom(const int id, const int size, char* values_str);
	virtual ~XDom();
private:
	void GenerateValues(char* values_str);
};

class XVar
{
public:
	int id;
	int dom_id;
	XVar(const int id, const int dom_id) :id(id), dom_id(dom_id) {}
	virtual ~XVar() {}
};

class XRel
{
public:
	int id;
	int arity;
	int size;
	int **tuples;
	Semantices sem;
	XRel(const int id, const int arity, const int size, const Semantices sem, char* tuples_str);
	virtual ~XRel();
private:
	void GenerateTuples(char *tuples_str);
};

class XCon
{
public:
	int id;
	int rel_id;
	int arity;
	int* scope;
	XCon(const int id, const int rel_id, const int arity, char* scope_str);
	~XCon();
};

struct XFeature
{
	int ds_size;
	int vs_size;
	int rs_size;
	int cs_size;
	u32 max_arity;
	u32 max_dom_size;
};

class XModel
{
public:
	XFeature feature;
	XDom* doms;
	XVar* vars;
	XRel* rels;
	XCon* cons;
	XModel() {};
	virtual ~XModel() {};
};
}


