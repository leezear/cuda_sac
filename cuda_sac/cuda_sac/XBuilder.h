#pragma once
#include <string>
#include "XModel.h"

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOMTreeWalker.hpp>
#include <iostream>
#include <string>
#include <stdio.h>
#include <map>
#include <vector>
#include "XModel.h"

using namespace xercesc_3_1;
typedef xercesc_3_1::DOMDocument DOMDOC;

namespace cudacp {

/**
*\brief some Types of bankmark file type
*/
enum BMFileType {
	BFT_TXT,			///<txt doc
	BFT_XCSP_2_0_INT,	///<xcsp v2.0 intension
	BFT_XCSP_2_1_INT,	///<xcsp v2.1 intension
	BFT_XCSP_3_0_INT,	///<xcsp v3.0 intension
	BFT_XCSP_2_0_EXT,	///<xcsp v2.0 extension
	BFT_XCSP_2_1_EXT,	///<xcsp v2.1 extension
	BFT_XCSP_3_0_EXT	///<xcsp v3.0 extension
};

enum XmlReaderType {
	XRT_BM_PATH,		///<banchmark path file
	XRT_BM				///<banchmark
};

enum DomType {
	standard,
	disperse
};

class DomMap {
public:
	int id;
	std::map<int, int> m;
	DomType dt = standard;

	DomMap() {};
	~DomMap() {};
	void MakeMap(XDom* d) {
		id = d->id;
		for (int i = 0; i < d->size; i++) {
			m[d->values[i]] = i;
			if (d->values[i] != i) {
				dt = disperse;
			}
		}
	}
};


class XBuilder {
public:

	/**
	* \brief Constructors and initialization
	* \param[in] i_file_name	The name of file.
	* \param[in] i_type		The type of file, path file or banckmark file.
	*/
	XBuilder(std::string file_name, XmlReaderType type);

	/**
	* \brief Constructors and initialization
	* \param[out] network	Pointer of model
	* \return model generate result
	*	-<em>false</em> fail
	*	-<em>true</em> succeed
	*/
	void GenerateModelFromXml(XModel *model);

	/**
	* \brief GetBMFile
	* \return Banchmark file path
	*/
	std::string GetBMFile() const;
	virtual ~XBuilder();
private:
	//XFeature feature_;		///<Feature of Network
	XModel *model_;				///<Constarint Network pointer
	std::string file_name_;			///<XML file name(path)
	XmlReaderType type_;			///<file type
	XercesDOMParser *parser_;
	DOMElement *root_;
	DOMDOC *document_;
	std::vector<DomMap> xds;

	bool initial();
	void getFeature();
	int getDomsNum();
	int getVarsNum();
	int getRelsNum();
	int getConsNum();
	int getMaxArity();

	void generateDomains();
	void generateVariables();
	void generateRelations();
	void generateConstraints();
	void modifyTuple();
};
}

