/*
 * XMLBuilder.cpp
 *
 *  Created on: 2016年6月21日
 *      Author: leezear
 */

#include "XMLBuilder.h"

namespace cp {
namespace parse {

XMLBuilder::XMLBuilder(std::string i_file_name, XmlReaderType i_type)
		: file_name_(i_file_name), type_(i_type) {
	if(initial()) {
		switch(type_) {
		case XRT_BM_PATH:
			break;
		case XRT_BM:
			break;
		default:
			break;
		}
	}
}

XMLBuilder::~XMLBuilder() {
	if(root_) {
		delete (parser_);
		parser_ = NULL;
		XMLPlatformUtils::Terminate();
	}
}

bool XMLBuilder::initial() {
	if(this->file_name_ == "") {
		return false;
	}

	try {
		XMLPlatformUtils::Initialize();
	}
	catch(const XMLException & toCatch) {
		// Do your failure processing here
		return false;
	}
	// Do your actual work with Xerces-C++ here.
	parser_ = new XercesDOMParser();
	parser_->setValidationScheme(XercesDOMParser::Val_Always);
	parser_->setDoNamespaces(true);

	parser_->parse(file_name_.c_str());
	//std::cout << file_name << std::endl;
	document_ = parser_->getDocument();
	root_ = document_->getDocumentElement();

	if(!root_) {
		delete (parser_);
		parser_ = NULL;
		return false;
	}
	return true;
}

std::string XMLBuilder::GetBMFile() const {
	std::string bm_file = "";

	if(true) {
		DOMNodeList *bmfile_list = root_->getElementsByTagName(XMLString::transcode("BMFile"));
		DOMNode *bmfile = bmfile_list->item(0);
		char* bm_name = XMLString::transcode(bmfile->getFirstChild()->getNodeValue());
		bm_file = bm_name;
		XMLString::release(&bm_name);
	}

	return bm_file;
}

u32 XMLBuilder::getConstraintsCount() {
	u32 num_cons = 0;
	DOMNode *cons_node = root_->getElementsByTagName(XMLString::transcode("constraints"))->item(0);
	num_cons = (u32) XMLString::parseInt(
			cons_node->getAttributes()->getNamedItem(XMLString::transcode("nbConstraints"))->getTextContent());
	return num_cons;
}

RelationType XMLBuilder::getRelationTpye() {
	RelationType type = REL_TABULAR;
	DOMNodeList * tab_nodes = root_->getElementsByTagName(XMLString::transcode("relation"));
	DOMNodeList * pre_nodes = root_->getElementsByTagName(XMLString::transcode("predicate"));

	if((tab_nodes->getLength()) && (!pre_nodes->getLength())) {
		type = REL_TABULAR;
	}
	else if((!tab_nodes->getLength()) && (pre_nodes->getLength())) {
		type = REL_PREDICATION;
	}

	return type;
}

void XMLBuilder::getNetworkFeature(XMLModel *i_model) {
	feature_.rel_type = getRelationTpye();
	feature_.cons_size = getConstraintsCount();
	i_model->features = feature_;
}

int XMLBuilder::getMaxArity()
{
	u32 max_arity = 0;
	DOMNode *node = root_->getElementsByTagName(XMLString::transcode("presentation"))->item(0);
	max_arity= (u32)XMLString::parseInt(
		node->getAttributes()->getNamedItem(XMLString::transcode("maxConstraintArity"))->getTextContent());
	return max_arity;
}

void XMLBuilder::generateDomains(XMLModel *i_model) {
	u32 max_d_s = 0;
	i_model->features.arity = getMaxArity();
	DOMNode *doms_nodes = root_->getElementsByTagName(XMLString::transcode("domains"))->item(0);
	u32 num_doms = (u32) XMLString::parseInt(
			doms_nodes->getAttributes()->getNamedItem(XMLString::transcode("nbDomains"))->getTextContent());
	DOMNodeList *dom_nodes = root_->getElementsByTagName(XMLString::transcode("domain"));
	u32 dom_id;
	u32 dom_size;
	char * dom_values;
	DOMNode * dom_node;

	CreateDomains(i_model, num_doms);

	for(u32 i = 0; i < num_doms; ++i) {
		dom_node = dom_nodes->item(i);
		dom_id = i;
		dom_size = (u32) XMLString::parseInt(
				dom_node->getAttributes()->getNamedItem(XMLString::transcode("nbValues"))->getTextContent());
		max_d_s = MAX(max_d_s, dom_size);
		dom_values = XMLString::transcode(dom_node->getFirstChild()->getNodeValue());
		i_model->domains[i] = CreateDomain(dom_id, dom_size, dom_values);
	}

	feature_.max_dom_size = max_d_s;
}

void XMLBuilder::generateVariables(XMLModel *i_model) {
	DOMNode *vars_node = root_->getElementsByTagName(XMLString::transcode("variables"))->item(0);
	u32 num_vars = (u32) XMLString::parseInt(
			vars_node->getAttributes()->getNamedItem(XMLString::transcode("nbVariables"))->getTextContent());
	DOMNodeList* var_nodes = root_->getElementsByTagName(XMLString::transcode("variable"));
	CreateVariables(i_model, num_vars);
	feature_.vars_size = num_vars;
	DOMNode* var_node;
	u32 var_id;
	u32 dom_id;
	char* dom_id_str;

	for(u32 i = 0; i < num_vars; ++i) {
		var_id = i;
		var_node = var_nodes->item(i);
		dom_id_str = XMLString::transcode(
				var_node->getAttributes()->getNamedItem(XMLString::transcode("domain"))->getTextContent());
		sscanf_s(dom_id_str, "D%d", &dom_id);
		i_model->variables[i] = CreateVariable(var_id, dom_id);
		XMLString::release(&dom_id_str);
	}
}

void XMLBuilder::generateTabulars(XMLModel * i_model) {
	DOMNode* relations_node = root_->getElementsByTagName(XMLString::transcode("relations"))->item(0);
	int relations_count = (int) XMLString::parseInt(
			relations_node->getAttributes()->getNamedItem(XMLString::transcode("nbRelations"))->getTextContent());
	DOMNodeList* relation_nodes = root_->getElementsByTagName(XMLString::transcode("relation"));
	CreateRelations(i_model, relations_count);
	DOMNode *relation_node;
	char* semantics;
	char* innertext;
	int relation_name;
	int tuple_arity;
	int tuples_count;
	Semantices sem;

	for(int i = 0; i < relations_count; ++i) {
		relation_name = i;
		relation_node = relation_nodes->item(i);
		tuple_arity = (int) XMLString::parseInt(
				relation_node->getAttributes()->getNamedItem(XMLString::transcode("arity"))->getTextContent());
		semantics = XMLString::transcode(
				relation_node->getAttributes()->getNamedItem(XMLString::transcode("semantics"))->getTextContent());
		tuples_count = XMLString::parseInt(
				relation_node->getAttributes()->getNamedItem(XMLString::transcode("nbTuples"))->getTextContent());
		//若属性semantics == supports则relation_type = SURPPOT
		sem = (strlen(semantics) == strlen("supports")) ? SEM_SUPPORT : SEM_CONFLICT;
		innertext = XMLString::transcode(relation_node->getFirstChild()->getNodeValue());
		i_model->relations[i] = CreateTabular(relation_name, tuples_count, tuple_arity, sem, innertext);
	}
}

void XMLBuilder::generatePredicates(XMLModel *i_model) {
	DOMNode *pres_node = root_->getElementsByTagName(XMLString::transcode("predicates"))->item(0);
	int num_pres = XMLString::parseInt(
			pres_node->getAttributes()->getNamedItem(XMLString::transcode("nbPredicates"))->getTextContent());
	DOMNodeList* pre_nodes = root_->getElementsByTagName(XMLString::transcode("predicate"));
	DOMNode* pre_node;
	u32 pre_id;
	CreatePredicates(i_model, num_pres);

	for(int i = 0; i < num_pres; ++i) {
		char* parameters_str;
		char* functional_str;
		pre_id = i;
		pre_node = pre_nodes->item(i);
		parameters_str = XMLString::transcode(pre_node->getChildNodes()->item(1)->getFirstChild()->getNodeValue());
		functional_str = XMLString::transcode(
				pre_node->getChildNodes()->item(3)->getChildNodes()->item(1)->getFirstChild()->getNodeValue());
		i_model->predicates[i] = CreatePredicate(pre_id, parameters_str, functional_str);
		XMLString::release(&parameters_str);
		XMLString::release(&functional_str);
	}
}

void XMLBuilder::generateRelations(XMLModel *i_model) {
	switch(i_model->features.rel_type) {
	case REL_TABULAR:
		generateTabulars(i_model);
		break;
	case REL_PREDICATION:
		generatePredicates(i_model);
		break;
	default:
		break;
	}
}

void XMLBuilder::generateConstraints(XMLModel *i_model) {
	int *con_scope;
	int pre_id;
	int con_id;
	char* pre_id_str;
	char* pars_str;
//	int* pars;
	DOMNode* con_node;
	char *con_scop_str;
	DOMNode* cons_node = root_->getElementsByTagName(XMLString::transcode("constraints"))->item(0);
	u32 num_cons = XMLString::parseInt(
			cons_node->getAttributes()->getNamedItem(XMLString::transcode("nbConstraints"))->getTextContent());
	DOMNodeList* con_nodes = root_->getElementsByTagName(XMLString::transcode("constraint"));
	CreateConstraints(i_model, num_cons);
	int pars_len;

	for(u32 i = 0; i < num_cons; ++i) {
		con_id = i;
		con_node = con_nodes->item(i);
		int con_arity = XMLString::parseInt(
				con_node->getAttributes()->getNamedItem(XMLString::transcode("arity"))->getTextContent());
		con_scop_str = XMLString::transcode(
				con_node->getAttributes()->getNamedItem(XMLString::transcode("scope"))->getTextContent());
		con_scope = CreateScope(con_scop_str, con_arity);
		pre_id_str = XMLString::transcode(
				con_node->getAttributes()->getNamedItem(XMLString::transcode("reference"))->getTextContent());
		switch(i_model->features.rel_type) {
		case REL_TABULAR:
			sscanf_s(pre_id_str, "R%d", &pre_id);
			i_model->constraints[con_id] = CreateConstraint(con_id, con_arity, con_scope, pre_id);
			break;
		case REL_PREDICATION: {
			sscanf_s(pre_id_str, "P%d", &pre_id);
			pars_str = XMLString::transcode(con_node->getChildNodes()->item(1)->getFirstChild()->getNodeValue());
			pars_len = i_model->predicates[pre_id].num_parameters;
			i_model->constraints[con_id] = CreateConstraint(con_id, con_arity, con_scope, pre_id, pars_str, pars_len);
			XMLString::release(&pars_str);
			XMLString::release(&pre_id_str);
			break;
		}
		default:
			break;
		}
	}
}

bool XMLBuilder::GenerateModelFromXml(XMLModel *i_model) {
//	getNetworkFeature();
	generateDomains(i_model);
	generateVariables(i_model);
	getNetworkFeature(i_model);
	generateRelations(i_model);
	generateConstraints(i_model);
	return true;
}

}
}
