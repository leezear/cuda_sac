#pragma once

#include <iostream>
#include <string>
#include "cudasac.cuh"
#undef DOMDocument

#include "XBuilder.h"

using namespace std;
using namespace cudacp;
const string X_PATH = "BMPath.xml";

int main()
{
	XBuilder path_builder(X_PATH, XRT_BM_PATH);
	string bm_path = path_builder.GetBMFile();
	cout << bm_path << endl;
	XBuilder builder(bm_path, XRT_BM);
	XModel* xmodel = new XModel();
	builder.GenerateModelFromXml(xmodel);
	//float build_time = BuidBitModel32bit(xmodel);
	//float exe_time = SACGPU();
	//DelGPUModel();

	delete xmodel;
	xmodel = NULL;
	//printf("Build time = %f\n", build_time);
	//printf("Execution time = %f\n", exe_time);
	printf("---end---\n");
	return 0;
}

