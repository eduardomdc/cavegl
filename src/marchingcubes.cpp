#include "marchingcubes.hpp"
#include "glm/ext/vector_float3.hpp"
#include "marchingtable.hpp"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/geometric.hpp"
#include "glm/trigonometric.hpp"
#include <iostream>
#include <vector>
#include <omp.h>

int mcIndex(float values[8], float level){
    int index = 0;
	if (values[0] > level) index |= 1;
	if (values[1] > level) index |= 2;
	if (values[2] > level) index |= 4;
	if (values[3] > level) index |= 8;
	if (values[4] > level) index |= 16;
	if (values[5] > level) index |= 32;
	if (values[6] > level) index |= 64;
	if (values[7] > level) index |= 128;
	return index;
}

std::vector<float> mcVertex(float vertex[8], float level){
	std::vector<float> verts;
	int idx = mcIndex(vertex, level);
	//std::cout<<"mcVertex:: idx="<<idx<<std::endl;
	const int* edges = triTable[idx];
	for (int i=0; edges[i] != -1; i++){
		const int* edgeVertices = edgeConnections[edges[i]];
		int vert1 = edgeVertices[0];
		glm::vec3 offset1 = cornerOffsets[vert1];
		int vert2 = edgeVertices[1];
		glm::vec3 offset2 = cornerOffsets[vert2];
		/*
		std::cout<<"offset1::"<<offset1[0]<<","<<offset1[1]<<","<<offset1[2]<<std::endl;
		std::cout<<"offset2::"<<offset2[0]<<","<<offset2[1]<<","<<offset2[2]<<std::endl;
		std::cout<<"value at 1::"<<vertex[vert1]<<std::endl;
		std::cout<<"value at 2::"<<vertex[vert2]<<std::endl;*/
		// now interpolate between two vertices offsets
		//glm::vec3 final = (offset1 + offset2)*0.5f;//((level-vertex[vert1]) * (offset2-offset1) / (vertex[vert2]-vertex[vert1]));
		glm::vec3 final = offset1 + (offset2-offset1)*(level-vertex[vert1])/(vertex[vert2]-vertex[vert1]);
		//std::cout<<"final::"<<final[0]<<","<<final[1]<<","<<final[2]<<std::endl;
		verts.push_back(final[0]);
		verts.push_back(final[1]);
		verts.push_back(final[2]);
	}
	return verts;
}

std::vector<float> marchingcubes_thread(std::vector<std::vector<std::vector<float>>> points, float level, int nthreads, int offset){
	std::vector<float> newVerts;
	for (int x=offset; x<points.size()-1; x+=nthreads){
		for(int y=0; y<points[x].size()-1; y++){
			for(int z=0; z<points[x][y].size()-1; z++){
				float values[8] = {
					points[x][y][z+1],
					points[x+1][y][z+1],
					points[x+1][y][z],
					points[x][y][z],
					points[x][y+1][z+1],
					points[x+1][y+1][z+1],
					points[x+1][y+1][z],
					points[x][y+1][z]
				};
				std::vector<float> mcVerts = mcVertex(values, level);
				// sum x,y,z to vertices offsets
				for (int i=0; i<mcVerts.size(); i+=3){
					mcVerts[i] += x;
					mcVerts[i+1] += y;
					mcVerts[i+2] += z;
				}
				//std::cout<<"mcVerts size="<<mcVerts.size()<<std::endl;
				//newVerts.insert(newVerts.end(), mcVerts.begin(), mcVerts.end());
				for (int i=0; i<mcVerts.size(); i+=3){
					newVerts.push_back(mcVerts[i]);
					newVerts.push_back(mcVerts[i+1]);
					newVerts.push_back(mcVerts[i+2]);
					if ((i/3)%3==0){
					newVerts.push_back(0.0);// texture coords
					newVerts.push_back(0.0);}
					else if ((i/3)%2==0){
					newVerts.push_back(0.5);// texture coords
					newVerts.push_back(1.0);}
					else {
					newVerts.push_back(1.0);// texture coords
					newVerts.push_back(0.0);}
				}
			}
		}
	}
	return newVerts;
}


std::vector<float> parallelmarchingcubes(std::vector<std::vector<std::vector<float>>> points, float level, int nthreads){
	std::vector<std::vector<float>> threadResults(nthreads); // To store results from each thread
    
    #pragma omp parallel num_threads(nthreads)
    {
        int thread_id = omp_get_thread_num();
        threadResults[thread_id] = marchingcubes_thread(points, level, nthreads, thread_id);
    }

    std::vector<float> newVerts;
    for (const auto& result : threadResults) {
        newVerts.insert(newVerts.end(), result.begin(), result.end());
    }

    return newVerts;
}