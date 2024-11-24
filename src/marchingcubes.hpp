#ifndef CUBES_HPP
#define CUBES_HPP

#include <vector>

int mcIndex(float values[8], float level); // returns marching cubes table index
std::vector<float> mcVertex(float vertex[8], float level);
std::vector<float> parallelmarchingcubes(std::vector<std::vector<std::vector<float>>> points, float level, int nthreads);


#endif