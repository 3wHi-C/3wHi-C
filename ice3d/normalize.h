#ifndef NORMALIZE_H_
#define NORMALIZE_H_
#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>

typedef struct
{
    int D0;
    int D1;
    int D2;
    float * value;
} cell;
void ice3DNormalization(std::string matrixFile, int iterNum, float diffExp);
void getMatrix(const std::string matrixFile, std::map<int,std::vector<float *> *> * matrix, std::vector<cell*> cells);
void writeNormalizedMatrix(std::string name,std::vector<cell*> * cells);
void writeBias(std::string name,std::map<int,float> *bias);
#endif