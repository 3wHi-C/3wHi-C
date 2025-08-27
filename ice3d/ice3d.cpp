#include "normalize.h"
#include <string>
using namespace std;
int main(int argc, char *argv[])
{
    string matrixFilename=argv[1];
    int iterNum=stoi(string(argv[2]));
    float diffExp=stof(string(argv[3]));
    ice3DNormalization(matrixFilename,iterNum,diffExp);
    return 0;
}
