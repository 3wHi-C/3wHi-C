#include "normalize.h"
using namespace std;
void getMatrix(const std::string matrixFile, std::map<int,std::vector<float *> *> * matrix, std::vector<cell*> * cells)
{
    ifstream f(matrixFile.c_str());
    char st[40];
    float sumall=0;
    memset(st,'\0',20);
    if (!f.is_open())
    {
        cout << 1 << endl;
        exit(EXIT_FAILURE);
    }
    
    f.peek();
    while (! f.eof())
    {
        cell * line = new cell;
        f.getline(st,20,'\t');
        line -> D0 = stoi(string(st));
        f.getline(st,20,'\t');
        line -> D1 = stoi(string(st));
        f.getline(st,20,'\t');
        line -> D2 = stoi(string(st));
        f.getline(st,20,'\n');
        float * value = new float;
        *value = stof(string(st));
        line -> value = value;
        sumall+=*value;
        cells->push_back(line);
        if (matrix -> count(line -> D0) == 0)
        {
            (*matrix)[line -> D0] = new vector<float *>;
        }
        (*matrix)[line -> D0]->push_back(line -> value);
        if (matrix -> count(line -> D1) == 0)
        {
            (*matrix)[line -> D1] = new vector<float *>;
        }
        (*matrix)[line -> D1]->push_back(line -> value);
        if (matrix -> count(line -> D2) == 0)
        {
            (*matrix)[line -> D2] = new vector<float *>;
        }
        (*matrix)[line -> D2]->push_back(line -> value);
        f.peek();
    }
    cout << sumall << endl;
}
void writeNormalizedMatrix(std::string name,std::vector<cell*> * cells)
{
    ofstream fout(("nor_"+name).c_str(),ios_base::trunc|ios_base::out);
    if (!fout.is_open())
    {
        cout << "File untouchable." << endl;
        exit(EXIT_FAILURE);
    }
    for (auto it:(*cells))
    {
        fout << it->D0 << "\t" << it->D1 << "\t" << it ->D2 <<"\t" << fixed << *(it->value) << endl;
    }
}
void writeBias(std::string name,std::map<int,float> *bias)
{
    ofstream fout(("nor_"+name+".bias").c_str(),ios_base::trunc|ios_base::out);
    if (!fout.is_open())
    {
        cout << "File untouchable." << endl;
        exit(EXIT_FAILURE);
    }
    for (auto it:(*bias))
    {
        fout << it.first << "\t" << fixed << it.second  << endl;
    }
}
void ice3DNormalization(std::string matrixFile, int iterNum, float diffExp )
{
    std::map<int,std::vector<float *> *> * matrix = new std::map<int,std::vector<float *> *>;
    std::vector<cell*> * cells = new std::vector<cell*>;
    std::map<int,float> *bias=new std::map<int,float>;
    std::map<int,float> *dbias=new std::map<int,float>;
    std::map<int,float> *oldBias=new std::map<int,float>;

    float sumV;
    float numV;
    float meanV;
    float sumReads;
    float biasDiff;
    int countCells;
    float ratio;
    float tempReadsSum;
    getMatrix(matrixFile,matrix,cells);
    for (auto it:*cells)
    {
        sumReads+=*(it->value);
    }
    for (auto it:(*matrix))
    {
        (*bias)[it.first]=1;
        (*dbias)[it.first]=0;
        (*oldBias)[it.first]=0;
    }

    for (int i=0; i<iterNum; ++i)
    {
        cout << "Iteration:" << i << "\t";
        tempReadsSum=0;
        biasDiff=0;
        for (auto it:(*matrix))
        {
            countCells=0;
            for (auto it1:(*it.second))
            {
                (*dbias)[it.first]+=*it1;
                countCells+=1;
            }
            (*dbias)[it.first]/=countCells;
        }

        numV=0;
        sumV=0;
        for (auto it:*dbias)
        {
            if (it.second>0)
            {
                sumV+=it.second;
                numV+=1;
            }
        }
        meanV=sumV/numV;
        for (auto it:(*dbias))
        {
            if (it.second>0)
            {
                (*dbias)[it.first]/=meanV;
                (*dbias)[it.first]=pow((*dbias)[it.first],1.0/3);
            }
            else
            {
                (*dbias)[it.first]=1;
            }
            (*bias)[it.first]*=(*dbias)[it.first];

        }

        for (auto it:*cells)
        {
            *(it->value)/=((*dbias)[it->D0]*(*dbias)[it->D1]*(*dbias)[it->D2]);
            tempReadsSum+=*(it->value);
        }


        for (auto it:*cells)
        {
            *(it->value)*=(sumReads/tempReadsSum);
        }
        float ratio=pow(tempReadsSum/sumReads,1.0/3);
        for (auto it:(*bias))
        {
            (*bias)[it.first]*=ratio;
        }

        if (i>0)
        {
            for (auto it:(*bias))
            {
                biasDiff+=abs(it.second-(*oldBias)[it.first]);
            }
            cout << "biasDiff" << biasDiff << "|" << diffExp << endl;
            if (biasDiff<diffExp)
            {
                cout << "Finished at Iteration " << i << '.' << endl;
                break;
            }
        } 
        else
        {
            cout << endl;
        }

        for (auto it:(*bias))
        {
            (*oldBias)[it.first]=(*bias)[it.first];
        }
    }
    writeNormalizedMatrix(matrixFile,cells);
    writeBias(matrixFile,bias);
}
