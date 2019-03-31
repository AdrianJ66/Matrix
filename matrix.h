#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
using namespace std;

void DeleteMatrix( double*** pTab, int nDim );
void InverseMatrix( double** pInv, double** pTab, int nSize, double det );
double Det( double** pTab, int nSize ); 
void LayoutEqu( double** pInv, double* pB, double* pRes, int nDim );
void PrintMatrix( double** pTab, int nDim );
int CreateMatrix( double*** pTab, int nDim );



  
void TransMatrix( double** pTab, int nDim ); 
//void TransMatrix( double** TransTab, double** pTab, int nDim );


#endif
