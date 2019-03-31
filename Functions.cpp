#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <cstdio>

using namespace std;

void Complement( double** pTabO, double** pTabI, int nRow, int nCol, int nDim );
void CompMatrix( double** pTabD, double** pTab, int nDim );

int CreateMatrix( double*** pTab, int nDim )
{
	double** p = *pTab = ( double** )malloc( sizeof( double* ) * nDim );
	for ( int i = 0; i < nDim; i++, p++ )
	{
		*p = ( double* )malloc( sizeof( double ) * nDim );
		memset( *p, 0, sizeof( double )*nDim );
	}
	return 1;
}


//---------------------------------------------------------------------------------------------------------

void DeleteMatrix( double*** pTab, int nDim )
{
	double** ptr = *pTab;
	for ( int i = 0; i < nDim; i++ )
		free( *ptr++ );
	
	free( *pTab );

	*pTab = NULL;
}

//---------------------------------------------------------------------------------------------------------

void PrintMatrix( double** pTab, int nDim )
{
	for ( int i = 0; i < nDim; i++ )
	{
		double* v = *pTab++;

		for ( int j = 0; j < nDim; j++ )
		{
			printf( "	%f	", *v++ );
		}
		printf( "\n" );
	}
}


//---------------------------------------------------------------------------------------------------------

double Det( double** pTab, int nDim )
{
	if ( nDim == 1 ) return **pTab;
	else if ( nDim == 2 ) return pTab[0][0] * pTab[1][1] - pTab[1][0] * pTab[0][1];
	else
	{
		double** pSTab = NULL;
		CreateMatrix( &pSTab, ( nDim - 1 ) );

		double res = 0;
		int m = 1;
		double* p = *pTab;
		for ( int c = 0; c < nDim; c++ )
		{
			Complement( pSTab, pTab, 0, c, nDim ); //wycinam z pTab wiersz zerowy

			res += *p++ * Det( pSTab, ( nDim - 1 ) ) * m;
			m = -m; //??
		}
		DeleteMatrix( &pSTab, ( nDim - 1 ) );
		return res;
	}
}

//---------------------------------------------------------------------------------------------------------

void CompMatrix( double** pTabD, double** pTab, int nDim )
{
	double** pCompl = NULL;
	CreateMatrix( &pCompl, nDim - 1 );

	double** pD = pTabD;

	for ( int r = 0; r < nDim; r++ )
	{
		double* v = *pD++;
		int m = ( r % 2 ) ? -1 : 1;
		
		for ( int c = 0; c < nDim; c++ )
		{
			Complement( pCompl, pTab, r, c, nDim );		//wycinam rzad i kolumne danego elementu
			*v++ = m * Det( pCompl, ( nDim - 1 ) );
			m = -m;
		}
	}
	DeleteMatrix( &pCompl, nDim - 1 );
}

//---------------------------------------------------------------------------------------------------------

void Complement( double** pTabO, double** pTabI, int nRow, int nCol, int nDim )
{
	for ( int r = 0; r < nDim; r++, pTabI++ )
	{
		double* vOut = *pTabO;
		double* vIn = *pTabI;

		if ( r == nRow ) continue;

		for ( int c = 0; c < nDim; c++, vIn++ )
		{
			if ( c == nCol )continue;
			*vOut++ = *vIn;
		}
		pTabO++;
	}
}

//---------------------------------------------------------------------------------------------------------

void TransMatrix( double** pTab, int nDim )
{
	double** v = pTab;

	for ( int i = 0; i < nDim-1; i++ )
	{
		double* p = *v++ + i + 1;

		for ( int j = i+1; j < nDim; j++, p++ )
		{
			double tmp = *p;
			*p = pTab[j][i];
			pTab[j][i] = tmp;
		}
	}
}

//---------------------------------------------------------------------------------------------------------

void InverseMatrix( double** pInv, double** pTab, int nDim, double det )
{
	
	CompMatrix( pInv, pTab, nDim );

	TransMatrix( pInv, nDim );

	//printf( "\n-----------------------------\n" );
	//printf( "Macierz dopelnien po transponowaniu: \n" );
	//PrintMatrix( pInv, nDim );
	//printf( "\n-----------------------------\n" );

	for ( int r = 0; r < nDim; r++, pInv++ )
	{
		double* pOut = *pInv;

		for ( int c = 0; c < nDim; c++ )
		{
			*pOut++ /= det;
		}

	}
}

//---------------------------------------------------------------------------------------------------------

void LayoutEqu( double** pInv, double* pB, double* pRes, int nDim )
{
	double* pr = pRes;

	for ( int r = 0; r < nDim; r++, pInv++, pr++ )
	{
		double* p = *pInv;
		double* b = pB;
		double sum = 0;

		for ( int c = 0; c < nDim; c++ )
		{
			sum += *p++ * *b++;
		}
		*pr = sum;
	}
}

//---------------------------------------------------------------------------------------------------------

void PrintToFile( FILE* fout, double *pRes, int nDim )
{
	//Wypisaywanie wynikow
	fprintf( fout, "Zestawienie wynikow: \n" );

	for ( int i = 0; i < nDim; i++ )
	{
		fprintf( fout, "%lf\n", *pRes++ );
	}
}

//---------------------------------------------------------------------------------------------------------
