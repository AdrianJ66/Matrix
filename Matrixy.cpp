// Matrixy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "matrix.h"
#include <iostream>
#include <conio.h>
#include <fstream>

using namespace std;

void ReadData( FILE* fin, double** pMatrix, double* b, int nDim );
void PrintToFile( FILE* fout, double *pRes, int nDim );
int CreateMatrix( double*** pTab, int nDim );
void PrintVector( double* pTab, int nDim );


#define _DEBUG_

int main( int arc, char* argv[] )
{
	if ( arc != 2 ) //program uruchamiwamy przez parametr ktory jest nazwa pliku z danymi
	{
		perror( "Start this program, using file containg matrix <input_file>");
		return 1;
	}

	//Otwieranie pliku z danymi do odczytu
	FILE* dane = fopen( argv[1], "r" );
	//FILE* dane = fopen( "dane.txt", "r" );
	if ( dane == NULL )
	{
		perror( "Failed to open the file." );
		return 1;
	}

	//Otwieranie pliku wynikowego do zapisu
	FILE* wyniki = fopen( "../Matrixy/wyniki.txt", "w" );
	if ( dane == NULL )
	{
		perror( "Failed to open the file." );
		return 1;
	}


	//Rozmiar macierzy
	int nDim;
	fscanf( dane, "%d", &nDim );

	//Macierz poczatkowa
	double** pMatrix = NULL;
	if ( !CreateMatrix( &pMatrix, nDim ) )
	{
		perror( "Matrix has not been successfully initialised." );
		return 1;
	}
	
	//Wektor wyrazow wolnych
	double* b = ( double* )malloc( sizeof( double ) * nDim );
	if ( b )
		memset( b, 0, sizeof( double ) );
	else
	{
		perror( "Matrix has not been successfully initialised." );
		return 1;
	}


	//MACIERZ DOPELNIEN
	double** pMatrixCompl = NULL;
	
	if ( !CreateMatrix( &pMatrixCompl, nDim ) )
	{
		perror( "Matrix has not been successfully initialised." );
		return 1;
	}



	//MACIERZ ODWROTNA
	double** pMatrixInv = NULL;

	if ( !CreateMatrix( &pMatrixInv, nDim ) )
	{
		perror( "Matrix has not been successfully initialised." );
		return 1;
	}

	//MACIERZ WYNIKOWA

	double* pRes = ( double* )malloc( sizeof( double ) * nDim );
	if ( pRes )
		memset( pRes, 0, sizeof( double ) );
	else
	{
		perror( "Matrix has not been successfully initialised." );
		return 1;
	}

	//Odczytanie danych

	ReadData( dane, pMatrix, b, nDim );


#ifdef _DEBUG_
	printf( "Macierz poczatkowa: \n" );
	PrintMatrix( pMatrix, nDim );

	printf( "\n-----------------------------\n" );

	printf( "Wektor wyrazow wolnych: \n" );
	PrintVector( b, nDim );

	printf( "\n-----------------------------\n" );

	printf( "Wyznacznik:	%f", Det(pMatrix, nDim));

	printf( "\n-----------------------------\n" );

#endif // _DEBUG_

	//Wywolanie funkcji

	InverseMatrix( pMatrixInv, pMatrix, nDim, Det( pMatrix, nDim ) );

	LayoutEqu( pMatrixInv, b, pRes, nDim );


#ifdef _DEBUG_
	
	printf( "Macierz dopelnien: \n" );
	PrintMatrix( pMatrixCompl, nDim );

	printf( "\n-----------------------------\n" );

	printf( "Macierz odwrotna: \n" );
	PrintMatrix( pMatrixInv, nDim );

	printf( "\n-----------------------------\n" );

	printf( "Macierz wynikowa: \n" );
	PrintVector( pRes, nDim );

	printf( "\n-----------------------------\n" );

	printf( "Wyniki zostaly takze zapisane do pliku wyniki.txt \n" );


#endif // _DEBUG_

	PrintToFile( wyniki, pRes, nDim );

	fclose( wyniki );
	fclose( dane );
	
	free( pMatrix );
	free( pMatrixInv );
	free( pMatrixCompl );
	free( pRes );
	free( b );

	do
	{
		printf( "\n Wcisnij ESC aby zamknac\n" );
	} while ( _getch() != '\x1B' );


	return 0;
}

void ReadData( FILE* fin, double** pMatrix, double* b, int nDim )
{
	double* pB = b;

	for ( int i = 0; i < nDim; i++ )
	{
		double* pM = *pMatrix++;
		for ( int j = 0; j < nDim ; j++ )
		{
			fscanf( fin, "%lf", pM++ );
		}
		fscanf( fin, "%lf", pB++ );
	}
}

//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------

void PrintVector( double* pTab, int nDim )
{
	for ( int i = 0; i < nDim; i++ )
	{
		printf( "	%lf\n", *pTab++ );
	}
}


//Rozw ukladu rownan( 6x6 ) : Ax = b = > x = A ^ -1 * b  //poczakowo zrob sobie dla cwiczen 3x3


//matrix.h---- - interface modulu matrix.cpp

//int CreateMatrix( double*** pTab, int nSize );
//void DeleteMatrix( double*** pTab, int nSize );
//void InverseMatrix( double** pInv, double** pTab, int nSize, double det );
//double Det( double** pTab, int nSize ); //rekurencja!! gdy zejdziesz do 2x2 to mozesz robic na krzyz , rob rozwiniecie wzgledem zerowego wiersza
//void LayoutEqu( double** pInv, double* pB, double* pRes, int nSize );
//void PrintMatrix( double** pTab, int nDim );

//implementacja modulu matrix.h

//void Complement( double** pTab, double** pTabI, int nRow, int nCol, int nDim ); //wycina z pTabI wiersz o numerze nRows  i wypisuje do pTabO, nie kopuj wierszy i kolumn, uzyj 2x continue
//void CompMatrix( double** pTabD, double** pTab, int nDim );  //liczy macierz dopelnien
//void TransMatrix( double** pTab, int nDim ); //transponowanie macierzy(chyba)

//Program glowny!!

//void ReadData( FILE* fin, double** pMatrix, double* b, int nDim ); //b to wektor wyrazow wolnych

//int main( int arc, char argv[] )
//{
//	if ( arc != 2 ) //program uruchamiwamy przez parametr ktory jest nazwa pliku z danymi
//	{
//		//wiadomosc Usage: macierz uklad rownan <input_file>;
//		return 1;
//	}
//}
//wszystkie wydruki maja byc na warunkowej kompilacji
//otworzyc plik do odczytu
//odczytaj rozmiar nDim
// double** pMx = NULL;
//wykreowac macierz ukladu nDim x nDim 
//wyreowac wektor wyrazow wolnych (nDim)
//wczytac dane (ReadData)
// obliczyc wyznacznik 
//obroc macierz jesli da
//obliczyc rozwiazanie ukladu
//zwolnienie pamieci (usuwanie macierzy)!!

//UWAGI
//pierwsza linia pliku wsadowego zawiera rozmiar tablicy 
//np 
//3
//1 2 3  1
//2 0 0  -1
//-1 -1 -1  1

//wykreowac dynamicznie maciekrz i wketor niewiadomych i wyrazow wolnych
//dostep do elementow tablic wszedzie(poza jednym wyjatkiem) przez wskazniki !!
//tylko standarowe we/wy
//  _CRT_NO_SECURE_WARNINGS  --- to trzeba wpisac w odpowiednie miejsce do projektu, on ma pokazac gdzie...