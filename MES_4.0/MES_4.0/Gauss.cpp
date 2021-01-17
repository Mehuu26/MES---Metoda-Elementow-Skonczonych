#pragma once
#include "pch.h"
#include "Header.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <math.h>


using namespace std;


double* gauss(int mN, double** HR, double* PR)
{
	//Tymaczasowa macierz aby scalic HG_temp z PG_temp
	double **a;
	a = new double*[mN];
	for (auto i = 0; i < mN; ++i) a[i] = new double[mN + 1];

	for (auto i = 0; i < mN; i++) {
		for (auto j = 0; j < mN + 1; j++) {
			if (j < mN) a[i][j] = HR[i][j];
			else a[i][j] = PR[i];
		}
	}

	//cout << "Macierz tymczasowa a, przed pivotyzacj¹" << endl;
	//for (int i = 0; i < mN; i++) {
	//	for (int j = 0; j < mN + 1; j++) {
	//		if (j > mN - 1) cout << "| ";
	//		//cout << a[i][j] << "   ";
	//	}
	//	cout << "\n";
	//}

	//zdefiniowany wektor, ktory pozniej zwrocimy
	double *x;
	x = new double[mN];

	//pivotyzacja
	for (int i = 0; i < mN; i++)
		for (int k = i; k < mN; k++)
			if (a[i][i] < a[k][i])
				for (int j = 0; j <= mN; j++) {
					double temp = a[i][j];
					a[i][j] = a[k][j];
					a[k][j] = temp;
				}

	/*cout << endl << "macierz a po pivotyzacji:" << endl;
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mN + 1; j++) {
			if (j > mN - 1) cout << "| ";
			cout << a[i][j] << "   ";
		}
		cout << "\n";
	}*/

	//Rozwiazywanie ukladu rownan
	for (int i = 0; i < mN - 1; i++) // -1 zeby nie przeszlo do ostatniego
		for (int k = i + 1; k < mN; k++)
		{
			double mult = a[k][i] / a[i][i];
			for (int j = 0; j <= mN; j++)
				a[k][j] = a[k][j] - mult * a[i][j]; //ponizej pivotyzacji rowne 0
		}

	/*cout << endl << "macierz a po eliminacji gausa" << endl;
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mN + 1; j++) {
			if (j > mN - 1) cout << "| ";
			cout << a[i][j] << "   ";
		}
		cout << "\n";
	}*/

	for (int i = mN - 1; i >= 0; i--)
	{
		x[i] = a[i][mN];				  //x[i] jest prawa strona rownania w lini [i]
		for (int j = i + 1; j < mN; j++)
			if (j != i)					 
				x[i] = x[i] - a[i][j] * x[j];
		x[i] = x[i] / a[i][i];            
	}

	/*cout << endl << "wektor rozwiazan x:" << endl;
	for (int i = 0; i < mN; i++) cout << x[i] << " ";
	cout << endl;*/


	return x;
}
