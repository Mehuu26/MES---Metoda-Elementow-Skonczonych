#include "pch.h"
#include "Header.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>


using namespace std;

int main()
{
	GlobalData global;
	int cp = 0;
	int ro = 0;
	int npc = 0;
	int alfa = 0;
	int numer = 1;	//numer linii
	double t0 = 0;
	int t_alfa = 0;
	int delta_tal = 0;
	int nt = 0;
	string y;		//zmienna pomocnicza
	fstream wczytaj;
	wczytaj.open("wczytaj.txt", ios::in);

	if (wczytaj.good() == false) {
		cout << "Nie udalo sie otworzyc pliku" << endl;
		exit(0);
	}

	while (getline(wczytaj, y)) {

		switch (numer) {
		case 1: global.H = atof(y.c_str()); break;
		case 2: global.W = atof(y.c_str()); break;
		case 3: global.mW = atof(y.c_str()); break;
		case 4: global.mH = atof(y.c_str()); break;
		case 5: npc = atoi(y.c_str()); break;
		case 6: cp = atoi(y.c_str()); break;
		case 7: ro = atoi(y.c_str()); break;
		case 8: alfa = atoi(y.c_str()); break;
		case 9: t0 = atoi(y.c_str()); break;
		case 10: t_alfa = atoi(y.c_str()); break;
		case 11: delta_tal = atoi(y.c_str()); break;
		case 12: nt = atoi(y.c_str()); break;
		}
		numer++;
	}

	cout << "H = " << global.H << endl;
	cout << "W = " << global.W << endl;
	cout << "liczba wezlow po wysokosci: mH = " << global.mH << endl;
	cout << "liczba wezlow po szerokosci: mW = " << global.mW << endl;
	cout << "NPC = " << npc << endl;
	cout << "cp = " << cp << endl;
	cout << "ro = " << ro << endl;
	cout << "afla = " << alfa << endl;
	cout << "t0 = " << t0 << endl;
	cout << "t_alfa = " << t_alfa << endl;
	cout << "delta_tal = " << delta_tal << endl;
	cout << "nt = " << nt << endl;

	int mN = global.mH * global.mW;				//wielkosc naszej siatki (jak prostokat), 3*5=15 elementów, liczba node'ów
	global.mN = mN;
	cout << "mN: " << mN << endl;
	
	int mE = (global.mH - 1)*(global.mW - 1);	//liczba elementow
	global.mE = mE;
	cout << "mE: " << mE << endl;

	global.c = cp;
	global.ro = ro;
	global.alfa = alfa;
	global.NPC = npc;

	int t_symulacji = nt / delta_tal;
	cout << endl << "Czas symulacji: " << t_symulacji << endl;

	element * tab; //wskaznik do tablicy
	tab = new element[mE];	//tablica dynamiczna o rozmiarze mE - liczba elementow

	element Elem;

	//element tab[30];

	//int *tab;
	int j = 1;


	//cout << "global mE = " << mE << endl;
	//cout << "global mN = " << mN << endl;
	for (int i = 0; i < mE; i++) {			//tworzenie siatki, kazdy element sklada sie z 4 wezlow, ktore tutaj tworzymy
		//if ((i + 1) == 6)j++;				//przeskakujemy na robienie wezlow obok
		//else if ((i + 1) == 11)j++;

		Elem.id[0] = j;
		Elem.id[1] = j + global.mH;
		Elem.id[2] = j + 1 + global.mH;
		Elem.id[3] = j + 1;

		tab[i] = Elem;

		j++;

		//tab[i].wyswietl();

		if ((i + 1) % (global.mH - 1) == 0) j++;

		//cout << "tab[" << i << "] = " << tab[i];
	}

	float deltaX, deltaY;

	deltaX = global.W / (global.mW - 1);
	deltaY = global.H / (global.mH - 1);


	//cout << "delta x " << deltaX << endl;
	//cout << "delta y " << deltaY << endl;

	node * tab1;	//tablica trzymajaca kazde x y dla kazdego wezla
	tab1 = new node[mN];	//tablica dynamiczna

	//node tab1[30];
	node Node;
	int h = 0;


	for (int i = 0; i < global.mW; i++) {	//wezly, dokladne ich polozenie (x,y)
		for (int p = 0; p < global.mH; p++) {
			Node.x = i * deltaX;
			Node.y = p * deltaY;

			tab1[h] = Node;

			//tab1[h].wyswietlnode();
			h++;
		}
	}

	oblicz();

	//cout << " global h " << global.H << endl;
	//cout << " global w " << global.W << endl;

	for (int i = 0; i < global.mN; i++) {	//przypisywanie flagi BC do kazdego node
	/*	cout << "x = " << tab1[i].x << endl;
		cout << "y = " << tab1[i].y << endl;
*//*
		cout << "tab[" << i << "].x = " << tab1[i].x << endl;
		cout << "tab[" << i << "].y = " << tab1[i].y << endl;
*/
		float xtemp = round(tab1[i].x * 10000) / 10000;
		float ytemp = round(tab1[i].y * 10000) / 10000;

		//if (tab1[i].y < global.H) cout << "y jest mniejszy niz global h";
		//if (tab1[i].y > global.H) cout << "y jest wiekszy niz global h";

		if (xtemp == 0 || xtemp == global.W || ytemp == 0 || ytemp == global.H) {
			tab1[i].BC = 1;
		}
		else
		{
			tab1[i].BC = 0;
		}
		//cout << "flaga bc dla " << i << " wezla node : " << tab1[i].BC << endl << endl;
	}

	int iteracja = 0; //dla kazdego elementu skonczonego, aby wpisalo dobrze wartosci do struktury

	double* t_min;
	t_min = new double[t_symulacji];

	double* t_max;
	t_max = new double[t_symulacji];

	int iteration_number = 0;

	double *t1;
	t1 = new double[mN];

	for (int i = 0; i < mN; i++) {
		t1[i] = t0;
	}

	for (int l = 0; l < t_symulacji; l++) {	//wywolanie t_symulacji razy
		for (int i = 0; i < mE; i++) {
			int w1 = 0, w2 = 0, w3 = 0, w4 = 0; //tymczasowe zmienne przechowujace numer wezla
			double tempx1, tempy1, tempx2, tempy2, tempx3, tempy3, tempx4, tempy4, tempbc1, tempbc2, tempbc3, tempbc4;
			
			w1 = tab[i].id[0] - 1;	//minus jeden, poniewaz pierwszy element tablicy to 0
			w2 = tab[i].id[1] - 1;
			w3 = tab[i].id[2] - 1;
			w4 = tab[i].id[3] - 1;

			tempx1 = tab1[w1].x;
			tempy1 = tab1[w1].y;
			tempbc1 = tab1[w1].BC;

			tempx2 = tab1[w2].x;
			tempy2 = tab1[w2].y;
			tempbc2 = tab1[w2].BC;

			tempx3 = tab1[w3].x;
			tempy3 = tab1[w3].y;
			tempbc3 = tab1[w3].BC;

			tempx4 = tab1[w4].x;
			tempy4 = tab1[w4].y;
			tempbc4 = tab1[w4].BC;

			//cout << endl << "-----------------------------------------------------------------------------------------------------------------------------------------------------";
			//cout << endl << "dla elementu skonczonego  " << i + 1 << endl;

			//t0 = t1[i];
			jacobi(tempx1, tempx2, tempx3, tempx4, tempy1, tempy2, tempy3, tempy4, tab, npc, cp, ro, tempbc1, tempbc2, tempbc3, tempbc4, alfa, deltaX, deltaY, t0, t_alfa, iteracja);
			iteracja++;
			//wypisywanie macierzy H dla danego elementu skonczonego
			//cout << endl << endl << "Macierz H dla  " << i + 1 << " elementu skonczonego = " << endl;
			/*for (int h = 0; h < 4; h++) {
				for (int j = 0; j < 4; j++) {
					cout << tab[i].matrix_h[h][j] << "  ";
				}
				cout << endl;
			}*/
		}

		//Global H
		//cout << "Globalna macierz H" << endl << endl;

		//deklaracja tablicy dynamicznej
		double** HG = new double*[mN];

		//drugi wymiar, dla kazdego elementu tworzymy tablcie mN elementowa (24)
		for (int i = 0; i < mN; i++) {
			HG[i] = new double[mN];
		}

		//inicjalizujemy na start zerami
		for (int i = 0; i < mN; i++) {
			for (int j = 0; j < mN; j++) {
				HG[i][j] = 0.0;
			}
		}

		for (int k = 0; k < mE; k++) {	//mE tyle ile elementow w siatce (15), 15 macierzy lokalnych H 
			for (int i = 0; i < 4; i++) {	//macierz h jest wymiaru 4x4
				for (int j = 0; j < 4; j++) {
					HG[tab[k].id[i] - 1][tab[k].id[j] - 1] += tab[k].matrix_h[i][j];
				}
			}
		}

		//wypisywanie, uzycie set, aby ladnie bylo widac w terminalu. 
		for (int i = 0; i < mN; i++) {
			for (int j = 0; j < mN; j++) {
				cout << setfill(' ') << setw(5) << setprecision(4) << HG[i][j] << "\t";
			}
			cout << endl;
		}

		//cout << endl << endl << "Macierz Globalna C" << endl;

		//deklaracja tablicy dynamicznej
		double** CG = new double*[mN];

		//drugi wymiar, dla kazdego elementu tworzymy tablcie mN elementowa (24)
		for (int i = 0; i < mN; i++) {
			CG[i] = new double[mN];
		}

		//inicjalizujemy na start zerami
		for (int i = 0; i < mN; i++) {
			for (int j = 0; j < mN; j++) {
				CG[i][j] = 0.0;
			}
		}


		for (int k = 0; k < mE; k++) {		//mE tyle ile elementow w siatce (15), 15 macierzy lokalnych H 
			for (int i = 0; i < 4; i++) {	//macierz h jest wymiaru 4x4
				for (int j = 0; j < 4; j++) {
					CG[tab[k].id[i] - 1][tab[k].id[j] - 1] += tab[k].macierz_c[i][j];
				}

			}
		}

		////wypisywanie, uzycie set, aby ladnie bylo widac w terminalu. 
		//for (int i = 0; i < mN; i++) {
		//	for (int j = 0; j < mN; j++) {
		//		cout << setfill(' ') << setw(5) << setprecision(4) << CG[i][j] << "\t";
		//	}
		//	cout << endl;
		//}

		//wektor globalny PG

		//cout << endl << endl << "Wektor globalny PG" << endl;

		//deklaracja tablicy dynamicznej
		double* PG = new double[mN];

		//inicjalizujemy na start zerami
		for (int i = 0; i <= mN; i++) {
			PG[i] = 0.0;
		}


		for (int k = 0; k < mE; k++) {		//mE tyle ile elementow w siatce (15), 15 macierzy lokalnych H 
			for (int i = 0; i < 4; i++) {	//macierz h jest wymiaru 4x4
				PG[tab[k].id[i]] += tab[k].macierz_p[0][i];
			}
		}

		////wypisywanie, uzycie set, aby ladnie bylo widac w terminalu. 
		//for (int i = 1; i <= mN; i++) {
		//	cout << fixed << setprecision(2) << PG[i] << endl;
		//	cout << endl;
		//}
		//cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		//cout << endl << "pg[nM]" << PG[mN] << endl;

		//rozwiazanie rownania 

		//tymczasowa macierz HG
		//cout << "HG Temp: " << endl;
		//deklaracja tablicy dynamicznej
		double** HG_temp = new double*[mN];

		//drugi wymiar, dla kazdego elementu tworzymy tablcie mN elementowa (24)
		for (int i = 0; i < mN; i++) {
			HG_temp[i] = new double[mN];
		}

		for (int i = 0; i < mN; i++) {
			for (int j = 0; j < mN; j++) {
				HG_temp[i][j] = HG[i][j] + (CG[i][j] / delta_tal);
			}
		}

		////wypisywanie, uzycie set, aby ladnie bylo widac w terminalu. 
		//for (int i = 0; i < mN; i++) {
		//	for (int j = 0; j < mN; j++) {
		//		cout << setfill(' ') << setw(5) << setprecision(4) << HG_temp[i][j] << "\t";
		//	}
		//	cout << endl;
		//}

		//cout << "PG Temp: " << endl;

		//deklaracja tablicy dynamicznej
		double* PG_temp = new double[mN];

		//inicjalizujemy na start zerami
		for (int i = 0; i < mN; i++) {
			PG_temp[i] = 0.0;
		}


		for (int i = 0; i < mN; i++) {
			PG_temp[i] = -PG[i + 1];
			for (int j = 0; j < mN; j++) {
				PG_temp[i] += ((CG[i][j] / delta_tal) * t1[j]);
			}
		}


		////wypisywanie, uzycie set, aby ladnie bylo widac w terminalu. 
		//for (int i = 0; i < mN; i++) {
		//	cout << fixed << setprecision(2) << PG_temp[i] << endl;
		//	cout << endl;
		//}

		t1 = gauss(mN, HG_temp, PG_temp);
		
		t_min[iteration_number] = t1[0]; //pierwszą liczbę przypisujemy do zmiennej min

		for (int i = 1; i < mN; i++) //przeszukanie pozostałych liczb
			if (t_min[iteration_number] > t1[i])
				t_min[iteration_number] = t1[i];
		cout << "Iteration number: " << iteration_number << "  ";
		cout << "Najmniejsza wczytana liczba to " << t_min[iteration_number];

		t_max[iteration_number] = t1[0]; //pierwszą liczbę przypisujemy do zmiennej min

		for (int i = 1; i < mN; i++) //przeszukanie pozostałych licz
			if (t_max[iteration_number] < t1[i])
				t_max[iteration_number] = t1[i];

		cout << "  Najwieksza wczytana liczba to " << t_max[iteration_number] << endl;



		iteracja = 0;	//iteracja dla kazdego elementu skonczonego

		iteration_number++;	//inkrementowanie dla t_min oraz t_max

		//zerowanie macierzy, przygotowanie do ponownego wykonania
		//for (int m = 0; m < mN; m++) {
		//	for (int n = 0; n < mN; n++) {
		//		HG[m][n] = { 0.0 };
		//		CG[m][n] = { 0.0 };


		//	}

		//	PG[m] = { 0.0 };
		//}

	}
}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
