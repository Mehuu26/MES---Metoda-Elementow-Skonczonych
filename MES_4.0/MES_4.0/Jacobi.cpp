#pragma once
#include "pch.h"
#include "Header.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <math.h>


using namespace std;

double jacobi(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, element H[], int npc, int cp, int ro, int tempbc1, int tempbc2, int tempbc3, int tempbc4, int alfa, float deltaX, float deltaY, double t0, int t_alfa, int iteracja) {
	elem4 X(npc);
	//static int iteracja = 0;

	int liczba_pkt = npc * npc;
	//int npc = npc;
	//cout << "liczba_pkt = " << liczba_pkt << endl;
	//cout << "-------------------------------------------------------------------------------------------------------" << endl;

	//cout << "ksi[0]: " << X.ksi[0];

	double ksi_pc[16][4];
	double eta_pc[16][4];
	/*
	double** ksi_pc;
	double** eta_pc;

	ksi_pc = new double*[liczba_pkt];
	eta_pc = new double*[liczba_pkt];

	for (int i = 0; i < liczba_pkt; i++)
	{
		ksi_pc[i] = new double[4]; //przydzielenie dla ka¿dego wiersza po 4 komorki
	}


	for (int i = 0; i < liczba_pkt; i++)
	{
		eta_pc[i] = new double[4]; //przydzielenie dla ka¿dego wiersza po 4 komorki
	}

	*/


	if (npc == 0) return 0;
	//X.eta[1] = -0.77;
	//X.eta[2] = 0.77;
	//X.eta[3] = 0.77;
	//X.eta[4] = -0.77;

	//X.ksi[1] = -0.77;
	//X.ksi[2] = -0.77;
	//X.ksi[3] = 0.77;
	//X.ksi[4] = 0.77;
/*
	for (int i = 0; i < 4; i++) {
		cout << "ksi[ " << i << "] = " << X.ksi[i] << endl;
		cout << "eta[ " << i << "] = " << X.eta[i] << endl;
	}
*/

/*
for (int i = 0; i < liczba_pkt; i++) {
	for (int j = 0; j < 4; j++) {
		if (j == 0)
			ksi_pc[i][j] = (-0.25*(1 - X.eta[i]));
		if (j == 1)
			ksi_pc[i][j] = (0.25*(1 - X.eta[i]));
		if (j == 2)
			ksi_pc[i][j] = (0.25*(1 + X.eta[i]));
		if (j == 3)
			ksi_pc[i][j] = (-0.25*(1 + X.eta[i]));
	}
}
*/
	for (int i = 0; i < liczba_pkt; i++) {
		ksi_pc[i][0] = (-0.25*(1 - X.eta[i]));
		ksi_pc[i][1] = (0.25*(1 - X.eta[i]));
		ksi_pc[i][2] = (0.25*(1 + X.eta[i]));
		ksi_pc[i][3] = (-0.25*(1 + X.eta[i]));
	}

	//cout << endl << endl;

	for (int i = 0; i < liczba_pkt; i++) {
		eta_pc[i][0] = (-0.25*(1 - X.ksi[i]));
		eta_pc[i][1] = (-0.25*(1 + X.ksi[i]));
		eta_pc[i][2] = (0.25*(1 + X.ksi[i]));
		eta_pc[i][3] = (0.25*(1 - X.ksi[i]));
	}

	//cout << "wypisywanie wartosci w punktach calkowania ksi" << endl;
	////wypisywanie wartosci na ekran ksi
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		if (j == 0)
	//			cout << "pc[" << i << "]:  ";
	//		cout << ksi_pc[i][j] << "  ";
	//	}
	//	cout << endl;

	//}
	//wypisywanie na ekran eta
	/*cout << endl << "wypisywanie wartosci w punktach calkowania eta" << endl;
	for (int i = 0; i < liczba_pkt; i++) {
		for (int j = 0; j < 4; j++) {
			if (j == 0)
				cout << "pc[" << i << "]:  ";
			cout << eta_pc[i][j] << "  ";
		}
		cout << endl;

	}*/

	double x[4] = { x1,x2,x3,x4 };
	double y[4] = { y1,y2,y4,y4 };

	//cout << "wypisywanie wartosci x oraz y" << endl << endl;

	//for (int i = 0; i < 4; i++) {
	//	cout << "x[" << i + 1 << "] = " << x[i] << endl;
	//	cout << "y[" << i + 1 << "] = " << y[i] << endl;

	//}

	double matrix[16][4];

	//macierz jacobina
	for (int i = 0; i < liczba_pkt; i++) {
		matrix[i][0] = pochodna(ksi_pc[i], x);
		matrix[i][1] = pochodna(ksi_pc[i], y);
		matrix[i][2] = pochodna(eta_pc[i], x);
		matrix[i][3] = pochodna(eta_pc[i], y);

	}
	/*cout << endl << "test tablicy matrix " << endl << endl;*/

	////sprawdzenie, wypisanie tablicy
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << matrix[i][j] << "  ";
	//	}
	//	cout << endl;

	//}

	double detJ[16];

	//wyznaczniki macierzy detJ
	for (int i = 0; i < liczba_pkt; i++) {
		detJ[i] = (matrix[i][0] * matrix[i][3]) - (matrix[i][1] * matrix[i][2]);
	}

	////wypisywanie wyznacznikow macierzy detJ
	//cout << endl << endl << "detJ: " << endl;

	//for (int j = 0; j < liczba_pkt; j++) {
	//	cout << "detJ[" << j << "]: " << detJ[j] << endl;
	//}


	////odwracanie macierzy
	//cout << endl << endl << "odwracanie macierzy: " << endl;

	double reverse_matrix[16][4];

	for (int i = 0; i < liczba_pkt; i++) {
		reverse_matrix[i][0] = matrix[i][3];
		reverse_matrix[i][3] = matrix[i][0];
		reverse_matrix[i][1] = matrix[i][1] * (-1);
		reverse_matrix[i][2] = matrix[i][2] * (-1);
	}

	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << reverse_matrix[i][j] << "  ";
	//	}
	//	cout << endl;

	//}

	//liczenie macierzy H

	double dn_dx[16][4];
	double dn_dy[16][4];

	//liczenie macierzy dn_dx, dla 9 punktow calkowania, 1 wiersz to 1 pc ... 4 wiersz to 4 pc. 
	for (int i = 0; i < liczba_pkt; i++) {
		for (int j = 0; j < 4; j++) {
			dn_dx[i][j] = 1 / detJ[i] * (ksi_pc[i][j] * reverse_matrix[i][0] + eta_pc[i][j] * (-reverse_matrix[i][1]));
		}
	}

	//cout << endl << "wypisanie macierzy dn/dx  " << endl;

	////sprawdzenie, wypisanie dn/dx
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << dn_dx[i][j] << "  ";
	//	}
	//	cout << endl;

	//}


	//liczenie macierzy dn_dy, dla 9 punktow calkowania, 1 wiersz to 1 pc ... 4 wiersz to 4 pc. 

	for (int i = 0; i < liczba_pkt; i++) {
		for (int j = 0; j < 4; j++) {
			dn_dy[i][j] = 1 / detJ[i] * (ksi_pc[i][j] * (-reverse_matrix[i][2]) + eta_pc[i][j] * reverse_matrix[i][3]);
		}
	}

	//cout << endl << "wypisanie macierzy dn/dy  " << endl;

	////sprawdzenie, wypisanie dn/dy
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << dn_dy[i][j] << "  ";
	//	}
	//	cout << endl;

	//}

	//mnozenie macierzy i macierzy transponowanej
	//double dn_dx_mn1[4][4];			//przemnozone tablice
	//double dn_dx_mn2[4][4];
	//double dn_dx_mn3[4][4];
	//double dn_dx_mn4[4][4];

	double dn_dx_mn[16][4][4];
	double dn_dy_mn[16][4][4];

	for (int k = 0; k < liczba_pkt; k++) {			//k = dla macierzy 
		for (int i = 0; i < 4; i++) {		//i = dla pc
			for (int j = 0; j < 4; j++) {	//j = dla n
				dn_dx_mn[k][i][j] = dn_dx[k][i] * dn_dx[k][j];
			}
		}
	}

	for (int k = 0; k < liczba_pkt; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				dn_dy_mn[k][i][j] = dn_dy[k][i] * dn_dy[k][j];
			}
		}
	}

	/*
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dn_dx_mn2[0][j] = dn_dx[0][0] * dn_dx[0][j];
			dn_dx_mn2[1][j] = dn_dx[1][0] * dn_dx[1][j];
			dn_dx_mn2[2][j] = dn_dx[2][0] * dn_dx[2][j];
			dn_dx_mn2[3][j] = dn_dx[3][0] * dn_dx[3][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dn_dx_mn3[0][j] = dn_dx[0][0] * dn_dx[0][j];
			dn_dx_mn3[1][j] = dn_dx[1][0] * dn_dx[1][j];
			dn_dx_mn3[2][j] = dn_dx[2][0] * dn_dx[2][j];
			dn_dx_mn3[3][j] = dn_dx[3][0] * dn_dx[3][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dn_dx_mn4[0][j] = dn_dx[0][0] * dn_dx[0][j];
			dn_dx_mn4[1][j] = dn_dx[1][0] * dn_dx[1][j];
			dn_dx_mn4[2][j] = dn_dx[2][0] * dn_dx[2][j];
			dn_dx_mn4[3][j] = dn_dx[3][0] * dn_dx[3][j];
		}
	}

	*/

	//sprawdzenie, wypisanie macierzy 4x4

	//for (int k = 0; k < liczba_pkt; k++) {
	//	cout << endl << endl << "iloczym macierzy dn/dx dla pc[" << k << "]: " << endl;
	//	for (int i = 0; i < 4; i++) {
	//		cout << endl;
	//		for (int j = 0; j < 4; j++) {
	//			cout << dn_dx_mn[k][i][j] << "  ";
	//		}
	//	}
	//}

	//for (int k = 0; k < liczba_pkt; k++) {
	//	cout << endl << endl << "iloczym macierzy dn/dy dla pc[" << k << "]: " << endl;
	//	for (int i = 0; i < 4; i++) {
	//		cout << endl;
	//		for (int j = 0; j < 4; j++) {
	//			cout << dn_dy_mn[k][i][j] << "  ";
	//		}
	//	}
	//}

	//sumowanie macierzy 4x4 dn/dx oraz dn/dy

	double dn_dx_suma[16][4];
	double dn_dy_suma[16][4];

	for (int i = 0; i < liczba_pkt; i++) {
		for (int j = 0; j < 4; j++) {
			dn_dx_suma[i][j] = dn_dx_mn[0][i][j] + dn_dx_mn[1][i][j] + dn_dx_mn[2][i][j] + dn_dx_mn[3][i][j];
			dn_dy_suma[i][j] = dn_dy_mn[0][i][j] + dn_dy_mn[1][i][j] + dn_dy_mn[2][i][j] + dn_dy_mn[3][i][j];
		}
	}


	////wypisanie zsumowanych macierzy
	//cout << endl << endl << "zsumowana macierz dn/dx dla" << endl;
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << dn_dx_suma[i][j] << "  ";
	//	}
	//	cout << endl;
	//}

	//cout << endl << endl << "zsumowana macierz dn/dy dla" << endl;
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << dn_dy_suma[i][j] << "  ";
	//	}
	//	cout << endl;
	//}

	//sumowanie dwoch zsumowanych tablic d{n}/dx oraz d{n}/dy

	double dn_suma[16][4];
	for (int i = 0; i < liczba_pkt; i++) {
		for (int j = 0; j < 4; j++) {
			dn_suma[i][j] = dn_dx_suma[i][j] + dn_dy_suma[i][j];
		}
	}

	////wypisanie zsumowanych macierzy
	//cout << endl << endl << "zsumowana macierz d{n}/dx oraz d{n}/dy" << endl;
	//for (int i = 0; i < liczba_pkt; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << dn_suma[i][j] << "  ";
	//	}
	//	cout << endl;
	//}

	//przypisywanie wartosci do macierzy H

	double matrix_H[4][4] = { 0.0 };
	int kt = 25;	//k(t) ze wzoru




	for (int k = 0; k < liczba_pkt; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				//matrix_H[i][j] += kt * dn_suma[k][j] * detJ[k] * X.w_eta[k] * X.w_ksi[k];
				matrix_H[i][j] += kt * (dn_dx_mn[k][i][j] + dn_dy_mn[k][i][j]) * detJ[k] * X.w_eta[k] * X.w_ksi[k];
			}
		}
	}
	////wypisywanie macierzy H
	//cout << endl << endl << "Macierz H = " << endl;
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << matrix_H[i][j] << "  ";
	//	}
	//	cout << endl;
	//}

	//zapisywanie macierzy h do struktury elementu skonczonego
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[iteracja].matrix_h[i][j] = matrix_H[i][j];
		}
		//cout << endl;
	}

	//tworzenie macierzy C lokalnej

	double N[16][4];

	//tworzenie N - funkcji ksztaltow
	for (int i = 0; i < liczba_pkt; i++) {
		N[i][0] = 0.25*(1 - X.ksi[i]) * (1 - X.eta[i]);
		N[i][1] = 0.25*(1 + X.ksi[i]) * (1 - X.eta[i]);
		N[i][2] = 0.25*(1 + X.ksi[i]) * (1 + X.eta[i]);
		N[i][3] = 0.25*(1 - X.ksi[i]) * (1 + X.eta[i]);
	}



	//przypisywanie wartosci do macierzy C

	double matrix_C[4][4] = { 0.0 };

	for (int k = 0; k < liczba_pkt; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				matrix_C[i][j] += N[k][i] * N[k][j] * detJ[k] * cp * ro * X.w_eta[k] * X.w_ksi[k]; //N * N transponowane
			}
		}
	}

	////wypisywanie macierzy C
	//cout << endl << endl << "Macierz C = " << endl;
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << matrix_C[i][j] << "  ";
	//	}
	//	cout << endl;
	//}

	//zapisywanie macierzy C do struktury elementu skonczonego
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[iteracja].macierz_c[i][j] = matrix_C[i][j];
		}
		//cout << endl;
	}

	//tworzenie macierzy P_l
	//liczenie macierzy [H]lbc

	double ksi[16];
	double eta[16];

	double waga[4];

	//ustawienie wartosci ksi i eta dla warunkow brzegowych
	if(npc==2){	
		waga[0] = waga[1] = 1;
		//ustawianie wartosci ksi i eta dla dolnego warunku brzegowego
		ksi[0] = -(1.0 / sqrt(3));
		ksi[1] = (1.0 / sqrt(3));
		eta[0] = -1;
		eta[1] = -1;

		//ustawianie wartosci ksi i eta dla prawego warunku brzegowego
		ksi[2] = 1;
		ksi[3] = 1;
		eta[2] = -(1.0 / sqrt(3));
		eta[3] = (1.0 / sqrt(3));

		//ustawianie wartosci ksi i eta dla gornego warunku brzegowego
		ksi[4] = (1.0 / sqrt(3));
		ksi[5] = -(1.0 / sqrt(3));
		eta[4] = 1;
		eta[5] = 1;

		//ustawianie wartosci ksi i eta dla lewego warunku brzegowego
		ksi[6] = -1;
		ksi[7] = -1;
		eta[6] = (1.0 / sqrt(3));
		eta[7] = -(1.0 / sqrt(3));
	}

	else if (npc == 3) {
		waga[0] = waga[2] = 5.0 / 9.0;
		waga[1] = 8.0 / 9.0;

		ksi[0] = ksi[6] = eta[3] = eta[9] = -1.0 * sqrt(3.0 / 5.0);
		ksi[1] = ksi[7] = eta[4] = eta[10] = 0;
		ksi[2] = ksi[8] = eta[5] = eta[11] = 1.0 * sqrt(3.0 / 5.0);
		
		for (auto i = 0; i < 3; i++)
		{
			eta[i] = ksi[9 + i] = -1;
			eta[6 + i] = ksi[3 + i] = 1;
		}
	}

	else if (npc == 4) {
		waga[0] = waga[3] = (18.0 - sqrt(30.0)) / 36;
		waga[1] = waga[2] = (18.0 + sqrt(30.0)) / 36;

		ksi[0] = eta[4] = ksi[11] = eta[15] = -1.0 * sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt((6.0 / 5.0)));
		ksi[1] = eta[5] = ksi[10] = eta[14] = -1.0 * sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt((6.0 / 5.0)));;
		ksi[2] = eta[6] = ksi[9] = eta[13] = 1.0 * sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt((6.0 / 5.0)));
		ksi[3] = eta[7] = ksi[8] = eta[12] = 1.0 * sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt((6.0 / 5.0)));

		for (auto i = 0; i < 4; i++)
		{
			eta[i] = ksi[12 + i] = -1;
			eta[8 + i] = ksi[4 + i] = 1;
		}
	}
	
	//else cout << "Blad wczytywania ksi i eta dla boundary condition" << endl << endl;

	//dla jednej plaszczyzny sa 2 punkty, w jednym sa 4 n i w drugim sa 4 n

	double N_bc[4][4];		

	double H_lbc[4][4] = { 0 };		//macierz lokalna H bc

	double detJ_bc_x = deltaX / 2;
	double detJ_bc_y = deltaY / 2;

	double detJ_temp = (deltaX*deltaX) + (deltaY*deltaY);
	double detJ_bc = sqrt(detJ_temp) / 2;

	double matrix_P_l[1][4] = { 0.0 };

	//cout << "detJ_bc_x= " << detJ_bc_x << endl;
	//cout << "detJ_bc_y= " << detJ_bc_y << endl;

	//cout << "t alfa= " << t_alfa << endl;
	//cout << "alfa =" << alfa << endl;

	for (int i = 0; i < 4; i++) {
		H[iteracja].macierz_p[0][i] = { 0.0 };
	}

	//int k = 0;
	//dolny warunek
	if (tempbc1 == 1 && tempbc2 == 1) {
		for (int i = 0; i < npc; i++) {		//funkcje ksztaltu dla bc, te same co dla macierzy H, ale ille ksi i eta
			N_bc[i][0] = 0.25*(1 - ksi[i]) * (1 - eta[i]);
			N_bc[i][1] = 0.25*(1 + ksi[i]) * (1 - eta[i]);
			N_bc[i][2] = 0.25*(1 + ksi[i]) * (1 + eta[i]);
			N_bc[i][3] = 0.25*(1 - ksi[i]) * (1 + eta[i]);
		}


		//cout << "Dolny warunek dla npc = " << npc << "------------------------------------------------------------" << endl;

		for (int k = 0; k < npc; k++) {
			for (int i = 0; i < 4; i++) {//4, bo macierz h_lbc 4x4
				for (int j = 0; j < 4; j++) {//4, bo macierz h_lbc 4x4
					H_lbc[i][j] += (N_bc[k][i] * N_bc[k][j])*alfa*detJ_bc_x*waga[k];
				}
				matrix_P_l[0][i] += (-1)*N_bc[k][i] * (alfa) * t_alfa * detJ_bc_x * waga[k];
			}
		}
	}

	//prawy warunek
	if (tempbc2 == 1 && tempbc3 == 1) {
		for (int i = 0; i < npc; i++) {		//dla jednej plaszczyzny sa 2 punkty, w jednym sa 4 n i w drugim sa 4 n
			N_bc[i][0] = 0.25*(1 - ksi[i + npc]) * (1 - eta[i + npc]);
			N_bc[i][1] = 0.25*(1 + ksi[i + npc]) * (1 - eta[i + npc]);
			N_bc[i][2] = 0.25*(1 + ksi[i + npc]) * (1 + eta[i + npc]);
			N_bc[i][3] = 0.25*(1 - ksi[i + npc]) * (1 + eta[i + npc]);
		}


		//cout << "Prawy warunek dla npc = " << npc << endl << "------------------------------------------------------------" << endl;

		for (int k = 0; k < npc; k++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					H_lbc[i][j] += (N_bc[k][i] * N_bc[k][j])*alfa*detJ_bc_y*waga[k];
				}
				matrix_P_l[0][i] += (-1)*N_bc[k][i] * (alfa)* t_alfa * detJ_bc_x * waga[k];
			}
		}
	}

	//gorny warunek
	if (tempbc3 == 1 && tempbc4 == 1) {
		for (int i = 0; i < npc; i++) {		//dla jednej plaszczyzny sa 2 punkty, w jednym sa 4 n i w drugim sa 4 n
			N_bc[i][0] = 0.25*(1 - ksi[i + (npc * 2)]) * (1 - eta[i + (npc * 2)]);
			N_bc[i][1] = 0.25*(1 + ksi[i + (npc * 2)]) * (1 - eta[i + (npc * 2)]);
			N_bc[i][2] = 0.25*(1 + ksi[i + (npc * 2)]) * (1 + eta[i + (npc * 2)]);
			N_bc[i][3] = 0.25*(1 - ksi[i + (npc * 2)]) * (1 + eta[i + (npc * 2)]);
		}
		//int temp_bc = 0;


		//cout << "Gorny warunek dla npc = " << npc << endl << "------------------------------------------------------------" << endl;

		//int k = 4;
		for (int k = 0; k < npc; k++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					H_lbc[i][j] += (N_bc[k][i] * N_bc[k][j])*alfa*detJ_bc_x*waga[k];
				}
				matrix_P_l[0][i] += (-1)*N_bc[k][i] * (alfa)* t_alfa * detJ_bc_x * waga[k];
			}
		}
	}

	//lewy warunek
	if (tempbc4 == 1 && tempbc1 == 1) {
		for (int i = 0; i < npc; i++) {		//dla jednej plaszczyzny sa 2 punkty, w jednym sa 4 n i w drugim sa 4 n
			N_bc[i][0] = 0.25*(1 - ksi[i + (npc * 3)]) * (1 - eta[i + (npc * 3)]);
			N_bc[i][1] = 0.25*(1 + ksi[i + (npc * 3)]) * (1 - eta[i + (npc * 3)]);
			N_bc[i][2] = 0.25*(1 + ksi[i + (npc * 3)]) * (1 + eta[i + (npc * 3)]);
			N_bc[i][3] = 0.25*(1 - ksi[i + (npc * 3)]) * (1 + eta[i + (npc * 3)]);
		}
	
		//cout << "Lewy warunek dla npc = " << npc << endl << "------------------------------------------------------------" << endl;

		//int k = 6;
		for (int k = 0; k < npc; k++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					H_lbc[i][j] += (N_bc[k][i] * N_bc[k][j])*alfa*detJ_bc_y*waga[k];
				}
				matrix_P_l[0][i] += (-1)*N_bc[k][i] * (alfa)* t_alfa * detJ_bc_x * waga[k];
			}
		}
	}

	//wypisywanie macierzy H_lbc
	//cout << endl << endl << "Macierz H_lbc = " << endl;
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << H_lbc[i][j] << "  ";
	//	}
	//	cout << endl;
	//}

	////wypisywanie wektora P_l
	//cout << endl << endl << "Wektora P_l = " << endl;
	//for (int i = 0; i < 4; i++) {
	//	cout << matrix_P_l[0][i] << "  ";
	//	cout << endl;
	//}


	//zapisywanie macierzy h + h_lbc do struktury elementu skonczonego
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[iteracja].matrix_h[i][j] += H_lbc[i][j];
		}
		//cout << endl;
	}

	//zapisywanie macierzy p lokalne do struktury elementu skoñczonego
	for (int i = 0; i < 4; i++) {
		H[iteracja].macierz_p[0][i] += matrix_P_l[0][i];
	}


	//iteracja++;
	return 0;

}

double pochodna(double x[], double y[]) {
	double suma = 0;					//zmienna lokalna
	int size = sizeof(x);				//sprawdzam wielkosc tablicy 

	for (int i = 0; i < size; i++) {	//do wielkosci tablicy x[] podanej jako argument
		suma += x[i] * y[i];			//wzor z cwiczen.
	}

	return suma;

}
