#pragma once
#include "pch.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <math.h>

using namespace std;

struct GlobalData {
	float H, W;
	int	mW, mH;
	int mN, mE;
	int NPC;
	int ro, c;
	int	alfa;
};

struct node {
	double x, y;
	float BC;	//boundary condition - 1 jezeli jest to wezel zewnetrzny (kontakt z otoczeniem), lub 0 jezeli jest w srodku (nie maja kontaktu z otoczeniem)
	//void wyswietlnode();
public:
	void wyswietlnode() {

		cout << "x = " << x << "  ";
		cout << "y = " << y << endl << endl;
	}
};

struct element {	//wspolrzedne elementu skonczonego
	int id[4];
	double matrix_h[4][4];

	double macierz_c[4][4];

	double macierz_p[1][4];
	//void wyswietl(); 
public:
	void wyswietl() {

		cout << "q = " << id[0] << "  ";
		cout << "w = " << id[1] << "  ";
		cout << "e = " << id[2] << "  ";
		cout << "r = " << id[3] << endl << endl;
	}
};

float funkcja(float x, float y);
//float funkcja2d();
float oblicz();

//float oblicz(float x);

/*struct elem4
{
	float q, w, e, r;

};
*/

double jacobi(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, element H[], int npc, int cp, int ro, int tempbc1, int tempbc2, int tempbc3, int tempbc4, int alfa, float deltaX, float deltaY, double t0, int t_alfa, int iteracja);

double* gauss(int mN, double** HG_temp, double* PG_temp);

/*struct jacobian
{
	float ksi[4];
	float eta[4];

	jacobian();
};
*/

/*jacobian::jacobian() {
	ksi[1] = -0.77;
	ksi[2] = 0.77;
	ksi[3] = 0.77;
	ksi[4] = -0.77;

	eta[1] = -0.77;
	eta[2] = -0.77;
	eta[3] = 0.77;
	eta[4] = 0.77;
}
*/

double pochodna(double x[], double y[]);

struct elem4
{
	int npc = 0;
	double eta[16];
	double ksi[16];

	double bc_ksi[16];
	double bc_eta[16];

	//double *ksi_pc;
	//double *eta_pc;

	double w_eta[16];
	double w_ksi[16];

	//double trzy_pkt_ksi_w[9] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	//double trzy_pkt_eta_w2[9] = { 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 };

public:
	elem4(int npc) {
		this->npc = npc;

		if (npc == 2) {
			//eta = new double[4];
			//ksi = new double[4];

			//waga
			//w_eta = new double[4];
			//w_ksi = new double[4];

			//eta
			eta[0] = -(1 / sqrt(3));
			eta[1] = -(1 / sqrt(3));
			eta[2] = (1 / sqrt(3));
			eta[3] = (1 / sqrt(3));

			//ksi
			ksi[0] = -(1 / sqrt(3));
			ksi[1] = (1 / sqrt(3));
			ksi[2] = (1 / sqrt(3));
			ksi[3] = -(1 / sqrt(3));

			//w_eta
			w_eta[0] = 1;
			w_eta[1] = 1;
			w_eta[2] = 1;
			w_eta[3] = 1;

			//w_ksi
			w_ksi[0] = 1;
			w_ksi[1] = 1;
			w_ksi[2] = 1;
			w_ksi[3] = 1;


			//bc - boundary coundition
			
		//	bc_ksi[0] = bc_eta[2] = bc_ksi[5] = bc_eta[7] = -1.0 / sqrt(3);
		//	bc_ksi[1] = bc_eta[3] = bc_ksi[4] = bc_eta[6] = 1.0 / sqrt(3);
		//	bc_ksi[2] = bc_ksi[3] = bc_eta[4] = bc_eta[5] = 1;
		//	bc_eta[0] = bc_eta[1] = bc_ksi[6] = bc_ksi[7] = -1;
		//	for (auto i = 0; i < 2; i++) multiplier[i] = 1;
		//	break;
		}

		else if (npc == 3) {
			//eta = new double[9];
			//ksi = new double[9];

			//w_eta = new double[9];
			//w_ksi = new double[9];

			//eta x = {0,0,0,1,1,1,2,2,2}
			eta[0] = -(sqrt(3.0 / 5.0));
			eta[1] = -(sqrt(3.0 / 5.0));
			eta[2] = -(sqrt(3.0 / 5.0));

			eta[3] = 0;
			eta[4] = 0;
			eta[5] = 0;

			eta[6] = (sqrt(3.0 / 5.0));
			eta[7] = (sqrt(3.0 / 5.0));
			eta[8] = (sqrt(3.0 / 5.0));

			//ksi x = {0,1,2,0,1,2,0,1,2}
			ksi[0] = -(sqrt(3.0 / 5.0));
			ksi[1] = 0;
			ksi[2] = (sqrt(3.0 / 5.0));

			ksi[3] = -(sqrt(3.0 / 5.0));
			ksi[4] = 0;
			ksi[5] = (sqrt(3.0 / 5.0));

			ksi[6] = -(sqrt(3.0 / 5.0));
			ksi[7] = 0;
			ksi[8] = (sqrt(3.0 / 5.0));


			//w_eta
			w_eta[0] = 5.0 / 9.0;
			w_eta[1] = 5.0 / 9.0;
			w_eta[2] = 5.0 / 9.0;

			w_eta[3] = 8.0 / 9.0;
			w_eta[4] = 8.0 / 9.0;
			w_eta[5] = 8.0 / 9.0;

			w_eta[6] = 5.0 / 9.0;
			w_eta[7] = 5.0 / 9.0;
			w_eta[8] = 5.0 / 9.0;

			//w_ksi
			w_ksi[0] = 5.0 / 9.0;
			w_ksi[1] = 8.0 / 9.0;
			w_ksi[2] = 5.0 / 9.0;

			w_ksi[3] = 5.0 / 9.0;
			w_ksi[4] = 8.0 / 9.0;
			w_ksi[5] = 5.0 / 9.0;

			w_ksi[6] = 5.0 / 9.0;
			w_ksi[7] = 8.0 / 9.0;
			w_ksi[8] = 5.0 / 9.0;
		}

		else if (npc == 4) {
			//wartosci dla eta i ksi, wziete z tabeli
			//eta = new double[16];
			//ksi = new double[16];

			//waga
			//w_eta = new double[16];
			//w_ksi = new double[16];

			//eta
			eta[0] = -0.861136;
			eta[1] = -0.861136;
			eta[2] = -0.861136;
			eta[3] = -0.861136;

			eta[4] = -0.339981;
			eta[5] = -0.339981;
			eta[6] = -0.339981;
			eta[7] = -0.339981;

			eta[8] = 0.339981;
			eta[9] = 0.339981;
			eta[10] = 0.339981;
			eta[11] = 0.339981;

			eta[12] = 0.861136;
			eta[13] = 0.861136;
			eta[14] = 0.861136;
			eta[15] = 0.861136;

			//ksi
			ksi[0] = -0.861136;
			ksi[1] = -0.339981;
			ksi[2] = 0.339981;
			ksi[3] = 0.861136;

			ksi[4] = -0.861136;
			ksi[5] = -0.339981;
			ksi[6] = 0.339981;
			ksi[7] = 0.861136;

			ksi[8] = -0.861136;
			ksi[9] = -0.339981;
			ksi[10] = 0.339981;
			ksi[11] = 0.861136;

			ksi[12] = -0.861136;
			ksi[13] = -0.339981;
			ksi[14] = 0.339981;
			ksi[15] = 0.861136;

			//w_eta
			w_eta[0] = 0.347855;
			w_eta[1] = 0.347855;
			w_eta[2] = 0.347855;
			w_eta[3] = 0.347855;

			w_eta[4] = 0.652145;
			w_eta[5] = 0.652145;
			w_eta[6] = 0.652145;
			w_eta[7] = 0.652145;

			w_eta[8] = 0.652145;
			w_eta[9] = 0.652145;
			w_eta[10] = 0.652145;
			w_eta[11] = 0.652145;

			w_eta[12] = 0.347855;
			w_eta[13] = 0.347855;
			w_eta[14] = 0.347855;
			w_eta[15] = 0.347855;

			//w_ksi
			w_ksi[0] = 0.347855;
			w_ksi[1] = 0.652145;
			w_ksi[2] = 0.652145;
			w_ksi[3] = 0.347855;

			w_ksi[4] = 0.347855;
			w_ksi[5] = 0.652145;
			w_ksi[6] = 0.652145;
			w_ksi[7] = 0.347855;

			w_ksi[8] = 0.347855;
			w_ksi[9] = 0.652145;
			w_ksi[10] = 0.652145;
			w_ksi[11] = 0.347855;

			w_ksi[12] = 0.347855;
			w_ksi[13] = 0.652145;
			w_ksi[14] = 0.652145;
			w_ksi[15] = 0.347855;
		}
	};
};