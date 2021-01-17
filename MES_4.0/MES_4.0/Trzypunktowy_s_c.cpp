#pragma once
#include "pch.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

float funkcja(float x, float y) {
	return -5 * (x*x)*y + 2 * x*y*y + 10;
}

float oblicz() {	//trzypunktowy schemat calkowania
	float pc[3] = { -0.77, 0, 0.77 };

	//cout << "pc1 = " << pc[0] << endl;
	//cout << "pc2 = " << pc[1] << endl;
	//cout << "pc3 = " << pc[2] << endl;

	//float pc1 = -0.77;
	//float pc2 = 0;
	//float pc3 = 0.77;

	float w[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

	//cout << "w1 = " << w[0] << endl;
	//cout << "w2 = " << w[1] << endl;
	//cout << "w3 = " << w[2] << endl;

	//float w1 = 5.0 / 9.0;
	//float w2 = 8.0 / 9.0;
	//float w3 = w1;

	float wynik = 0;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)

			wynik += funkcja(pc[i], pc[j]) * w[i] * w[j];
		//cout << "wynik = " << wynik << endl;
	}

	cout << "wynik = " << wynik << endl;
	return wynik;
}