#pragma once
#include <iostream>
#include <string>
#include <fstream>

using namespace std;


void covarAux(double *Ei, double* E, int* combs, int card, int nmed, int indc) 
{

	int m = 0, n = 0;
	//int nTmed = nmed + 1;
	for (int i = 0; i < nmed; i++) {
		if (combs[i + nmed * indc] == 1) {

			Ei[m * card + m] = E[i * nmed + i];
			
			n = m + 1;
			for (int j = i + 1; j < nmed; j++) {
				if (combs[j + nmed * indc] == 1) {
					Ei[m * card + n] = E[i * nmed + j];

					Ei[n * card + m] = E[j * nmed + i];

					n++;
				}
			}

			m++;
		}
	}
}

bool isinvertible(double mat[], int m) {
	//monta inv eye(m)
	bool inv;

	double pivo = 0.;
	for (int i = 0; i < m; i++) {

		//li = 0.;


		//pivotiamento 
		int indmaior = i;
		double maior = mat[i * m + i];
		for (int j = i; j < m; j++)
		{
			if (abs(maior) < abs(mat[j * m + i]))
			{
				maior = mat[j * m + i];
				indmaior = j;
			}
		}
		for (int j = 0; j < m; j++)
		{
			double swap = mat[i * m + j];
			mat[i * m + j] = mat[indmaior * m + j];
			mat[indmaior * m + j] = swap;
		}
		/*for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++)
				printf("%f ", mat[i * m + j]);
			printf("\n");
		}
		printf("\n");*/

		pivo = mat[i * m + i];
		if (abs(pivo) < 1e-10) {
			inv = 0;
			//printf("matriz nao possui inversa\n");

			return inv;
		}


		for (int j = 0; j < m; j++) {
			mat[i * m + j] = mat[i * m + j] / pivo;
			//Inv[i*m + j] = Inv[i*m + j] / pivo;

		}


		for (int j = i; j < m; j++) {
			if (j != i) {
				pivo = mat[j * m + i];
				for (int k = 0; k < m; k++) {
					mat[j * m + k] = mat[j * m + k] - pivo * mat[i * m + k];

				}
			}
		}



	}

	inv = 1;
	return inv;
}

void lerE(double* E, string Efilename , int nmed)
{
	
	ifstream Efile(Efilename.c_str());
	for (int i = 0; i < nmed; i++)
	{
		for (int j = 0; j < nmed; j++)
		{
			Efile >> E[i * nmed + j];

		}

	}
	Efile.close();

}
