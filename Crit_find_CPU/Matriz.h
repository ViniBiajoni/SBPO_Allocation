#pragma once
#include<iostream>


void printmat(double *mat, int m, int n) {
	for (int i = 0; i < m; i++) {
		printf("[");
		for (int j = 0; j < n; j++) {
			printf("%3.4lf ", mat[i * n + j]);
		}
		printf("]\n");
	}
}
void fprintmat(double *mat, int m, int n, FILE* fp) {
	for (int i = 0; i < m; i++) {
		printf("[");
		for (int j = 0; j < n; j++) {
			fprintf(fp, "%3.4lf & ", mat[i * n + j]);
		}
		fprintf(fp, "\\\\ \n");
	}
}
bool isinvertible(double mat[], int m) {
	//monta inv eye(m)
	bool inv;
	
	double pivo = 0.;
	for (int i = 0; i < m; i++) {
		pivo = mat[i * m + i];
		//li = 0.;
		for (int k = i; k < m; k++) {
			//li = li + mat[i*m + k];
		}
		//if (m == 5) 
			//printf("pivo %f= \n", pivo);
		if (abs(pivo) < 0.0001) {//|| abs(li) < 0.0001) {
			inv = 0;
			//printf("matriz nao possui inversa\n");
			//free(Inv);
			return inv;
		}

		for (int j = 0; j < m; j++) {
			mat[i * m + j] = mat[i * m + j] / pivo;
			//Inv[i*m + j] = Inv[i*m + j] / pivo;

		}


		for (int j = 0; j < m; j++) {
			if (j != i) {
				pivo = mat[j * m + i];
				for (int k = 0; k < m; k++) {
					mat[j * m + k] = mat[j * m + k] - pivo * mat[i * m + k];
					//Inv[j*m + k] = Inv[j*m + k] - pivo * Inv[i*m + k];

				}
			}
		}



	}

	inv = 1;
	//free(Inv);
	return inv;
}