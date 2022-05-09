#pragma once
#include<iostream>
#include<string>
#include <fstream>

using namespace std;
//---Combinatórios


unsigned long long choose(unsigned long long n, unsigned long long k);




unsigned long long choose(unsigned long long n, unsigned long long k) {
	if (k > n) {
		return 0;
	}
	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d) {
		r *= n--;
		r /= d;
	}
	return r;
}

void lerCn(unsigned int* Cn, string Cnfilename,int nmed, int kmax)
{
	ifstream Efile(Cnfilename.c_str());
	for (int i = 0; i < nmed+1; i++)
	{
		for (int j = 0; j < kmax+1; j++)
		{
			Efile >> Cn[i * (kmax +1) + j];

		}

	}
	Efile.close();
}

