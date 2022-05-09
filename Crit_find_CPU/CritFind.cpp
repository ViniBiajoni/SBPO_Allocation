#include <stdio.h>
#include <iostream>

#include<omp.h>
#include<math.h>


#include "AnalisCrit.h"
#include "Combinatoric.h"

#pragma warning(disable:4996)

#include<string>
#include <fstream>
#include<time.h>
using namespace std;


#define min(a, b) (((a) < (b)) ? (a) : (b))

const int nbar =30;
//const int nmed = 33;
const int kmax =3;
const int solSize = 1000;


//Estrututa do Conjunto Solucao
struct SolStack
{
    int nSols;

    int *Sols;

};




//// Etapas da Analise de criticalidade
void enumerar(int* combs, int card, int combsIdini, int combsInWave, unsigned int* Cn,int  nmed );
//
void procurar(int* iscrit, int combscard, int card, int* combs, double* E, int  nmed);
//
void confirmar(struct SolStack conjSol, int* combs, int combscard, int* iscrit, int  nmed);
//
void atualizaConjSol(SolStack* conjSol, int* iscrit, int* combs, int combscard, int  nmed);


int main()
{
    int nmed;
    ifstream nmedFile("nmed.txt");
    
    nmedFile >> nmed;

    nmedFile.close();

    // Variaveis iniciais
    const unsigned int wave_size = (int)pow(2, 20);

    // Combinacoes pre calculadas (CpC)
    unsigned int Cn[(400 + 1) * (4 + 1)] = { 0 };
    string nomeCnfile = "combs400em4.txt";
    lerCn(Cn, nomeCnfile, nmed, 4);

    // Matriz Covariancia E
    string nomeEfile = "Covariancia.txt";
    double* E;
    E = (double*)malloc(nmed * (size_t)nmed * sizeof(double));
    lerE(E, nomeEfile, nmed);

    // Conjunto Solucao
    SolStack conjSol;
    conjSol.nSols = -1;
    conjSol.Sols = (int*)malloc(nmed * (size_t)solSize * sizeof(int));

    int* combs; //Combinacoes enumeradas por ondas
    combs = (int*)malloc((size_t)wave_size * nmed * sizeof(int));

    int* iscrit;// Vetor booleano  1: Combinacao critica 0: Combinacao nao 
    iscrit = (int*)malloc(wave_size * sizeof(int));
   
    for (int card = 1; card <= kmax; card++)
    {

        //cout << "cardinalidae " << card;
        //cout << ": iniciado...";
        unsigned int totalwaves = 0; // Combinacoes vizitadas em todas as ondas
        //int combsId = 0;    // Identificador de combinacoes :  (1 = 0011) (2 = 0101) ... (6 = 1100) 

        while (totalwaves < Cn[nmed * (4 + 1) + card])
        {

            int combsInWave = min(wave_size, Cn[nmed * (4 + 1) + card] - totalwaves); // Combinacoes Visitada em onda

            //1-Enumeracao--------------------------- 
            
            enumerar(combs, card, totalwaves + 1, combsInWave, Cn,nmed);
            
           
            //2-Procura-------------------------------
          
            procurar(iscrit, combsInWave, card, combs, E,nmed);
            

            //3-Confirmacao---------------------------
           
            confirmar(conjSol, combs, combsInWave, iscrit,nmed);
           

            //4-Atualizacao do Conjunto Solucao-------
            atualizaConjSol(&conjSol, iscrit, combs, combsInWave,nmed);


            totalwaves += wave_size;

        }
        //cout << "finalizado." << endl;

    }
    
    free(combs);
    free(iscrit);

    

    

    //Print Resultados ----------------------------------------------------------------------------------------------
    string file_name;
    //asctime(localtime(&timetoday))+
    file_name = "Criticalidades.txt";
    FILE* Output_file;
    Output_file = fopen(file_name.c_str(), "w");

    //printf("RESULTADOS:\n");
    //printf("Conjunto solucao:\n");
    int ncrits[kmax];
    for (int i = 0; i < kmax; i++)
        ncrits[i] = 0;
    int card;
    for (int i = 0; i < conjSol.nSols + 1; i++) {
        card = 0;
        //fprintf(Output_file,"[");
        for (int j = 0; j < nmed; j++) {
            if (conjSol.Sols[i * nmed + j] == 1) {
                card++;
                //fprintf(Output_file," %i-%i ", medida[j*3], medida[j * 3 +1]);
            }
        }//fprintf(Output_file,"]\n");
        ncrits[card - 1] += 1;
    }
    //printf("Cardinalidade -> numero de tuplas criticas\n");
    for (int i = 0; i < kmax; i++){
        //printf("%i -> %i\n", i + 1, ncrits[i]);
        fprintf(Output_file, "%i\n", ncrits[i]);
    }
    //printf(" total de tuplas criticas : %i\n", conjSol.nSols + 1);

    //Print Tempos e CSV----------------------------------------------------------------------------------------------



    fclose(Output_file);

}

void enumerar(int* combs, int card, int combsIdinicial, int combsInWave, unsigned int* Cn,int nmed)
{
    #pragma omp parallel
    {
    #pragma omp for
        for (int linha = 0; linha < combsInWave; linha++)
        {
            //int test = Cn[nmed * (kmax + 1) + card];
            int n = linha + combsIdinicial;

            int freq[2] = { 0 };

            freq[0] = nmed - card;
            freq[1] = card;

            // iterate till sum equals n 
            int sum = 0;
            int k = 0;

            while (sum != n) {

                sum = 0;

                for (int i = 0; i < 2; i++) {
                    if (freq[i] == 0)
                        continue;

                    // Remove character 
                    freq[i]--;


                    unsigned int xsum;
                    int num = nmed - 1 - k;
                    if (num == 0) xsum = 1;
                    else {
                        xsum = Cn[num * (4 + 1) + min(freq[0], freq[1])];
                    }
                    sum += xsum;

                    if (sum >= n) {

                        combs[linha * nmed + k] = i;
                        k++;
                        n -= (sum - xsum);
                        break;
                    }
                    if (sum < n)
                        freq[i]++;
                }

            }

            for (int i = 2 - 1; k < nmed && i >= 0; i--)
                if (freq[i]) {

                    combs[linha * nmed + k] = i;
                    k++;

                    //cout << out[k - 1] << endl;
                    freq[i++]--;
                }

            // append string termination 
            // character and print result 


        }
    }
}

void procurar(int* iscrit, int combscard, int card, int* combs, double* E, int nmed)
{
    #pragma omp parallel
    {
        #pragma omp for
        for (int combi = 0; combi < combscard; combi++)
        {
            
            //monta matriz de combinacao auxiliar
            double Ei[kmax * kmax];
            covarAux(Ei, E, combs, card, nmed, combi);

            

            if (!(isinvertible(Ei, card)))
            {
                
                iscrit[combi] = 1;
            }
            else {
                iscrit[combi] = 0;
            }
        }
    }
}

void confirmar(struct SolStack conjSol, int* combs, int combscard, int* iscrit,int nmed)
{


    for (int sol = 0; sol < conjSol.nSols + 1; sol++)
    {
    #pragma omp parallel
        {
            #pragma omp for
            for (int crit = 0; crit < combscard; crit++)
            {
                if (iscrit[crit] == 1)
                {
                    int is = 0;
                    for (int i = 0; i < nmed; i++)
                    {
                        int result = combs[crit * nmed + i] - conjSol.Sols[sol * nmed + i];

                        if (result == -1) is = 1;
                    }
                    iscrit[crit] = is;

                }

            }
        }
    }

}

void atualizaConjSol(SolStack* conjSol, int* iscrit, int* combs, int combscard, int nmed)
{


    for (int i = 0; i < combscard; i++)
    {
        if (iscrit[i] == 1)
        {
            conjSol->nSols += 1;
            for (int j = 0; j < nmed; j++)
                conjSol->Sols[conjSol->nSols * nmed + j] = combs[i * nmed + j];
        }
        //addSol(conjSol, &combs[i * nmed]);
    }
}