#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include "funzioni.h"
#include "random.h"


using namespace std;

const int M = 100000;  // Total number of throws
const int N = 100;     // Number of blocks
const int L = M / N;    // Number of throws in each block, ensure M is a multiple of N

 
int main (int argc, char *argv[]){

    // GENERATION OF PSEUDO-RANDOM NUMBERS:
    Random rnd;
    int seed[4];
    int p1, p2;

    ifstream Primes("Primes");
        if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();


    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);   // Fixing random seed for reproducibility
            }
        } 
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Generate and print random numbers
    for(int i=0; i<20; i++){
        cout << rnd.Rannyu() << endl;
    }
    rnd.SaveSeed();  

    // DEFINITION OF VARIABLES
    vector<double> r(M);
    vector<int> x(N);
    vector<double> ave(N, 0.0);
    vector<double> av2(N, 0.0);
    vector<double> ave_var(N, 0.0);
    vector<double> av2_var(N, 0.0);
    vector<double> sum_prog(N, 0.0);
    vector<double> su2_prog(N, 0.0);
    vector<double> err_prog(N, 0.0);
    vector<double> sum_prog_var(N, 0.0);
    vector<double> su2_prog_var(N, 0.0);
    vector<double> err_prog_var(N, 0.0);

    // U[0,1) uniform distribution
    for (int i = 0; i < M; ++i) {
        r[i] = rnd.Rannyu();  
    }
    // [0,1,2,...,N-1]
    for (int i = 0; i < N; ++i) {
        x[i] = i;  
    }

    // AVERAGE AND VARIANCE CALCULATION
    for (int i = 0; i < N; ++i) {
        double sum1 = 0.0;
        // double sum2 = 0.0; // PROVA
        double sum3 =0.0;
        for (int j = 0; j < L; ++j) {
            int k = j + i * L;

            //Average
            sum1 += r[k];
            //sum2 += r[k] * r[k]; // PROVA

            // Variance
            sum3 += pow(r[k]-0.5 ,2 ); 
        }
        // Average
        ave[i] = sum1 / L; // r_i 
        // av2[i] = sum2 / L; // (r_i)^2 // PROVA
        av2[i] = pow(ave[i], 2); // (r_i)^2 

        // Variance
        ave_var[i] = sum3 / L;              
        av2_var[i] = pow(ave_var[i], 2);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {

            //Average
            sum_prog[i] += ave[j]; // SUM_{j=0,i} r_j
            su2_prog[i] += av2[j]; // SUM_{j=0,i} (r_j)^2

            //Variance
            sum_prog_var[i] += ave_var[j]; 
            su2_prog_var[i] += av2_var[j];
        }
        //Average
        sum_prog[i] /= (i + 1); // Cumulative average
        su2_prog[i] /= (i + 1); // Cumulative square average
        err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty

        //Variance
        sum_prog_var[i] /= (i + 1); 
        su2_prog_var[i] /= (i + 1); 
        err_prog_var[i] = error(sum_prog_var, su2_prog_var, i);
    }

    // x *= L; // Number of throws = block * (Number of throws in each block)

    // CREATION OF OUTPUT FILE
    string nomeFileOutput = "output_dati.txt";
    scriviSuFile(nomeFileOutput, sum_prog, err_prog, sum_prog_var, err_prog_var);


return 0;
}