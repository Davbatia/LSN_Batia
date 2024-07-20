#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include "funzioni.h"
#include "random.h"

using namespace std;

// DEFINITION OF VARIABLES
const int M = 100000;  // Total number of throws
const int N = 100;     // Number of blocks
const int n = M / N;   // Number of throws in each block, ensure M is a multiple of N
const double d = 3.0;   // Distance between two straight lines
const double L = 2.0;   // needle's length

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
    vector<double> Pi(M, 0.0);          //Pi calculated
    vector<double> Random(M);            //Random numbers UD [0,d); yi is the inizial values of height + L*sin(Theta), Theta uniformly distributed beetween (0,PI) (without knowing Pi)
    vector<double> ave(N, 0.0);         //Average
    vector<double> av2(N, 0.0);         //Square average
    vector<double> sum_prog(N, 0.0);    //Cumulative average
    vector<double> su2_prog(N, 0.0);    //Cumulative square average
    vector<double> err_prog(N, 0.0);    // Statistical uncertainty


    //CALCULATION OF PI
    int N_hit=0;
    double hit=0;
    for(int i=0; i<M; i++){
        N_hit=0;
        for(int j=0; j<M; j++){
            hit= rnd.Rannyu(0,d) + L*sin(rnd.UnifTheta());
            if(hit >= d){
                N_hit +=1;
                //cout << "ok" << endl;
            }else{}
        }
        //cout << "ok" << endl;
        Pi[i]=2*L/d * M/N_hit;
    }
    cout << "ok" << endl;

    // DATA BLOCKING
    double sum1 = 0;
    for (int i = 0; i < N; ++i) {
        sum1 = 0.0;
        for (int j = 0; j < n; ++j) {
            int k = j + i * n;
            sum1 += Pi[k];
        }
        ave[i] = sum1 / n;          // Pi_i
        av2[i] = pow(ave[i], 2);    // (Pi_i)^2 
    }
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j <= k; ++j) {
                sum_prog[k] += ave[j]; // SUM_{j=0,i} Pi_j
                su2_prog[k] += av2[j]; // SUM_{j=0,i} (Pi_j)^2
        }
        sum_prog[k] /= (k + 1); // Cumulative average
        su2_prog[k] /= (k + 1); // Cumulative square average
        err_prog[k] = error(sum_prog, su2_prog, k); // Statistical uncertainty
    }
    //cout << "ok" << endl;

    // CREATION OF OUTPUT FILE
    string nomeFileOutput = "output_dati.txt";
    scriviSuFile(nomeFileOutput, sum_prog, err_prog);

    return 0;
}