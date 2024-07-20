#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include "funzioni.h"
#include "random.h"


using namespace std;

const int M = 100;  // number of sub-interval   
const int n = 10000;     // Number of throws

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
    vector<double> r(n);
    vector<double> chi2(10000, 0.0);
    int n_i=0;                   //Numbers observed in each sub-interval

    //CHI-SQUARED CALCULATION
    for(int k=0;k<10000;k++){   
            // U[0,1) uniform distribution
            for (int l = 0; l < n; ++l) {
                r[l] = rnd.Rannyu();  
            }

            for(int i=0; i<M; i++) {
                for(int j=0 ; j<n; j++){
                    if(0.01*i <= r[j] && r[j] < 0.01*i+0.01){    // 0.01 because the number of sub-intervalls is M=100
                        n_i +=1 ;
                    }
                }   

                chi2[k] += pow((n_i-n/M),2)/(n/M) ;
                n_i=0 ;                                         // Resetting the counter
            }
    }
    // CREATION OF OUTPUT FILE
    string nomeFileOutput = "output_dati.txt";
    scriviSuFile(nomeFileOutput, chi2);






    return 0;
}