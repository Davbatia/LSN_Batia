#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include "funzioni.h"
#include "random.h"

using namespace std;

// DEFINITION OF VARIABLES
const int M = 10000;     // Number of realizations
const double lambda = 1;    // Decay factor
const double Gamma = 1;     // Lorentzian factor (assuming mu = 0)

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
    vector<double> N(4);

    // Assign values to the elements
    N[0] = 1;
    N[1] = 2;
    N[2] = 10;
    N[3] = 100;

    //Uniform sampled
    vector<double> Suni1(M); 
    vector<double> Suni2(M); 
    vector<double> Suni10(M); 
    vector<double> Suni100(M); 

    //Exponential sampled
    vector<double> Sexp1(M); 
    vector<double> Sexp2(M); 
    vector<double> Sexp10(M); 
    vector<double> Sexp100(M); 

    //Lorentzian sampled
    vector<double> Slor1(M); 
    vector<double> Slor2(M); 
    vector<double> Slor10(M); 
    vector<double> Slor100(M); 

    //Ausiliar variables to compute S_N
    double sum_uni_1=0;
    double sum_exp_1=0;
    double sum_lor_1=0;
    double sum_uni_2=0;
    double sum_exp_2=0;
    double sum_lor_2=0;
    double sum_uni_10=0;
    double sum_exp_10=0;
    double sum_lor_10=0;
    double sum_uni_100=0;
    double sum_exp_100=0;
    double sum_lor_100=0;

    //S_N CALCULATION
    for(int i=0; i<M; i++) {
        sum_uni_1=0;
        sum_exp_1=0;
        sum_lor_1=0;
        sum_uni_2=0;
        sum_exp_2=0;
        sum_lor_2=0;
        sum_uni_10=0;
        sum_exp_10=0;
        sum_lor_10=0;
        sum_uni_100=0;
        sum_exp_100=0;
        sum_lor_100=0;
        for(int k=0; k<N[0]; k++) {
                sum_uni_1 += rnd.Rannyu()/N[0];
                sum_exp_1 += rnd.Exponential(lambda)/N[0];
                sum_lor_1 += rnd.Lorentz(Gamma)/N[0];
        }
        for(int k=0; k<N[1]; k++) {
                sum_uni_2 += rnd.Rannyu()/N[1];
                sum_exp_2 += rnd.Exponential(lambda)/N[1];
                sum_lor_2 += rnd.Lorentz(Gamma)/N[1];
        }
        for(int k=0; k<N[2]; k++) {
                sum_uni_10 += rnd.Rannyu()/N[2];
                sum_exp_10 += rnd.Exponential(lambda)/N[2];
                sum_lor_10 += rnd.Lorentz(Gamma)/N[2];
        }
        for(int k=0; k<N[3]; k++) {
                sum_uni_100 += rnd.Rannyu()/N[3];
                sum_exp_100 += rnd.Exponential(lambda)/N[3];
                sum_lor_100 += rnd.Lorentz(Gamma)/N[3];
        }
        //Uniform
        Suni1[i]=sum_uni_1; 
        Suni2[i]=sum_uni_2;
        Suni10[i]=sum_uni_10; 
        Suni100[i]=sum_uni_100;
        //Exponential
        Sexp1[i]=sum_exp_1; 
        Sexp2[i]=sum_exp_2;
        Sexp10[i]=sum_exp_10; 
        Sexp100[i]=sum_exp_100;
        //Lorentzian
        Slor1[i]=sum_lor_1; 
        Slor2[i]=sum_lor_2;
        Slor10[i]=sum_lor_10; 
        Slor100[i]=sum_lor_100;
    }


    // CREATION OF OUTPUT FILE
    string nomeFileOutput = "output_dati.txt";
    scriviSuFile(nomeFileOutput, Suni1, Suni2, Suni10, Suni100, Sexp1, Sexp2, Sexp10, Sexp100, Slor1, Slor2, Slor10, Slor100);

    return 0;
}