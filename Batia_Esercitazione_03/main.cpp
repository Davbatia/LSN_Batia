#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include "funzioni.h"
#include "random.h"


double max(double numero1, double numero2) {
    if (numero1 > numero2) {
        return numero1;
    } else {
        return numero2;
    }
}


using namespace std;

const int M = 100000;  // Total number of throws
const int N = 100;     // Number of blocks
const int L = M / N;    // Number of throws in each block, ensure M is a multiple of N

//const double t = 0. ;
const double T = 1 ;
const double S0 = 100 ;
const double K = 100  ;
const double r = 0.1 ;
const double sigma = 0.25 ;

 
int main (int argc, char *argv[]){

    //Generatore numeri casuali
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

    for(int i=0; i<20; i++){
        cout << rnd.Rannyu() << endl;
    }
    rnd.SaveSeed();  

   
    //Asset price (Directly)
    vector<double> S(M);
    for (int i=0 ; i<M; i++){
        S[i]=S0 * exp((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0,T)); // Normal distribution
    }

    //Asset price (Discretized)
    vector<double> s(M);
    for(int j = 0; j<M ; j++){
        for(int i=0; i<100 ; i++) {
            if(i==0){
                s[j]=S0;
            }
            // s[j]=s[j-1]* exp((r-0.5*sigma*sigma)*1+sigma*rnd.Gauss(0,1)*sqrt(1));  //Ho tenuto indicate le i ma sono superflue
            s[j]=s[j]* exp((r-0.5*sigma*sigma)*0.01+sigma*rnd.Gauss(0,1)*sqrt(0.01));              
        }
    }

    //Call price in (S0,0) (Directly)
    vector<double> C(M);
    for (int i=0 ; i<M; i++){
        C[i]= exp(-r*T)*max(0 , S[i]-K) ;
    }

    //Put price in (S0,0) (Directly)
    vector<double> P(M);
    for (int i=0 ; i<M; i++){
        P[i]= exp(-r*T)*max(0 , K-S[i]) ;
    }

    //Call price in (S0,0) (Discretized)
    vector<double> c(M);
    for (int i=0 ; i<M; i++){
        c[i]= exp(-r*T)*max(0 , s[i]-K) ;
    }

    //Put price in (S0,0) (Discretized)
    vector<double> p(M);
    for (int i=0 ; i<M; i++){
        p[i]= exp(-r*T)*max(0 , K-s[i]) ;
    }

    //Call price (Directly)
    vector<double> aveC(N, 0.0);
    vector<double> av2C(N, 0.0);
    vector<double> sum_progC(N, 0.0);
    vector<double> su2_progC(N, 0.0);
    vector<double> err_progC(N, 0.0);

    //Put price (Directly)
    vector<double> aveP(N, 0.0);
    vector<double> av2P(N, 0.0);
    vector<double> sum_progP(N, 0.0);
    vector<double> su2_progP(N, 0.0);
    vector<double> err_progP(N, 0.0);

    //Call price (Discretized)
    vector<double> avec(N, 0.0);
    vector<double> av2c(N, 0.0);
    vector<double> sum_progc(N, 0.0);
    vector<double> su2_progc(N, 0.0);
    vector<double> err_progc(N, 0.0);

    //Put price (Discretized)
    vector<double> avep(N, 0.0);
    vector<double> av2p(N, 0.0);
    vector<double> sum_progp(N, 0.0);
    vector<double> su2_progp(N, 0.0);
    vector<double> err_progp(N, 0.0);

     for (int i = 0; i < N; ++i) {
        double sumC = 0.0;
        double sumP =0.0;
        double sumc = 0.0;
        double sump =0.0;
        for (int j = 0; j < L; ++j) {
            int k = j + i * L;
            sumC += C[k];
            sumP += P[k];
            sumc += c[k];
            sump += p[k];
          
             }
        aveC[i] = sumC / L; // C_i 
        av2C[i] = pow(aveC[i], 2); // (C_i)^2 

        aveP[i] = sumP / L; // P_i 
        av2P[i] = pow(aveP[i], 2); // (P_i)^2 

        avec[i] = sumc / L; // C_i 
        av2c[i] = pow(avec[i], 2); // (C_i)^2 

        avep[i] = sump / L; // P_i 
        av2p[i] = pow(avep[i], 2); // (P_i)^2 
       
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            sum_progC[i] += aveC[j]; // SUM_{j=0,i} C_j
            su2_progC[i] += av2C[j]; // SUM_{j=0,i} (C_j)^2

            sum_progP[i] += aveP[j]; // SUM_{j=0,i} P_j
            su2_progP[i] += av2P[j]; // SUM_{j=0,i} (P_j)^2

            sum_progc[i] += avec[j]; // SUM_{j=0,i} C_j
            su2_progc[i] += av2c[j]; // SUM_{j=0,i} (C_j)^2

            sum_progp[i] += avep[j]; // SUM_{j=0,i} P_j
            su2_progp[i] += av2p[j]; // SUM_{j=0,i} (P_j)^2
        }

        sum_progC[i] /= (i + 1); // Cumulative average
        su2_progC[i] /= (i + 1); // Cumulative square average
        err_progC[i] = error(sum_progC, su2_progC, i); // Statistical uncertainty

        sum_progP[i] /= (i + 1); // Cumulative average
        su2_progP[i] /= (i + 1); // Cumulative square average
        err_progP[i] = error(sum_progP, su2_progP, i); // Statistical uncertainty

        sum_progc[i] /= (i + 1); // Cumulative average
        su2_progc[i] /= (i + 1); // Cumulative square average
        err_progc[i] = error(sum_progc, su2_progc, i); // Statistical uncertainty
        
        sum_progp[i] /= (i + 1); // Cumulative average
        su2_progp[i] /= (i + 1); // Cumulative square average
        err_progp[i] = error(sum_progp, su2_progp, i); // Statistical uncertainty
    }



    string nomeFileOutput = "output_dati.txt";
    scriviSuFile(nomeFileOutput, sum_progC, err_progC, sum_progP, err_progP, sum_progc, err_progc, sum_progp, err_progp);


return 0;
}