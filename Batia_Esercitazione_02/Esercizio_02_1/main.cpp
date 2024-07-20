#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include <utility>  // per std::pair
#include "funzioni.h"
#include "random.h"

using namespace std;

// DEFINITION OF VARIABLES
const int M = 10000;  // Total number of throws
const int N = 100;     // Number of blocks
// const int n = M / N;   // Number of throws in each block, ensure M is a multiple of N


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

    //MC INTEGRAL WITH UNIFORM DISTRIBUTION

    //definition of variables
    vector<double> x(M);
    vector<double> f(M);             

    // U[0,1) uniform distribution
    for (int i = 0; i < M; ++i) {
        x[i] = rnd.Rannyu();  
    }

    // Function f(x)
    for (int i = 0; i < M; ++i) {
        f[i] = M_PI / 2 * cos(M_PI / 2 * x[i]); 
    }

    //Average calculation
    auto result_ave_unif = datablocking_ave(f, M, N);
    std::vector<double> ave_unif = result_ave_unif.first;
    std::vector<double> err_ave_unif = result_ave_unif.second;

    //Variance calculation
    auto result_var_unif = datablocking_var(f, M, N, ave_unif[N-1]);
    std::vector<double> var_unif = result_var_unif.first;
    std::vector<double> err_var_unif = result_var_unif.second;


    // CREATION OF OUTPUT FILE
    string Uniform_distribution = "Uniform_distribution.txt";
    scriviSuFile(Uniform_distribution, ave_unif, err_ave_unif, var_unif, err_var_unif);


    //MC INTEGRAL WITH p(r)=2(1-r)

    //definition of variables
    vector<double> r(M);
    vector<double> g(M);

    //Inverso della comulativa:
    for (int i = 0; i < M; ++i) {
        r[i] = 1 - sqrt( 1 - rnd.Rannyu() );
    }

    // Function g(x)
    for (int i = 0; i < M; ++i) {
        g[i] = M_PI / 2 * cos(M_PI / 2 * r[i])/(2 * (1 - r[i])); 
    }

    //Average calculation
    auto result_ave_imp = datablocking_ave(g, M, N);
    std::vector<double> ave_imp = result_ave_imp.first;
    std::vector<double> err_ave_imp = result_ave_imp.second;

    //Variance calculation
    auto result_var_imp = datablocking_var(g, M, N, ave_imp[N-1]);
    std::vector<double> var_imp = result_var_imp.first;
    std::vector<double> err_var_imp = result_var_imp.second;


    // CREATION OF OUTPUT FILE
    string Importance_sampling = "Importance_sampling.txt";
    scriviSuFile(Importance_sampling, ave_imp, err_ave_imp, var_imp, err_var_imp);

    return 0;
}