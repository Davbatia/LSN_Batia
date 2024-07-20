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
const int M = pow(10, 4);  // Total number of throws
const int N = pow(10,2);     // Number of blocks
const int n = M / N;   // Number of throws in each block, ensure M is a multiple of N
const int steps =100;

//RANDOM WALK 

//Random Walk discreto: prendiamo un numero pseudocasuale tra [0,1], 
//lo moltiplichiamo per 6 e ne prendiamo la parte intera;
//in base al valore di quest'ultima facciamo un passo lungo x, y o z, in avanti o indietro.
vector<double> DRW(vector<double> &v, double y){
    y = 6. * y;
    int direction = static_cast<int>(floor(y));
    int axis = direction / 2;
    int step = (direction % 2 == 0) ? 1 : -1;
    v[axis] += step;
    return v;
}

vector<double> CRW(vector<double> &v, double theta, double phi){
	v[0] += sin(theta)*cos(phi);
	v[1] += sin(theta)*sin(phi);
	v[2] += cos(theta);
	return v;
};

// MAIN

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

    //RANDOM WALK

    //Initialize vectors to store positions

    //Discrete
    vector<double> v_discreto = {0., 0., 0.};
    double norm_discreta = 0.;
    vector<double> discreto_average_sqrt(steps, 0.0);
    vector<double> discreto_block(steps, 0.0);
    vector<double> discreto_block_sqrt(steps, 0.0);
    vector<double> discreto_average(steps, 0.0);
    vector<double> discreto_average2(steps, 0.0);
    vector<double> discreto_error (steps,0.0);

    //Continuum
    vector<double> v_continuum = {0., 0., 0.};
    double norm_continuum = 0.;
    vector<double> continuum_average_sqrt(steps, 0.0);
    vector<double> continuum_block(steps, 0.0);
    vector<double> continuum_block_sqrt(steps, 0.0);
    vector<double> continuum_average(steps, 0.0);
    vector<double> continuum_average2(steps, 0.0);
    vector<double> continuum_error (steps,0.0);

    // Loop over blocks
    for (int j = 0; j < N; j++) {

        // Reset the block accumulators

        fill(discreto_block.begin(), discreto_block.end(), 0.0);
        fill(discreto_block_sqrt.begin(), discreto_block_sqrt.end(), 0.0);

        fill(continuum_block.begin(), continuum_block.end(), 0.0);
        fill(continuum_block_sqrt.begin(), continuum_block_sqrt.end(), 0.0);

        // Loop over throws in block
        for (int k = 0; k < n; k++) {
            v_discreto = {0., 0., 0.};
            v_continuum = {0., 0., 0.};

            // Loop over random walk steps
            for (int i = 0; i < steps; i++) {

                DRW(v_discreto, rnd.Rannyu());
                norm_discreta = pow(v_discreto[0], 2) + pow(v_discreto[1], 2) + pow(v_discreto[2], 2);
                discreto_block[i] += norm_discreta;
                discreto_block_sqrt[i] += sqrt(norm_discreta);

                CRW(v_continuum, acos( 2.*rnd.Rannyu() - 1.), rnd.Rannyu(0., 2 * M_PI));
                norm_continuum = pow(v_continuum[0], 2) + pow(v_continuum[1], 2) + pow(v_continuum[2], 2);
                continuum_block[i] += norm_continuum;
                continuum_block_sqrt[i] += sqrt(norm_continuum);
            }
        }

        // Compute the block average for each step
        for (int i = 0; i < steps; i++) {

            discreto_average[i] += discreto_block[i] / (double)n;
            discreto_average_sqrt[i] += discreto_block_sqrt[i] / (double)n;
            discreto_average2[i] += pow(discreto_block_sqrt[i] / (double)n, 2);

            continuum_average[i] += continuum_block[i] / (double)n;
            continuum_average_sqrt[i] += continuum_block_sqrt[i] / (double)n;
            continuum_average2[i] += pow(continuum_block_sqrt[i] / (double)n, 2);
        }
        
    }

    // Finalize the averages and compute the errors
    for (int i = 0; i < steps; i++) {

        discreto_average[i] /= N;
        discreto_average_sqrt[i] /= N;
        discreto_average2[i] /= N;
        discreto_error[i] = sqrt((discreto_average2[i] - pow(discreto_average_sqrt[i], 2)) / (N - 1));

        continuum_average[i] /= N;
        continuum_average_sqrt[i] /= N;
        continuum_average2[i] /= N;
        continuum_error[i] = sqrt((continuum_average2[i] - pow(continuum_average_sqrt[i], 2)) / (N - 1));
    }
    
    // CREATION OF OUTPUT FILE

    string DiscreteRW = "DiscreteRW.txt";
    scriviSuFile(DiscreteRW, discreto_average_sqrt, discreto_error);

    string ContinuousRW = "ContinuousRW.txt";
    scriviSuFile(ContinuousRW, continuum_average_sqrt, continuum_error);



    return 0;
}