#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib> // For the exit function
#include <string>  // For std::stoi
#include "Metropolis.h"


using namespace std;

// Directory where results will be saved
const string ResultsDirectory = "Results/";

// Quantum numbers for the wave function
const int n = 2;
const int l = 1;
const int m = 0;

// Function that returns the probability density of the wave function for given coordinates (x, y, z)
double wavef_prob(double x, double y, double z){
    if (n == 1 and l == 0 and m == 0){
        return (1. / M_PI) * exp(-2 * sqrt(x * x + y * y + z * z));
    }
    else if (n == 2 and l == 1 and m == 0){
        return (1. / 32.) * (1. / M_PI) * exp(-sqrt(x * x + y * y + z * z)) * z * z;
    }
    else { 
        return 0; 
    }
}

// Function that returns the radial distance from the origin for given coordinates (x, y, z)
double r(double x, double y, double z){
    return sqrt(x * x + y * y + z * z);
}

int main(){


    //METROPOLIS
    Metropolis metro;  // Create a Metropolis object
    metro.initialize("metropolis.dat");  // Initialize the Metropolis object with parameters from an input file

    // Retrieve the total number of steps, number of blocks, and steps per block
    const int M = metro.get_nsteps();
    const int N = metro.get_nblocks();
    const int L = metro.get_blksteps();

    string resultsf;  // String to hold the results file name

    // Determine the results file name based on the quantum numbers and distribution type
    if (n == 1 and l == 0 and m == 0 and metro.get_distype() == 0) {
        resultsf = "1s-uniform.out";
    } else if (n == 1 and l == 0 and m == 0 and metro.get_distype() == 1) {
        resultsf = "1s-gaussian.out";
    } else if (n == 2 and l == 1 and m == 0 and metro.get_distype() == 0) {
        resultsf = "2p-uniform.out";
    } else if (n == 2 and l == 1 and m == 0 and metro.get_distype() == 1) {
        resultsf = "2p-gaussian.out";
    }

    vec ave(N);  // Vector to store the averages of the radial distances

    double appo;  // Variable to accumulate the radial distances
    bool appo_;  // Variable to store the result of the Metropolis step

    // Open a file to write the coordinates of accepted steps
    ofstream print_steps;
    print_steps.open(ResultsDirectory + "steps-" + resultsf);

    // Outer loop over the number of blocks
    for (int i = 0; i < N; i++){
        appo = 0.;  // Reset the accumulator

        // Inner loop over the steps in each block
        for (int j = i * L; j < L * (i + 1); j++){
            appo_ = metro.step(wavef_prob);  // Perform a Metropolis step
            appo += r(metro.get_coord(0), metro.get_coord(1), metro.get_coord(2));  // Accumulate the radial distance

            // If the step was accepted, write the coordinates to the output file
            if (appo_){
                print_steps << metro.get_coord(0) << " " << metro.get_coord(1) << " " << metro.get_coord(2) << endl;
            }
        }

        ave[i] = appo / L;  // Store the average radial distance for the current block
    }

    // Close the file that recorded the coordinates of accepted steps
    print_steps.close();

    metro.print_results(ave, ResultsDirectory+resultsf);

    return 0;  
}