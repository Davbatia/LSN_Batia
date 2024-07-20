#include "Metropolis.h"  // Include the header file for the Metropolis class
#include <cmath>
using namespace std;  // Use the standard namespace

// Sets the coordinate value at index i to x
void Metropolis :: setcoord(int i, double x){
    _coords(i) = x;
}

// Returns the coordinate value at index i
double Metropolis :: get_coord(int i){  
    return _coords(i);
}

// Increments the step index by 1
void Metropolis :: increase_step(){
    _step_index++;
}

// Initializes the Metropolis object with parameters from an input file
void Metropolis :: initialize(string inputf){

    double coords_temp;
    ifstream input(inputf);  // Open the input file

    // Read the input file until the end of the file
    while (!input.eof()){
        input >> property;
        if (property == "DIMENSION"){
            input >> _ndim;
            _coords.resize(_ndim);  // Resize the coordinate vector based on the number of dimensions
        }
        else if (property == "STARTING_POSITION"){
            for (int i = 0; i < _ndim; i++){
                input >> coords_temp;
                this->setcoord(i, coords_temp);  // Set the starting positions
            }
        }
        else if (property == "DISTRIBUTION"){
            input >> _dis_type;  // Set the distribution type
        }
        else if (property == "STEPS"){
            input >> _nsteps;  // Set the number of steps
        }
        else if (property == "BLOCKS"){
            input >> _nblocks;  // Set the number of blocks
        }
        else if (property == "ENDINPUT"){
            cout << "Reading input completed!" << endl;
            break;  // End of input
        }
    }
    input.close();
    _blksteps = _nsteps / _nblocks;  // Calculate the number of steps per block
    _acceptance.resize(_nblocks);  // Resize the acceptance vector
}

// Performs a single Metropolis step
bool Metropolis :: step(double (*fun)(double, double, double)){
    double a = 0.;
    double r = 0.;
    vec proposed_coord;
    vec attempted_step;

    proposed_coord.resize(_ndim);
    attempted_step.resize(_ndim);

    if (_dis_type == 0){
        for (int i=0; i<3; i++)
            attempted_step(i) = _rnd.Rannyu(-3,3);   //Uniform: _rnd.Rannyu(-1.3,1.3) ok for [n,l,m]=[1,0,0], _rnd.Rannyu(-3,3) ok for [n,l,m]=[2,1,0]
    }                                                // ok => acceptance ≈ 50%
    else if (_dis_type == 1){
        for (int i=0; i<3; i++)
            attempted_step(i) = _rnd.Gauss(0,1.8);   //Gaussian: _rnd.Gauss(0,0.8) ok for [n,l,m]=[1,0,0], _rnd.Gauss(0,1.8) ok for [n,l,m]=[2,1,0]
    }                                                // ok => acceptance ≈ 50%

    proposed_coord = _coords + attempted_step;

    // Calculate the acceptance probability
    a = min(1., fun(proposed_coord(0), proposed_coord(1), proposed_coord(2)) / fun(_coords(0), _coords(1), _coords(2)));
    r = _rnd.Rannyu();

    this->increase_step();  // Increment the step index

    // Update acceptance statistics at the end of each block
    if (_step_index % _blksteps == 0) {
        _acceptance(_blk_index) = static_cast<double>(_acc_step) / _blksteps;
        _blk_index++;
        _acc_step = 0;  // Reset the acceptance count for the new block
    }

    // Accept or reject the proposed step
    if (r <= a) {
        _coords = proposed_coord;  // Update coordinates
        _acc_step++;  // Increment the acceptance count
        return true;
    }

    return false;  // Step was rejected
}

// Returns the total number of steps
int Metropolis :: get_nsteps(){
    return _nsteps;
}

// Returns the total number of blocks
int Metropolis :: get_nblocks(){
    return _nblocks;
}

// Returns the number of steps per block
int Metropolis :: get_blksteps(){
    return _blksteps;
}

//Error
double Metropolis :: error_value(double sum, double sum2, int n) {
    if (n == 0) {
        return 0; // Per evitare la divisione per zero
    } else {
        return sqrt((sum2 - sum * sum) / n);
    }
}

// Definition of the min function
double Metropolis::min(double a, double b) {
    return (a < b) ? a : b;
}

// Writes the results of the simulation to an output file
void Metropolis :: print_results(vec ave, string outputf){

    ofstream WriteResults;
    WriteResults.open(outputf);

    vec ave2 = ave % ave;

    vec prog_sum(_nblocks);
    vec prog_sum2(_nblocks);

    WriteResults << "# BLOCK" << setw(10) << "OBS_AVE" << setw(15) << "ERROR" << setw(15) << "ACCEPTANCE" << endl;

    // Calculate and write progressive averages and errors for each block
    for (int i = 0; i < _nblocks; i++){
        for(int j = 0; j <= i; j++){
            prog_sum[i] += ave[j];
            prog_sum2[i] += ave2[j];
        }

        prog_sum[i] /= (i + 1);
        prog_sum2[i] /= (i + 1);

       WriteResults << setw(7) << i << setw(10) << prog_sum[i] << setw(15) << Metropolis :: error_value(prog_sum[i], prog_sum2[i], i) << setw(15) << _acceptance[i] << endl;
    }

    WriteResults.close();
}

// Returns the type of distribution
int Metropolis :: get_distype(){
    return _dis_type;
}