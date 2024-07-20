#include "funzioni.h"
#include <utility>  // per std::pair

using namespace std;

//FUNCTION'S LIBRARY

//AVERAGE CALCULATION
double calcolaMedia(const std::string& nomeFile) {
    std::ifstream fileInput(nomeFile);

    if (!fileInput.is_open()) {
        std::cerr << "Impossibile aprire il file: " << nomeFile << std::endl;
        return 0.0; // Ritorna 0 in caso di errore
    }

    std::vector<double> dati;
    double valore;
    
    // Reading from file
    while (fileInput >> valore) {
        dati.push_back(valore);
    }

    fileInput.close();

    if (dati.empty()) {
        std::cerr << "Il file non contiene dati validi." << std::endl;
        return 0.0; // Ritorna 0 se non ci sono dati validi nel file
    }

    // Average calculation
    double somma = 0.0;
    for (const auto& dato : dati) {
        somma += dato;
    }

    double media = somma / dati.size();

    return media;
}

//ERROR CALCULATION
double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n) {
    if (n == 0) {
        return 0.0;
    } else {
        return std::sqrt((AV2[n] - std::pow(AV[n], 2)) / n);
    }
}

// CREATION OF OUTPUT FILE
void scriviSuFile(const string& nomeFile, const vector<double>& dati1, const vector<double>& dati2) {
    std::ofstream fileOutput(nomeFile);
    if (!fileOutput.is_open()) {
        std::cerr << "Impossibile creare il file di output: " << nomeFile << std::endl;
        return;
    }

    // Writing
    for (size_t i = 0; i < dati1.size(); ++i) {
        fileOutput << dati1[i] << "\t" << dati2[i] << std::endl;
    }
    fileOutput.close();

    std::cout << "I risultati sono stati scritti su: " << nomeFile << std::endl;
}

//DATA BLOCKING AVERAGE
// Function that calculates the average of a function using the data blocking method. 
// Takes as input a function f, the number of throws M, and the number of blocks N.
std::pair<std::vector<double>, std::vector<double>> datablocking_ave(const std::vector<double>& f, int M, int N) {

    //Definition of variables
    const int n = M / N;                    // Number of throws in each block, ensure M is a multiple of N
    vector<double> ave(N, 0.0);             // Average
    vector<double> av2(N, 0.0);             // Square Average
    vector<double> sum_prog(N, 0.0);        // Progessive sum
    vector<double> su2_prog(N, 0.0);        // Progressive square sum
    vector<double> err_prog(N, 0.0);        // Statistical uncertainty
    double sum1 = 0.0;

    // Average calculation
    for (int i = 0; i < N; ++i) {
        sum1 = 0.0;
        for (int j = 0; j < n; ++j) {
            int k = j + i * n;
            sum1 += f[k]; 
        }
        ave[i] = sum1 / n; // f(x_i) 
        av2[i] = pow(ave[i], 2); // f(x_i)^2 
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            sum_prog[i] += ave[j]; // SUM_{j=0,i} f(x_i)
            su2_prog[i] += av2[j]; // SUM_{j=0,i} f(x_i)^2
        }
        sum_prog[i] /= (i + 1); // Cumulative average
        su2_prog[i] /= (i + 1); // Cumulative square average
        err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
    }
    return std::make_pair(sum_prog, err_prog);
}

//DATA BLOCKING VARIANCE
// Function that calculates the variance using the data blocking method. 
// Takes as input a function f, the number of throws M, and the number of blocks N and the expected value I.
std::pair<std::vector<double>, std::vector<double>> datablocking_var(const std::vector<double>& f, int M, int N, double I) {

    //Definition of variables
    const int n = M / N;                    // Number of throws in each block, ensure M is a multiple of N
    vector<double> ave(N, 0.0);             // Average
    vector<double> av2(N, 0.0);             // Square Average
    vector<double> sum_prog(N, 0.0);        // Progessive sum
    vector<double> su2_prog(N, 0.0);        // Progressive square sum
    vector<double> sum_progeps(N, 0.0);        // Progessive sum
    vector<double> su2_progeps(N, 0.0);        // Progressive square sum
    vector<double> err_prog(N, 0.0);        // Statistical uncertainty
    double sum1 = 0.0;

    // Average calculation
    for (int i = 0; i < N; ++i) {
        sum1 = 0.0;
        for (int j = 0; j < n; ++j) {
            int k = j + i * n;
            sum1 += pow(f[k] ,2 ) - pow(I,2); 
        }
        ave[i] = sum1 / n ; //- pow(I,2); // f(x_i) 
        av2[i] = pow(ave[i], 2); // f(x_i)^2 
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            sum_prog[i] += ave[j]; // SUM_{j=0,i} f(x_i)
            su2_prog[i] += av2[j]; // SUM_{j=0,i} f(x_i)^2
        }
        sum_prog[i] /= (i + 1); // Cumulative average
        su2_prog[i] /= (i + 1); // Cumulative square average
    //   sum_progeps[i] = sqrt(sum_prog[i]/(i+1)); // Cumulative average
    //  su2_progeps[i] = sqrt(su2_prog[i]/(i+1)); // Cumulative square average
        err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
    }
    return std::make_pair(sum_prog, err_prog);
}
