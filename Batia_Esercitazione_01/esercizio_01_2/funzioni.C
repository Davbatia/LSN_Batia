#include "funzioni.h"

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
void scriviSuFile(const string& nomeFile, const vector<double>& dati1, const vector<double>& dati2, const vector<double>& dati3, const vector<double>& dati4, const vector<double>& dati5, const vector<double>& dati6, const vector<double>& dati7, const vector<double>& dati8, const vector<double>& dati9, const vector<double>& dati10, const vector<double>& dati11, const vector<double>& dati12) {
    std::ofstream fileOutput(nomeFile);
    if (!fileOutput.is_open()) {
        std::cerr << "Impossibile creare il file di output: " << nomeFile << std::endl;
        return;
    }

    // Writing
    for (size_t i = 0; i < dati1.size(); ++i) {
        fileOutput << dati1[i] << "\t" << dati2[i] << "\t" << dati3[i] << "\t" << dati4[i] << "\t" << dati5[i] << "\t" << dati6[i] << "\t" << dati7[i] << "\t" << dati8[i] << "\t" << dati9[i] << "\t" << dati10[i] << "\t" << dati11[i] << "\t" << dati12[i] << std::endl;
    }
    fileOutput.close();

    std::cout << "I risultati sono stati scritti su: " << nomeFile << std::endl;
}
