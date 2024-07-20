#include "funzioni.h"

using namespace std;

//Libreria generale di funzioni 



//Prende dati da un file e ne calcola la media
double calcolaMedia(const std::string& nomeFile) {
    std::ifstream fileInput(nomeFile);

    if (!fileInput.is_open()) {
        std::cerr << "Impossibile aprire il file: " << nomeFile << std::endl;
        return 0.0; // Ritorna 0 in caso di errore
    }

    std::vector<double> dati;
    double valore;
    
    // Leggi i dati dal file
    while (fileInput >> valore) {
        dati.push_back(valore);
    }

    fileInput.close();

    if (dati.empty()) {
        std::cerr << "Il file non contiene dati validi." << std::endl;
        return 0.0; // Ritorna 0 se non ci sono dati validi nel file
    }

    // Calcola la media
    double somma = 0.0;
    for (const auto& dato : dati) {
        somma += dato;
    }

    double media = somma / dati.size();

    return media;
}


//Calcola l'errore prendendo come input media, media^2 e n blocchi
double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n) {
    if (n == 0) {
        return 0.0;
    } else {
        return std::sqrt((AV2[n] - std::pow(AV[n], 2)) / n);
    }
}

//Scrive su un file di output i risultati di un calcolo. In questo caso per ogni riga inserisce il risultato di due set di dati (media e errore)
void scriviSuFile(const string& nomeFile, const vector<double>& dati1) {
    std::ofstream fileOutput(nomeFile);

    if (!fileOutput.is_open()) {
        std::cerr << "Impossibile creare il file di output: " << nomeFile << std::endl;
        return;
    }

    // Scrivi i risultati sul file di output
    for (size_t i = 0; i < dati1.size(); ++i) {
        fileOutput << dati1[i] << std::endl;
    }
    fileOutput.close();

    std::cout << "I risultati sono stati scritti su: " << nomeFile << std::endl;
}
