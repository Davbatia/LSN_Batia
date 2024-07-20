
#ifndef FUNZIONI_H
#define FUNZIONI_H


#include <iostream>                      // #ifndef __header_h__
#include <iomanip>                       // #define __header_h__
#include <cmath>                         // #include<cstdlib>
#include <string>                        // In fondo devo mettere: #endif
#include <fstream>                       // Sono le istruzioni per le Turbo Funzioni

#include <vector>
#include <cstdlib>
#include <ctime>


using namespace std;

double calcolaMedia(const std::string& nomeFile); //Prende dati da un file e ne calcola la media

double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n); //Calcola l'errore prendendo come input media, media^2 e n blocchi

void scriviSuFile(const string& nomeFile, const vector<double>& dati1, const vector<double>& dati2, const vector<double>& dati3, const vector<double>& dati4); //Scrive su un file di output i risultati di un calcolo

#endif // FUNZIONI_H
