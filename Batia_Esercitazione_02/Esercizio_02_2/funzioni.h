
#ifndef FUNZIONI_H
#define FUNZIONI_H


#include <iostream>                      // #ifndef __header_h__
#include <iomanip>                       // #define __header_h__
#include <cmath>                         // #include<cstdlib>
#include <string>                        // In fondo devo mettere: #endif
#include <fstream>                      
#include <vector>
#include <cstdlib>
#include <ctime>


using namespace std;

//Prende dati da un file e ne calcola la media
double calcolaMedia(const std::string& nomeFile); 

//Calcola l'errore prendendo come input media, media^2 e n blocchi
double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n); 

//Calcola l'errore 
double error_value(double sum, double sum2, int n);

//Scrive su un file di output i risultati di un calcolo
void scriviSuFile(const string& nomeFile, const vector<double>& dati1, const vector<double>& dati2); 

// Function that calculates the average of a function using the data blocking method
std::pair<std::vector<double>, std::vector<double>> datablocking_ave(const std::vector<double>& f, int M, int N);

// Function that calculates the variance using the data blocking method
std::pair<std::vector<double>, std::vector<double>> datablocking_var(const std::vector<double>& f, int M, int N, double I);

#endif // FUNZIONI_H
