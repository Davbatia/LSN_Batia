#include<iostream>

#include"Genetic_Algorithm.h"
#include"random.h"
#include<fstream>
#include<iomanip>

using namespace std;



void initRand(Random &rnd){
}

int main(int argc, char *argv[]){

    // Variabili per i parametri
    unsigned int N; // numero citta
    unsigned int a; // disposizione, circolare se a=0, quadrata altrimenti
    unsigned int M; // numero cromosomi
    unsigned int nsteps; // numero di passi
    double pcrossover; // probabilita di riproduzione
    double pmutation; // probabilita di mutare

    // Lettura dei dati dal file input.dat
    ifstream inputFile("input.dat");
    if (inputFile.is_open()) {
        string line;
        
        while (getline(inputFile, line)) {
            if (line.find("N_Città") != string::npos) {
                N = stoi(line.substr(line.find_last_of(" \t") + 1));
            } else if (line.find("Disposizione") != string::npos) {
                a = stoi(line.substr(line.find_last_of(" \t") + 1));
            } else if (line.find("N_Cromosomi") != string::npos) {
                M = stoi(line.substr(line.find_last_of(" \t") + 1));
            } else if (line.find("N_steps") != string::npos) {
                nsteps = stoi(line.substr(line.find_last_of(" \t") + 1));
            } else if (line.find("p_crossover") != string::npos) {
                pcrossover = stod(line.substr(line.find_last_of(" \t") + 1));
            } else if (line.find("p_mutation") != string::npos) {
                pmutation = stod(line.substr(line.find_last_of(" \t") + 1));
            }
        }
        inputFile.close();
    } else {
        cerr << "Errore nell'apertura del file input.dat" << endl;
        return 1;
    }

    initRand(rnd);



    vector<Posizione> map=Mappa(N,a);
    // la stampo su un file
    ofstream MappaData;
    string type;
    if(a==0)
        type="circle";
        
    else
        type="square";
        

    MappaData.open("Results/map_"+type+".data");
    for(unsigned int i=0; i<map.size();i++){
        MappaData<<map[i].GetX()<<setw(12)<<map[i].GetY()<<endl;
    }

    Population pop(map,M,pcrossover,pmutation);

    ofstream bestpathLength, bestpath,input;

    if(pcrossover==0){
        bestpathLength.open("Results/bestlength_"+type+"_onlymutation.data");
        bestpath.open("Results/best_onlymutation_"+type+".data"); 
        //input.open("input_onlymutation_"+type+".data");  
    }else{
        bestpathLength.open("Results/bestlength_"+type+".data");
        bestpath.open("Results/best_"+type+".data");
        //input.open("input_"+type+".data");   
    }

    

    for( unsigned int i=0; i<nsteps; i++){
        pop.Move();
        cout<<"Mossa fatta"<<endl;
        for(unsigned int i=0; i< M; i++)
            bestpathLength<<setw(12)<<pop.m_pop[i].L2;
            
        bestpathLength<<endl;
        
    }

    
    
    for(unsigned int i=0; i<N; i++){
        bestpath<<setw(12)<<pop.m_pop[0].m_p[i];
    }

    

//    input<<N<<endl<<M<<endl<<nsteps<<endl<<pcrossover<<endl<<pmutation<<endl;

    
    
    
    




    return 0;
}