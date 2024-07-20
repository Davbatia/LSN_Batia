/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  if (argc != 2){
        cout << "Usage: " << argv[0] << " <SOLID/LIQUID/GAS>" << endl;
        exit(-1);
    }

  
  string modello_input = string(argv[1])+".dat"; // Costruisce il nome del file di input
  string modello_output = string(argv[1]);       //Per il file di output
  int nconf = 1;
  System SYS;
  //SYS.set_temp(stod(argv[2]));
  SYS.initialize(modello_input, modello_output);
  SYS.initialize_properties(modello_output);
  SYS.block_reset(0, modello_output);

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(j%10 == 0){
//        SYS.write_XYZ(nconf, modello_output); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
      // SYS.print_istant(outputd);
    }
    SYS.averages(i+1, modello_output);
    SYS.block_reset(i+1, modello_output);
  }
  SYS.finalize(modello_output);

  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
