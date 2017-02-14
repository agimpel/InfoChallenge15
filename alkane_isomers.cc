/* Generation and enumeration of all alkane constitutional isomers from methane to icosane   */
/* ========================================================================================= */
/* using an altered canonical isomer code and Morgan's Algorithm                             */
/*                                                                                           */
/* by Andreas Gimpel, agimpel@student.ethz.ch                                                */
/* student of the Department of Chemistry and Applied Biosciences, ETH Zürich                */
/*                                                                                           */
/* entry for the challenge issued by Prof. Philippe H. Hünenberger                           */
/* in his lecture "Informatik I for Chemists" on 29.10.2015                                  */
/*                                                                                           */
/*                                                                                           */
/*                                                                                           */
/* Using the g++ compiler, the program was compiled with                                     */
/* g++ -O3 alkane_isomers.cc                                                                 */
/*                                                                                           */
/* For an explanation of the algorithms and general structure used, refer to the .pdf.       */
/* There are many comments in the code which should allow for a basic understanding of the   */
/* methodology as well as the concept this algorithm employs.                                */
/*                                                                                           */
/* NOTE: In order for the file output to work, an (empty) directory isomer/ has to be        */
/* present in the directory the program is compiled to. It is included in this .zip file.    */


#include <iostream>                   //iostream for general output and input
#include <fstream>                    //fstream to create output files with all isomer codes, if user wishes to do so
#include <vector>                     //vector library to create an array-like structure capable of storing 22x500000x22 entries

using namespace std;



// ==============================================================================================================================================================================================================
// GLOBAL VARIABLES
// ==============================================================================================================================================================================================================

const int firstCarbon = 1;            //carbon amount of first alkane, cannot be changed as algorithm is recursive and depends on methane as its base
const int lastCarbon  = 20;           //amount of carbons the last alkane generated has. Can be changed for debugging purposes. Not tested for >20 (maxIsomers has to be changed for >20)
const int maxIsomers  = 500000;       //maximum amount of individual isomers that can be stored. This does not include isomers that are ignored due to not being unique



// ==============================================================================================================================================================================================================
// FUNCTION PROTOTYPES (function descriptions can be found at the function declaration)
// ==============================================================================================================================================================================================================

//generation function group
void generate_Isomers(vector<vector<vector<int> > > &isomers, int CarbonAmount, int maxIsomers, vector<vector<int> > &morgans);
int Isomer_digit_validity_check(vector<vector<vector<int> > > &isomers, int previousC, int C, int I, int CarbonAmount); 

//examination function group
int check_Isomers(vector<vector<vector<int> > > &isomers, int isomer, int CarbonAmount, vector<vector<int> > &morgans);
void morgans_splicing(vector<vector<vector<int> > > &isomers, int isomer, int CarbonAmount, vector<vector<int> > &morgans, vector<vector<int> > &connections);
void morgans_algorithm(int isomer, int CarbonAmount, vector<vector<int> > &morgans, vector<vector<int> > &connections);
void morgans_sort(int isomer, int CarbonAmount, vector<vector<int> > &morgans);
int check_morgan_uniqueness(int isomer, int CarbonAmount, vector<vector<int> > &morgans);

//ui function group
bool print_intro();
void print_structure(vector<vector<vector<int> > > &isomers);
void print_Isomers(vector<vector<vector<int> > > &isomers, int CarbonAmount, string filenames[], vector<vector<int> > &morgans, bool generate_files);



// ==============================================================================================================================================================================================================
// MAIN
// ==============================================================================================================================================================================================================

int main(){
  int CarbonCount = lastCarbon - firstCarbon + 1;                                                                       //amount of alkanes generated in total
  int arrayCount = CarbonCount + 2;                                                                                     //for usability, every vector uses vector[1...CarbonCount+2] with vector[0] reserved for
                                                                                                                        //additional information (ie length of code, amount of isomers

  vector<vector<vector<int> > > isomers(arrayCount, vector<vector<int> >(maxIsomers, vector<int>(arrayCount, -1)));     //main vector for isomer storage is initialized with a dimension of 22x500000x22,
                                                                                                                        //where the first dimension is the carbon amount of current alkane, second dimension
                                                                                                                        //represents the current isomer and the third dimension representing the digit of code


  string filenames[21] = {"isomer/0.isomers",  "isomer/1.isomers",                                                      //string array to store the file names for output files
                          "isomer/2.isomers",  "isomer/3.isomers", 
                          "isomer/4.isomers",  "isomer/5.isomers", 
                          "isomer/6.isomers",  "isomer/7.isomers",
                          "isomer/8.isomers",  "isomer/9.isomers", 
                          "isomer/10.isomers", "isomer/11.isomers",
                          "isomer/12.isomers", "isomer/13.isomers",
                          "isomer/14.isomers", "isomer/15.isomers",
                          "isomer/16.isomers", "isomer/17.isomers", 
                          "isomer/18.isomers", "isomer/19.isomers", 
                          "isomer/20.isomers"};


  bool generate_files=print_intro();  //user introduction and input on whether output files will be generated




 //METHANE DECLARATION
 isomers[1][1][1]=0;    //declare methane in isomer code  -> methane is [0]
 isomers[1][1][0]=0;
 isomers[1][0][0]=1;    //methane has 1 isomer

 print_structure(isomers);       //print table structure and methane



 //MAIN ALKANE LOOP
 for(int CarbonAmount=2; CarbonAmount<=lastCarbon; CarbonAmount++){                    //for all amounts of carbon atoms between 2 and lastCarbon (including)
   vector<vector<int> > morgans(maxIsomers+2, vector<int>(CarbonAmount+2, 0));         //initialize a vector for this alkane isomer's morgans codes
   generate_Isomers(isomers, CarbonAmount, maxIsomers, morgans);                       //generate all possible isomers using the altered canonical representation
   print_Isomers(isomers, CarbonAmount, filenames, morgans, generate_files);           //output all information
 }
}



// ==============================================================================================================================================================================================================
// FUNCTION DECLARATIONS
// ==============================================================================================================================================================================================================


//GENERATION FUNCTION GROUP

//generate_Isomers will use the previous alkane's isomer codes to generate all possible isomers of the current alkane by incrementing each digit once. This is equivalent to generating isomers by attaching a
//carbon atom to an existing isomer structure to obtain an isomer of the next alkane.
//The function loops over every isomer and every digit of the previous alkane and checks if the chosen digit can be incremented (see function Isomer_digit_validity_check) and does so if it is allowed.
//It copies all digits up to the chosen digit to the current alkane isomer code, increments the current digit and adds a 0 behind it. Then the rest of the isomer code is copied and the whole isomer code
//is passed to function check_isomer to examine the isomer's uniqueness. If it is unique, the amount of valid isomers (stored at isomers[CarbonAmount][0][0]) is incremented and a new isomer is generated.
//Otherwise, the isomer will be dropped and the next isomer will be generated in its place.

void generate_Isomers(vector<vector<vector<int> > > &isomers, int CarbonAmount, int maxIsomers, vector<vector<int> > &morgans){
  int  previousC = CarbonAmount - 1;                             //previous alkane has previousC carbon atoms
  int  previousIsomers = isomers[previousC][0][0];               //previous alkane has previousIsomers different isomers
  int  isomer = 1;                                               //current alkane starts with isomer 1
  isomers[CarbonAmount][0][0] = 0;                               //so far 0 valid isomers of new alkane were found


  for(int I=1; I<=previousIsomers; I++){                                                    //for all isomers of previous alkane

    for(int C=1; C<=previousC; C++){                                                        //for every carbon atom of this isomer

      if(Isomer_digit_validity_check(isomers, previousC, C, I, CarbonAmount)){              //check validity of chosen isomer digit if incremented and if it is allowed:

	for(int transfer=1; transfer<C; transfer++){                                        //transfer all characters up to C into new alkane isomer code 
	  isomers[CarbonAmount][isomer][transfer] = isomers[previousC][I][transfer]; 
        }

	isomers[CarbonAmount][isomer][C] = isomers[previousC][I][C]+1;                      //increment chosen carbon by 1 and transfer new value to current alkane isomer
	isomers[CarbonAmount][isomer][C+1]=0;                                               //add 0 after it (new branch is not connected any further)

	for(int transfer=C+1; transfer<=CarbonAmount; transfer++){                          //transfer all digits after current C into new alkane isomer at C+1
	  isomers[CarbonAmount][isomer][transfer+1] = isomers[previousC][I][transfer];
	}

	if(check_Isomers(isomers, isomer, CarbonAmount, morgans)){                          //check new generated isomer for validity and uniqueness
	  isomer++;                                                                         //if it is unique, the next isomer will be generated on the next position and the amount of valid codes is increased
          isomers[CarbonAmount][0][0]++;                                                    //otherwise, current isomer will be overwritten with the next isomer

        }
      }
    }
  }
}



//Isomer_digit_validity_check is necessary to limit the possibilities of isomer generation based on the valence properties and the canonical code representation. As every carbon atom may only have 4 bonds,
//the root may not be higher than 4 (-> 4 forward bonds) and every non-root may not be higher than 3 (-> 3 forward bonds and 1 backward bond = 4 bonds). In addition, the non-roots cannot have more bonds than
//the root, as the root is always one of the atoms with the highest amount of bonds. This implementation simply checks for these conditions and thus prevents the entire generation and examination of an isomer
//that violates these rules. This check is vital to limit the amount of isomers generated and reduce computation time based on simple conditions.

int Isomer_digit_validity_check(vector<vector<vector<int> > > &isomers, int previousC, int C, int I, int CarbonAmount){

  bool is_root_less_4 = (C==1 && isomers[previousC][I][1]<4);                                                                           //root cannot be greater than 4
  bool not_root_less_3_less_root = (!(C==1) && isomers[previousC][I][C]<3 && ((isomers[previousC][I][C]+1)<isomers[previousC][I][1]));  //non-roots cannot be greater than root and 3

  if(is_root_less_4 || not_root_less_3_less_root){                                                                                      //if digit meets condition, return true, else false
    return 1;
  } else {
    return 0;
  }
}



//EXAMINATION FUNCTION GROUP

//check_isomer will return TRUE or FALSE depending on the uniqueness of the generated isomer. In order to judge the uniqueness, the canonical isomer code is translated into a sequence of morgans's algorithm
//codes ordered by value which is compared to all existing morgan's codes. If no other isomer with this code exists, the isomer is unique and TRUE is returned.
//The translation is done by morgans_splicing, the algorithm is executed by morgans_algorithm and the values are sorted by morgans_sort using insertion sort. The comparison to all generated isomers is
//done by check_morgan_uniqueness.

int check_Isomers(vector<vector<vector<int> > > &isomers, int isomer, int CarbonAmount,  vector<vector<int> > &morgans){
 
  vector<vector<int> > connections(CarbonAmount+2, vector<int>(CarbonAmount+2, 0));   //initialize a vector for this isomer code containing information on the connectivity between carbon atoms in this isomer

  morgans_splicing(isomers, isomer, CarbonAmount, morgans, connections);              //split the isomer code into morgan's code and determine the connectivity between carbon atoms
  morgans_algorithm(isomer, CarbonAmount, morgans, connections);                      //use the morgan's algorithm for CarbonAmount/2 iterations to generate canonical AND comparable isomer code
  morgans_sort(isomer, CarbonAmount, morgans);                                        //sort the values of the morgan's code by value, to make the comparison easier


  return check_morgan_uniqueness(isomer, CarbonAmount, morgans);                      //check if the generated morgan's code is unique and thus represents a valid new isomer
}



//morgans_splicing is responsible to create the connectivity table based on the canonical isomer representation of generate_isomers. In addition, it copies the initial values for the Morgan's values
//into the morgans vector.
//The canonical representation of generate_isomers contains connectivity information. However, to facilitate easy generation of isomers, the representation has to be analyzed to obtain the information.
//This is done by a simple concept: Each branch of the isomer ends with a 0. Every atom with a 1 has a connection to the next atom in the code, and the previous one. For >1, the idea of "foreign connections" 
//is used. If for example the code [...2100...] represents a split branch with a methyl and ethyl-subsituent, the carbon atom represented by 2 is connected to the 1 and the second 0, as the first 0 is only
//connected to the 1 and belongs to the end of the ethyl-substituent. So, the first 0 is a "foreign connection", an atom that is not connected to the carbon atom we are looking at, but is found in the code
//before the second (or third, ...) connection. So while there are still foreign connections unfilled, the atom we are looking at does not have a connection to the current atom. This is true for all values
//>1.
//This implementation first creates the initial Morgan's code by incrementing all values by 1 but the root, and copying it into the morgans vector. This code is later used for Morgan's Algorithm.
//Then, the amount of connections of the current atom is determined as its value in the canonical representation. If it has at least one connection, it is connected to the next atom in the code.
//While there are still connections open, the algorithm will add the value of the current atom as the amount of foreign connections (ie. a value of 1 means 1 additional foreign connection, see example).
//If there are no foreign connections open, the current atom is connected to the one we are examining. This means, that for both atoms, the connectivity to the other atom is saved in the connections table.
//In addition, the amount of connections of each atom in this bond saved at connection[digit][0] is incremented.
//If there are foreign connections open, the current atom's value is added to the amount of foreign connections and the amount of foreign connections is decreased by 1.
//In the end, the connections vector contains the amount of bonds of each atom at n=connections[digit][0] and the atoms it is connected to at connections[digit][1..n].
//As the algorithm only checks for forward connectivity, only connections to atoms that appear later in the code are recognized. This is why the canonical representation of generate_isomers is important, as
//it decreases each digit but the root by 1 and thus accounts for the backward connectivity of each atom (the root is not backwards connected at it is the first digit in the representation).
//Because both atoms recieve the connection information if a forward bond is found, the backwards connection is still present in the connections table for the other atom. 

void morgans_splicing(vector<vector<vector<int> > > &isomers, int isomer, int CarbonAmount, vector<vector<int> > &morgans, vector<vector<int> > &connections) {

  for(int main_digit=1; main_digit<=CarbonAmount; main_digit++) {                           //cycle through all digits of the current canonical isomer code

    int connections_n=isomers[CarbonAmount][isomer][main_digit], foreign_connections_n=0;   //the amount of (forward) connections is the same as the value in the isomer code

    if(main_digit==1) {                                                                     //TRANSLATION INTO MORGAN'S CODE OF 0th ITERATION
      morgans[isomer][main_digit]=isomers[CarbonAmount][isomer][main_digit];                //If digit is root, copy same value into morgans vector
    } else {
      morgans[isomer][main_digit]=isomers[CarbonAmount][isomer][main_digit]+1;              //if digit is not root, copy the value+1 into the morgans vector, as there is a backwards bond
    }

    if(!(connections_n==0)) {                                                               //if there is at least one forward connection (code=[1,2,3,4])
      connections[main_digit][0]++;                                                         //active atom has one additional connection in connections table
      connections[main_digit][connections[main_digit][0]]=main_digit+1;                     //this connections is to the next atom in the code
      connections[main_digit+1][0]++;                                                       //the next atom in the code has a backwards connection to the active atom, so increment connections amount
      connections[main_digit+1][connections[main_digit+1][0]]=main_digit;                   //the connection is to the current atom
      connections_n--;                                                                      //one connection was found, decrease remaining connections by 1
      foreign_connections_n+=isomers[CarbonAmount][isomer][main_digit+1];                   //the value of the next atom represents the amount of foreign connections 
    }

    for(int digit=main_digit+2; digit<=CarbonAmount && connections_n>0; digit++) {          //for every digit in the code following active atom+2 and while not all connections have been found:

      if(foreign_connections_n==0) {                                                        //if there are no open foreign connections, the digit atom is connected to the main_digit atom
	connections[main_digit][0]++;                                                       //increment both the connections amount of the active_digit and digit
	connections[main_digit][connections[main_digit][0]]=digit;                          //and insert the connection to the other atom in the connections table
        connections[digit][0]++;
        connections[digit][connections[digit][0]]=main_digit;
        connections_n--;                                                                    //one additional connection was found, decrease remanining connections by 1
        foreign_connections_n+=isomers[CarbonAmount][isomer][digit];                        //the amount of foreign connections is the amount of forwards connections of digit atom 
 
      } else {                                                                              //if there are foreign connections open:
        foreign_connections_n+=isomers[CarbonAmount][isomer][digit];                        //the amount of remaining foreign connections is increased by the value of digit atom
	foreign_connections_n--;                                                            //as one foreign connection was found and ignored, decrement the amount of foreign connections
      }
    } 
  }
}



//morgans_algorithm will use Morgan's Algorithm to generate a canonical value for each of the isomer's carbon atoms that contains information of the isomers complete structure. Instead of terminating when no 
//change in the amount of different values is found, as Morgan's Algorithm suggests, the algorithm will be applied a total of (CarbonAtom/2)+1 times. This way, the generated values can be used without 
//translation into a ranking. With (CarbonAtom/3)+1 iterations, every value has information of the connectivity of at least 50% of the entire structure, which is more than enough to ensure that the maximum 
//diversity for the Morgan's values has been found.
//This implementation cycles a total of (CarbonAmount/3)+1 times. Each time, the previous representation is copied to a temporary array and each carbon atom is assigned the sum of the values of each carbon
//atom it is connected to. This is done using the connectivity table generated by morgans_splicing. The new value is stored in the morgans vector.

void morgans_algorithm(int isomer, int CarbonAmount, vector<vector<int> > &morgans, vector<vector<int> > &connections){

  int iteration=0;                                                                              //initialize the iteration counter as 0
  int temp[CarbonAmount+1];                                                                     //initialize the temporary array to store previous values

  while(iteration<=(CarbonAmount/3)){                                                           //while there were less than (CarbonAmount/2)+1 iterations

    for(int digit=1; digit<=CarbonAmount; digit++){                                             //copy all previous values to temporary array
      temp[digit]=morgans[isomer][digit];
    }

    for(int digit=1; digit<=CarbonAmount; digit++){                                             //for every carbon atom:
      int connections_n=connections[digit][0];                                                  //check the amount of connections it has
      morgans[isomer][digit]=0;                                                                 //reset the current carbon atoms Morgan's value to 0

      for(int current_connection=1; current_connection<=connections_n; current_connection++){   //for every connection this carbon atom has
	morgans[isomer][digit]+=temp[connections[digit][current_connection]];                   //add the value of the connected carbon atom to the current atoms Morgan's value
      }
    }
    iteration++;                                                                                //increment the amount of iterations done
  }
}



//morgans_sort uses insertion sort to reorder the values stored for the current isomers morgan's code in the morgans vector. morgans[isomer][0] acts as sentinel. After this function is called, the code at
//morgans[isomer][1...CarbonAmount+1] is sorted by value, beginning with the highest value

void morgans_sort(int isomer, int CarbonAmount, vector<vector<int> > &morgans){

  int main, compare;                                               //two positions on the code are initialized, the digit that moves, main, and the digits that main is compared to, compare

  for(main=2; main<=CarbonAmount; main++) {                        //cycle through every digit of the code, beginning at 2 (first digit is "sorted" as it is the only digit up to current main)
    morgans[isomer][0]=morgans[isomer][main];                      //copy current value to position 0 to act as sentinel
    compare=main;                                                  //initialize compare digit to be at the current main's position

    while(morgans[isomer][0]>morgans[isomer][compare-1]){          //while the values in front of main are less than the main's value
      morgans[isomer][compare]=morgans[isomer][compare-1];         //transfer the value in front to current position
      compare--;                                                   //check for the next (previous) digit
    }
    morgans[isomer][compare]=morgans[isomer][0];                   //if the previous digit is greater, insert digit at current position
  }
}



//check_morgan_uniqueness will go through the list of all Morgan's codes of the current alkane to search for an identical code. The implementation uses the Morgan's code it is checking for as a sentinel
//at the end of the morgans vector to ensure termination. The implementation stops at the first occurence of the code.

int check_morgan_uniqueness(int isomer, int CarbonAmount, vector<vector<int> > &morgans){
  int twin_found=0;                                                                                           //twin_found is the exit condition for the while loop, it is triggered if the same isomer is found
  int current_isomer=0;                                                                                       //current_isomer keeps track of the isomer that is compared to

  while(!twin_found){                                                                                         //as long as there is no identical isomer already in the vector found, keep searching
    int current_digit=1;                                                                                      //current_digit keeps track of the digit in the current isomer that is compared to
    current_isomer++;                                                                                         //if previous isomer was not identical, the next one will be checked

    while((current_digit <= CarbonAmount) &&                                                                  //check if digit in both isomers are identical and the isomer end has not been reached
          (morgans[isomer][current_digit]==morgans[current_isomer][current_digit])){                          //(shortcut property of && used to prevent possible segmentation faults)
      current_digit++;                                                                                        //if the digits are identical, check for the next digit
    }

    if(current_digit>CarbonAmount){                                                                           //if the while loop reaches the end the isomer code, the isomers are identical
      twin_found=1;                                                                                           //exit condition is set to true, escaping the loop over every isomer
    }                                                                                                         //because there is a sentinel at the end of the vector, the wile loop always terminates
  }

  if(current_isomer==isomer){                                                                                 //check position of the identical code found in the vector
    return 1;                                                                                                 //if the only identical isomer found is the sentinel, the isomer itself, the isomer is unique
  } else {                                                                                                    //otherwise an identical isomer has already been generated and this one is to be ignored
    return 0;
  }
}



//UI FUNCTION GROUP

//print_intro will print the introduction and ask whether output files should be generated

bool print_intro(){
  char answer;
  bool boolean;

  cout << endl;
  cout << "Generation and enumeration of all alkane constitutional isomers up to icosane   " << endl;
  cout << "================================================================================" << endl;
  cout << "by Andreas Gimpel, agimpel@student.ethz.ch                                      " << endl;
  cout << "an entry to Prof. Philippe H. Hünenberger's challenge HS15                      " << endl;
  cout << endl;
  cout << "This program will enumerate the constitutional isomers from methane to icosane  " << endl;
  cout << "using a canonical representation of isomers and a modified Morgan's algorithm.  " << endl;
  cout << "Depending on the system, computation time may exceed 1h. Termination is possible" << endl;
  cout << "with CTRL-C. Documentation available in the attached .pdf file and the code.    " << endl;
  cout << endl;
  cout << "Should all valid isomer codes be output to files to allow for reconstruction    " << endl;
  cout << "(This will require the directory '/isomer' to be present)?  [y/n]               " << endl;  
  cin  >> answer;

  if(answer=='y') {
    boolean=1;
    cout << "Output will be generated in /isomer as [Carbon Atoms in alkane].isomers. " << endl;
  } else if(answer=='n') {
    boolean=0;
    cout << "No output files will be generated." << endl;
  } else {
    cout << "No valid input. Defaulting to no output." << endl;
    boolean=0;
  }
  return boolean;
}



//print_structure simply prints table captions and the predefined isomer code for methane

void print_structure(vector<vector<vector<int> > > &isomers){
  cout << endl;
  cout << "n \t" << "#isomers" << endl;
  cout << "____________________________________" << endl;
  cout << "1 \t" << isomers[1][0][0] << endl;
}



//print_Isomers will ouput the amount of isomers found after each main alkane cycle as well as create the output file, if it was enabled

void print_Isomers(vector<vector<vector<int> > > &isomers, int CarbonAmount, string filenames[], vector<vector<int> > &morgans, bool generate_files){

  cout << CarbonAmount << " \t" << isomers[CarbonAmount][0][0] << endl;

  if(isomers[1][1][0]==0 && generate_files) {
    ofstream file(filenames[1].c_str());
    file << "# Carbon atoms in this alkane: " << 1 << endl << "# Amount of isomers found for this alkane: " << isomers[1][0][0] << endl;
    for(int z=1;z<=isomers[1][0][0]; z++){
      for(int y=1; y<=1; y++){ 
	file << isomers[1][z][y];
      }
      file << endl;
    }
    file.close();
    isomers[1][1][0]=1;
  }

  if(generate_files) {
    ofstream file(filenames[CarbonAmount].c_str());
    file << "# Carbon atoms in this alkane: " << CarbonAmount << endl << "# Amount of isomers found for this alkane: " << isomers[CarbonAmount][0][0] << endl;
    for(int z=1;z<=isomers[CarbonAmount][0][0]; z++){
      for(int y=1; y<=CarbonAmount; y++){ 
	file << isomers[CarbonAmount][z][y];
      }
      file << endl;
    }
    file.close();
  }
}
