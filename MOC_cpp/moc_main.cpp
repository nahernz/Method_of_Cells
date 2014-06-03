/* This is the main structure for the Method of Cells
   code. It takes in input files, and writes to output files
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <assert.h>
#include "struct.h"


using namespace std;

// Takes in an arbitrary number of input files:
int main(int argc, char *argv[]) {
// --------------------------------------------------------
	// Open first file and check to see if successful:
	// exit otherwise
	ifstream input1;
	input1.open(argv[1]);
	if(!input1.is_open()) {
		cerr << "Opening of file " << argv[1] << " failed" << endl;
		exit(1);
	}

	// Read in the variables here:



	input1.close();

// --------------------------------------------------------
	{ Quad ruc; }

	return 0;
}