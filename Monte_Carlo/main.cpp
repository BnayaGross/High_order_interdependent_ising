#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>
#include <string.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <bitset>
#include <stdexcept>
#include <exception>
#include <map>
#include <cassert>
#include <iomanip>
#include <limits>
#include <cmath>
#include <omp.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <set>
#include <unistd.h>
#include <sys/resource.h>



// Classes
#include "ising_HyperGraph_3_q.hpp"
#include "ising_HyperGraph_2_q.hpp"
#include "ising_interdependent_3_layers.hpp"
#include "ising_interdependent_2_layers.hpp"



using std::pair;
using std::make_pair;
using std::vector;
using std::logic_error;
using std::map;
using std::cout;
using std::cerr;
using std::setprecision;
using namespace std;

int main(int argc, char *argv[]){

	
	//*****************************************************
	// Annealing interdependent ising model with two layers -- thermal, high order
	
	double L = 100; // N = LXL
	double k = 4;
	double q = 1;
	double beta = 1;
	string where_to = "/home/bnaya/Desktop/ising code/finished/test";

	ising_interdependent_2_layers I(L, k, q);
	//I.scan_T_thermal(0.01,0.6*k,0.01, where_to, I.G1.gen(), 100,500*L*L);
	I.scan_T_high_order(0.01,0.6*k, 0.01, where_to, I.G1.gen(), 100,L*L);
	//*****************************************************
	
	
	
	/*
	//*****************************************************
	// Annealing interdependent ising model with three layers -- thermal, high order
	
	double L = 100; // N = LXL
	double k = 4;
	double q = 1;
	double beta = 1;
	string where_to = "/home/bnaya/Desktop/ising code/finished/test";

	ising_interdependent_3_layers I(L, k, q);
	//I.scan_T_thermal( 0.01,0.6*k,0.01, where_to, I.G1.gen(), 100,500*L*L);
	I.scan_T_high_order(0.01,0.6*k, 0.01, where_to, I.G1.gen(), 100,L*L);
	//*****************************************************
	*/
	
	
	/*
	//*****************************************************
	// Annealing  ising Hypergraph 2+q model
	
	double L = 100; // N = LXL
	double k = 4;
	double q = 1;
	double beta = 1;
	string where_to = "/home/bnaya/Desktop/ising code/finished/test";

	ising_HyperGraph_2_q I(L, k, q);
	I.scan_T(0.01,0.6*k, 0.01, where_to, I.G.gen(), 100,L*L);
	//*****************************************************
	*/
	
	
	/*
	//*****************************************************
	// Annealing  ising Hypergraph 3+q model
	
	double L = 100; // N = LXL
	double k = 4;
	double q = 1;
	double beta = 1;
	string where_to = "/home/bnaya/Desktop/ising code/finished/test";

	ising_HyperGraph_3_q I(L, k, q);
	I.scan_T(0.01,0.6*k, 0.01, where_to, I.G.gen(), 100,L*L);
	//*****************************************************
	*/
	
	
}
