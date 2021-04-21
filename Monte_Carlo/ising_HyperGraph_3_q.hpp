
#ifndef SRC_ising_HyperGraph_3_q_HPP_
#define SRC_ising_HyperGraph_3_q_HPP_

#include "HyperGraph_3_q.hpp"
#include <math.h>
#include <string>
#include <numeric>



class ising_HyperGraph_3_q{
public:


	double k, q;
	double L,N;  // size, N=LXL
	HyperGraph_3_q G; // Hypergraph 3+q
	vector<int> spins;  // The state of the spins in each network. [-1,1].
	vector<double> local_m;

	double M;
	double E;

	double giant;

	ising_HyperGraph_3_q(double L_, double k_, double q);
	void initialize_spins_disordered();
	void initialize_spins_ordered_up();
	void initialize_spins_ordered_down();
	void initializing_local_m();
	void flip_spin(int spin);

	void revive_ordered_up();
	void revive_ordered_down();
	void revive_disordered();

	void scan_T(double T_i,double T_f,double dT, string where_to, int TestNum, double MCS, double MCF);
};

ising_HyperGraph_3_q::ising_HyperGraph_3_q(double L_, double k_, double q_):L(L_),N(L_*L_),k(k_),q(q_){


	cout<<"start constructor"<<endl;
	G.setHyperGraph_3_q(N);
	G.createHyperGraphER_3_q(k,q) ;
	giant = G.testConnectivity(); 
				    
	//initialize_spins_disordered();
	initialize_spins_ordered_up();
	initializing_local_m();



	M = 0;
	for(int i=0;i<N;i++){
		if(G.connected[i])
			M+=spins[i];
	}


/*
	E = 0;
	
	for(int i = 0; i < N; i++){
		if(G.connected[i]){
			for(int j = 0; j < G.v[i].size();j++){
				E += -spins[i]*spins[G.v[i][j].first]*spins[G.v[i][j].second];
			}
		}
	}
	*/
	
	cout<<"finish constructor"<<endl;

}

void ising_HyperGraph_3_q::initialize_spins_disordered(){
	double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

	for(int i=0;i<N;i++){ //network1
		rand_num = (double)G.gen()/G.gen.max();
		if(rand_num < 0.5){
			spins.push_back(1);
		}
		else{
			spins.push_back(-1);
		}
	}
}

void ising_HyperGraph_3_q::initialize_spins_ordered_up(){
	spins = vector<int>(N,1); //all the spins are up in both layers
}

void ising_HyperGraph_3_q::initialize_spins_ordered_down(){
	spins = vector<int>(N,-1); //all the spins are down in both layers
}

void ising_HyperGraph_3_q::initializing_local_m(){

	for(int i=0;i<N;i++){
		double sum = 0;
		for(int j=0;j<G.v[i].size();j++){
			if(G.v[i][j].size() == 2)
				sum+=spins[G.v[i][j][0]]*spins[G.v[i][j][1]];
			else
				sum+=spins[G.v[i][j][0]]*spins[G.v[i][j][1]]*spins[G.v[i][j][2]];
		}
		if(G.connected[i])
			local_m.push_back((double)sum);
		else{
			local_m.push_back(0);
		}
	}
}

void ising_HyperGraph_3_q::flip_spin(int spin){

	if(G.connected[spin]){
		spins[spin]*=-1;

		for(int i=0;i<G.v[spin].size();i++){
			
			//triangle
			if(G.v[spin][i].size() == 2){
				int j = G.v[spin][i][0], u = G.v[spin][i][1];
				local_m[j]+= (double)2*spins[spin]*spins[u];
				local_m[u]+= (double)2*spins[spin]*spins[j];

				//E += -6*spins[spin]*spins[j]*spins[u];
			}
			//square
			else{
				int j = G.v[spin][i][0], u = G.v[spin][i][1], w = G.v[spin][i][2] ;
				local_m[j]+= (double)2*spins[spin]*spins[u]*spins[w];
				local_m[u]+= (double)2*spins[spin]*spins[j]*spins[w];
				local_m[w]+= (double)2*spins[spin]*spins[j]*spins[u];
				//E += -4*spins[spin]*spins[j]*spins[u];
			}
		}
		
		M+=2*spins[spin];
	}


}
void ising_HyperGraph_3_q::scan_T(double T_i,double T_f,double dT, string where_to, int TestNum, double MCS, double MCF){

	vector<double> T_vec, M_vec, E_vec;
	double epsilon = 0.1;

	double rand_num;
	int u;
	double pi;
	for(double T = T_i;T<=T_f;T+=dT){
		vector<double> E_dist;
		T_vec.push_back(T);
		double beta = 1/T;
		for(int j=0;j<MCS;j++){

			for(int i=0;i<MCF;i++){
				u = G.gen()%int(N);
				pi = (0.5)*(1-spins[u]*tanh(beta*local_m[u]));
				
				rand_num = (double)G.gen()/G.gen.max();
				if(rand_num < pi){
					flip_spin(u);
				}
			}
		}

		M_vec.push_back(M/N);
		E_vec.push_back(E/(N));
		cout<<"M = "<<M/giant<<",    T="<<T<<", L = "<<L<<", Tf - T = "<<T_f- T<<", MCS = "<<MCS<<endl;
		
	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = where_to + "/spins";
	ss<<TestNum;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<T_vec.size();i++)
	{
		dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<endl;
	}

	dataInfo.close();
}

void ising_HyperGraph_3_q::revive_ordered_up(){
	for(int i=0;i<N;i++){   
		if(spins[i] == -1)
			flip_spin(i);
	}
}
void ising_HyperGraph_3_q::revive_ordered_down(){
	for(int i=0;i<N;i++){    
		if(spins[i] == 1)
			flip_spin(i);
	}
}

void ising_HyperGraph_3_q::revive_disordered(){
	double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

	for(int i=0;i<N;i++){ //network
		rand_num = (double)G.gen()/G.gen.max();
		if(rand_num < 0.5){
			flip_spin(i);
		}
	}
}


#endif /* SRC_ising_HyperGraph_3_q_HPP_ */
