#ifndef SRC_ising_interdependent_2_layers_HPP_
#define SRC_ising_interdependent_2_layers_HPP_


#include "Graph.hpp"
#include <math.h>
#include <string>
#include <numeric>



class ising_interdependent_2_layers{
public:

	double L,N;  // size, N=LXL
	double q;    // fraction of interdependent nodes
	double k;    // Network average degree

	Graph G1, G2; // The networks
	vector<bool> inter_connected; //vector which states if the node is interdependent or not
	vector<int> spins1, spins2;  // The state of the spins in each network. [-1,1].
	vector<double> local_m1, local_m2; // The local magnetization of each node

	double M1, M2; //Global magnetization of each layer
	double E1, E2, E; // THe energy

	double giant1, giant2; // The giant component of each layer

	ising_interdependent_2_layers(double L_, int k_, double q_); //Constructor

	void create_interlinks();
	void calc_energy();
	void initialize_spins_disordered();
	void initialize_spins_ordered_up();
	void initialize_spins_ordered_down();
	void initializing_local_m();
	void flip_spin(int spin, int net_idx);
	void revive_ordered_up();
	void revive_ordered_down();
	void revive_disordered();
	
	void scan_T_thermal(double T_i,double T_f,double dT, string where_to, int TestNum, double NOI, double NOF);
	void scan_T_high_order(double T_i,double T_f,double dT, string where_to, int TestNum, double NOI, double NOF);

};

ising_interdependent_2_layers::ising_interdependent_2_layers(double L_, int k_, double q_):L(L_),N(L_*L_),k(k_),q(q_){

	G1.setGraph(N); // set network 1 size to N
	G1.createGraphER(k); // Create ER network
	giant1 = G1.testConnectivity(); 	
	
	for(double ii = 0; ii<3000000000; ii++) //pause  before creating another network
		int dd = 5;
	
	G2.setGraph(N); //same as G1
	G2.createGraphER(k);
	giant2 = G2.testConnectivity();


	create_interlinks();
	
	initialize_spins_ordered_up();
	//initialize_spins_disordered();
	
	initializing_local_m();



	M1 = 0;
	M2 = 0;
	for(int i=0;i<N;i++){
		if(G1.active[i]){
			M1+=spins1[i];
			M2+=spins2[i];
		}
	}

	cout<<"finish constructor"<<endl;

}

void ising_interdependent_2_layers::create_interlinks(){
	inter_connected = vector<bool>(N, 0); //create vector of size N full with zeroes

	double rand_num;
	for(int i=0;i<N;i++){ //for each node
		rand_num = (double)G1.gen()/G1.gen.max(); //random number between 0 to 1
		//if the random number is smaller then q and both nodes are connected, flag the node as interconnected
		//else, leave it as zero
		//if(rand_num < q and G1.connected[i] and G2.connected[i])
		if(rand_num < q)
			inter_connected[i] = 1;
	}
}

void ising_interdependent_2_layers::calc_energy(){
	
	E1 = 0;
	E2 = 0;
	E = 0;
	
	for(int i = 0; i < N; i++){
		for(int j = 0 ;j < G1.v[i].size();j++)
			E1 -= spins1[i]*spins1[G1.v[i][j]];
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0 ;j < G2.v[i].size();j++)
			E2 -= spins2[i]*spins2[G2.v[i][j]];
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0 ;j < G1.v[i].size();j++)
			E -= spins1[i]*spins1[G1.v[i][j]]*spins2[i];
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0 ;j < G2.v[i].size();j++)
			E -= spins2[i]*spins2[G1.v[i][j]]*spins1[i];
	}

}
void ising_interdependent_2_layers::initialize_spins_disordered(){
	double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

	for(int i=0;i<N;i++){ //network1
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			spins1.push_back(1);
		}
		else{
			spins1.push_back(-1);
		}
	}

	for(int i=0;i<N;i++){ //network2
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			spins2.push_back(1);
		}
		else{
			spins2.push_back(-1);
		}
	}
}

void ising_interdependent_2_layers::initialize_spins_ordered_up(){
	spins1 = vector<int>(N,1); //all the spins are up in both layers
	spins2 = vector<int>(N,1);

}

void ising_interdependent_2_layers::initialize_spins_ordered_down(){
	spins1 = vector<int>(N,-1); //all the spins are down in both layers
	spins2 = vector<int>(N,-1);

}

void ising_interdependent_2_layers::initializing_local_m(){

	for(int i=0;i<N;i++){
		if(!G1.connected[i]){
			local_m1.push_back(0);
			continue;
		}
		
		double sum = 0, count = 0;
		for(int j=0;j<G1.v[i].size();j++){
			if(G1.connected[G1.v[i][j]]){
				count++;
				sum+=spins1[G1.v[i][j]];
			}
		}

		local_m1.push_back((double)sum/count);

	}

	for(int i=0;i<N;i++){
		if(!G2.connected[i]){
			local_m2.push_back(0);
			continue;
		}
		double sum = 0, count = 0;;
		for(int j=0;j<G2.v[i].size();j++){
			if(G2.connected[G2.v[i][j]]){
				sum+=spins2[G2.v[i][j]];
				count++;
			}
		}

		local_m2.push_back((double)sum/count);

	}
}

void ising_interdependent_2_layers::flip_spin(int spin, int net_idx){

	if(net_idx == 1 and G1.connected[spin]){
		spins1[spin]*=-1;

		for(int i=0;i<G1.v[spin].size();i++){
			int j = G1.v[spin][i];
			if(!G1.connected[j])
				continue;
			
			double count = 0;
			for(int jj = 0; jj < G1.v[j].size(); jj++)
				if(G1.connected[G1.v[j][jj]])
					count++;
			
			local_m1[j]+= (double)2*spins1[spin]/count;

			//E1 += -4*spins1[spin]*spins1[j];

		}
		M1+=2*spins1[spin];
	}
	else if(net_idx == 2 and G2.connected[spin]){
		spins2[spin]*=-1;

		for(int i=0;i<G2.v[spin].size();i++){
			int j = G2.v[spin][i];
			if(!G2.connected[j])
				continue;
			
			double count = 0;
			for(int jj = 0; jj < G2.v[j].size(); jj++)
				if(G2.connected[G2.v[j][jj]])
					count++;
				
			local_m2[j]+= (double)2*spins2[spin]/count;
			//E2 += -4*spins2[spin]*spins2[j];
		}

		M2+=2*spins2[spin];
	}

}

void ising_interdependent_2_layers::scan_T_thermal(double T_i,double T_f,double dT, string where_to, int TestNum, double NOI, double NOF){

	vector<double> T_vec, M_vec, E_vec, E1_vec, E2_vec;


	double rand_num;
	int u;
	double pi;
	for(double T = T_i;T<=T_f;T+=dT){
		cout<<"T = "<<T<<", M = "<<M1/N<<endl;
		T_vec.push_back(T);
		double beta = 1/T;
		for(int j=0;j<NOI;j++){

			for(int i=0;i<NOF;i++){
				u = G1.gen()%int(N);
				if(!G1.connected[u])
					continue;
				
				if(!inter_connected[u])
					pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));
				else
					pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
				// To reduce fluctuations use (M2/N) instead of local_m2[u]
				
				

					
				rand_num = (double)G1.gen()/G1.gen.max();
				if(rand_num < pi){
					flip_spin(u,1);
				}
			}
		
			for(int i=0;i<NOF;i++){
					u = G2.gen()%int(N);
					if(!G2.connected[u])
						continue;
					
					if(!inter_connected[u])
						pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m1[u]*G2.v[u].size()));
					// To reduce fluctuations use (M1/N) instead of local_m1[u]

					rand_num = (double)G2.gen()/G2.gen.max();
					if(rand_num < pi){
						flip_spin(u,2);
					}
				}
		}
		
	
		calc_energy();
		M_vec.push_back(M1/N);
		E_vec.push_back(E/N);
		E1_vec.push_back(E1/N);
		E2_vec.push_back(E2/N);
		cout<<"M = "<<M1/N<<", E1 = "<<E1/N<<", E2 = "<<E2/N<<", E = "<<E/N<<",    T="<<T<<endl;

	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = where_to + "/giant";
	ss<<TestNum;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<T_vec.size();i++)
	{
		dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<"	"<<(double) E1_vec[i]<<"	"<<(double) E2_vec[i]<<"	"<<(double) E_vec[i]<<endl;
	}

	dataInfo.close();
	
	
}

void ising_interdependent_2_layers::scan_T_high_order(double T_i,double T_f,double dT, string where_to, int TestNum, double NOI, double NOF){

	vector<double> T_vec, M_vec, E_vec, E1_vec, E2_vec;


	double rand_num;
	int u;
	double pi;
	for(double T = T_i;T<=T_f;T+=dT){
		T_vec.push_back(T);
		double beta = 1/T;
		for(int j=0;j<NOI;j++){

			for(int i=0;i<2*NOF;i++){
				int ppi = G1.gen()%2;
				if(ppi == 0){
					u = G1.gen()%int(N);
					if(!G1.connected[u])
						continue;
					
					if(!inter_connected[u])
						pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
					// To reduce fluctuations use (M2/N) instead of local_m2[u]
					//pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*G1.v[u].size()));
					

						
					rand_num = (double)G1.gen()/G1.gen.max();
					if(rand_num < pi){
						flip_spin(u,1);
					}
				}
			else{
				
					u = G2.gen()%int(N);
					if(!G2.connected[u])
						continue;
					
					if(!inter_connected[u])
						pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m1[u]*G2.v[u].size()));
					// To reduce fluctuations use (M1/N) instead of local_m1[u]
					// pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(M1/N)*G2.v[u].size()));

					rand_num = (double)G2.gen()/G2.gen.max();
					if(rand_num < pi){
						flip_spin(u,2);
					}
				}
			}
		}
	
		calc_energy();
		M_vec.push_back(M1/N);
		E_vec.push_back(E/N);
		E1_vec.push_back(E1/N);
		E2_vec.push_back(E2/N);
		cout<<"M = "<<M1/N<<", E1 = "<<E1/N<<", E2 = "<<E2/N<<", E = "<<E/N<<",    T="<<T<<endl;

	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = where_to + "/giant";
	ss<<TestNum;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<T_vec.size();i++)
	{
		dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<"	"<<(double) E1_vec[i]<<"	"<<(double) E2_vec[i]<<"	"<<(double) E_vec[i]<<endl;
	}

	dataInfo.close();
	
	
}

void ising_interdependent_2_layers::revive_ordered_up(){
	for(int i=0;i<N;i++){    //flip all the spins back up
		if(spins1[i] == -1)
			flip_spin(i,1);
		if(spins2[i] == -1)
			flip_spin(i,2);
	}
}

void ising_interdependent_2_layers::revive_ordered_down(){
	for(int i=0;i<N;i++){    //flip all the spins back up
		if(spins1[i] == 1)
			flip_spin(i,1);
		if(spins2[i] == 1)
			flip_spin(i,2);
	}
}

void ising_interdependent_2_layers::revive_disordered(){
	double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

	for(int i=0;i<N;i++){ //network1
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			flip_spin(i,1);
		}
	}

	for(int i=0;i<N;i++){ //network2
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			flip_spin(i,2);
		}
	}
}


#endif /* SRC_ising_interdependent_2_layers_HPP_ */
