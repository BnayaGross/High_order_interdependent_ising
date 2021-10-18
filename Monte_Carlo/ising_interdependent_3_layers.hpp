

#ifndef SRC_ising_interdependent_3_layers_HPP_
#define SRC_ising_interdependent_3_layers_HPP_


#include "Graph.hpp"
#include <math.h>
#include <string>
#include <numeric>



class ising_interdependent_3_layers{
public:


	double k; //Networks properties
	double L,N;  // size, N=LXL
	double q;    // fraction of interdependent nodes

	Graph G1, G2, G3; // The networks
	vector<bool> inter_connected1, inter_connected2, inter_connected3; //vector which states  the dependence of node to the other layers
	vector<int> spins1, spins2, spins3;  // The state of the spins in each network. [-1,1].
	vector<double> local_m1, local_m2, local_m3;
	
	vector<double> local_m1_thermal_avg, local_m2_thermal_avg, local_m3_thermal_avg; //Used for thermal averaging
	vector<double> local_m1_counter, local_m2_counter, local_m3_counter;

	double M1, M2, M3;
	double E1, E2, E3, E;

	double giant1, giant2, giant3;

	ising_interdependent_3_layers(double L_, int k_, double q_); //constructor
	void create_interlinks();
	void calc_energy();
	void initialize_spins_disordered();
	void initialize_spins_ordered_up();
	void initialize_spins_ordered_down();
	void initializing_local_m();
	void flip_spin(int spin, int net_idx, int step);
	void revive_ordered_up();
	void revive_ordered_down();
	void revive_disordered();

	void scan_T_high_order(double T_i,double T_f,double dT,string where_to, int TestNum, double NOI, double NOF);
	void scan_T_thermal(double T_i,double T_f,double dT,string where_to, int TestNum, double NOI, double NOF);
};

ising_interdependent_3_layers::ising_interdependent_3_layers(double L_, int k_, double q_):L(L_),N(L_*L_),k(k_),q(q_){

	G1.setGraph(N); // set network 1 size to N
	G1.createGraphER(k); // Create ER network
	giant1 = G1.testConnectivity(); 	
	
	for(double ii = 0; ii<3000000000; ii++) //pause  before creating another network
		int dd = 5;
	
	
	G2.setGraph(N); //same as G1
	G2.createGraphER(k);
	giant2 = G2.testConnectivity();
	
	for(double ii = 0; ii<3000000000; ii++) //pause  before creating another network
		int dd = 5;
	
	G3.setGraph(N); //same as G1
	G3.createGraphER(k);
	giant3 = G3.testConnectivity();


	create_interlinks();
	initialize_spins_ordered_up();
	//initialize_spins_disordered();

	initializing_local_m();



	M1 = 0;
	M2 = 0;
	M3 = 0;
	for(int i=0;i<N;i++){
		if(G1.active[i]){
			M1+=spins1[i];
			M2+=spins2[i];
			M3+=spins3[i];
		}
	}

	cout<<"finish constructor"<<endl;

}

void ising_interdependent_3_layers::create_interlinks(){
	inter_connected1 = vector<bool>(N, 0); //create vector of size N full with zeroes
	inter_connected3 = vector<bool>(N, 0); //create vector of size N full with zeroes
	inter_connected2 = vector<bool>(N, 0); //create vector of size N full with zeroes

	double rand_num;
	for(int i=0;i<N;i++){ //for each node
		
	
		rand_num = (double)G1.gen()/G1.gen.max(); 
		if(rand_num < q){
			inter_connected2[i] = 1;
			inter_connected3[i] = 1;
			inter_connected1[i] = 1;
		}
		
			
	}
}

void ising_interdependent_3_layers::calc_energy(){
	
	// Define better the total energy E
	E1 = 0;
	E2 = 0;
	E3 = 0;
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
		for(int j = 0 ;j < G3.v[i].size();j++)
			E3 -= spins3[i]*spins3[G3.v[i][j]];
	}
	
	/*for(int i = 0; i < N; i++){
		for(int j = 0 ;j < G1.v[i].size();j++)
			E -= spins1[i]*spins1[G1.v[i][j]]*spins2[i];
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0 ;j < G2.v[i].size();j++)
			E -= spins2[i]*spins2[G1.v[i][j]]*spins1[i];
	}*/

}
void ising_interdependent_3_layers::initialize_spins_disordered(){
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
	
	for(int i=0;i<N;i++){ //network3
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			spins3.push_back(1);
		}
		else{
			spins3.push_back(-1);
		}
	}
}

void ising_interdependent_3_layers::initialize_spins_ordered_up(){
	spins1 = vector<int>(N,1); //all the spins are up in both layers
	spins2 = vector<int>(N,1);
	spins3 = vector<int>(N,1);

}

void ising_interdependent_3_layers::initialize_spins_ordered_down(){
	spins1 = vector<int>(N,-1); //all the spins are down in both layers
	spins2 = vector<int>(N,-1);
	spins3 = vector<int>(N,-1);

}


void ising_interdependent_3_layers::initializing_local_m(){

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
	
	for(int i=0;i<N;i++){
		if(!G3.connected[i]){
			local_m3.push_back(0);
			continue;
		}
		double sum = 0, count = 0;;
		for(int j=0;j<G3.v[i].size();j++){
			if(G3.connected[G3.v[i][j]]){
				sum+=spins3[G3.v[i][j]];
				count++;
			}
		}

		local_m3.push_back((double)sum/count);

	}
}


void ising_interdependent_3_layers::flip_spin(int spin, int net_idx, int step){

	if(net_idx == 1){
		spins1[spin]*=-1;

		for(int i=0;i<G1.v[spin].size();i++){
			int j = G1.v[spin][i];
			
			double count = 0;
			for(int jj = 0; jj < G1.v[j].size(); jj++)
					count++;
			
			local_m1[j]+= (double)2*spins1[spin]/count;
			
			if(step > 10000000 and step%10000 == 0){
				local_m1_thermal_avg[j] += local_m1[j];
				local_m1_counter[j]++;
			}

			//E1 += -4*spins1[spin]*spins1[j];

		}
		M1+=2*spins1[spin];
	}
	else if(net_idx == 2){
		spins2[spin]*=-1;

		for(int i=0;i<G2.v[spin].size();i++){
			int j = G2.v[spin][i];
			
			double count = 0;
			for(int jj = 0; jj < G2.v[j].size(); jj++)
					count++;
				
			local_m2[j]+= (double)2*spins2[spin]/count;
			
			if(step > 10000000 and step%10000 == 0){
				local_m2_thermal_avg[j] += local_m2[j];
				local_m2_counter[j]++;
			}
			//E2 += -4*spins2[spin]*spins2[j];
		}

		M2+=2*spins2[spin];
	}	
	else if(net_idx == 3){
		spins3[spin]*=-1;

		for(int i=0;i<G3.v[spin].size();i++){
			int j = G3.v[spin][i];
			
			double count = 0;
			for(int jj = 0; jj < G3.v[j].size(); jj++)
					count++;
				
			local_m3[j]+= (double)2*spins3[spin]/count;
			
			if(step > 10000000 and step%10000 == 0){
				local_m3_thermal_avg[j] += local_m3[j];
				local_m3_counter[j]++;
			}
			//E2 += -4*spins2[spin]*spins2[j];
		}

		M3+=2*spins3[spin];
	}

}

void ising_interdependent_3_layers::scan_T_high_order(double T_i,double T_f,double dT,string where_to, int TestNum, double NOI, double NOF){

	vector<double> T_vec, M_vec, E_vec, E1_vec, E2_vec, E3_vec;

	double rand_num;
	int u;
	double pi;
	for(double T = T_i;T<=T_f;T+=dT){
		cout<<"T = "<<T<<", M = "<<M1/N<<endl;
		T_vec.push_back(T);
		double beta = 1/T;
		for(int j=0;j<NOI;j++){
			//cout<<"NOI = "<<j<<"  ,M1="<<M1/giant1<<"  ,M2="<<M2/giant2<<",  beta="<<beta<<endl;
			//change the local magnetization to global magnetization to reduce fluctuations
			for(int i=0;i<NOF;i++){
				
				int ppi = G1.gen()%3;
				if(ppi == 0){
					u = G1.gen()%int(N);

					if(inter_connected1[u])
						pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*local_m3[u]*G1.v[u].size()));
						//pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*(M3/N)*G1.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
						//pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*G1.v[u].size()));
					
					rand_num = (double)G1.gen()/G1.gen.max();
					if(rand_num < pi){
						flip_spin(u,1,i);
					}
				}
				else if (ppi == 1){
					u = G2.gen()%int(N);

					if(inter_connected2[u])
						pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
						//pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*(M3/N)*G2.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
						//pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(M3/N)*G2.v[u].size()));
					
					rand_num = (double)G2.gen()/G2.gen.max();
					if(rand_num < pi){
						flip_spin(u,2,i);
					}
				}
				else if (ppi == 2){
					u = G3.gen()%int(N);

					if(inter_connected3[u])
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
						//pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*(M2/N)*local_m3[u]*G3.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*local_m1[u]*G3.v[u].size()));
						//pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(M1/N)*G3.v[u].size()));
					
					rand_num = (double)G3.gen()/G3.gen.max();
					if(rand_num < pi){
						flip_spin(u,3,i);
					}
				}
				
				
			}
		}
		
		calc_energy();
		M_vec.push_back(M1/N);
		E_vec.push_back(E/N);
		E1_vec.push_back(E1/N);
		E2_vec.push_back(E2/N);
		E3_vec.push_back(E3/N);
		cout<<"M = "<<M1/N<<", E1 = "<<E1/N<<", E2 = "<<E2/N<<", E3 = "<<E3/N<<", E = "<<E/N<<",    T="<<T<<endl;


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
			dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<"	"<<(double) E1_vec[i]<<"	"<<(double) E2_vec[i]<<"	"<<(double) E3_vec[i]<<"	"<<(double) E_vec[i]<<endl;
	}

	dataInfo.close();
}

void ising_interdependent_3_layers::scan_T_thermal(double T_i,double T_f,double dT,string where_to, int TestNum, double NOI, double NOF){

	vector<double> T_vec, M_vec, E_vec, E1_vec, E2_vec, E3_vec;

	double rand_num;
	int u;
	double pi;
	for(double T = T_i;T<=T_f;T+=dT){
		cout<<"T = "<<T<<", M = "<<M1/N<<endl;
		T_vec.push_back(T);
		double beta = 1/T;
		for(int j=0;j<NOI;j++){
			//cout<<"NOI = "<<j<<"  ,M1="<<M1/giant1<<"  ,M2="<<M2/giant2<<",  beta="<<beta<<endl;
			
			//change the local magnetization to global magnetization to reduce fluctuations
			
			//local magnetization
			for(int i=0;i<NOF;i++){
		
					u = G1.gen()%int(N);

					if(inter_connected1[u]){
						if(local_m2_counter[u] != 0 and local_m3_counter[u] != 0)
							pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
						else if(local_m2_counter[u] != 0 and local_m3_counter[u] == 0)
							pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*local_m3[u]*G1.v[u].size()));
						else if(local_m2_counter[u] == 0 and local_m3_counter[u] != 0)
							pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
						else
							pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*local_m3[u]*G1.v[u].size()));
					}
					else{
						if(local_m2_counter[u] != 0)
							pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G1.v[u].size()));
						else
							pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
					}

					
					rand_num = (double)G1.gen()/G1.gen.max();
					if(rand_num < pi){
						flip_spin(u,1,i);
					}
					
			}
			
			for(int i = 0;i < N; i++){
				local_m2_thermal_avg[i] = 0;
				local_m2_counter[i] = 0;
			}
			
			for(int i=0;i<NOF;i++){
				
					u = G2.gen()%int(N);

					if(inter_connected2[u]){
						if(local_m1_counter[u] != 0 and local_m3_counter[u] != 0)
							pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
						else if(local_m1_counter[u] != 0 and local_m3_counter[u] == 0)
							pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*local_m3[u]*G2.v[u].size()));
						else if(local_m1_counter[u] == 0 and local_m3_counter[u] != 0)
							pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
						else
							pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
					}
					else{
						if(local_m3_counter[u] != 0)
							pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
						else
							pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
					}

					
					
					rand_num = (double)G2.gen()/G2.gen.max();
					if(rand_num < pi){
						flip_spin(u,2,i);
					}
				}
				
			for(int i = 0;i < N; i++){
				local_m3_thermal_avg[i] = 0;
				local_m3_counter[i] = 0;
			}
			
			for(int i=0;i<NOF;i++){
				
				u = G3.gen()%int(N);

				if(inter_connected3[u]){
					if(local_m1_counter[u] != 0 and local_m2_counter[u] != 0)
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
					else if(local_m1_counter[u] != 0 and local_m2_counter[u] == 0)
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*local_m2[u]*G3.v[u].size()));
					else if(local_m1_counter[u] == 0 and local_m2_counter[u] != 0)
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
				}
				else{
					if(local_m1_counter[u] != 0)
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*G3.v[u].size()));
					else
						pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*G3.v[u].size()));
				}

				
				rand_num = (double)G3.gen()/G3.gen.max();
				if(rand_num < pi){
					flip_spin(u,3,i);
				}
			}
			
			for(int i = 0;i < N; i++){
				local_m1_thermal_avg[i] = 0;
				local_m1_counter[i] = 0;
			}
				
				
			//global magnetization
			/*for(int i=0;i<NOF;i++){
	
				u = G1.gen()%int(N);
				
				if(inter_connected1[u])
					pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*(M3/N)*G1.v[u].size()));
				else
					pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*G1.v[u].size()));
		
				
				rand_num = (double)G1.gen()/G1.gen.max();
				if(rand_num < pi){
					flip_spin(u,1,i);
				}
					
			}
			for(int i=0;i<NOF;i++){
				
				u = G2.gen()%int(N);
				
				if(inter_connected2[u])
					pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*(M3/N)*G2.v[u].size()));
				else
					pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(M3/N)*G2.v[u].size()));

				
				rand_num = (double)G2.gen()/G2.gen.max();
				if(rand_num < pi){
					flip_spin(u,2,i);
				}
			}
			for(int i=0;i<NOF;i++){
				
				u = G3.gen()%int(N);
				
				if(inter_connected3[u])
					pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*(M2/N)*local_m3[u]*G3.v[u].size()));
				else
					pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*local_m3[u]*G3.v[u].size()));

				
				rand_num = (double)G3.gen()/G3.gen.max();
				if(rand_num < pi){
					flip_spin(u,3,i);
				}
			}*/
			}
		
		calc_energy();
		M_vec.push_back(M1/N);
		E_vec.push_back(E/N);
		E1_vec.push_back(E1/N);
		E2_vec.push_back(E2/N);
		E3_vec.push_back(E3/N);
		cout<<"M = "<<M1/N<<", E1 = "<<E1/N<<", E2 = "<<E2/N<<", E3 = "<<E3/N<<", E = "<<E/N<<",    T="<<T<<endl;


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
			dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<"	"<<(double) E1_vec[i]<<"	"<<(double) E2_vec[i]<<"	"<<(double) E3_vec[i]<<"	"<<(double) E_vec[i]<<endl;
	}

	dataInfo.close();

}

void ising_interdependent_3_layers::revive_ordered_up(){
	for(int i=0;i<N;i++){    //flip all the spins back up
		if(spins1[i] == -1)
			flip_spin(i,1,i);
		if(spins2[i] == -1)
			flip_spin(i,2,i);
		if(spins3[i] == -1)
			flip_spin(i,3,i);
	}
}

void ising_interdependent_3_layers::revive_ordered_down(){
	for(int i=0;i<N;i++){    //flip all the spins back up
		if(spins1[i] == 1)
			flip_spin(i,1,i);
		if(spins2[i] == 1)
			flip_spin(i,2,i);
		if(spins2[i] == 1)
			flip_spin(i,3,i);
	}
}

void ising_interdependent_3_layers::revive_disordered(){
	double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

	for(int i=0;i<N;i++){ //network1
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			flip_spin(i,1,i);
		}
	}

	for(int i=0;i<N;i++){ //network2
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			flip_spin(i,2,i);
		}
	}
	
	for(int i=0;i<N;i++){ //network3
		rand_num = (double)G1.gen()/G1.gen.max();
		if(rand_num < 0.5){
			flip_spin(i,3,i);
		}
	}
}


#endif /* SRC_ising_interdependent_3_layers_HPP_ */
