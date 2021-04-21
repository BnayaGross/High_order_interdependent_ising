
using namespace std;
#include <queue>
#include <omp.h>
#include <pthread.h>
typedef std::mersenne_twister_engine< uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;

class HyperGraph_3_q {
public:
        int n;
        vector<int> components; //used for testconnectivity
        vector<vector<vector<int> > > v;
        vector<int> active;
        vector <int> connected;	//connected to the giant
        mt11213b gen;

        HyperGraph_3_q();
        ~HyperGraph_3_q();
		void setHyperGraph_3_q(int size);
		void createHyperGraphER_3_q(double k, double q);
		int testConnectivity();
};
HyperGraph_3_q::HyperGraph_3_q() {n=0;}
void HyperGraph_3_q::setHyperGraph_3_q(int size) {
      if (size < 2) n = 2;
      else n = size;
      v.resize(n);
      connected.resize(n);
      components.resize(n);
      active = vector<int>(n,1);
      gen.seed(time(0));
}
HyperGraph_3_q::~HyperGraph_3_q() {
  for(int i=0;i<v.size();i++)
	  v[i].clear();
}


void HyperGraph_3_q::createHyperGraphER_3_q(double k, double q) {
	//k is average degree
	
	srand(time(0));
	int s,t,u,w,i=0;
	
	double E4 = (n*k*q)/4; //number of squares
	double E3 = (n*k*(1-q))/3; // number of triangles
	
	while(i<E4)
	{
		s = gen()%n;
		t = gen()%n;
		u = gen()%n;
		w = gen()%n;
		if(s!=t and t!=u and s!=u and s!=w and t!=w and u!=w)
		{
			vector<int> vec {s,t,w};
			v[u].push_back(vec);
			vector<int> vec1 {u,t,w};
			v[s].push_back(vec1);
			vector<int> vec2 {s,u,w};
			v[t].push_back(vec2);
			vector<int> vec3 {s,t,u};
			v[w].push_back(vec3);
			i++;
		}
	}
	
	i = 0;
	while(i<E3)
	{
		s = gen()%n;
		t = gen()%n;
		u = gen()%n;
		if(s!=t and t!=u and s!=u)
		{
			vector<int> vec {s,t};
			v[u].push_back(vec);
			vector<int> vec1 {u,t};
			v[s].push_back(vec1);
			vector<int> vec2 {s,u};
			v[t].push_back(vec2);
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/3<<" edges."<<endl;
}


int HyperGraph_3_q::testConnectivity()
{
	queue<int> Q;
	int j, size,expCount=0, cluster_id=0, giant_size=0, giant_id;
	int *explored = new int[n];//0-white, 1-grey, 2-black
	vector<int>* sources = new vector<int>;
	bool flag = true;

	for (int i = 0; i < n; ++i) //initialization explored
	{
		    if (active[i]==0)
		    {
		    	explored[i] = 2;
		    	expCount++;
		    	components[i]=cluster_id++;
		    }
		    else
		    	explored[i] = 0;
	}
	j=0;
	cluster_id--;
	while(expCount<n &&flag)
	{
	while(explored[j]!=0) //find a source for BFS algorithm
	{
	    j++;
	}
	cluster_id++;
	Q.push(j);
	components[j]=cluster_id;
	sources->push_back(j);
	size=1;
	explored[j] = 1;
	expCount++;

	while (!Q.empty()) {

	    int u = Q.front();
	    Q.pop();
	    explored[u] = 2;

	for(int i=0;i<v[u].size();i++)
	{
		for(int jjjj = 0; jjjj < v[u][i].size(); jjjj++){
			
			if(explored[v[u][i][jjjj]]==0){
				size++;
				Q.push(v[u][i][jjjj]);
				components[v[u][i][jjjj]] = cluster_id;
				explored[v[u][i][jjjj]] = 1;
				expCount++;
			}
			
		}
	}

	}
	if(size>giant_size)
	{
		giant_size=size;
		giant_id=cluster_id;
	}
	}
	for(int k=0;k<n;k++)
	{
		if(components[k] == giant_id)
			connected[k] = 1;
		else
			connected[k] = 0;
	}

	delete [] explored;
	return giant_size;
}
