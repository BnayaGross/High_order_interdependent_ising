
using namespace std;
#include <queue>
#include <omp.h>
#include <pthread.h>
typedef std::mersenne_twister_engine< uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;

//HyperGrpah 2+q

class HyperGraph_2_q {
public:
        int n;
        vector<int> components; //used for testconnectivity
        vector<vector<pair<int,int> > > v;
        vector<int> active;
        vector <int> connected;	//connected to the giant
        mt11213b gen;

        HyperGraph_2_q();
        ~HyperGraph_2_q();
		void setHyperGraph_2_q(int size);
		void createHyperGraphER_2_q(double k, double q);
		bool isConnected(int u1, int u2, int u3);
		bool isConnected2(int u1, int u2);
		int testConnectivity();
};
HyperGraph_2_q::HyperGraph_2_q() {n=0;}
void HyperGraph_2_q::setHyperGraph_2_q(int size) {
      if (size < 2) n = 2;
      else n = size;
      v.resize(n);
      connected.resize(n);
      components.resize(n);
      active = vector<int>(n,1);
      gen.seed(time(0));
}
HyperGraph_2_q::~HyperGraph_2_q() {
  for(int i=0;i<v.size();i++)
	  v[i].clear();
}

void HyperGraph_2_q::createHyperGraphER_2_q(double k, double q) {
	//k is average degree
	srand(time(0));
	int s,t,u,i=0;
	
	double E3 = (n*k*q)/3; //number of triangles
	double E2 = (n*k*(1-q))/2; // number of pairs
	
	while(i<E3)
	{
		s = gen()%n;
		t = gen()%n;
		u = gen()%n;
		if(s!=t and t!=u and s!=u && !isConnected(s,t,u))
		{
			pair <int, int> p = make_pair(s,t);
			v[u].push_back(p);
			p = make_pair(s,u);
			v[t].push_back(p);
			p = make_pair(u,t);
			v[s].push_back(p);
			i++;
		}
	}
	
	i = 0;
	while(i<E2)
	{
		s = gen()%n;
		t = gen()%n;
		if(s!=t && !isConnected2(s,t))
		{
			pair <int, int> p = make_pair(s,-1);
			v[t].push_back(p);
			p = make_pair(t,-1);
			v[s].push_back(p);
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/3<<" edges."<<endl;
}

bool HyperGraph_2_q::isConnected(int u1, int u2, int u3) {
    if(active[u1]==0||active[u2]==0||active[u3]==0)
    	return false;
    for(int i=0;i<v[u1].size();i++)
    	if((v[u1][i].first==u2 and v[u1][i].second == u3) or (v[u1][i].first==u3 and v[u1][i].second == u2))
    		return true;
    return false;
}
bool HyperGraph_2_q::isConnected2(int u1, int u2) {
    if(active[u1]==0||active[u2]==0)
    	return false;
    for(int i=0;i<v[u1].size();i++)
    	if(v[u1][i].first==u2 and v[u1][i].second == -1)
    		return true;
    return false;
}


int HyperGraph_2_q::testConnectivity()
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
	    if(explored[v[u][i].first]==0)
	    {
	    	size++;
	    	Q.push(v[u][i].first);
	    	components[v[u][i].first] = cluster_id;
	    	explored[v[u][i].first] = 1;
	    	expCount++;
	    }
		
		if(v[u][i].second != -1){
			if(explored[v[u][i].second]==0)
			{
				size++;
				Q.push(v[u][i].second);
				components[v[u][i].second] = cluster_id;
				explored[v[u][i].second] = 1;
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
	//cout<<"unconnected: "<< n - giant_size<<endl;
	//cout<<giant_size<<endl;
	delete [] explored;
	return giant_size;
}
