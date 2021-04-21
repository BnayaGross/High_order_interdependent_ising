

#ifndef SRC_MYGRAPH_HPP_
#define SRC_MYGRAPH_HPP_


using namespace std;
#include <queue>
#include <omp.h>
#include <pthread.h>
typedef std::mersenne_twister_engine< uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;

class Graph {
public:
        int n;
        vector<int> components; //used for testconnectivity
        vector<vector<int> > v;
        vector<int> active;
        vector <int> connected;	//connected to the giant
        mt11213b gen;

        Graph();
        ~Graph();
        void setGraph(int size);
        bool isConnected(int, int);
        void addEdge(int , int );
        int getDegree(int);
        void createGraphER(double);
        int testConnectivity();
	
};
Graph::Graph() {n=0;}
void Graph::setGraph(int size) {
      if (size < 2) n = 2;
      else n = size;
      v.resize(n);
      connected.resize(n);
      components.resize(n);
      active = vector<int>(n,1);
      gen.seed(time(0));
}
Graph::~Graph() {
  for(int i=0;i<v.size();i++)
	  v[i].clear();
}

bool Graph::isConnected(int u1, int u2) {
    if(active[u1]==0||active[u2]==0)
    	return false;
    for(int i=0;i<v[u1].size();i++)
    	if(v[u1][i]==u2)
    		return true;
    return false;
}
inline void Graph::addEdge(int u1, int u2) {
    v[u1].push_back(u2);
    v[u2].push_back(u1);
}

int Graph::getDegree(int u)
{
	return v[u].size();

}

void Graph::createGraphER(double k) {
	//k is average degree
	srand(time(0));
	int s,t,i=0;
	while(i<(n*k)/2)
	{
		s = gen()%n;
		t = gen()%n;
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/2<<" edges."<<endl;
}

int Graph::testConnectivity()
{// we are going to change the BFS algorithem a little bit
	// we will remember the nodes of the giant and then we will make connected
	queue<int> Q;
	int j, size,expCount=0, cluster_id=0, giant_size=0, giant_id;
	int *explored = new int[n];//0-white, 1-grey, 2-black
	vector<int>* sources = new vector<int>;
	//vector<int>* giant_nodes;
	bool flag = true;

	for (int i = 0; i < n; ++i) //initialization explored
	{
		    if (active[i]==0)
		    {
		    	explored[i] = 2;
		    	expCount++;
		    	components[i]=cluster_id++;
		    }
		    else if(getDegree(i)==0)
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
	    if(explored[v[u][i]]==0)
	    {
	    	size++;
	    	Q.push(v[u][i]);
	    	components[v[u][i]] = cluster_id;
	    	explored[v[u][i]] = 1;
	    	expCount++;
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

#endif
