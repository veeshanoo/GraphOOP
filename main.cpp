#include <bits/stdc++.h>

using namespace std;

// Generic function that swaps elements
template <class T>
void swap(T *a, T *b) {
	T aux = *a;
	*a = *b;
	*b = aux;
}

// Generic function that sorts an array
template <class T> 
void Sort(T a[], int N) {
	for (int i = 1; i < N; i++)
		for (int j = i; j > 0; j--)
			if (a[j] < a[j - 1])
				swap<T>(&a[j], &a[j - 1]);
}

template <class T>
class Node {
private:
	int id;
	T cost;
	Node<T> *nxt;
public:
    Nod() {
        id = 0;
        cost = 0;
        nxt = NULL;
    }
	Node(const int& v, const T& c) {
		id = v;
		cost = c;
		nxt = NULL;
	}
	~Node() {}
	int getId() {
		return id;
	}
	int getCost() {
		return cost;
	}
	Node *getNxt() {
		return nxt;
	}
	template <class U>
	friend class Graph;
}; 

template <class T>
class Edge {
private:
	int x;
	int y;
	T cost;
public:
	Edge() {
		x = 0;
		y = 0;
		cost = 0;
	}
	Edge(const int& a, const int& b, const T& c) {
		x = a;
		y = b;
		cost = c;
	}
	~Edge() {}
	int getX() {
		return x;
	}
	int getY() {
		return y;
	}
	T getCost() {
		return cost;
	}
	bool operator < (const Edge& other) const {
		return cost < other.cost;
	}
	Edge& operator = (const Edge& other) {
		x = other.x;
		y = other.y;
		cost = other.cost;
		return *this;
	}
	bool operator == (const Edge& other) const {
		return (x == other.x && y == other.y) || (x == other.y && y == other.x);
	}
	bool operator != (const Edge& other) const {
		!(*this == other);
	}
	template <class U> friend class Graph;
};

template <class T>
class Graph {
private:
	int nrNodes;
	int nrEdges;
	int nrComp;
	Node<T> **v;
	Node<T> **comp;
	Edge<T> *edges;
	int *seen;
	int *comp_id;
	int *comp_sz;
public:
	// Constructor
	Graph() {
		nrNodes = 0;
		nrEdges = 0;	
		nrComp = 0;
	}
	// Copy constructor
	Graph(const Graph<T>& other) {
		nrNodes = other.nrNodes;
		nrEdges = other.nrEdges;
		nrComp = 0;
		Reserve();
		for (int i = 0; i < nrEdges; i++)
			edges[i] = other.edges[i];
		for (int i = 0; i < nrEdges; i++) {
			addEdge(&v[edges[i].x], edges[i].y, edges[i].cost);
			addEdge(&v[edges[i].y], edges[i].x, edges[i].cost);
		}
	}
	// = operator
	Graph<T>& operator = (const Graph& other) {
		if (this != &other) {
			this -> Free();
			nrNodes = other.nrNodes;
			nrEdges = other.nrEdges;
			nrComp = 0;
			Reserve();
			for (int i = 0; i < nrEdges; i++)
				edges[i] = other.edges[i];
			for (int i = 0; i < nrEdges; i++) {
				addEdge(&v[edges[i].x], edges[i].y, edges[i].cost);
				addEdge(&v[edges[i].y], edges[i].x, edges[i].cost);
			}
		}
	}
	// Destructor
	~Graph() {
		Free();
	}
	// Function declarations
	template <class U>
	friend Graph<U> operator * (Graph<U> lhs, const Graph<U>& rhs);
	int getNrNodes() {
		return nrNodes;
	}
	int getNrEdges() {
		return nrEdges;
	}
	void buildGraph();
	void Reserve();
	void Free();
    void addEdge(Node<T>**, const int&, const int&);
	int **RoyFloydMatrix();
	void Dfs(const int&, const int&);
	int *DijkstraPath(const int&, const int&);
	void getConnectedComponents();
	int isConnected();
	Edge<T> *Apm(const int&);
	// I/O
	template <class U> friend istream& operator >>(istream&, Graph<U>&);
	template <class U> friend ostream& operator <<(ostream&, const Graph<U>&);
};

// Graph input
template <class T>
istream& operator >>(istream& in, Graph<T>& g) {
	in >> g.nrNodes;
	in >> g.nrEdges;
	g.Reserve();
	for (int i = 0; i < g.nrEdges; i++) {
		int x, y; T c;
		in >> x >> y >> c;
		g.addEdge(&g.v[x], y, c);
		g.addEdge(&g.v[y], x, c);
		g.edges[i] = Edge<T>(x, y, c);
	}
	return in;
}

// Graph output
template <class T>
ostream& operator <<(ostream& out, const Graph<T>& g) {
	out << g.nrNodes << '\n';
	out << g.nrEdges << '\n';
	for (int i = 0; i < g.nrNodes; i++) {
		out << i << " : ";
		for (Node<T> *x = g.v[i]; x != NULL; x = x -> getNxt())
			out << "(" << x -> getId() << ", " << x -> getCost() << ") ";
		out << '\n';
	}
	return out;
}

// Memory alloc to store the edges of the graph
template <class T>
void Graph<T>::Reserve() {
	v = new Node<T>*[nrNodes];
	comp = new Node<T>*[nrNodes];
	seen = new int[nrNodes];
	comp_id = new int[nrNodes];
	comp_sz = new int[nrNodes];
    for (int i = 0; i < nrNodes; i++) {
        v[i] = comp[i] = NULL;
        seen[i] = comp_id[i] = comp_sz[i] = 0;
    }
	if (nrEdges > 0) 
		edges = new Edge<T>[nrEdges];
}

template <class T>
void Graph<T>::Free() {
    for (int i = 0; i < nrNodes; i++) {
	    Node<T> *aux;
		while (v[i] != NULL) {
			aux = v[i] -> nxt;
			delete v[i];
			v[i] = aux;
		}
	}
	for (int i = 0; i < nrComp; i++) {
		Node<T> *aux;
		while (comp[i] != NULL) {
			aux = comp[i] -> nxt;
			delete comp[i];
			comp[i] = aux;
		}
	}
	if (nrNodes != 0) {
		delete[] v;
		delete[] comp;
		delete[] seen;
		delete[] comp_id;
		delete[] comp_sz;
	}
	if (nrEdges != 0) {
		delete[] edges;
	}
	nrNodes = nrEdges = nrComp = 0;		    
}

// Funtion that adds an edge to the graph
template <class T>
void Graph<T>::addEdge(Node<T> **x, const int& y, const int& c) {
	if (*x == NULL) {
		*x = new Node<T>(y, c);
		return;
	}
	Node<T> *aux = new Node<T>(y, c);
	aux -> nxt = *x;
	*x = aux;
}

// Roy-Floyd algorithm
template <class T>
int **Graph<T>::RoyFloydMatrix() {
	int **dp;
	dp = new int*[nrNodes];
	
	for (int i = 0; i < nrNodes; i++)
		dp[i] = new int[nrNodes];
	for (int i = 0; i < nrNodes; i++)
		for (int j = 0; j < nrNodes; j++)
			dp[i][j] = 0;

	for (int i = 0; i < nrNodes; i++) {
		int aux = i;
		for (Node<T> *x = v[i]; x != NULL; x = x -> nxt) {
			if (dp[aux][x -> id] == 0 || x -> cost < dp[aux][x -> id]) 
				dp[aux][x -> id] = x -> cost;
		}
	}
	
	for (int k = 0; k < nrNodes; k++)
		for (int i = 0; i < nrNodes; i++)
			for (int j = 0; j < nrNodes; j++) 
				if (dp[i][k] != 0 && dp[k][j] != 0 && i != j && (dp[i][j] == 0 || dp[i][k] + dp[k][j] < dp[i][j]))
					dp[i][j] = dp[i][k] + dp[k][j];
	
	cout << "Roy-Floyd matrix:\n";
	for (int i = 0; i < nrNodes; i++, cout << '\n')
		for (int j = 0; j < nrNodes; j++)
			cout << dp[i][j] << ' ';

	return dp;
}

// Dijkstra algorithm (output[0] will store the size of the path)
template <class T>
int *Graph<T>::DijkstraPath(const int& from, const int& to) {
	int dp[nrNodes];
	int prv[nrNodes];
	int q[1 << 13];
	int q2[1 << 13];
	int tp = 0, tp2 = 1;

	for (int i = 0; i < nrNodes; i++)
		dp[i] = (int) 2e9, prv[i] = -1;
	
	dp[from] = 0;
	q[tp2] = 0;
	q2[tp2] = from;
	
	while (tp != tp2) {
		++tp;
		int c = q[tp];
		int nod = q2[tp];
		if (c > dp[nod])
			continue;
		for (Node<T> *x = v[nod]; x != NULL; x = x -> nxt)
			if (c + x -> cost < dp[x -> id]) {
				dp[x -> id] = c + x -> cost;
				prv[x -> id] = nod;
				++tp2;
				q[tp2] = dp[x -> id];
				q2[tp2] = x -> id;	
			}
	}

	if (dp[to] == (int) 2e9) {
		cout << "There is no path from " << from << " to " << to << '\n';
		return NULL;
	}
	
	int *path;
	int p = to;
	int cnt = 0;
	
	while (p != -1) {
		cnt++;
		p = prv[p];
	}
	
	path = new int[cnt + 1];
	path[0] = cnt;
	p = to;
	
	while (p != -1) {
		path[cnt--] = p;
		p = prv[p];
	}

	cout << "Shortest path from " << from << " to " << to << " is:\n";
	for (int i = 1; i <= path[0]; i++)
		cout << path[i] << ' ';
	cout << '\n';

	return path;
}

// Funtion that finds a connected component
template <class T>
void Graph<T>::Dfs(const int& from, const int& component) {
	seen[from] = 1;
	comp_id[from] = component;
	comp_sz[component]++;
	addEdge(&comp[component], from, 0);
	for (Node<T> *x = v[from]; x != NULL; x = x -> nxt)
		if (!seen[x -> id])
			Dfs(x -> id, component);
}

// Function that builds connected components
template <class T>
void Graph<T>::getConnectedComponents() {
	for (int i = 0; i < nrNodes; i++)
		seen[i] = comp_id[i] = 0;
	
	nrComp = 0;
	for (int i = 0; i < nrNodes; i++)
		if (seen[i] == 0) {
			Dfs(i, nrComp);
			nrComp++;
		}
	
	cout << "Connected components:\n";
	for (int i = 0; i < nrComp; i++, cout << '\n') {
		cout << i << " : ";
		for (Node<T> *j = comp[i]; j != NULL; j = j -> nxt)
			cout << j -> id << ' ';
	}
}

// Funtion that checks if a given graph is connected
template <class T>
int Graph<T>::isConnected() {
	return nrComp == 1;
}

// Funtion that returns minimum spanning tree of a connected component
template <class T>
Edge<T> *Graph<T>::Apm(const int& x) {
	Sort(edges, nrEdges);
	int id = comp_id[x];
	int dad[nrNodes];
	int cnt = 0;
	int ans = 0;
	Edge<T> *apm = new Edge<T>[comp_sz[id] - 1];
	
	for (int i = 0; i < nrNodes; i++)
		dad[i] = i;
	
	for (int i = 0; i < nrEdges; i++) {
		if (comp_id[edges[i].getX()] != id || comp_id[edges[i].getY()] != id)
			continue;
		int X = edges[i].getX();
		int Y = edges[i].getY();
		while (dad[X] != X) 
			X = dad[X];
		while (dad[Y] != Y) 
			Y = dad[Y];
		if (X == Y)
			continue;
		ans += edges[i].cost;
		apm[cnt++] = edges[i];
		dad[X] = Y;
	}
	cout << "MST containing vertex " << x << " has a cost of " << ans << ":\n";
	for (int i = 0; i < cnt; i++) {
		cout << apm[i].x << ' ' << apm[i].y << ' ' << apm[i].cost << '\n';
	}
	return apm;
}


// Returns edge intersection of two graphs with the same set of vertices
template <class T>
Graph<T> operator * (Graph<T> lhs, const Graph<T>& rhs) {
	Edge<T> aux[lhs.nrEdges];
	int cnt = 0;
	
	for (int i = 0; i < lhs.nrEdges; i++)
		for (int j = 0; j < rhs.nrEdges; j++)
			if (lhs.edges[i] == rhs.edges[j]) 
				if (lhs.edges[i] < rhs.edges[j])
					aux[cnt++] = lhs.edges[i];
				else aux[cnt++] = rhs.edges[j];
	
	Sort(aux, cnt);
	int p = 0, pp, cnt2 = 0;;
	
	while (p < cnt) {
		pp = p;
		while (pp < cnt && aux[pp] == aux[p])
			pp++;
		aux[cnt2++] = aux[p];
		p = pp;
	}
	
	lhs.~Graph();
	lhs.nrNodes = rhs.nrNodes;
	lhs.nrEdges = cnt2;
	lhs.Reserve();
	
	for (int i = 0; i < cnt2; i++) {
		lhs.edges[i] = aux[i];
		lhs.addEdge(&lhs.v[aux[i].getX()], aux[i].getY(), aux[i].getCost());
		lhs.addEdge(&lhs.v[aux[i].getY()], aux[i].getX(), aux[i].getCost());
	}
	return lhs;
}

int main() {
	ifstream in("tst.in");
	ofstream out("tst.out");
	Graph<int> a, b;
	in >> a >> b;
	out << a;
	a.getConnectedComponents();
	a.Apm(0);
	a.DijkstraPath(0, 2);
	a.RoyFloydMatrix();
	cout << a.isConnected() << '\n';
	Graph<int> c;
	c = a * b;
	out << c;
    Graph<int> d = a * b * c;
    return 0;
}