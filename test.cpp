#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <memory>

using namespace std;

class Node {
  private:
    pair<int,int> loc;
    vector<shared_ptr<Node> > adj;
  public:

    Node(int i, int j) {
      loc = make_pair(i,j);
    }

    void makeEdge(shared_ptr<Node> n) {
      adj.push_back(n);
    }

    const pair<int, int> getLoc(){
      return loc;
    }

    friend ostream& operator<<(ostream& os, const Node& n);
};

ostream& operator<<(ostream& os, const Node& n) {
      os << n.loc.first << ' ' << n.loc.second << ": ";
      for (auto a: n.adj) {
        auto loc = a->getLoc();
        os << loc.first << " " << loc.second << "  ";
      }
      return os;
    }

class Graph{
  private:
    vector<shared_ptr<Node> > nodes;

  public:
    // Initialize Graph as grid of nodes
    Graph(int h, int w){
      // Create and store nodes in a temp grid
      vector<vector<shared_ptr<Node> > > grid;
      for (int i = 0; i < h; i++) {
        grid.push_back(vector<shared_ptr<Node> >());
        for (int j = 0; j < w; j++) {
          auto node_ptr = make_shared<Node>(i,j);
          grid[i].push_back(node_ptr);
        }
      }
      for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
          if (i > 0) {
            grid[i][j]->makeEdge(grid[i-1][j]);
          }
          if (i < h-1) {
            grid[i][j]->makeEdge(grid[i+1][j]);
          }
          if (j > 0) {
            grid[i][j]->makeEdge(grid[i][j-1]);
          }
          if (j < h-1) {
            grid[i][j]->makeEdge(grid[i][j+1]);
          }
        }
      }
      for (auto row: grid) {
        nodes.insert(nodes.end(), row.begin(), row.end());
      }
    }

    void printGraph() {
      for (auto n: nodes) {
        cout << *n << endl;
      }
    }
};

float h(pair<int, int> p, pair<int, int> g) {
  return sqrt(pow(p.first-g.first, 2) + pow(p.second-g.second,2));
}

int main() {
  Graph g(4,4);
  g.printGraph();

  pair<int,int> start (0,0);
  // priority_queue<pair<int,int> > q;
  // set<pair<int,int> > opened;
  // unordered_map<pair<int, int>, float> g;
  // g.setVal(3,3,1);

  // s.push(start);
  // while (!s.empty()) {
  //   pair<int,int> curr = s.top();
  //   s.pop();
    
  //   if (opened.count(curr) == 0) {
  //     g.printPair(curr);
  //     opened.insert(curr);
  //     if (g.getVal(curr) == 1) {
  //     cout << "Found!" << endl;
  //     break;
  //     }
      
  //     vector<pair<int,int> > adjlist = g.getAdj(curr);
  //     for (pair<int,int> p: adjlist) {
  //       s.push(p);
  //     }
  //   }
    
    
  }
  


