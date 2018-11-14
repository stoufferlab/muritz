
#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

// Node structure
struct Node {
    string name;
    int idx;
    int side;
    vector<Node *> predators;
    vector<Node *> prey;
    map<unsigned int, set<Node *> > neighbors;
};

// Position structure
typedef struct {
    string name;
    int size;
    double frequency;
} Position;

// Role structure
typedef struct {
    string name;
    vector<Position> f;
} Role;

// Network structure
typedef struct {
    string name;
    bool bipartite;
    vector<Node *> nodes;
    map<string,int> node_i;
    vector<Role> roles;
} Network;

// Alignment information
typedef struct {
    bool doneflag;
    vector<pair<int,int> > matches;
    vector<pair<int,int> > set_pairs;
    vector<int> fixed_pairs;
    vector<int> unfixed_pairs_A;
    vector<int> unfixed_pairs_B;
    vector<int> match1;
    vector<int> match2;
    double (*dfunc)(Role*,Role*);
    
    // The matches contributing towards the total energy, and the number of them.
    // Only used if degree != 0.
    map<pair<int, int>, int> matchesContributing;
    
    int p1, p2;// The proposed pair to switch.
    
    double energy;
    double proposedEnergy;
    map<pair<int, int>, int> proposedContributionDeltas;
    set<pair<int, int> > proposedContributionWipes;
    
    vector<pair<int,int> > proposedMatches;
    vector<int> proposedMatch1;
    vector<int> proposedMatch2;
    
    unsigned int degree;
    //add a set pairs vector<pair<int, int>>
} Alignment;

#endif
