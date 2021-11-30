
#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

// The different distance measures and their flag values.
// Make sure to keep these matching the flag values!
enum class cost_func_flag_t {euc=0, cor=1, chisq=2, maha=3};


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
    
    // The matches contributing towards the total energy, and the number of them.
    // Only used if degree != 0.
    map<pair<int, int>, int> matchesContributing;
    
    int p1, p2;// The proposed pair to switch.
    
    double energy;
    double proposedEnergy;
    
    // Proposed changes to matchesContributinhg, also only used if degree != 0.
    map<pair<int, int>, int> proposedContributionDeltas;
    set<pair<int, int> > proposedContributionWipes;
    
    unsigned int degree;
    //add a set pairs vector<pair<int, int>>
} Alignment;

//The type of a pointer to a distance function.
typedef double (*dist_func_t)(Role*,Role*);

#endif
