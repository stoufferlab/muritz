
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
    vector<Node *> predators;
    vector<Node *> prey;
    map<unsigned int, set<Node *> > neighbors;
};

// Position structure
typedef struct {
        string name;
        int size;
        int frequency;
} Position;

// Role structure
typedef struct {
        string name;
        vector<Position> f;
} Role;

// Network structure
typedef struct {
        string name;
        vector<Node *> nodes;
        map<string,int> node_i;
        vector<Role> roles;
} Network;

// Alignment information
typedef struct {
        vector<pair<int,int> > matches;
        vector<int> match1;
        vector<int> match2;
        double (*dfunc)(Role*,Role*);

        int iters_fixed_T;
        double t_initial;
        double mu_t;
        double t_min;
        unsigned int degree;
} Alignment;

#endif
