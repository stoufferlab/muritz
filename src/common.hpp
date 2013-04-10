#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <string>
#include <vector>
using namespace std;

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
        vector<Role> roles;
} Network;

// Alignment information
typedef struct {
	Network *net1;
	Network *net2;

        vector<pair<int,int> > matches;
        vector<int> match1;
} Alignment;

#endif