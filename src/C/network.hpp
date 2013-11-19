#ifndef __NETWORK_HPP
#define __NETWORK_HPP

#include <iostream>
#include <set>
#include <common.hpp>

Network read_network(char *filename, char separator);
Network read_network(char separator);
void read_networks(char separator,Network& A, Network& B);
void print_network(Network N);

set<Node *> neighbors(Network *N, Node *n, unsigned int degree, int direction);

#endif
