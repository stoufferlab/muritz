#ifndef __NETWORK_HPP
#define __NETWORK_HPP

#include <iostream>
#include <common.hpp>

Network read_network(char *filename, char separator);
Network read_network(char separator);
void print_network(Network N);

#endif