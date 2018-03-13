#ifndef __ROLES_HPP
#define __ROLES_HPP

#include "common.hpp"

vector<Role> read_roles(char *filename, char separator, Network N);
vector<Role> read_roles(char separator, Network N);
void print_roles(Network N);

#endif
