
#ifndef __ALIGNMENT_HPP
#define __ALIGNMENT_HPP

#include <common.hpp>

void read_alignment_data(char separator,Network& A,Network& B);
Alignment * setup_alignment();
Alignment * alignment_alloc(size_t n);
void alignment_free(Alignment *a);

#endif
