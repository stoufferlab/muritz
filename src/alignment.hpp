
#ifndef __ALIGNMENT_HPP
#define __ALIGNMENT_HPP

#include <common.hpp>

Alignment * setup_alignment();

Alignment * alignment_alloc(size_t n);
void alignment_free(Alignment *a);

#endif
