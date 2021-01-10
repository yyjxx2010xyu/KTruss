#ifndef EDGE_H
#define EDGE_H
#include <iostream>
#include <vector>
#include "util.h"

void ReadBaseLine(const char * path,  std::vector<Edge> & Edges);
void LoadEdge(const char * const bptr, const uint32_t len, std::vector<Edge> & Edges);
uint32_t mysscanf(const char * const eptr, const char *& cptr);
void sscanfSkip(const char * const eptr, const char *& cptr);

#endif
