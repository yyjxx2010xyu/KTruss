#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <set>
#include <tuple>
#include "edge.h"
#include "tools.h"
#include "util.h"
#include "libpopcnt.h"
#include "fixedvector.h"
#include "cachebuffer.h"

class Graph
{
private:
    uint32_t nodeNum;
    uint32_t edgeNum;
    uint32_t supportNum;

    uint32_t *supports;

    uint32_t *edgeList;
    uint32_t *edgeId;
    uint32_t *edgeRev;
    uint32_t *edgeListBegPtr;
    uint32_t *edgeListEndPtr;
    std::vector<Edge> *ptrEdges;

    bool *curStage;
    uint32_t k;
    uint32_t realDel;
    double tcInitTime;
    double tcBiasTime;
    double triTime;

    uint32_t *LastToCur;
    uint32_t *CurToLast;

    uint32_t Triangles;


public:
    std::vector<std::vector<int>> bitIdList;
    std::vector<std::vector<WordType>> bitList;

public:
    Graph(int _nodeNum, int _edgeNum)
    {

        tcInitTime = 0.0;
        tcBiasTime = 0.0;
        triTime = 0.0;
        Triangles = 0;
        nodeNum = _nodeNum;
        edgeNum = _edgeNum;

        bitIdList.resize(nodeNum + 1);
        bitList.resize(nodeNum + 1);

        uint32_t mallocEdgeNum = edgeNum << 1 | 1;

        edgeList = (uint32_t *)malloc(sizeof(uint32_t) * mallocEdgeNum);
        edgeId = (uint32_t *)malloc(sizeof(uint32_t) * mallocEdgeNum);
        edgeRev = (uint32_t *)malloc(sizeof(uint32_t) * mallocEdgeNum);
        edgeListBegPtr = (uint32_t *)malloc(sizeof(uint32_t) * (2 * nodeNum + 2));
        edgeListEndPtr = (uint32_t *)malloc(sizeof(uint32_t) * (2 * nodeNum + 2));
        supports = (uint32_t *)malloc(sizeof(uint32_t) * mallocEdgeNum);
        curStage = (bool *)malloc(sizeof(bool) * mallocEdgeNum);

        LastToCur = (uint32_t *)malloc(sizeof(uint32_t) * mallocEdgeNum);
        CurToLast = (uint32_t *)malloc(sizeof(uint32_t) * mallocEdgeNum);

        memset(edgeList, 0, sizeof(uint32_t) * mallocEdgeNum);
        memset(edgeListBegPtr, 0, sizeof(uint32_t) * (2 * nodeNum + 2));
        memset(edgeListEndPtr, 0, sizeof(uint32_t) * (2 * nodeNum + 2));
        memset(supports, 0, sizeof(uint32_t) * mallocEdgeNum);
        memset(curStage, false, sizeof(bool) * mallocEdgeNum);

    }
    ~Graph()
    {
        free(edgeList);
        free(edgeId);
        free(edgeRev);
        free(edgeListBegPtr);
        free(edgeListEndPtr);
        free(supports);
        free(curStage);
        free(CurToLast);
        free(LastToCur);

    }

    void LoadEdge(std::vector<Edge> &Edges);
    void DeleteDup();
    void SortEdge();
    void CountTriangle();

    uint32_t BinarySearch(uint32_t l, uint32_t r, const uint32_t value);
    uint32_t GetRevId(const uint32_t u, const uint32_t v);

    void RemoveEdge(const uint32_t u, const uint32_t v, const uint32_t uvId);

    uint32_t CountInter(const uint32_t u, const uint32_t v);
    void FindInter(const uint32_t u, const uint32_t v, std::vector<std::pair<uint32_t, uint32_t>> &AffectEdgeList);
    void TrussDecomp();
    void FindKMax();
    void PeelTriangle(const uint32_t uvId, const uint32_t uwId, const uint32_t vwId, CacheBuffer &Buffer);
    void PackVertex(uint32_t u);
    void Update(const uint32_t eId, CacheBuffer &Buffer) noexcept;
    void DeleteAll(FixedVector &QCur);
    void PeelAll(FixedVector &QCur, FixedVector &QNext);
    void CompactCSR();
    void ShrinkCSR();
    void ReCountTriangle(FixedVector &QCur, FixedVector &QNext);
    void NewTC();
    uint32_t SearchFirst(uint32_t u);

    void MordernTruss();
    void InitRev();
};

#endif