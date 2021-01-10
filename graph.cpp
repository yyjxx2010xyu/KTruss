#include "graph.h"
#include "mybitset.h"
#include "timer.h"
#include "hashfilter.h"
#include "fixedvector.h"
#include "cachebuffer.h"

/*
    Input:      没有被处理过的边信息。
    Function:   构建CSR数组。
    CSR初始化：
        edgeList[uPtr] = v;     //  CSR存放与之连接的顶点编号v
        edgeId[uPtr] = i;       //  当前顶点所在原来边的位置i
        edgeRev[uPtr] = vPtr;   //  双向边，相反的边在CSR数组中的Id

    #define USELOCK 在构造的时候可以使用锁，不过无锁操作__sync_bool_compare_and_swap更快一些。
*/
void Graph::LoadEdge(std::vector<Edge> &Edges)
{

    ptrEdges = &Edges;
#ifdef TIMER_ON
    Timer CurTimer;
    CurTimer.StartTimer();
#endif

//  统计每个顶点的度数，e.u + 1的 意义是为了能够在直接后面用于构造CSR数组的BegPtr.
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < edgeNum; i++)
    {
        struct Edge e = Edges[i];
#pragma omp atomic
        edgeListBegPtr[e.u + 1]++;
#pragma omp atomic
        edgeListBegPtr[e.v + 1]++;
    }

#ifdef TIMER_ON
    CurTimer.GapTimer("AddPtr");
#endif

    // Prefix Sum，将每一个顶点度数做前缀和就是CSR数组中的BegPtr。
    for (uint32_t i = 2; i <= nodeNum; i++)
        edgeListBegPtr[i + 1] += edgeListBegPtr[i];
    edgeListBegPtr[1] = 0;

#ifdef TIMER_ON
    CurTimer.GapTimer("PrefixSum");
#endif

#ifdef USELOCK
    omp_lock_t *locks = new omp_lock_t[nodeNum + 2];
    for (uint32_t i = 0; i < nodeNum + 2; i++)
        omp_init_lock(&locks[i]);
#endif

        //  CSR 数组初始化
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < edgeNum; i++)
    {
        struct Edge e = Edges[i];

        uint32_t uPtr;
        do
        {
            uPtr = edgeListBegPtr[e.u];
        } while (!__sync_bool_compare_and_swap(&edgeListBegPtr[e.u], uPtr, uPtr + 1));

        uint32_t vPtr;
        do
        {
            vPtr = edgeListBegPtr[e.v];
        } while (!__sync_bool_compare_and_swap(&edgeListBegPtr[e.v], vPtr, vPtr + 1));

        edgeList[uPtr] = e.v;
        edgeId[uPtr] = i;
#ifdef REV_ON
        edgeRev[uPtr] = vPtr;
#endif
        edgeList[vPtr] = e.u;
        edgeId[vPtr] = i;
#ifdef REV_ON
        edgeRev[vPtr] = uPtr;
#endif

#ifdef USELOCK
        omp_set_lock(&locks[e.u]);
        uint32_t uPtr = edgeListBegPtr[e.u];
        edgeListBegPtr[e.u]++;
        omp_unset_lock(&locks[e.u]);
        edgeList[uPtr] = e.v;
        edgeList[uPtr] = i;

        omp_set_lock(&locks[e.v]);
        uint32_t vPtr = edgeListBegPtr[e.v];
        edgeListBegPtr[e.v]++;
        omp_unset_lock(&locks[e.v]);
        edgeList[vPtr] = e.u;
        edgeList[vPtr] = i;
#endif
    }
#ifdef USELOCK
    free(locks);
#endif

#ifdef TIMER_ON
    CurTimer.GapTimer("LoadCSR");
#endif

    //  Shift CSR PTR to Origin Status，由于在初始化CSR之后，需要恢复BegPtr的位置。
    for (uint32_t i = nodeNum; i >= 2; i--)
        edgeListBegPtr[i] = edgeListBegPtr[i - 1];
    edgeListBegPtr[1] = 0;
    edgeListBegPtr[nodeNum + 1] = edgeNum * 2;

#ifdef TIMER_ON
    CurTimer.GapTimer("ShiftCSR");
#endif
#ifdef DEBUG_ON
    std::cout << "Edges Loaded" << std::endl;
#endif
}

/*
    Input:      
    Function:    删除重边

    重点在于维护edgeRev数组，后面CompactCSR也同理：
        原本需要通过critical 操作，实际上要对以双向的链表维护时需要做一个原子操作，这样会非常慢
        #pragama omp critical
        {
            edgeRev[edgeRev[cureId]] = lasteId;
            edgeRev[lasteId] = cureId;
        }

        需要先与处理一下
        CurToLast[cureId] = lasteId;    // 原来的位置到现在的位置的映射
        LastToCur[lasteId] = cureId;    // 现在的位置到原来的位置的映射
        最后在通过并行的维护：edgeRev[eId] = CurToLast[edgeRev[LastToCur[eId]]];
    做这个目的其实就是为了能够在Truss 的时候能够在O(1)的得到edgeRev
*/
void Graph::DeleteDup()
{

#ifdef TIMER_ON
    Timer CurTimer;
    CurTimer.StartTimer();
#endif

#ifdef REV_ON
//  预处理edgeRev，因为即使是相同的边edgeRev只是指向自己的，预处理使得相同的边的Rev是统一指向最前面的
#pragma omp parallel for schedule(static)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        uint32_t minRev = edgeRev[edgeListBegPtr[u]], firsteId = edgeListBegPtr[u];
        for (uint32_t cureId = edgeListBegPtr[u] + 1; cureId < edgeListBegPtr[u + 1]; cureId++)
            if (edgeList[cureId] == edgeList[cureId - 1])
                minRev = std::min(minRev, edgeRev[cureId]);
            else
            {
                edgeRev[firsteId] = minRev;
                minRev = edgeRev[cureId];
                firsteId = cureId;
            }
        edgeRev[firsteId] = minRev;
    }
#endif

#ifdef TIMER_ON
    CurTimer.GapTimer("InitRev");
#endif

    //  OpenMP 负优化！
    //  去重的主要操作，同时用tot统计真实的边的故事
    //#pragma omp parallel for schedule(dynamic)
    uint32_t tot = 0;
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        uint32_t lasteId = edgeListBegPtr[u];
        for (uint32_t cureId = edgeListBegPtr[u] + 1; cureId < edgeListBegPtr[u + 1]; cureId++)
            if (edgeList[lasteId] != edgeList[cureId])
            {
                edgeList[++lasteId] = edgeList[cureId];
                edgeId[lasteId] = edgeId[cureId];

#ifdef REV_ON
                edgeRev[lasteId] = edgeRev[cureId];
                edgeRev[edgeRev[cureId]] = lasteId;
#endif
            }

        if (edgeListBegPtr[u] == edgeListBegPtr[u + 1])
            edgeListEndPtr[u] = lasteId;
        else
        {
            edgeListEndPtr[u] = lasteId + 1;
            tot += edgeListEndPtr[u] - edgeListBegPtr[u];
        }
    }

#ifdef TIMER_ON
    CurTimer.GapTimer("DelDup");
#endif

    edgeNum = tot >> 1;

#ifdef DEBUG_ON
    std::cout << "tot : " << tot << std::endl;
    std::cout << "Real Edge Num : " << edgeNum << std::endl;

    std::cout << "-----------------------------" << std::endl;
    std::cout << "CSR preview: " << std::endl;

    for (int i = 0; i < 10; i++)
        std::cout << edgeList[i] << " ";
    std::cout << std::endl;
    std::cout << "BegPtr preview: " << std::endl;
    for (int i = 1; i <= 10; i++)
        std::cout << edgeListBegPtr[i] << " ";
    std::cout << std::endl;
    std::cout << "EndPtr preview: " << std::endl;
    for (int i = 1; i <= 10; i++)
        std::cout << edgeListEndPtr[i] << " ";
    std::cout << std::endl;
#if 0
    for (int u = 1; u <= 3; u++)
    {
        std::cout << "u: " << u << " v: ";
        //return;
        for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId += 1)
        {
            std::cout << edgeList[eId] << "-" << edgeId[eId] << " ";
        }
        std::cout << std::endl;
    }
#endif
    std::cout << "-----------------------------" << std::endl;
#endif
}

uint32_t *mapTo;
#ifdef REV_ON
uint32_t *revmapTo;
#endif
uint32_t *globalList;
void Graph::SortEdge()
{
#ifdef TIMER_ON
    Timer CurTimer;
    CurTimer.StartTimer();
#endif

    mapTo = (uint32_t *)malloc(sizeof(uint32_t) * ((edgeNum << 1) + 2));
#ifdef REV_ON
    revmapTo = (uint32_t *)malloc(sizeof(uint32_t) * ((edgeNum << 1) + 2));
#endif
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < (edgeNum << 1 | 1); i++)
        mapTo[i] = i;
#ifdef TIMER_ON
    CurTimer.GapTimer("InitMap");
#endif
    globalList = edgeList;
#pragma omp parallel for schedule(dynamic)
    for (uint32_t i = 1; i <= nodeNum; i++)
    {
        uint32_t *beginPos = &edgeList[edgeListBegPtr[i]];
        uint32_t Len = edgeListBegPtr[i + 1] - edgeListBegPtr[i];
        if (Len <= 1)
            continue;

        std::sort(mapTo + (beginPos - edgeList), mapTo + (beginPos - edgeList) + Len, [](uint32_t a, uint32_t b) {
            return (globalList[a] < globalList[b]);
        });
    }

#ifdef TIMER_ON
    CurTimer.GapTimer("SortMapto");
#endif

#ifdef REV_ON
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < (edgeNum << 1); i++)
        revmapTo[mapTo[i]] = i;
#endif
    uint32_t *newEdgeList = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
    uint32_t *newEdgeId = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
#ifdef REV_ON
    uint32_t *newEdgeRev = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
#endif
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < (edgeNum << 1); i++)
    {
        newEdgeList[i] = edgeList[mapTo[i]];
        newEdgeId[i] = edgeId[mapTo[i]];
#ifdef REV_ON
        newEdgeRev[i] = revmapTo[edgeRev[mapTo[i]]];
#endif
    }

#ifdef TIMER_ON
    CurTimer.GapTimer("ApplyMapTo");
#endif
    free(edgeList);
    free(edgeId);
#ifdef REV_ON
    free(edgeRev);
#endif

    edgeList = newEdgeList;
    edgeId = newEdgeId;
#ifdef REV_ON
    edgeRev = newEdgeRev;
#endif

    free(mapTo);
#ifdef REV_ON
    free(revmapTo);
#endif
#ifdef DEBUG_ON
    std::cout << "Edges Sorted" << std::endl;
#endif
}

void Graph::PackVertex(uint32_t u)
{
    if (edgeListBegPtr[u] >= edgeListEndPtr[u] || ((edgeListEndPtr[u] - edgeListBegPtr[u]) <= 512))
        return;

    uint32_t lastBlockId = MINUS;
    uint32_t blockNum = 0;
    for (uint32_t j = edgeListBegPtr[u]; j < edgeListEndPtr[u]; j++)
    {
        uint32_t v = edgeList[j];
        uint32_t curBlockId = v >> ShiftBits;
        if (curBlockId != lastBlockId)
        {
            lastBlockId = curBlockId;
            blockNum++;
        }
    }
    lastBlockId = MINUS;

    //  Packing Threshold (wpt)
    if ((edgeListEndPtr[u] - edgeListBegPtr[u]) >= 512 || (edgeListEndPtr[u] - edgeListBegPtr[u]) / blockNum > 2)
    {
        bitIdList[u].clear();
        bitList[u].clear();
        for (uint32_t j = edgeListBegPtr[u]; j < edgeListEndPtr[u]; j++)
        {
            uint32_t v = edgeList[j];
            uint32_t curBlockId = v >> ShiftBits;
            if (curBlockId != lastBlockId)
            {
                lastBlockId = curBlockId;
                blockNum++;

                bitIdList[u].emplace_back(curBlockId);
                bitList[u].emplace_back(0);
            }
            //  origin v % WordBits
            bitList[u].back() |= (WordType)(1) << (v & ModBits);
        }
    }
}

uint32_t Graph::BinarySearch(uint32_t l, uint32_t r, const uint32_t value)
{
    uint32_t End = r;
    r = r - 1;
    while (l <= r)
    {
        uint32_t mid = (l + r) >> 1;
        if (edgeList[mid] < value)
            l = mid + 1;
        else
            r = mid;
        if (r - l <= 2)
            break;
    }
    while (edgeList[l] < value && l < End)
        l++;
    return l;
}

uint32_t Graph::SearchFirst(uint32_t u)
{
    uint32_t eId = edgeListBegPtr[u], endeId = edgeListEndPtr[u];

    if (edgeList[eId] > u)
        return eId;
    if ((eId + 1) < endeId && edgeList[eId + 1] > u)
        return eId + 1;
    if ((eId + 2) < endeId && edgeList[eId + 2] > u)
        return eId + 2;
    if (eId + 2 >= endeId)
        return eId;

    //  Binary Search would be slightly slower.
    uint32_t skip = 4;
    while (true)
    {
        uint32_t curId = eId + skip;
        if (curId >= endeId)
            return eId + (skip >> 1) + 1;
        //return BinarySearch(eId + (skip >> 1) + 1, endeId, u);

        if (edgeList[curId] < u)
            skip <<= 1;
        else
            return eId + (skip >> 1) + 1;
        //return BinarySearch(eId + (skip >> 1) + 1, eId + skip + 1, u);
    }
}

void Graph::CountTriangle()
{
    Triangles = 0;
#ifdef TRI_ON
    uint32_t ClrCnt = 0;

#endif
#ifdef TIMER_ON
    // double InterTime = 0.0;
    Timer CurTimer;
    CurTimer.StartTimer();
#endif

    Timer InitTimer;
    InitTimer.StartTimer();
#pragma omp parallel for schedule(dynamic)
    for (uint32_t u = 1; u <= nodeNum; u++)
        PackVertex(u);
    InitTimer.EndTimer();
    tcInitTime = InitTimer.QueryTimer();

#ifdef TIMER_ON
    CurTimer.GapTimer("PackVertex");
#endif

    uint32_t ThreadNum = omp_get_num_procs();
    MyBitSet Arr[ThreadNum];
    for (uint32_t i = 0; i < ThreadNum; i++)
        Arr[i].Init((nodeNum + 2) + 2);

#ifdef FILTER_ON
    HashFilter Filter[ThreadNum];
    for (uint32_t i = 0; i < ThreadNum; i++)
        Filter[i].Init(nodeNum + 2);
#endif

#ifdef TIMER_ON
    CurTimer.GapTimer("BitsetInit");
#endif

#pragma omp parallel for schedule(dynamic, 4096) reduction(+ \
                                                           : Triangles)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        if (edgeListEndPtr[u] == edgeListBegPtr[u])
            continue;

        uint32_t tId = omp_get_thread_num();
#ifdef FILTER_ON
        Filter[tId].Construct(edgeListEndPtr[u] - edgeListBegPtr[u]);
#endif
        for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
#ifdef ASSERT_ON
            assert(v >= 1 && v < (nodeNum + 2));
            assert(tId >= 0 && tId < ThreadNum);
#endif

            Arr[tId].Set(v);
#ifdef FILTER_ON
            Filter[tId].Set(v);
#endif
        }

#ifdef TIMER_ON
        Timer InterTimer;
        InterTimer.StartTimer();
#endif

        uint32_t eSupport = 0;

        uint32_t starteId = SearchFirst(u);
        for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            if (v < u)
                continue;

            eSupport = 0;
            if (__unlikely(!bitIdList[v].empty()))
            {
                for (uint32_t bId = 0; bId < bitIdList[v].size(); bId++)
                {
                    WordType res = Arr[tId].GetWord(bitIdList[v][bId]) & bitList[v][bId];
                    //  似乎比__builtin_popcountll 快
                    eSupport += __builtin_popcount(res & 0xFFFFFFFF) + __builtin_popcount(res >> 32);
                }
            }
            else
            {
                for (uint32_t vId = edgeListBegPtr[v]; vId < edgeListEndPtr[v]; vId++)
                {
                    uint32_t w = edgeList[vId];

#ifdef FILTER_ON
                    if (Filter[tId].Try(w))
#endif
                        if (Arr[tId].Get(w))
                            eSupport++;
                }
            }

            supports[eId] = eSupport;
            Triangles += eSupport;

#ifdef TRI_ON
            //#pragma omp atomic
            // Triangles += eSupport;
#endif
        }

#ifdef TIMER_ON
        InterTimer.EndTimer();
        // InterTime += InterTimer.QueryTimer();
#endif

//  在开头直接就可以进行初始化..
#if 1
#if 1
        //  针对稠密图的优化，对于SE数据有bug，bug fixed
        if (__unlikely(nodeNum < ((edgeListEndPtr[u] - edgeListBegPtr[u]) << 6)))
        {
#ifdef TRI_ON
            ClrCnt++;
#endif
            Arr[tId].Clr();
        }
        else
#endif
        {
            {
                uint32_t lastBlock = MINUS;
                for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
                {
                    uint32_t v = edgeList[eId];
                    uint32_t curBlock = v >> ShiftBits;
                    if (lastBlock != curBlock)
                    {
                        lastBlock = curBlock;
                        Arr[tId].SetWord(curBlock, 0);
                    }
                }
            }
        }
#endif
    }

    Triangles = Triangles / 3;
#ifdef TRI_ON
    // std::cout << "Inter Time : " << InterTime / (double)ThreadNum << std::endl;
    std::cout << "Bitset Clr Cnt : " << ClrCnt << std::endl;
    std::cout << "Triangle Found : " << Triangles << std::endl;
#endif
}

__always_inline void Graph::RemoveEdge(const uint32_t u, const uint32_t v, const uint32_t uvId)
{
    edgeList[uvId] = EDGEDELETE;
    uint32_t vuId = edgeRev[uvId];
    edgeList[vuId] = EDGEDELETE;
}

void Graph::FindInter(const uint32_t u, const uint32_t v, std::vector<std::pair<uint32_t, uint32_t>> &AffectEdgeList)
{
    uint32_t uPtr = edgeListBegPtr[u];
    uint32_t vPtr = edgeListBegPtr[v];
    AffectEdgeList.clear();
    while (uPtr < edgeListEndPtr[u] && vPtr < edgeListEndPtr[v])
    {
        uint32_t w1 = edgeList[uPtr];
        if (w1 == EDGEDELETE)
        {
            uPtr++;
            continue;
        }
        uint32_t w2 = edgeList[vPtr];
        if (w2 == EDGEDELETE)
        {
            vPtr++;
            continue;
        }
        if (w1 == w2)
        {
            uint32_t uwId = uPtr;
            uint32_t vwId = vPtr;
            if (w1 < u)
                uwId = edgeRev[uwId];
            if (w1 < v)
                vwId = edgeRev[vwId];
            AffectEdgeList.push_back(std::move(std::make_pair(uwId, vwId)));
            uPtr++;
            vPtr++;
            continue;
        }
        if (w1 < w2)
            uPtr++;
        else
            vPtr++;
    }
}

uint32_t Graph::CountInter(const uint32_t u, const uint32_t v)
{
    uint32_t cnt = 0;
    uint32_t uPtr = edgeListBegPtr[u];
    uint32_t vPtr = edgeListBegPtr[v];
    while (uPtr < edgeListEndPtr[u] && vPtr < edgeListEndPtr[v])
    {
        uint32_t w1 = edgeList[uPtr];
        uint32_t w2 = edgeList[vPtr];
        if (w1 == EDGEDELETE)
        {
            uPtr++;
            continue;
        }
        if (w2 == EDGEDELETE)
        {
            vPtr++;
            continue;
        }
        if (w1 == w2)
        {
            cnt++;
            uPtr++;
            vPtr++;
            continue;
        }
        if (w1 < w2)
            uPtr++;
        else
            vPtr++;
    }
    return cnt;
}

__always_inline void Graph::Update(const uint32_t eId, CacheBuffer &Buffer) noexcept
{
    uint32_t *addr = &supports[eId];
    if (*addr <= k - 2)
        return;
#ifdef ASSERT_ON
    assert(curStage[eId] == false);
#endif

    if (__sync_fetch_and_sub(addr, 1) == k - 1)
        Buffer.push(eId);
}
void Graph::PeelTriangle(const uint32_t uvId, const uint32_t uwId, const uint32_t vwId, CacheBuffer &Buffer)
{
    register bool curuw = curStage[uwId];
    register bool curvw = curStage[vwId];
    if (curuw && curvw)
        return;

    if (edgeList[uwId] != EDGEDELETE && edgeList[vwId] != EDGEDELETE)
    {
        if (!curuw && !curvw)
        {
            Update(uwId, Buffer);
            Update(vwId, Buffer);
        }
        else if (curuw && uvId < uwId)
            Update(vwId, Buffer);
        else if (uvId < vwId)
            Update(uwId, Buffer);
    }
}

void Graph::DeleteAll(FixedVector &QCur)
{
#pragma omp parallel for schedule(static) reduction(+ \
                                                    : realDel)
    for (uint32_t i = 0; i < QCur.size(); i++)
    {
        uint32_t uvId = QCur[i];
        if (edgeList[uvId] == EDGEDELETE)
            continue;
        //#pragma omp atomic
        realDel++;
        uint32_t u = (*ptrEdges)[edgeId[uvId]].u;
        uint32_t v = (*ptrEdges)[edgeId[uvId]].v;
        RemoveEdge(u, v, uvId);
    }

#ifdef DEBUG_ON
    // std::cout << "Delete Edge Num For " << k << " Real Delete: " << realDel << std::endl;
#endif
}
void Graph::PeelAll(FixedVector &QCur, FixedVector &QNext)
{

    uint32_t ThreadNum = omp_get_num_procs();
    CacheBuffer Buffer[ThreadNum];
    for (uint32_t i = 0; i < ThreadNum; i++)
    {
        Buffer[i].Init(BUFFER_SIZE);
        Buffer[i].Set(QNext.GetArrPtr(), QNext.GetSizePtr());
    }

#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < QCur.size(); i++)
        curStage[QCur[i]] = true;

#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < QCur.size(); i++)
    {
        uint32_t uvId = QCur[i];
        uint32_t u = (*ptrEdges)[edgeId[uvId]].u;
        uint32_t v = (*ptrEdges)[edgeId[uvId]].v;

        static thread_local std::vector<std::pair<uint32_t, uint32_t>> AffectEdgeList;
        FindInter(u, v, AffectEdgeList);

        uint32_t tId = omp_get_thread_num();
        for (uint32_t j = 0; j < AffectEdgeList.size(); j++)
        {
            uint32_t uwId, vwId;
            std::tie(uwId, vwId) = AffectEdgeList[j];
            PeelTriangle(uvId, uwId, vwId, Buffer[tId]);
        }
    }

    for (uint32_t i = 0; i < ThreadNum; i++)
        Buffer[i].pushBufferAnyway();

    DeleteAll(QCur);
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < QCur.size(); i++)
    {
        curStage[QCur[i]] = false;
    }
    QCur.swap(QNext);
    QNext.clear();
}

void Graph::CompactCSR()
{
#pragma omp parallel for schedule(static, 4096)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        if (edgeListBegPtr[u] == edgeListEndPtr[u])
            continue;
        uint32_t lasteId = edgeListBegPtr[u];
        for (uint32_t cureId = edgeListBegPtr[u]; cureId < edgeListEndPtr[u]; cureId++)
            if (edgeList[cureId] != EDGEDELETE)
            {
                if (lasteId != cureId)
                {
                    edgeList[lasteId] = edgeList[cureId];
                    edgeId[lasteId] = edgeId[cureId];
                    supports[lasteId] = supports[cureId];
                }
                LastToCur[lasteId] = cureId;
                CurToLast[cureId] = lasteId;

                lasteId++;

                // 我似乎没法解决 快速的 维护edgeRev数组，critical 太慢了，解决了..
                /*
#pragma omp critical
                {
                    edgeRev[edgeRev[cureId]] = lasteId;
                    edgeRev[lasteId] = edgeRev[cureId];
                }
                */
            }
        edgeListEndPtr[u] = lasteId;
    }

#ifdef REV_ON

#if 1
// Cache Coherence
#pragma omp parallel for schedule(static)
    for (uint32_t u = 1; u <= nodeNum; u++)
        for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
            edgeRev[eId] = CurToLast[edgeRev[LastToCur[eId]]];
#endif
#endif
}
void Graph::ShrinkCSR()
{
    return;

    uint32_t *compactsupports = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
    uint32_t *compactList = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
    uint32_t *compactId = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
    uint32_t *compactRev = (uint32_t *)malloc(sizeof(uint32_t) * (edgeNum << 1 | 1));
    uint32_t *compactBegPtr = (uint32_t *)malloc(sizeof(uint32_t) * (2 * nodeNum + 2));

#ifdef DEBUG_ON
    if (compactsupports == NULL || compactList == NULL || compactId == NULL || compactRev == NULL || compactBegPtr == NULL)
    {
        std::cout << "Malloc Error " << std::endl;
        exit(-1);
    }
#endif

    memset(compactBegPtr, 0, sizeof(uint32_t) * (2 * nodeNum + 2));
    for (uint32_t u = 1; u <= nodeNum; u++)
        compactBegPtr[u + 1] = edgeListEndPtr[u] - edgeListBegPtr[u];
    for (uint32_t u = 2; u <= nodeNum; u++)
        compactBegPtr[u + 1] += compactBegPtr[u];
    compactBegPtr[1] = 0;

    // #pragma omp parallel for schedule(static)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        uint32_t len = edgeListEndPtr[u] - edgeListBegPtr[u];
        if (len == 0)
            continue;
        if (compactBegPtr[u] > (edgeNum << 1 | 1))
        {
            puts("Error memcpy");
            exit(-1);
        }
        memcpy(compactList + compactBegPtr[u], edgeList + edgeListBegPtr[u], sizeof(uint32_t) * len);
        memcpy(compactId + compactBegPtr[u], edgeId + edgeListBegPtr[u], sizeof(uint32_t) * len);
        memcpy(compactRev + compactBegPtr[u], edgeRev + edgeListBegPtr[u], sizeof(uint32_t) * len);
        memcpy(compactsupports + compactBegPtr[u], supports + edgeListBegPtr[u], sizeof(uint32_t) * len);
        edgeListEndPtr[u] = compactBegPtr[u + 1];
    }
#if 1
    // #pragma omp parallel for schedule(static)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        std::cout << "u : " << u << std::endl;
        if (compactBegPtr[u] == edgeListBegPtr[u])
            continue;

        for (uint32_t eId = compactBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = compactList[eId];
            if (v == EDGEDELETE)
                continue;
            uint32_t cureId = edgeListBegPtr[u] + eId - compactBegPtr[u];
            compactRev[eId] = edgeRev[cureId] - (edgeListBegPtr[v] - compactBegPtr[v]);
            // uint32_t origin = compactRev[eId];
            //compactRev[eId] -= (edgeListBegPtr[v] - compactBegPtr[v]);
            if (!(compactRev[eId] <= (edgeNum << 1 | 1)))
            {
                std::cout << edgeListBegPtr[v] << " " << compactBegPtr[v] << std::endl;
                uint32_t cureId = edgeListBegPtr[u] + eId - compactBegPtr[u];

                std::cout << cureId << " " << edgeRev[cureId] << " " << edgeList[cureId] << " " << edgeList[edgeRev[cureId]] << std::endl;
                // std::cout << u << " " << v << " " << origin << " " << compactList[origin] << " " << edgeListBegPtr[v] - compactBegPtr[v] << std::endl;

                std::cout << cureId << " " << edgeRev[cureId] << " " << edgeList[cureId] << " " << edgeList[edgeRev[cureId]] << std::endl;
                //std::cout << u << " " << v << " " << compactRev[eId] << " " << compactList[origin] << " " << edgeListBegPtr[v] - compactBegPtr[v] << std::endl;

                puts("Error !!");
                exit(-1);
            }
        }
    }
#endif
    puts("!!!");
    for (uint32_t u = 1; u <= nodeNum; u++)
        for (uint32_t eId = compactBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = compactList[eId];
            if (compactList[compactRev[eId]] != u)
            {
                std::cout << u << " " << v << " " << compactList[compactRev[eId]] << std::endl;
                puts("Wrong Shrink");
                exit(-1);
            }
        }

    free(edgeList);
    free(edgeId);
    free(edgeRev);
    free(supports);
    free(edgeListBegPtr);

    edgeList = compactList;
    edgeId = compactId;
    edgeRev = compactRev;
    supports = compactsupports;
    edgeListBegPtr = compactBegPtr;
}

void Graph::ReCountTriangle(FixedVector &QCur, FixedVector &QNext)
{
    DeleteAll(QCur);
    CompactCSR();

    Triangles = 0;
#pragma omp parallel for schedule(dynamic)
    for (uint32_t u = 1; u <= nodeNum; u++)
        PackVertex(u);

    uint32_t ThreadNum = omp_get_num_procs();
    CacheBuffer Buffer[ThreadNum];
    for (uint32_t i = 0; i < ThreadNum; i++)
    {
        Buffer[i].Init(BUFFER_SIZE);
        Buffer[i].Set(QNext.GetArrPtr(), QNext.GetSizePtr());
    }

    MyBitSet Arr[ThreadNum];
    for (uint32_t i = 0; i < ThreadNum; i++)
        Arr[i].Init(nodeNum + 2);

#pragma omp parallel for schedule(dynamic) reduction(+ \
                                                     : Triangles)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        uint32_t tId = omp_get_thread_num();
        for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            Arr[tId].Set(v);
        }

        register uint32_t eSupport = 0;
        uint32_t starteId = SearchFirst(u);
        for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            if (v < u)
                continue;

            eSupport = 0;

            if (__unlikely(!bitIdList[v].empty()))
            {
                for (uint32_t bId = 0; bId < bitIdList[v].size(); bId++)
                {
                    WordType res = Arr[tId].GetWord(bitIdList[v][bId]) & bitList[v][bId];
                    eSupport += __builtin_popcount(res & 0xFFFFFFFF) + __builtin_popcount(res >> 32);
                }
            }
            else
            {
                for (uint32_t vId = edgeListBegPtr[v]; vId < edgeListEndPtr[v]; vId++)
                {
                    uint32_t w = edgeList[vId];
                    if (Arr[tId].Get(w))
                        eSupport++;
                }
            }

            supports[eId] = eSupport;

            Triangles = Triangles + eSupport;
            if (eSupport <= k - 2)
            {
                Buffer[tId].push(eId);
                // QNext.push_back(eId);
            }
        }

        {
            uint32_t lastBlock = MINUS;
            for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
            {
                uint32_t v = edgeList[eId];
                uint32_t curBlock = v >> ShiftBits;
                if (lastBlock != curBlock)
                {
                    lastBlock = curBlock;
                    Arr[tId].SetWord(curBlock, 0);
                }
            }
        }
    }

    for (uint32_t i = 0; i < ThreadNum; i++)
        Buffer[i].pushBufferAnyway();

    Triangles = Triangles / 3;

    QCur.swap(QNext);
    QNext.clear();
}
void Graph::NewTC()
{
    Triangles = 0;
#pragma omp parallel for schedule(dynamic)
    for (uint32_t u = 1; u <= nodeNum; u++)
        PackVertex(u);

    uint32_t ThreadNum = omp_get_num_procs();

    MyBitSet Arr[ThreadNum];
    for (uint32_t i = 0; i < ThreadNum; i++)
        Arr[i].Init(nodeNum + 2);

#pragma omp parallel for schedule(dynamic) reduction(+ \
                                                     : Triangles)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        uint32_t tId = omp_get_thread_num();
        for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            Arr[tId].Set(v);
        }

        register uint32_t eSupport = 0;
        uint32_t starteId = SearchFirst(u);
        for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            if (v < u)
                continue;

            eSupport = 0;

            if (__unlikely(!bitIdList[v].empty()))
            {
                for (uint32_t bId = 0; bId < bitIdList[v].size(); bId++)
                {
                    WordType res = Arr[tId].GetWord(bitIdList[v][bId]) & bitList[v][bId];
                    eSupport += __builtin_popcount(res & 0xFFFFFFFF) + __builtin_popcount(res >> 32);
                }
            }
            else
            {
                for (uint32_t vId = edgeListBegPtr[v]; vId < edgeListEndPtr[v]; vId++)
                {
                    uint32_t w = edgeList[vId];
                    if (Arr[tId].Get(w))
                        eSupport++;
                }
            }

            supports[eId] = eSupport;

            Triangles = Triangles + eSupport;
        }

        uint32_t lastBlock = MINUS;
        for (uint32_t eId = edgeListBegPtr[u]; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            uint32_t curBlock = v >> ShiftBits;
            if (lastBlock != curBlock)
            {
                lastBlock = curBlock;
                Arr[tId].SetWord(curBlock, 0);
            }
        }
    }

    Triangles = Triangles / 3;
}
uint32_t Graph::GetRevId(const uint32_t u, const uint32_t v)
{
    uint32_t l = edgeListBegPtr[v], r = edgeListEndPtr[v] - 1;
    if (edgeList[l] == u)
        return l;
    if ((l + 1) <= r && edgeList[l + 1] == u)
        return l + 1;
    if ((l + 2) <= r && edgeList[l + 2] == u)
        return l + 2;

    uint32_t skip = 4;
    while (true)
    {
        uint32_t curId = l + skip;
        if (curId > v)
        {
            l = l + (skip >> 1) + 1;
            break;
        }

        if (edgeList[curId] == u)
            return curId;
        if (edgeList[curId] < u)
            skip <<= 1;
        else
        {

            l = l + (skip >> 1) + 1;
            r = curId;
        }
    }
    while (l <= r)
    {
        uint32_t mid = (l + r) >> 1;
        if (edgeList[mid] == u)
            return mid;

        if (edgeList[mid] < u)
            l = mid + 1;
        else
            r = mid - 1;
    }
    return l;
}
void Graph::InitRev()
{
#pragma omp parallel for schedule(static)
    for (uint32_t u = 1; u <= nodeNum; u++)
    {
        uint32_t starteId = SearchFirst(u);
        for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
        {
            uint32_t v = edgeList[eId];
            if (v < u)
                continue;
            edgeRev[eId] = GetRevId(u, v);
            edgeRev[edgeRev[eId]] = eId;
        }
    }
}
void Graph::TrussDecomp()
{
    k = 1;
    FixedVector QCur(edgeNum << 1 | 1);
    FixedVector QNext(edgeNum << 1 | 1);

    realDel = 0;

    while (true)
    {
        k++;
        QCur.clear();
        QNext.clear();
        // InitRev();

#pragma omp parallel for schedule(static)
        for (uint32_t u = 1; u <= nodeNum; u++)
        {
            //   Assume Compact
            uint32_t starteId = SearchFirst(u);
            for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
            {
                uint32_t v = edgeList[eId];
                if (v < u)
                    continue;
                if (supports[eId] <= k - 2)
                {
                    QCur.push_back(eId);
#ifdef ASSERT_ON
                    assert(curStage[eId] == false);
#endif
                    curStage[eId] = true;
                }
            }
        }

        if (k == 2)
        {
            DeleteAll(QCur);
#pragma omp parallel for schedule(static)
            for (uint32_t i = 0; i < QCur.size(); i++)
            {
                curStage[QCur[i]] = false;
            }

            CompactCSR();
            continue;
        }

        while (QCur.size() > 0)
        {

            // double tcTime = (double) (edgeNum - realDel) / edgeNum * tcInitTime + tcBiasTime;
            // double peelTime = (double)(QCur.size() + 1) * (k + 1) / (double)2e7;
            PeelAll(QCur, QNext);
        }

        if (realDel == edgeNum)
            break;

        CompactCSR();
    }

    std::cout << k << std::endl;
}

void Graph::FindKMax()
{
#ifdef TIMER_ON
    Timer CurTimer;
    double CompactTime = 0.0;
    double PeelTime = 0.0;
    CurTimer.StartTimer();
#endif

    k = Triangles / edgeNum + 2 - 1;

#ifdef DEBUG_ON
    std::cout << "Initial K : " << k << std::endl;
#endif
    FixedVector QCur(edgeNum + 2);
    FixedVector QNext(edgeNum + 2);

    uint32_t ThreadNum = omp_get_num_procs();
    CacheBuffer Buffer[ThreadNum];

    for (uint32_t i = 0; i < ThreadNum; i++)
        Buffer[i].Init(BUFFER_SIZE);

    uint32_t lastDel = 0;
    realDel = 0;

    bool Shrinked = false;

    while (true)
    {
    REST_TRUSS:
        k++;
        QCur.clear();
        QNext.clear();

        for (uint32_t i = 0; i < ThreadNum; i++)
            Buffer[i].Set(QCur.GetArrPtr(), QCur.GetSizePtr());

#pragma omp parallel for schedule(static)
        for (uint32_t u = 1; u <= nodeNum; u++)
        {
            uint32_t tId = omp_get_thread_num();
            //   Assume Compact
            uint32_t starteId = SearchFirst(u);
            for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
            {
                uint32_t v = edgeList[eId];
                if (v < u)
                    continue;
                if (supports[eId] <= k - 2)
                {
                    Buffer[tId].push(eId);
                    // QCur.push_back(eId);
                }
            }
        }

        for (uint32_t i = 0; i < ThreadNum; i++)
            Buffer[i].pushBufferAnyway();

        if (k == 2)
        {
            DeleteAll(QCur);
#pragma omp parallel for schedule(static)
            for (uint32_t i = 0; i < QCur.size(); i++)
            {
                curStage[QCur[i]] = false;
            }
            CompactCSR();
            continue;
        }

        while (QCur.size() > 0)
        {

#if 1
            while (__unlikely(QCur.size() > (edgeNum >> 6) && QCur.size() > 100024))
            {
                ReCountTriangle(QCur, QNext);
                if (k > Triangles / (edgeNum - realDel - QCur.size()) + 2)
                {
                    k = Triangles / (edgeNum - realDel - QCur.size()) + 2;
                    // DeleteAll(QCur);
#ifdef TIMER_ON
                    CurTimer.GapTimer("ReCount");
#endif
                    goto REST_TRUSS;
                }

#ifdef TIMER_ON
                CurTimer.GapTimer("ReCount");
#endif
                goto REST_TRUSS;
            }
#endif
            {
#ifdef TIMER_ON
                Timer PeelTimer;
                PeelTimer.StartTimer();
#endif
                PeelAll(QCur, QNext);
#ifdef TIMER_ON
                PeelTimer.EndTimer();
                PeelTime += PeelTimer.QueryTimer();
#endif
            }
        }

        if (__unlikely(realDel == edgeNum))
            break;
        lastDel = realDel;
#ifdef TIMER_ON
        Timer CompactTimer;
        CompactTimer.StartTimer();
#endif
        CompactCSR();
#ifdef TIMER_ON
        CompactTimer.EndTimer();
        CompactTime += CompactTimer.QueryTimer();
#endif

//  Shrink还有Bug不知道怎么解决
#if 0
        if (!Shrinked && realDel > edgeNum >> 1)
        {
#ifdef TIMER_ON
            Timer ShrinkTimer;
            ShrinkTimer.StartTimer();
#endif
            puts("Shrink !!!");
            ShrinkCSR();
            Shrinked = true;
#ifdef TIMER_ON
            ShrinkTimer.EndTimer();
            double ShrinkTime = ShrinkTimer.QueryTimer();
            std::cout << "Shrink Time : " << ShrinkTime << std::endl;
#endif
        }
#endif
    }
#ifdef TIMER_ON
    std::cout << "Compact Time : " << CompactTime << std::endl;
    std::cout << "Peel Time : " << PeelTime << std::endl;
#endif
    printf("kmax = %d, Edges in kmax-truss = %d.\n", k, realDel - lastDel);
    // std::cout << k << " " << realDel - lastDel << std::endl;
}

void Graph::MordernTruss()
{
#ifdef TIMER_ON
    Timer CurTimer;
    CurTimer.StartTimer();
#endif
    uint32_t lastDel = 0;
    realDel = 0;

    FixedVector QCur(edgeNum + 2);
    FixedVector QNext(edgeNum + 2);
    QCur.clear();
    QNext.clear();

    uint32_t ThreadNum = omp_get_num_procs();
    CacheBuffer Buffer[ThreadNum];

    for (uint32_t i = 0; i < ThreadNum; i++)
        Buffer[i].Init(BUFFER_SIZE);

    k = 2;
#if 0
    while (true)
    {
        std::cout << "k : " << k << " " << realDel << std::endl;
        NewTC();
        // ReCountTriangle(QCur, QNext);
        uint32_t bound = Triangles / (edgeNum - realDel - QCur.size()) + 2;
        if (k < bound)
            k = bound;
        else
            break;

        for (uint32_t i = 0; i < ThreadNum; i++)
            Buffer[i].Set(QCur.GetArrPtr(), QCur.GetSizePtr());
#pragma omp parallel for schedule(static)
        for (uint32_t u = 1; u <= nodeNum; u++)
        {
            uint32_t tId = omp_get_thread_num();
            //   Assume Compact
            uint32_t starteId = SearchFirst(u);
            for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
            {
                uint32_t v = edgeList[eId];
                if (v < u)
                    continue;
                if (supports[eId] <= k - 2)
                    Buffer[tId].push(eId);
            }
        }
        for (uint32_t i = 0; i < ThreadNum; i++)
            Buffer[i].pushBufferAnyway();
        DeleteAll(QCur);
        CompactCSR();
    }

#ifdef TIMER_ON
    CurTimer.GapTimer("TC Reduce");
    std::cout << "TC k : " << k << std::endl;
#endif
#endif

    while (true)
    {
        // std::cout << k << " " << realDel << " " << edgeNum << std::endl;
        k++;
        QCur.clear();
        QNext.clear();

        for (uint32_t i = 0; i < ThreadNum; i++)
            Buffer[i].Set(QCur.GetArrPtr(), QCur.GetSizePtr());
#pragma omp parallel for schedule(static)
        for (uint32_t u = 1; u <= nodeNum; u++)
        {
            uint32_t tId = omp_get_thread_num();
            //   Assume Compact
            uint32_t starteId = SearchFirst(u);
            for (uint32_t eId = starteId; eId < edgeListEndPtr[u]; eId++)
            {
                uint32_t v = edgeList[eId];
                if (v < u)
                    continue;
                if (supports[eId] <= k - 2)
                    Buffer[tId].push(eId);
            }
        }
        for (uint32_t i = 0; i < ThreadNum; i++)
            Buffer[i].pushBufferAnyway();

        while (QCur.size() > 0)
            PeelAll(QCur, QNext);
        if (__unlikely(realDel == edgeNum))
            break;
        lastDel = realDel;
        CompactCSR();
    }

#ifdef TIMER_ON
    CurTimer.GapTimer("Peel Reduce");
#endif
    printf("kmax = %d, Edges in kmax-truss = %d.\n", k, realDel - lastDel);
}
