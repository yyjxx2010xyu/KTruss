#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <map>
#include <algorithm>
#include "timer.h"
#include "fileloader.h"
#include "edge.h"
#include "graph.h"

void par_memcpy(
    void *const dest,
    void const *const src,
    size_t const bytes)
{
#pragma omp parallel
    {
        int const nthreads = omp_get_num_threads();
        int const tid = omp_get_thread_num();

        /* block distribution */
        size_t const n_per_thread = (bytes + nthreads - 1) / nthreads;
        size_t const n_begin = std::min(n_per_thread * tid, bytes);
        size_t const n_end = std::min(n_begin + n_per_thread, bytes);

        char *const __restrict__ dest_data = (char *)dest;
        char *const __restrict__ src_data = (char *)src;
        memcpy(dest_data + n_begin, src_data + n_begin, n_end - n_begin);
    }
}

uint32_t CntDup(std::vector<Edge> &Edges)
{
    std::map<uint64_t, int> M;
    M.clear();
    for (uint32_t i = 0; i < Edges.size(); i++)
    {
        uint64_t e = ((uint64_t)Edges[i].u << 32) | (Edges[i].v);
        if (Edges[i].u > Edges[i].v)
            e = ((uint64_t)Edges[i].v << 32) | (Edges[i].u);
        M[e] = 1;
    }
    return M.size();
}

int main(int argc, char *argv[])
{
    if (argc != 3 || strcmp(argv[1], "-f") != 0)
    {
        printf("Usage: -f [data_file_path]\n");
        exit(1);
    }
    const char *filepath = argv[2];

    omp_set_num_threads(ThreadNum);

#if 0
    Timer CurTimer;
    CurTimer.StartTimer();
#endif
#ifdef TIMER_ON
    Timer CurTimer;
    CurTimer.StartTimer();
#endif

    // std::cout << "Core Num: " << omp_get_num_procs() << std::endl;
#ifdef DEBUG_ON
    std::cout << "Core Num: " << omp_get_num_procs() << std::endl;
#endif

#if 0
    freopen(filepath, "rb", stdin);
    struct Edge curEdge;
    uint32_t w;
    while (true)
    {
        
        
        Edges.push_back(std::move(curEdge));
        if (feof(fp))
            break;
    }
    fclose(fp);
#endif
#if 1
    fileLoader loader;

    // loader.loadFile("./data/truss.tsv");
    // loader.loadFile("./data/smalldata.tsv");
    // loader.loadFile("./data/tri.tsv");
    // loader.loadFile("./data/s18.e16.rmat.edgelist.tsv");
    // loader.loadFile("./data/s19.e16.rmat.edgelist.tsv");
    // loader.loadFile("./data/cit-Patents.tsv");
    // loader.loadFile("./data/soc-LiveJournal.tsv");

    loader.loadFile(filepath);

    std::vector<struct Edge> Edges;
    const char *const addr = (const char *const)loader.getAddr();
    const int len = loader.getLen();
    LoadEdge(addr, len, Edges);
    loader.release();
#endif
#if 0
    std::vector<struct Edge> Edges;
    ReadBaseLine(filepath, Edges);
#endif
#if 0
        fileLoader loader;
        // loader.loadFile("./data/tri.tsv");
        // loader.loadFile("./data/s18.e16.rmat.edgelist.tsv");
        //loader.loadFile("./data/s19.e16.rmat.edgelist.tsv");
        // loader.loadFile("./data/cit-Patents.tsv");
        loader.loadFile("./data/soc-LiveJournal.tsv");

        std::vector<struct Edge> Edges;
        const char *const addr = (const char *const)loader.getAddr();
        const int len = loader.getLen();


        Edge * edgeList = (Edge *) (new char[len]);
        par_memcpy(edgeList, addr, len);
        free(edgeList);
#endif
#ifdef TIMER_ON
    CurTimer.GapTimer("Load File");
#endif
    uint64_t edgeNum = Edges.size();

#ifdef DEBUG_ON
    std::cout << "-----------------------------" << std::endl;
    std::cout << "data preview: " << std::endl;
    std::cout << std::setw(8) << "u" << std::setw(8) << "v" << std::endl;
    for (int i = 0; i < 5; i++)
        std::cout << std::setw(8) << Edges[i].u << std::setw(8) << Edges[i].v << std::endl;
    std::cout << "-----------------------------" << std::endl;
#endif

    uint32_t nodeNum = 0;

#ifdef TIMER_ON
    CurTimer.UpdateLastTimer();
#endif

#pragma omp parallel for schedule(static) reduction(max \
                                                    : nodeNum)
    for (uint32_t i = 0; i < edgeNum; i++)
    {
        struct Edge &e = Edges[i];
        nodeNum = std::max(nodeNum, e.u);
        nodeNum = std::max(nodeNum, e.v);
        if (e.v < e.u)
            std::swap(e.u, e.v);
    }
#ifdef DEBUG_ON
    std::cout << "node num: " << nodeNum << std::endl;
    std::cout << "edge num: " << edgeNum << std::endl;
#endif
#ifdef ASSERT_ON
    // std::cout << "Dup Edge Num: " << CntDup(Edges) << std::endl;
#endif

#ifdef TIMER_ON
    CurTimer.GapTimer("FindMax");
#endif

    Graph grf(nodeNum, edgeNum);
#ifdef TIMER_ON
    CurTimer.GapTimer("MallocGraph");
#endif

    grf.LoadEdge(Edges);

#ifdef TIMER_ON
    CurTimer.GapTimer("LoadEdge");
#endif

    grf.SortEdge();

#ifdef TIMER_ON
    CurTimer.GapTimer("SortEdge");
#endif

    grf.DeleteDup();

#ifdef TIMER_ON
    CurTimer.GapTimer("DeleteDup");
#endif

    grf.CountTriangle();
#ifdef TIMER_ON
    CurTimer.GapTimer("CountTriangle");
#endif

    // return 0;
#if 0
    grf.TrussDecomp();
#ifdef TIMER_ON
    CurTimer.GapTimer("TrussDecomp");
#endif
#endif

#if 0
    grf.FindKMax();
#ifdef TIMER_ON
    CurTimer.GapTimer("FindKMax");
#endif
#endif

#if 1
    grf.MordernTruss();
#ifdef TIMER_ON
    CurTimer.GapTimer("MordernTruss");
#endif
#endif

#if 0
    CurTimer.EndTimer();
    std::cout << "Total Time Cost: " << CurTimer.QueryTimer() << std::endl;
#endif

#ifdef TIMER_ON
    CurTimer.EndTimer();
    std::cout << "Total Time Cost: " << CurTimer.QueryTimer() << std::endl;
#endif
#
    return 0;
}