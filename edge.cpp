#include "edge.h"
#include "tools.h"
#include <cassert>

static bool IS_EOF = false;
static char buf[1 << 24], *fs, *fe;
__always_inline char getc()
{

    return __unlikely((fs == fe && (fe = (fs = buf) + fread(buf, 1, 1 << 24, stdin)), fs == fe)) ? IS_EOF = true : *fs++;
}

__always_inline uint32_t Get_Int()
{
    uint32_t x = 0;
    register char c;
    c = getc();
    while ((c < '0' || c > '9') && !IS_EOF)
        c = getc();
    while (c >= '0' && c <= '9')
    {
        x = x * 10 + c - '0';
        c = getc();
    }
    return x;
}

__always_inline uint32_t Skip_W()
{
    getc();
    getc();
}
void ReadBaseLine(const char *path, std::vector<Edge> &Edges)
{
    FILE *_ = freopen(path, "rb", stdin);

    struct Edge curEdge;

    while (true)
    {
        curEdge.u = Get_Int();
        curEdge.v = Get_Int();
        Skip_W();
        if (__unlikely(IS_EOF))
            break;
        Edges.push_back(std::move(curEdge));
    }

    fclose(stdin);
#if 0
    FILE *fp;
    fp = fopen(path, "r");
    Edges.clear();

    while (true)
    {
        uint32_t _ = fscanf(fp, "%u %u %u", &curEdge.u, &curEdge.v, &w);
        Edges.push_back(std::move(curEdge));
        if (feof(fp))
            break;
    }
    fclose(fp);
#endif
}

void LoadEdge(const char *const bptr, const uint32_t len, std::vector<Edge> &Edges)
{
    const char *cptr = bptr;
    const char *const eptr = bptr + len;
    Edges.clear();

    while (cptr < eptr)
    {

        struct Edge curEdge;
        curEdge.u = mysscanf(eptr, cptr);
        curEdge.v = mysscanf(eptr, cptr);
        if (cptr >= eptr)
            break;

        //  自环判断
        if (curEdge.u == curEdge.v)
            continue;

        Edges.push_back(curEdge);

        sscanfSkip(eptr, cptr);
    }
}

uint32_t mysscanf(const char *const eptr, const char *&cptr)
{
    uint32_t ret = 0;
    bool readF = false;
    while (cptr < eptr)
    {
        const char ch = *cptr;
        if (ch >= '0' && ch <= '9')
        {
            ret = ret * 10 + (ch - '0');
            readF = true;
        }
        else
        {
            if (readF)
                return ret;
        }
        cptr++;
    }
    if (readF)
        return ret;
}

void sscanfSkip(const char *const eptr, const char *&cptr)
{
    while (*cptr != '1' && cptr < eptr)
        cptr++;
    cptr++;
}
