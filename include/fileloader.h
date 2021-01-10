#ifndef FILELOADER
#define FILELOADER
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>  
#include <sys/socket.h>  
#include <netinet/in.h>  
#include <arpa/inet.h>  
#include <unistd.h>  
class fileLoader
{
private:
    void *addr;
    int fd;
    int len;

public:
    fileLoader()
    {
        fd = -1;
        len = 0;
        addr = NULL;
    }
    ~fileLoader()
    {
    }
    void release()
    {
        if (fd != -1)
        {
            close(fd);
            munmap(addr, len);
        }
    }
    void loadFile(const char *path)
    {
        struct stat statbuf;
        fd = open(path, O_RDONLY);
        if (fd < 0)
        {
            puts("can not open file");
            exit(-1);
        }
        int ret = fstat(fd, &statbuf);
        if (ret < 0)
        {
            puts("can not fstat");
            exit(-1);
        }

        len = statbuf.st_size;
        addr = mmap(0, len, PROT_READ, MAP_SHARED, fd, 0);
        if (addr == (void *)-1)
        {
            puts("can not mmap");
            exit(-1);
        }
    }

    int getLen()
    {
        return len;
    }

    void *getAddr()
    {
        return addr;
    }
};

#endif