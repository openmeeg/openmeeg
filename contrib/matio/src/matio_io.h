#ifndef MATIO_IO_H
#define MATIO_IO_H

#include <stdio.h>

#if defined (WIN32)
    #define inline
#endif

typedef void*  MatioOpenFunction(const char* path,const char* mode);
typedef int    MatioSeekFunction(void* stream,long offset,int whence);
typedef size_t MatioReadFunction(void* ptr,size_t size,size_t nmemb, void* stream);
typedef size_t MatioWriteFunction(const void* ptr,size_t size,size_t nmemb,void* stream);
typedef int    MatioCloseFunction(void* fp);

extern MatioOpenFunction*  default_matio_open;
extern MatioSeekFunction*  default_matio_seek;
extern MatioReadFunction*  default_matio_read;
extern MatioWriteFunction* default_matio_write;
extern MatioCloseFunction* default_matio_close;

extern inline void* matio_open(const char* path,const char* mode);
extern inline int matio_seek(void* stream,long offset,int whence);
extern inline size_t matio_read(void* ptr,size_t size,size_t nmemb,void* stream);
extern inline size_t matio_write(const void* ptr,size_t size,size_t nmemb,void* stream);
extern inline int matio_close(void* fp);

#endif  //  ! MATIO_IO_H
