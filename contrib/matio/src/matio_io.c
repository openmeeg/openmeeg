#include "matio_io.h"

MatioOpenFunction*  default_matio_open  = (MatioOpenFunction*) fopen;
MatioSeekFunction*  default_matio_seek  = (MatioSeekFunction*) fseek;
MatioReadFunction*  default_matio_read  = (MatioReadFunction*) fread;
MatioWriteFunction* default_matio_write = (MatioWriteFunction*) fwrite;
MatioCloseFunction* default_matio_close = (MatioCloseFunction*) fclose;

inline void* matio_open(const char* path,const char* mode) {
    return default_matio_open(path,mode);
}

inline int matio_seek(void* stream,long offset,int whence) {
    return default_matio_seek(stream,offset,whence);
}

inline size_t matio_read(void* ptr,size_t size,size_t nmemb,void* stream) {
    return default_matio_read(ptr,size,nmemb,stream);
}

inline size_t matio_write(const void* ptr,size_t size,size_t nmemb,void* stream) {
    return default_matio_write(ptr,size,nmemb,stream);
}

inline int matio_close(void* fp) {
    return default_matio_close(fp);
}

