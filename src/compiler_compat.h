#pragma once

#ifndef _MSC_VER
#include <cstring>
#include <stdio.h>
#define _stricmp strcasecmp
#define _strnicmp strncasecmp
#endif

#ifndef _WIN32
#include <stdlib.h>
#define __forceinline inline
#ifndef __cdecl
#define __cdecl
#endif
#define _InterlockedCompareExchangePointer(a,b,c) __sync_val_compare_and_swap(a,c,b)

static inline void* _aligned_malloc(size_t size, size_t alignment)
{
    void *tmp;
    if (posix_memalign(&tmp, alignment, size))
    {
        tmp = 0;
    }
    return tmp;
}
#define _aligned_free free
#else
#include <intrin.h>
    // ICL complains about unresolved external symbol
    #if __INTEL_COMPILER && !_WIN64
    __forceinline void* _InterlockedCompareExchangePointer(
        void* volatile *Destination, void* Exchange, void* Comperand) {
    return (void*) _InterlockedCompareExchange((long volatile *) Destination, (long) Exchange, (long) Comperand);
    }
    #endif
#endif


#define ALIGNED_ARRAY(type, decl, alignment) alignas(alignment) type decl
