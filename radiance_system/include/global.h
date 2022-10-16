/* global.h: defining global compile settings */
/* should be included by all "base*.h" files. */
#pragma once
#define NULL 0

/* define debug */
#if defined(_MSC_VER) 
#ifdef _DEBUG
#define DEBUG
#endif
/* other platforms */
/* ... */
#endif

/* * * * * * * * * * * * * * * * * * * * */
/*        GLOBAL COMPILE SETTINGS        */
/* * * * * * * * * * * * * * * * * * * * */
/* use single/double precision floats */
#define USE_SINGLE_PRECISION // USE_DOUBLE_PRECISION
#define USE_SAFE_DIVISION



