#ifndef ENHANCE_H_
#define ENHANCE_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "typedef.h"
#include "basic_op.h"
#include "oper_32b.h"

typedef short Word16;
typedef long  Word32;
typedef unsigned short UWord16;
typedef unsigned long  UWord32;
#define MUL_16_16(a, b) \
    ((int32_t) (((int16_t)(a)) * ((int16_t)(b))))
#define MUL_16_16_RSFT(a, b, c) \
    (MUL_16_16(a, b) >> (c))
#define DIV(a, b) \
    ((int32_t) ((int32_t)(a) / (int32_t)(b)))
#define UDIV(a, b) \
    ((uint32_t) ((uint32_t)(a) / (uint32_t)(b)))

typedef signed char         int8_t;
typedef signed short        int16_t;
typedef signed int          int32_t;
typedef __int64             int64_t;
typedef unsigned char       uint8_t;
typedef unsigned short      uint16_t;
typedef unsigned int        uint32_t;
typedef unsigned __int64    uint64_t;

#define	TRUE			1
#define	FALSE			0
#define	PI				3.1415926535897932384626433832795
#define FS				8000
#define FRM_LEN			160
#define	NUM_STAGE		7
#define	FFT_LEN			(2<<NUM_STAGE)
#define	DELAY			24
#define	EMP_FAC_Q11	    1638


#endif
