#ifndef SCREAM_GOCART_PREPROCESSOR_DEFINES_H
#define SCREAM_GOCART_PREPROCESSOR_DEFINES_H

/*-------------------------------------------------------------------------------------------------
 * This Header is for preprocessor #define directives used by the stubs from the 
 * {E3SM_ROOT}/share/util directory
 * so that they do not need to be individually set
 *-----------------------------------------------------------------------------------------------*/

/*  preprocessor variables may need to be set manually  */
#define HAVE_IEEE_ARITHMETIC
 
/* CPRGNU is needed for nan testing if HAVE_IEEE_ARITHMETIC is NOT defined */
#define CPRGNU 

#endif //SCREAM_GOCART_PREPROCESSOR_DEFINES_H