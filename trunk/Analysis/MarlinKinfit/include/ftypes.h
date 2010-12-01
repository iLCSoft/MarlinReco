#ifndef H_FTYPES
#define H_FTYPES
// ftypes.h
// typedef's to be used for interfacing FORTRAN routines to C++
typedef   short       FInteger2;      // INTEGER*2
typedef   int         FInteger;       // INTEGER
typedef   long long   FInteger8;      // INTEGER*8
typedef   float       FReal;          // REAL
typedef   double      FReal8;         // REAL*8
typedef   long double FReal16;        // REAL*16
typedef   double      FDouble;        // DOUBLE PRECISION
typedef   char        FCharacter;     // CHARACTER
typedef   int	      FLogical;	      // LOGICAL
#define LOGICAL_IS_INTEGER
typedef   short	      FLogical2;      // LOGICAL*2
#define LOGICAL2_IS_INTEGER2
typedef   bool	      FLogical1;      // LOGICAL*1
typedef   int         FStrLength;     // Type for passing string lengths

#endif   
