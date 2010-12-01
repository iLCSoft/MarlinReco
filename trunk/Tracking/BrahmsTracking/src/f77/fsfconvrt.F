*  
* $Id: fsfconvrt.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFCONVRT(ITE)
********************************************************************
* conversion of BRAHMS detector IDs (of the order 10^2..10^3) into
* IDs usable for DELSEARCH (bit-encoding -> need IDs <32)
*
* INPUT:  ITE : TE array
* OUTPUT: ITE : TE array with different module ID
*
*                                         written by Kristian Harder
********************************************************************
      IMPLICIT NONE
#include "include/detnums.inc"
#include "include/fsfparams.inc"
      INTEGER ID,ITE(*),DELDET

      ID=DELDET(ITE(1))
      ITE(1)=ID

      RETURN
      END

