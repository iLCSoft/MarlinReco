*  
* $Id: fxdump.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
       SUBROUTINE FXDUMP(IFXDBG)
*********************************************************************
*      Debug facilities
*********************************************************************
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fxapar.inc"
*
       IPRNT = 6
       DO 100 I=1,NTKR
       WRITE(IPRNT,1)I,IDTKR(I),NTERUS(I),MODUSE(I)
       IF(IFXDBG.GT.1)THEN
         WRITE(IPRNT,5)ITKDAT(2,I),ITKDAT(3,I),ITKDAT(4,I),ITKDAT(6,I),
     /                 (TKDAT(K,I),K=17,22)
         WRITE(IPRNT,2)(IDTERU(K,I),K=1,NTERUS(I))
c         DO J=1,MXEXLI
c            WRITE(IPRNT,3)(IDEXTE(J,K,I),K=1,NTERUS(I))
c         ENDDO
         WRITE(IPRNT,4)(IMODTE(K,I),K=1,NTERUS(I))
         WRITE(IPRNT,6)(ILABTE(K,I),K=1,NTERUS(I))
c         DO J=1,MXEXLI
c           WRITE(IPRNT,7)(IDLITE(J,K,I),K=1,NTERUS(I))
c         ENDDO
       ENDIF
  100  CONTINUE
    1  FORMAT(1X,I5,' TKR: ',1X,I8,I9,1X,I8)
    2  FORMAT(6X,'TER ',12(1X,I8)/(10X,12(1X,I8)))
    3  FORMAT(6X,'Excl',12(1X,I8)/(10X,12(1X,I8)))
    4  FORMAT(6X,'Dets',12I9/(10X,12I9))
    5  FORMAT(6X,'DATA',1X,I8,1X,I8,1X,I8,I9,6F9.3)
    6  FORMAT(6X,'Labl',12I9/(10X,12I9))
    7  FORMAT(6X,'Link',12(1X,I8)/(10X,12(1X,I8)))
       RETURN
       END

