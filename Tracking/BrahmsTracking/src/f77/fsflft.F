*  
* $Id: fsflft.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFLFT(X,Y,HIN,RA,RHOUT,ICASE,NPOINT,RCHI2)
*******************************************************************************
*
*
*
*    subroutine to fit a straight line allowing for errors in both directons
*
*
*
*******************************************************************************
*
      implicit none
      integer npoint, icase
      REAL X(NPOINT),Y(NPOINT),HIN(3,NPOINT),RHOUT(3),RA(2),RCHI2
      double precision GX(3), HOUT(3)
      double precision A(2),AGC(2),C,GB,GBI,chi2,chinew
      integer i, ier

*
*     first take x error free
*
      do i=1,3
         gx(i) = 0.
      enddo
      do i=1,2
         agc(i) = 0.
      enddo
      DO I=1,NPOINT
        IF (HIN(3,I) .LE. 0.) GOTO 5
        GX(1)  = GX(1)+DBLE(X(I))**2/DBLE(HIN(3,I))
        GX(2)  = GX(2)+DBLE(X(I))/DBLE(HIN(3,I))
        GX(3)  = GX(3)+1./DBLE(HIN(3,I))
        AGC(1) = AGC(1)+DBLE(X(I))*DBLE(Y(I))/DBLE(HIN(3,I))
        AGC(2) = AGC(2)+DBLE(Y(I))/DBLE(HIN(3,I))
 5    enddo
      CALL FSFIC2(GX,HOUT,IER)
      IF (IER .NE. 0) THEN
C -KH-  CHI2 = 99999.
        RCHI2 = 99999.
        RETURN
      ENDIF
      A(1)     = HOUT(1)*AGC(1)+HOUT(2)*AGC(2)
      A(2)     = HOUT(2)*AGC(1)+HOUT(3)*AGC(2)
      IF (ICASE.EQ.1)GOTO 10
*
* *** now improve once using full error matrix
*
      do i=1,3
         gx(i) = 0.
      enddo
      do i=1,2
         agc(i) = 0.
      enddo
      DO I=1,NPOINT
        IF (HIN(3,I) .LE. 0.) GOTO 7
        GB     = A(1)*A(1)*HIN(1,I)+2.*A(1)*HIN(2,I)+HIN(3,I)
        GBI    = 1./(A(1)*A(1)*HIN(1,I)+2.*A(1)*HIN(2,I)+HIN(3,I))
        GX(1)  = GX(1)+DBLE(X(I))**2*GBI
        GX(2)  = GX(2)+DBLE(X(I))*GBI
        GX(3)  = GX(3)+GBI
        C      = DBLE(Y(I))-A(1)*DBLE(X(I))-A(2)
        AGC(1) = AGC(1)+DBLE(X(I))*C*GBI
        AGC(2) = AGC(2)+C*GBI
 7    enddo
      CALL FSFIC2(GX,HOUT,IER)
      IF (IER .NE. 0) GOTO 10
      A(1)     = A(1)+HOUT(1)*AGC(1)+HOUT(2)*AGC(2)
      A(2)     = A(2)+HOUT(2)*AGC(1)+HOUT(3)*AGC(2)
*
* *** calculate chi**2
*
   10 CHI2     = 0.
      DO I=1,NPOINT
        IF (HIN(3,I) .LE. 0.) GOTO 20
        GB     = A(1)*A(1)*HIN(1,I)+2.*A(1)*HIN(2,I)+HIN(3,I)
        CHINEW = (Y(I)-A(1)*X(I)-A(2))**2/GB
        CHI2   = CHI2+CHINEW
 20   enddo
      RA(1) = A(1)
      RA(2) = A(2)
      RHOUT(1) = HOUT(1)
      RHOUT(2) = HOUT(2)
      RHOUT(3) = HOUT(3)
      rchi2    = chi2
      RETURN
      END
