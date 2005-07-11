*  
* $Id: fkwgmn.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKWGMN(KRANK,D1,W1,ND1,CH1,IRNK1,D2,W2,ND2,CH2,IRNK2,
     & DM,WM,NDM,CHM,IERR)
*                                                                      *
*   weighted mean of measurements and variation of chi2                *
*                                                                      *
*   Input :  KRANK : rank of combination                               *
*            D1  : deviations from reference )                         *
*            W1  : weight matrix             )                         *
*            ND1 : nb of deg. of freedom     )   of first measurement  *
*            CH1 : chi2                      )                         *
*            IRNK1 : rank                    )                         *
*            D2  : deviations from reference )                         *
*            W2  : weight matrix             )                         *
*            ND2 : nb of deg. of freedom     )   of second measurement *
*            CH2 : chi2                      )                         *
*            IRNK2 : rank                    )                         *
*                                                                      *
*   Output : DM  : deviations from reference )                         *
*            WM  : weight matrix             )                         *
*            NDM : nb of deg. of freedom     )   of combination        *
*            CHM : chi2                      )                         *
*            IERR :  error flag : 0 if OK, 2 if weight matrix not pos. *
*                                                                      *
*     DM and WM may overwrite D1 and W1, or D2 and W2                  *
*     NDM and CHM     "       ND1 and CH1, or ND2 and CH2              *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL CH1,CH2,CHM
*
      DIMENSION D1(5),D2(5),DM(5),V(5),W1(15),W2(15),WM(15),CV(15)
     & ,DS1(5),DS2(5),WS1(15),WS2(15)
*
*  save initial quantities to compute variation of chi2
*
      CALL FKDPCO(D1,DS1,5)
      CALL FKDPCO(D2,DS2,5)
      CALL FKDPCO(W1,WS1,15)
      CALL FKDPCO(W2,WS2,15)
*
*  "raw" chi2
      CHR=FKDCH2(W1,D1)+FKDCH2(W2,D2)
*
*  compute linear system
*
      V(1)=
     & W1( 1)*D1(1)+W1( 2)*D1(2)+W1( 4)*D1(3)+W1( 7)*D1(4)+W1(11)*D1(5)+
     & W2( 1)*D2(1)+W2( 2)*D2(2)+W2( 4)*D2(3)+W2( 7)*D2(4)+W2(11)*D2(5)
      V(2)=
     & W1( 2)*D1(1)+W1( 3)*D1(2)+W1( 5)*D1(3)+W1( 8)*D1(4)+W1(12)*D1(5)+
     & W2( 2)*D2(1)+W2( 3)*D2(2)+W2( 5)*D2(3)+W2( 8)*D2(4)+W2(12)*D2(5)
      V(3)=
     & W1( 4)*D1(1)+W1( 5)*D1(2)+W1( 6)*D1(3)+W1( 9)*D1(4)+W1(13)*D1(5)+
     & W2( 4)*D2(1)+W2( 5)*D2(2)+W2( 6)*D2(3)+W2( 9)*D2(4)+W2(13)*D2(5)
      V(4)=
     & W1( 7)*D1(1)+W1( 8)*D1(2)+W1( 9)*D1(3)+W1(10)*D1(4)+W1(14)*D1(5)+
     & W2( 7)*D2(1)+W2( 8)*D2(2)+W2( 9)*D2(3)+W2(10)*D2(4)+W2(14)*D2(5)
      V(5)=
     & W1(11)*D1(1)+W1(12)*D1(2)+W1(13)*D1(3)+W1(14)*D1(4)+W1(15)*D1(5)+
     & W2(11)*D2(1)+W2(12)*D2(2)+W2(13)*D2(3)+W2(14)*D2(4)+W2(15)*D2(5)
*
      CALL FKDPVA(W1,W2,WM,15)
*
*  check positivity of total weight matrix and solve the system
*
      CALL FKMI5D(WM,CV,IERR)
      IF(IERR.NE.0) RETURN
*
      DM(1)=CV( 1)*V(1)+CV( 2)*V(2)+CV( 4)*V(3)+CV( 7)*V(4)+CV(11)*V(5)
      DM(2)=CV( 2)*V(1)+CV( 3)*V(2)+CV( 5)*V(3)+CV( 8)*V(4)+CV(12)*V(5)
      DM(3)=CV( 4)*V(1)+CV( 5)*V(2)+CV( 6)*V(3)+CV( 9)*V(4)+CV(13)*V(5)
      DM(4)=CV( 7)*V(1)+CV( 8)*V(2)+CV( 9)*V(3)+CV(10)*V(4)+CV(14)*V(5)
      DM(5)=CV(11)*V(1)+CV(12)*V(2)+CV(13)*V(3)+CV(14)*V(4)+CV(15)*V(5)
*
*   variation of chi2
*
      NDM=ND1+IRNK1+ND2+IRNK2-KRANK
      DO 4 I=1,5
      DS1(I)=DS1(I)-DM(I)
    4 DS2(I)=DS2(I)-DM(I)
*
      IF (NDM.GT.0) THEN
        CHM=CH1+CH2+FKDCH2(WS1,DS1)+FKDCH2(WS2,DS2)
      ELSE
        CHM=0
      END IF
*
      END
