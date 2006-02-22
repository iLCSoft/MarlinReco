*  
* $Id: fkexmn.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKEXMN(D1,W1,D2,W2,DM,WM,CV,IERR)
*                                                                      *
*   weighted mean of measurements and variation of chi2                *
*                                                                      *
*   Input :  D1  : deviations from reference )                         *
*            W1  : weight matrix             )   of first measurement  *
*            D2  : deviations from reference )                         *
*            W2  : weight matrix             )   of second measurement *
*                                                                      *
*   Output : DM  : deviations from reference )                         *
*            WM  : weight matrix             )   of combination        *
*            CV  : covariance matrix         )                         *
*            IERR :  error flag : 0 if OK                              *
*                                 2 if weight matrix not pos.          *
*                                                                      *
*     DM and WM may overwrite D1 and W1, or D2 and W2                  *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D1(5),D2(5),DM(5),V(5),W1(15),W2(15),WM(15),CV(15)
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
      END
