CDECK  ID>, PXLUC5.
      SUBROUTINE PXLUC5 (N,NRLUDM,P,MSTU46,MSTU47,
     +                   PARU44,PARU45,NJET,IJMUL,IPASS)
*.*********************************************************
*. ------
*. PXLUC5
*. ------
*. An "in-house" version of the Jetset jet-finding algorithm
*. which works entirely through an argument list rather than
*. through e.g. the Jetset common blocks. Its operation is
*. therefore entirely decoupled from the the operation of Jetset
*. (i.e. the values of MST or MSTJ etc. in Jetset common do not
*. affect this routine whatsoever).
*. The main purpose of an in-house version of the
*. Jetset jetfinding algorithm is to have a version
*. which is compatible with both Jetset6.3 and Jetset7.1 etc.
*. (because of the change in the Jetset common blocks between
*. these two versions, the version of this algorithm in the
*. Jetset library is version specific).
*. The input arguments MSTU46, MSTU47, PARU44, PARU45 correspond
*. to the parameters MSTU(46), MSTU(47), PARU(44), PARU(45) of
*. Jetset7.1, see "A manual to ... Jetset7.1," T.Sjostrand
*. (filename "Jetset7.1 MANUAL A" on the Opal generator disk).
*.
*.      INTEGER  NRLUDM,MXJET
*.      PARAMETER (NRLUDM=1000.or.so,MXJET=10.or.so)
*.      INTEGER  IPASS (NRLUDM),IJMUL (MXJET)
*.      INTEGER  NTRAK,NJET
*.      REAL PLUND (NRLUDM,5)
*.      REAL  MSTU46,MSTU47,PARU44,PARU45
*.
*.      (define NTRAK, fill PLUND)
*.      CALL PXLUC5 (NTRAK,NRLUDM,PLUND,MSTU46,MSTU47,
*.     +             PARU44,PARU45,NJET,IJMUL,IPASS)
*.
*. INPUT     : NTRAK    Number of tracks
*. INPUT     : NRLUDM   First dimension of PLUND
*. IN/OUTPUT : PLUND    4-momenta in Jetset format
*. INPUT     : MSTU46   same as MSTU(46) in Jetset7.1:jet-finder mode
*. INPUT     : MSTU47   same as MSTU(47) in Jetset7.1
*. INPUT     : PARU44   same as PARU(44) in Jetset7.1
*. INPUT     : PARU45   same as PARU(45) in Jetset7.1
*. OUTPUT    : NJET      Number of jets found
*. OUTPUT    : IJMUL(i)  Jet i contains IJMUL(i) particles
*. OUTPUT    : IPASS(k)  Particle k belongs to jet number IPASS(k)
*.
*. CALLS     : PXANXY,PXPLU3,PXRMX3,PXROF3,PXROB3
*. CALLED    : PXLTH4
*.
*. AUTHOR    : Modified from LUCLUS (T.Sjostrand) by J.W.Gary
*. CREATED   : 31-Jan-89
*. LAST MOD  : 31-Jan-89
*.
*. Modification Log.
*. 05-May-97 D. Chrisman, remove declaration of unused variables
*.             I1, I2 and PMAS.
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  LOCDIM
      PARAMETER  (LOCDIM=2000)
      REAL PARU48,PIMAS,PARU43
      PARAMETER (PARU48=0.0001,PIMAS=0.13957,PARU43=0.25)
      INTEGER  MSTU42,MSTU48,MSTU43,MSTU46,MSTU47,ITRY1,ITRY2,
     +         I,J,IDEL,IREC,IEMP,IJET,IMIN,IMAX,NJET,NREM,
     +         INEW,IORI,NPRE,ISPL,NSAV,ITRY,NP,IMIN1,IMIN2,NRLUDM,
     +         N,NLOOP,IP,IJ
      INTEGER  K (LOCDIM,5),IJMUL (*),IPASS (*)
      REAL  PARU44,PARU45,R2MAX,PEMAX,RINIT,PXRR2M,PXRR2T,PSS,
     +      PMAX,R2,TSAV,PSJT,R2ACC,R2MIN,PXMAS
      REAL  P(NRLUDM,*),V(LOCDIM,5),PS(5)
      DATA  MSTU42 / 2 /,MSTU48 / 0 /, MSTU43 / 1 /

CC...If first time, reset. If reentering, skip preliminaries.
**JWG      IF(MSTU48.LE.0) THEN
        NP=0
        DO 100 J=1,5
  100   PS(J)=0.
        PSS=0.
**JWG      ELSE
**JWG        NJET=NSAV
**JWG        IF(MSTU43.GE.2) N=N-NJET
**JWG        DO 110 I=N+1,N+NJET
**JWG  110   P(I,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
**JWG        IF(MSTU46.LE.3) R2ACC=PARU44**2
**JWG        IF(MSTU46.GE.4) R2ACC=PARU45*PS(5)**2
**JWG        NLOOP=0
**JWG        GOTO 290
**JWG      ENDIF
C...Find which particles are to be considered in cluster search.
      DO 140 I=1,N
      IF(N+2*NP.GE.LOCDIM-5) THEN
          WRITE (6,FMT='('' PXLUC5: Error, not enough buffer'',
     +           ''space for jet-finer calculation'')')
          NJET = -1
          GO TO 990
      ENDIF
C...Take copy of these particles, with space left for jets later on.
      NP=NP+1
      K(N+NP,3)=I
      DO 120 J=1,5
  120 P(N+NP,J)=P(I,J)
**JWG      IF(MSTU42.EQ.0) P(N+NP,5)=0.
      PXMAS = P (I,5)
**JWG      IF(MSTU42.EQ.1.AND.PXMAS.NE.0) P(N+NP,5)=PIMAS
      P(N+NP,4)=SQRT(P(N+NP,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
      P(N+NP,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
      DO 130 J=1,4
  130 PS(J)=PS(J)+P(N+NP,J)
      PSS=PSS+P(N+NP,5)
  140 CONTINUE
      DO 150 I=N+1,N+NP
      K(I+NP,3)=K(I,3)
      DO 150 J=1,5
  150 P(I+NP,J)=P(I,J)
      PS(5)=SQRT(MAX(0.,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2))
C...Very low multiplicities not considered.
      IF(NP.LT.MSTU47) THEN
        NJET=-1
        RETURN
      ENDIF
C...Find precluster configuration. If too few jets, make harder cuts.
      NLOOP=0
      IF(MSTU46.LE.3) R2ACC=PARU44**2
      IF(MSTU46.GE.4) R2ACC=PARU45*PS(5)**2
      RINIT=1.25*PARU43
      IF(NP.LE.MSTU47+2) RINIT=0.
  160 RINIT=0.8*RINIT
      NPRE=0
      NREM=NP
      DO 170 I=N+NP+1,N+2*NP
  170 K(I,4)=0
C...Sum up small momentum region. Jet if enough absolute momentum.
      IF(MSTU46.LE.2) THEN
        DO 180 J=1,4
  180   P(N+1,J)=0.
        DO 200 I=N+NP+1,N+2*NP
        IF(P(I,5).GT.2.*RINIT) GOTO 200
        NREM=NREM-1
        K(I,4)=1
        DO 190 J=1,4
  190   P(N+1,J)=P(N+1,J)+P(I,J)
  200   CONTINUE
        P(N+1,5)=SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2)
        IF(P(N+1,5).GT.2.*RINIT) NPRE=1
        IF(RINIT.GE.0.2*PARU43.AND.NPRE+NREM.LT.MSTU47) GOTO 160
      ENDIF
C...Find fastest remaining particle.
  210 NPRE=NPRE+1
      PMAX=0.
      DO 220 I=N+NP+1,N+2*NP
      IF(K(I,4).NE.0.OR.P(I,5).LE.PMAX) GOTO 220
      IMAX=I
      PMAX=P(I,5)
  220 CONTINUE
      DO 230 J=1,5
  230 P(N+NPRE,J)=P(IMAX,J)
      NREM=NREM-1
      K(IMAX,4)=NPRE
C...Sum up precluster around it according to pT separation.
      IF(MSTU46.LE.2) THEN
        DO 250 I=N+NP+1,N+2*NP
        IF(K(I,4).NE.0) GOTO 250
        R2=PXRR2T(NRLUDM,P,I,IMAX)
        IF(R2.GT.RINIT**2) GOTO 250
        NREM=NREM-1
        K(I,4)=NPRE
        DO 240 J=1,4
  240   P(N+NPRE,J)=P(N+NPRE,J)+P(I,J)
  250   CONTINUE
        P(N+NPRE,5)=SQRT(P(N+NPRE,1)**2+P(N+NPRE,2)**2+P(N+NPRE,3)**2)
C...Sum up precluster around it according to mass separation.
      ELSE
  260   IMIN=0
        R2MIN=RINIT**2
        DO 270 I=N+NP+1,N+2*NP
        IF(K(I,4).NE.0) GOTO 270
        R2=PXRR2M(NRLUDM,P,I,N+NPRE)
        IF(R2.GE.R2MIN) GOTO 270
        IMIN=I
        R2MIN=R2
  270   CONTINUE
        IF(IMIN.NE.0) THEN
          DO 280 J=1,4
  280     P(N+NPRE,J)=P(N+NPRE,J)+P(IMIN,J)
          P(N+NPRE,5)=SQRT(P(N+NPRE,1)**2+P(N+NPRE,2)**2+P(N+NPRE,3)**2)
          NREM=NREM-1
          K(IMIN,4)=NPRE
          GOTO 260
        ENDIF
      ENDIF
C...Check if more preclusters to be found. Start over if too few.
      IF(RINIT.GE.0.2*PARU43.AND.NPRE+NREM.LT.MSTU47) GOTO 160
      IF(NREM.GT.0) GOTO 210
      NJET=NPRE
C...Reassign all particles to nearest jet. Sum up new jet momenta.
  290 TSAV=0.
      PSJT=0.
  300 IF(MSTU46.LE.1) THEN
        DO 310 I=N+1,N+NJET
        DO 310 J=1,4
  310   V(I,J)=0.
        DO 340 I=N+NP+1,N+2*NP
        R2MIN=PSS**2
        DO 320 IJET=N+1,N+NJET
        IF(P(IJET,5).LT.RINIT) GOTO 320
        R2=PXRR2T(NRLUDM,P,I,IJET)
        IF(R2.GE.R2MIN) GOTO 320
        IMIN=IJET
        R2MIN=R2
  320   CONTINUE
        K(I,4)=IMIN-N
        DO 330 J=1,4
  330   V(IMIN,J)=V(IMIN,J)+P(I,J)
  340   CONTINUE
        PSJT=0.
        DO 360 I=N+1,N+NJET
        DO 350 J=1,4
  350   P(I,J)=V(I,J)
        P(I,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
  360   PSJT=PSJT+P(I,5)
      ENDIF
C...Find two closest jets.
      R2MIN=2.*R2ACC
      DO 370 ITRY1=N+1,N+NJET-1
      DO 370 ITRY2=ITRY1+1,N+NJET
      IF(MSTU46.LE.2) R2=PXRR2T(NRLUDM,P,ITRY1,ITRY2)
      IF(MSTU46.GE.3) R2=PXRR2M(NRLUDM,P,ITRY1,ITRY2)
      IF(R2.GE.R2MIN) GOTO 370
      IMIN1=ITRY1
      IMIN2=ITRY2
      R2MIN=R2
  370 CONTINUE
C...If allowed, join two closest jets and start over.
      IF(NJET.GT.MSTU47.AND.R2MIN.LT.R2ACC) THEN
        IREC=MIN(IMIN1,IMIN2)
        IDEL=MAX(IMIN1,IMIN2)
        DO 380 J=1,4
  380   P(IREC,J)=P(IMIN1,J)+P(IMIN2,J)
        P(IREC,5)=SQRT(P(IREC,1)**2+P(IREC,2)**2+P(IREC,3)**2)
        DO 390 I=IDEL+1,N+NJET
        DO 390 J=1,5
  390   P(I-1,J)=P(I,J)
        IF(MSTU46.GE.2) THEN
          DO 400 I=N+NP+1,N+2*NP
          IORI=N+K(I,4)
          IF(IORI.EQ.IDEL) K(I,4)=IREC-N
  400     IF(IORI.GT.IDEL) K(I,4)=K(I,4)-1
        ENDIF
        NJET=NJET-1
        GOTO 290
C...Divide up broad jet if empty cluster in list of final ones.
      ELSEIF(NJET.EQ.MSTU47.AND.MSTU46.LE.1.AND.NLOOP.LE.2) THEN
        DO 410 I=N+1,N+NJET
  410   K(I,5)=0
        DO 420 I=N+NP+1,N+2*NP
  420   K(N+K(I,4),5)=K(N+K(I,4),5)+1
        IEMP=0
        DO 430 I=N+1,N+NJET
  430   IF(K(I,5).EQ.0) IEMP=I
        IF(IEMP.NE.0) THEN
          NLOOP=NLOOP+1
          ISPL=0
          R2MAX=0.
          DO 440 I=N+NP+1,N+2*NP
          IF(K(N+K(I,4),5).LE.1.OR.P(I,5).LT.RINIT) GOTO 440
          IJET=N+K(I,4)
          R2=PXRR2T(NRLUDM,P,I,IJET)
          IF(R2.LE.R2MAX) GOTO 440
          ISPL=I
          R2MAX=R2
  440     CONTINUE
          IF(ISPL.NE.0) THEN
            IJET=N+K(ISPL,4)
            DO 450 J=1,4
            P(IEMP,J)=P(ISPL,J)
  450       P(IJET,J)=P(IJET,J)-P(ISPL,J)
            P(IEMP,5)=P(ISPL,5)
            P(IJET,5)=SQRT(P(IJET,1)**2+P(IJET,2)**2+P(IJET,3)**2)
            IF(NLOOP.LE.2) GOTO 290
          ENDIF
        ENDIF
      ENDIF
C...If generalized thrust has not yet converged, continue iteration.
      IF(MSTU46.LE.1.AND.NLOOP.LE.2.AND.PSJT/PSS.GT.TSAV+PARU48)
     +THEN
        TSAV=PSJT/PSS
        GOTO 300
      ENDIF
C...Reorder jets according to energy.
      DO 460 I=N+1,N+NJET
      DO 460 J=1,5
  460 V(I,J)=P(I,J)
      DO 490 INEW=N+1,N+NJET
      PEMAX=0.
      DO 470 ITRY=N+1,N+NJET
      IF(V(ITRY,4).LE.PEMAX) GOTO 470
      IMAX=ITRY
      PEMAX=V(ITRY,4)
  470 CONTINUE
      K(INEW,1)=31
      K(INEW,2)=97
      K(INEW,3)=INEW-N
      K(INEW,4)=0
      DO 480 J=1,5
  480 P(INEW,J)=V(IMAX,J)
      V(IMAX,4)=-1.
  490 K(IMAX,5)=INEW
C...Clean up particle-jet assignments and jet information.
      DO 500 I=N+NP+1,N+2*NP
      IORI=K(N+K(I,4),5)
      K(I,4)=IORI-N
      IF(K(K(I,3),1).NE.3) K(K(I,3),4)=IORI-N
      K(IORI,4)=K(IORI,4)+1
  500 CONTINUE
      IEMP=0
      PSJT=0.
      DO 520 I=N+1,N+NJET
      K(I,5)=0
      PSJT=PSJT+P(I,5)
      P(I,5)=SQRT(MAX(P(I,4)**2-P(I,5)**2,0.))
      DO 510 J=1,5
  510 V(I,J)=0.
  520 IF(K(I,4).EQ.0) IEMP=I
C...Select storing option. Output variables. Check for failure.
      IF(IEMP.NE.0) THEN
        NJET=-1
      ENDIF
      NSAV=NJET
      DO 560  IJ = 1,NJET
          IJMUL (IJ) = K (N+IJ,4)
 560  CONTINUE
      DO 580  IP = 1,NP
          IPASS (IP) = K (IP,4)
 580  CONTINUE

 990  RETURN
      END

C...Functions: distance measure in pT or (pseudo)mass.
      FUNCTION PXRR2T (NRLUDM,P,I1,I2)
      IMPLICIT NONE
      INTEGER  NRLUDM,I1,I2
      REAL  P (NRLUDM,*),PXRR2T
      PXRR2T = (P(I1,5)*P(I2,5)-P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)-
     +P(I1,3)*P(I2,3))*2.*P(I1,5)*P(I2,5)/(0.0001+P(I1,5)+P(I2,5))**2
      RETURN
      END

      FUNCTION PXRR2M (NRLUDM,P,I1,I2)
      IMPLICIT NONE
      INTEGER  NRLUDM,I1,I2
      REAL  P (NRLUDM,*),PXRR2M
      PXRR2M = 2.*P(I1,4)*P(I2,4)*(1.-(P(I1,1)*P(I2,1)+P(I1,2)*
     +P(I2,2)+P(I1,3)*P(I2,3))/(P(I1,5)*P(I2,5)))
      RETURN
      END
