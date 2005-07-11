CDECK  ID>, PXPRNT.
      SUBROUTINE PXPRNT (IFST,ILST,ISZE,IPRT,ARRY,COPT)
*.*********************************************************
*. ------
*. PXPRNT
*. ------
*. SOURCE: J.W.Gary
*. Print-out of vectors from an array of vectors
*. Usage     :
*.
*.      INTEGER  ITKDM,MXTRAK
*.      PARAMETER  (ITKDM=1.or.more,MXTRAK=1.or.more)
*.      INTEGER IFST,ILST,ILEN
*.      REAL  ARRY (ITKDM,MXTRAK)
*.      CHARACTER*1  COPT
*.
*.      ... fill vectors PTRAK
*.
*.      IFST = 1.to.MXTRAK
*.      ILST = IFST.to.MXTRAK
*.      ILEN = 1.to.ITKDM
*.      COPT = '4'
*.      CALL PXPRNT (IFST,ILST,ITKDM,ILEN,ARRY,COPT)
*.
*.      ... example printout (with IFST=1,ILST=3,COPT='4')
*.
*.      Vector    Px          Py          Pz           E
*.         1 -0.1431E+00  0.2408E+00  0.2636E+00  0.4092E+00
*.         2  0.4027E+00  0.3949E+00 -0.4255E+00  0.7202E+00
*.         3 -0.1195E+00 -0.1654E+00  0.7137E-01  0.2573E+00
*.
*. INPUT     : IFST    First vector to print
*. INPUT     : ILST    Last vector to print
*. INPUT     : ITKDM   The length of each vector in the array
*. INPUT     : ILEN    The number of elements in each vector
*.                     to print (1 to ILEN)
*. INPUT     : ARRY    The array containing the vectors
*. INPUT     : COPT    Printing option
*.                       ='4' for 4-momenta
*.                       ='E' for exponential (E12.4) formatting
*.                       =' ' for "standard" (F8.2) formatting
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  IOL
      PARAMETER  (IOL=6)
      INTEGER  IFST,ILST,ISZE,IX,JX,IPRT
      REAL ARRY (ISZE,*)
      CHARACTER*1  COPT

      IF (IPRT.LT.1.OR.IPRT.GT.10) THEN
          WRITE (IOL,FMT='('' PXPRNT: Error,'',
     +      ''Array size must be between 1 and 10'')')
          RETURN
      END IF
      IF (COPT.EQ.'4') THEN
          WRITE (IOL,FMT='('' Vector'',4X,''Px'',10X,
     +           ''Py'',10X,''Pz'',10X,'' E'')')
      ELSE IF (COPT.EQ.'E') THEN
          WRITE (IOL,FMT='(:'' Vector'',I6,9(I12))')
     +          (IX,IX=1,IPRT)
      ELSE
          WRITE (IOL,FMT='(:'' Vector'',I6,9(I8))')
     +          (IX,IX=1,IPRT)
      END IF
      DO 160  JX = IFST,ILST
          IF (COPT.EQ.'4') THEN
              WRITE (IOL,FMT='(1X,I4,4E12.4)')
     +               JX,(ARRY (IX,JX),IX=1,4)
          ELSE IF (COPT.EQ.'E') THEN
              WRITE (IOL,FMT='(:I5,1X,10(E12.4))')
     +           JX,(ARRY (IX,JX),IX=1,IPRT)
          ELSE
              WRITE (IOL,FMT='(:I5,1X,10(F8.2))')
     +           JX,(ARRY (IX,JX),IX=1,IPRT)
          END IF
 160  CONTINUE
      RETURN
      END
