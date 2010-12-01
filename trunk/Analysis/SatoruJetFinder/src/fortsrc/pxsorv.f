CDECK  ID>, PXSORV.
      SUBROUTINE PXSORV (ISZ,ARY,KIX,COPT)
*.*********************************************************
*. ------
*. PXSORV
*. ------
*. SOURCE: HERWIG (B.Webber,G.Marchesini)
*. Sort a real array into assending order based on
*. the magnitude of its elements; provide an
*. integer "index array" which specifies the ordering
*. of the array.
*. Usage     :
*.
*.      PARAMETER  (NDIM=1.or.more)
*.      REAL  ARY (NDIM)
*.      INTEGER  KIX (NDIM)
*.      INTEGER  ISIZ
*.      CHARACTER*1  COPT
*.
*.      ISIZ = 1.to.NDIM
*.      COPT = 'I'
*.      CALL PXSORV (ISIZ,ARY,KIX,COPT)
*.
*. INPUT     : ISIZ  The dimension of the input array
*. INPUT     : ARY   The input array
*. OUTPUT    : KIX   The index array
*. CONTROL   : COPT  Control of output vector ARY
*.              = ' ' : return sorted ARY and index array KIX
*.              = 'I' : return index array KIX only, don't
*.                      modify ARY
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  MXSZ
      PARAMETER  (MXSZ=500)
      INTEGER  ISZ,IX,JX
      INTEGER  KIX (*),IL (MXSZ),IR (MXSZ)
      REAL  ARY (*),BRY (MXSZ)
      CHARACTER*1  COPT
      IF (ISZ.GT.MXSZ) THEN
          WRITE (6,FMT='('' PXSORT: Error,'',
     +           '' Max array size ='',I6)') MXSZ
          KIX (1) = -1
          GO TO 990
      END IF
      IL (1) = 0
      IR (1) = 0
      DO 10 IX = 2,ISZ
          IL (IX) = 0
          IR (IX) = 0
          JX = 1
   2      IF (ARY (IX).GT.ARY (JX)) GO TO 5
   3      IF (IL (JX).EQ.0) GO TO 4
          JX = IL (JX)
          GO TO 2
   4      IR (IX) = -JX
          IL (JX) =  IX
          GO TO 10
   5      IF (IR (JX).LE.0) GO TO 6
          JX = IR (JX)
          GO TO 2
   6      IR (IX) = IR (JX)
          IR (JX) = IX
  10  CONTINUE
      IX = 1
      JX = 1
      GO TO 8
  20  JX = IL (JX)
   8  IF (IL (JX).GT.0) GO TO 20
   9  KIX (IX) = JX
      BRY (IX) = ARY (JX)
      IX = IX + 1
      IF (IR (JX)) 12,30,13
  13  JX = IR (JX)
      GO TO 8
  12  JX = -IR (JX)
      GO TO 9
  30  IF (COPT.EQ.'I') RETURN
      DO 31 IX = 1,ISZ
  31  ARY (IX) = BRY (IX)
 990  RETURN
      END
