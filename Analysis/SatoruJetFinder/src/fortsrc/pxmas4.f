CDECK  ID>, PXMAS4.
      SUBROUTINE PXMAS4 (PMOM,XMAS)
*.*********************************************************
*. ------
*. PXMAS4
*. ------
*. SOURCE: HERWIG (B.Webber,G.Marchesini)
*. Calculate the invariant mass of a 4-vector
*. (negative if spacelike)
*. Usage     :
*.
*.      REAL  PMOM (4.or.more)
*.      REAL  XMAS
*.
*.      CALL PXMAS4 (PMOM,XMAS)
*.
*. INPUT     : PMOM   Particle 4-momentum (Px,Py,Pz,E)
*. OUTPUT    : XMAS   The invariant mass
*.
*.*********************************************************
      IMPLICIT NONE
      REAL  PMOM (*)
      REAL  XMAS2,XMAS
      XMAS2 = ((PMOM (4) + PMOM (3)) * (PMOM (4) - PMOM (3))
     +        - PMOM (1)**2 - PMOM(2)**2)
      CALL PXROOT (XMAS2,XMAS)
      RETURN
      END
