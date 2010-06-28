      real pin(3), q, vin(3),B,EB,zbcal,pout(3),vout(3)
      character*500 pxst,pyst,pzst,qst,vxst,vyst,vzst
      B=3.5
      EB=250.
      zbcal=3000.
      Nargs = IARGC ()
      if (Nargs .eq. 7 ) then
          CALL GETARG(1 ,pxst )
          CALL GETARG(2 ,pyst )
          CALL GETARG(3 ,pzst )
          CALL GETARG(4 ,qst )
          CALL GETARG(5 ,vxst )
          CALL GETARG(6 ,vyst )
          CALL GETARG(7 ,vzst )
          read(pxst,*) pin(1)
          read(pyst,*) pin(2)
          read(pzst,*) pin(3)
          read(qst,*) q
          read(vxst,*) vin(1)
          read(vyst,*) vin(2)
          read(vzst,*) vin(3)
      else
        print *, ' please give seven arguments ! '
        stop
      endif   
 
      call bcalhit(pin,q,vin,B,EB,zbcal,pout,vout)
      print *, (pout(i),i=1,3)
      print *, (vout(i),i=1,3)
      end

