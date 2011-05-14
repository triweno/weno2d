program tri
  include 'triLRIM.comm.inc'
  
  character (len = 50) fname
  


  real    (kind = 8) :: atan, la, lb, lc, alpha, qre, qrb, atmp, rnx, rny
  integer (kind = 4) :: i,j,k
  REAL (kind = 8), EXTERNAL :: BO1,BO2

  GAM=1.4d0
  AGAM=GAM-1.D0
  BGAM=2.D0*SQRT(GAM/AGAM)
  CGAM=1.D0/GAM
  DGAM=2.D0/AGAM
  EGAM=AGAM/(GAM+1.D0)
  GGAM=DSQRT(GAM*AGAM)
  HGAM=AGAM/2.D0
  FGAM=3.D0*GAM-1.D0
  OGAM=AGAM/(2.D0*GAM)
  QGAM=GAM+1.D0
  PGAM=QGAM/(2.D0*GAM)
  RGAM=4.D0*GAM
  SGAM=GAM*AGAM
  TGAM=QGAM/2.D0

  call ReadFile('tube.1')
  
  
  ! calculation start
  allocate(r(ntri))
  allocate(u(ntri))
  allocate(v(ntri))
  allocate(p(ntri))
  allocate(ro_1(ntri))
  allocate(ru_1(ntri))
  allocate(rv_1(ntri))
  allocate(e(ntri))
  allocate(en_1(ntri))
  allocate(ro(ntri))
  allocate(ru(ntri))
  allocate(rv(ntri))
  allocate(en(ntri))
  
!  do i = 1, ntri
!    r(i) = 1.4d0
!    p(i) = 1.0d0
!    u(i) = 3.0d0
!    v(i) = 0.0d0
!  enddo

!  do i = 1, ntri
!  	if (x(2,tri_v(1,i))<=(y_-epsilon) .or. x(2,tri_v(2,i))<=(y_-epsilon) .or. x(2,tri_v(3,i))<=(y_-epsilon)) then
!    	r(i) = RO1_
!    	u(i) = U1_
!    	p(i) = P1_
!    	v(i) = V1_
!    else
!     if (x(2,tri_v(1,i))<=-x(1,tri_v(1,i))/5-epsilon .or. x(2,tri_v(2,i))<=-x(1,tri_v(2,i))/5-epsilon .or.   &
!         x(2,tri_v(3,i))<=-x(1,tri_v(3,i))/5-epsilon) then
!      r(i) = RO2_
!      u(i) = U2_
!      p(i) = P2_
!      v(i) = V2_
!      else
!      r(i) = RO3_
!      u(i) = U3_
!      p(i) = P3_
!      v(i) = V3_
!     endif
!   endif
!  enddo

 do i = 1, ntri
 	if (cm(2,i)<=(y_-epsilon)) then
   	r(i) = RO1_
   	u(i) = U1_
   	p(i) = P1_
   	v(i) = V1_
   else
    if (cm(2,i)<=-cm(1,i)/5-epsilon .and. cm(2,i)>(y_-epsilon)) then
     r(i) = RO2_
     u(i) = U2_
     p(i) = P2_
     v(i) = V2_
     else
     r(i) = RO3_
     u(i) = U3_
     p(i) = P3_
     v(i) = V3_
    endif
  endif
 enddo

!   do i = 1, ntri
! 	  IF (cm(2,i).LT.BO1(cm(1,i))) THEN
! 	    R(I) = RO1_
! 	    U(I) = U1_
! 	    V(I) = V1_
! 	    P(I) = P1_
! 	  ENDIF
! 	  IF ((cm(2,i).LT.BO2(cm(1,i))).AND.(cm(2,i).GE.BO1(cm(1,i)))) THEN
! 	    R(I) = RO2_
! 	    U(I) = U2_
! 	    V(I) = V2_
! 	    P(I) = P2_
! 	  ENDIF
! 	  IF (cm(2,i).GE.BO2(cm(1,i))) THEN
! 	    R(I) = RO3_
! 	    U(I) = U3_
! 	    V(I) = V3_
! 	    P(I) = P3_
! 	  ENDIF
!   enddo



  t = 0
  step = 0
  do while (t<tmax)
    do i=1, ntri 
      ro_1(i)=r(i)
      ru_1(i)=u(i)*r(i)
      rv_1(i)=v(i)*r(i)
      e(i)=p(i)/(AGAM*r(i))
      en_1(i)=r(i)*(e(i)+0.5d0*((u(i))**2+(v(i))**2))
      
      ro(i) = ro_1(i)
      ru(i) = ru_1(i)
      rv(i) = rv_1(i)
      en(i) = en_1(i) 
      
      if ((r(i) <= 0.d0).or.(p(i) <= 0.d0)) then
        print *, 'E R R O R ! ! !'
        print *, 't = ', t
        print *, 'i = ', i
        print *, 'bnd = ', xbound(i)  
        print *, 'xci = ', centers(1,i)  
        print *, 'yci = ', centers(2,i)
        print *, 'ri = ', r(i)  
        print *, 'pi = ', p(i)  
        print *, 'ui = ', u(i)  
        print *, 'vi = ', v(i)  
        print *, 'ei = ', e(i) 
        stop
      endif
    enddo
    ! step 1/2
    do i = 1, ntri
      rtmp = 0
      utmp = 0
      vtmp = 0
      etmp = 0
      is_bound = .false.
      do j =0, 2
        call ENO(i, j, rnx, rny, la)
!         call RIM(RI,EI,PI,UIn,VIn,WIn,CI)
!            
!         UI = rnx*UIn-rny*VIn
!         VI = rny*UIn+rnx*VIn
!         WI = WIn
!                 
!         FR= UI*RI
!         FU= FR*UI + PI
!         FV= FR*VI
!         FE= FR*(EI+0.5d0*(UI**2+VI**2))+ PI*UI
!          
!         GR= VI*RI
!         GU= GR*UI
!         GV= GR*VI + PI
!         GE= GR*(EI+0.5d0*(UI**2+VI**2))+ PI*VI

        alpha = max(abs(-rny*ub+rnx*vb) + dsqrt(GAM*pb/rb), & 
                    abs(-rny*ue+rnx*ve) + dsqrt(GAM*pe/re), & 
                    abs(rnx*ub+rny*vb)  + dsqrt(GAM*pb/rb), & 
                    abs(rnx*ue+rny*ve)  + dsqrt(GAM*pe/re) )
				EB = PB/(AGAM*RB)        
				EE = PE/(AGAM*RE)        
				        

        FR= (UB*RB+UE*RE              -alpha*(RE-RB))*5.d-1
        FU= (RB*UB*UB+PB+RE*UE*UE+PE  -alpha*(RE*UE-RB*UB))*5.d-1
        FV= (RB*UB*VB+PB+RE*UE+VE+PE  -alpha*(RE*VE-RB*VB))*5.d-1
        FE= (RE*UE*(EE+0.5d0*((UE)**2+(VE)**2))+ PE*UE + &
             RB*UB*(EB+0.5d0*((UB)**2+(VB)**2))+ PB*UB - &
             alpha*( &
             RE*(EE+0.5d0*((UE)**2+(VE)**2)) -  &
             RB*(EB+0.5d0*((UB)**2+(VB)**2))   )   )*5.d-1
         
        GR= (VB*RB+VE*RE             -alpha*(RE-RB))*5.d-1
        GU= (RB*VB*UB+RE*VE*UE       -alpha*(RE*UE-RB*UB))*5.d-1
        GV= (RB*VB*VB+PB+RE*VE*VE+PE -alpha*(RE*VE-RB*VB))*5.d-1
        GE= (RE*VE*(EE+0.5d0*((UE)**2+(VE)**2))+ PE*VE + &
             RB*VB*(EB+0.5d0*((UB)**2+(VB)**2))+ PB*VB - &
             alpha*( &
             RE*(EE+0.5d0*((UE)**2+(VE)**2)) -  &
             RB*(EB+0.5d0*((UB)**2+(VB)**2))   )   )*5.d-1
        
        rtmp = rtmp + (FR*rnx + GR*rny)*la
        utmp = utmp + (FU*rnx + GU*rny)*la
        vtmp = vtmp + (FV*rnx + GV*rny)*la
        etmp = etmp + (FE*rnx + GE*rny)*la
      enddo
      atmp = tau/s(i)

      ro_1(i) = ro_1(i)-rtmp*atmp
      ru_1(i) = ru_1(i)-utmp*atmp
      rv_1(i) = rv_1(i)-vtmp*atmp
      en_1(i) = en_1(i)-etmp*atmp
    enddo
    do i = 1, ntri
      r(i)= ro_1(i)
      u(i)= ru_1(i)/ro_1(i)
      v(i)= rv_1(i)/ro_1(i)
      e(i)= en_1(i)/ro_1(i)-0.5d0*(u(i)*u(i)+v(i)*v(i))
      p(i)= AGAM*e(i)*ro_1(i)

!------------------------------
      if ((r(i) <= 0.d0).or.(p(i) <= 0.d0)) then
        print *, 'E R R O R ! ! !    (step 1/2)'
        print *, 't = ', t
        print *, 'i = ', i
        print *, 'xci = ', centers(1,i)  
        print *, 'yci = ', centers(2,i)
        print *, 'ri = ', r(i)  
        print *, 'pi = ', p(i)  
        print *, 'ui = ', u(i)  
        print *, 'vi = ', v(i)  
        print *, 'ei = ', e(i) 
        print *, '---'
        print *, 'rui = ', ru_1(i)  
        print *, 'rvi = ', rv_1(i)  
        print *, 'eni = ', en_1(i) 
        print *, '---'
        print *, 'triangle vertices:'
        do j=1,3
          xa = x(1, tri_v(j, i))
          ya = x(2, tri_v(j, i))
          print *, j, ': ', xa, ya, xbound(tri_v(j, i))
        enddo
        print *, 'Press ENTER...'
        read(*,*)
      endif
!--------------------------------------------

    enddo
    do i = 1, ntri
      ro_1(i) = (ro_1(i)+ro(i))/2.d0
      ru_1(i) = (ru_1(i)+ru(i))/2.d0
      rv_1(i) = (rv_1(i)+rv(i))/2.d0
      en_1(i) = (en_1(i)+en(i))/2.d0
    enddo
    ! step 2/2
    do i = 1, ntri
      rtmp = 0
      utmp = 0
      vtmp = 0
      etmp = 0
      is_bound = .false.
      do j =0, 2
        call ENO(i, j, rnx, rny, la)
!         call RIM(RI,EI,PI,UIn,VIn,WIn,CI)
! 
!       if ((RI.ne.RI).or.(PI.ne.PI).or.(UI.ne.UI).or.(VI.ne.VI)) then
!         print *, 'E R R O R ! ! !    (step 2/2, RIM)'
!         print *, 't = ', t
!         print *, 'i = ', i
!         print *, 'bnd = ', xbound(i)  
!         print *, 'xci = ', centers(1,i)  
!         print *, 'yci = ', centers(2,i)
!         print *, 'ri = ', RI  
!         print *, 'pi = ', pi  
!         print *, 'ui = ', ui  
!         print *, 'vi = ', vi  
!         print *, 'ei = ', ei 
!         print *, '---'
!         print *, 'RB = ', RB
!         print *, 'PB = ', PB
!         print *, 'UB = ', UB
!         print *, 'VB = ', VB
!         print *, 'RE = ', RE
!         print *, 'PE = ', PE
!         print *, 'UE = ', UE
!         print *, 'VE = ', VE
!         stop
!       endif
!         UI = rnx*UIn-rny*VIn
!         VI = rny*UIn+rnx*VIn
!         WI = WIn
!                 
!         FR= UI*RI
!         FU= FR*UI + PI
!         FV= FR*VI
!         FE= FR*(EI+0.5d0*(UI**2+VI**2))+ PI*UI
!          
!         GR= VI*RI
!         GU= GR*UI
!         GV= GR*VI + PI
!         GE= GR*(EI+0.5d0*(UI**2+VI**2))+ PI*VI

        alpha = max(abs(-rny*ub+rnx*vb) + dsqrt(GAM*pb/rb), & 
                    abs(-rny*ue+rnx*ve) + dsqrt(GAM*pe/re), & 
                    abs(rnx*ub+rny*vb)  + dsqrt(GAM*pb/rb), & 
                    abs(rnx*ue+rny*ve)  + dsqrt(GAM*pe/re) )
				EB = PB/(AGAM*RB)        
				EE = PE/(AGAM*RE)        
				        

        FR= (UB*RB+UE*RE              -alpha*(RE-RB))*5.d-1
        FU= (RB*UB*UB+PB+RE*UE*UE+PE  -alpha*(RE*UE-RB*UB))*5.d-1
        FV= (RB*UB*VB+PB+RE*UE+VE+PE  -alpha*(RE*VE-RB*VB))*5.d-1
        FE= (RE*UE*(EE+0.5d0*((UE)**2+(VE)**2))+ PE*UE + &
             RB*UB*(EB+0.5d0*((UB)**2+(VB)**2))+ PB*UB - &
             alpha*( &
             RE*(EE+0.5d0*((UE)**2+(VE)**2)) -  &
             RB*(EB+0.5d0*((UB)**2+(VB)**2))   )   )*5.d-1
         
        GR= (VB*RB+VE*RE             -alpha*(RE-RB))*5.d-1
        GU= (RB*VB*UB+RE*VE*UE       -alpha*(RE*UE-RB*UB))*5.d-1
        GV= (RB*VB*VB+PB+RE*VE*VE+PE -alpha*(RE*VE-RB*VB))*5.d-1
        GE= (RE*VE*(EE+0.5d0*((UE)**2+(VE)**2))+ PE*VE + &
             RB*VB*(EB+0.5d0*((UB)**2+(VB)**2))+ PB*VB - &
             alpha*( &
             RE*(EE+0.5d0*((UE)**2+(VE)**2)) -  &
             RB*(EB+0.5d0*((UB)**2+(VB)**2))   )   )*5.d-1
        
        rtmp = rtmp + (FR*rnx + GR*rny)*la
        utmp = utmp + (FU*rnx + GU*rny)*la
        vtmp = vtmp + (FV*rnx + GV*rny)*la
        etmp = etmp + (FE*rnx + GE*rny)*la
      enddo
      atmp = 0.5d0*tau/s(i)
      ro_1(i) = ro_1(i)-rtmp*atmp
      ru_1(i) = ru_1(i)-utmp*atmp
      rv_1(i) = rv_1(i)-vtmp*atmp
      en_1(i) = en_1(i)-etmp*atmp
    enddo
    do i = 1, ntri
      r(i)= ro_1(i)
      u(i)= ru_1(i)/ro_1(i)
      v(i)= rv_1(i)/ro_1(i)
      e(i)= en_1(i)/ro_1(i)-0.5d0*(u(i)*u(i)+v(i)*v(i))
      p(i)= AGAM*e(i)*ro_1(i)

!-------------------------------------------------------------
      if ((r(i).ne.r(i)).or.(p(i).ne.p(i)).or.(u(i).ne.u(i)).or.(v(i).ne.v(i))) then
        print *, 'E R R O R ! ! !    (step 2/2)'
        print *, 't = ', t
        print *, 'i = ', i
        print *, 'bnd = ', xbound(i)  
        print *, 'xci = ', centers(1,i)  
        print *, 'yci = ', centers(2,i)
        print *, 'ri = ', r(i)  
        print *, 'pi = ', p(i)  
        print *, 'ui = ', u(i)  
        print *, 'vi = ', v(i)  
        print *, 'ei = ', e(i) 
        print *, '---'
        print *, 'rui = ', ru_1(i)  
        print *, 'rvi = ', rv_1(i)  
        print *, 'eni = ', en_1(i) 
        print *, '---'
        print *, 'triangle vertices:'
        do j=1,3
          xa = x(1, tri_v(j, i))
          ya = x(2, tri_v(j, i))
          print *, j, ': ', xa, ya, xbound(tri_v(j, i))
        enddo
        print *, 'Press ENTER...'
        read(*,*)
      endif
!-------------------------------------------------------------

    enddo
    t = t+tau
    step = step + 1

!   +==================================================+
!   |   Вывод данных в файл и на экран                 |
!   +--------------------------------------------------+
    if (mod(step, SAVE_STEP) == 0) then
      write(fname,30) step
      open(1, file=fname)
      open(2, file='c_'//fname)
      open(3, file='me_r_'//fname)
      write(2,*) ntri+301*21
      write(3,*) '{'
      do i = 1, ntri
        do j = 1, 3
          write(1,20) x(1, tri_v(j,i)), x(2, tri_v(j,i)), r(i),u(i),v(i),p(i)
        enddo
        write(1,20) x(1, tri_v(1,i)), x(2, tri_v(1,i)), r(i),u(i),v(i),p(i)
        write(1,*)
        write(2,20) centers(1, i), centers(2,i), r(i),u(i),v(i),p(i)
        write(3,21) centers(1, i), centers(2,i), r(i)
      enddo
!      do i=0,300
!        do j=0,20
!          write(2,20) 0.3d0+i*(3.d0-0.3d0)/300.d0, j*0.2d0/20, 0.d0,0.d0,0.d0,0.d0
!        enddo
!      enddo
      write(3,*) '}'
      close(1)
      close(2)
      close(3)
      write(*,10) fname
    endif
    if (mod(step, PRINT_STEP) == 0) then
      write(*,11) t, step
    endif
!   +==================================================+

  enddo

  
  deallocate(x)
  deallocate(xbound)
  deallocate(tri_v)
  deallocate(tri_neigh)
  deallocate(r)
  deallocate(u)
  deallocate(v)
  deallocate(p)
  deallocate(ro_1)
  deallocate(ru_1)
  deallocate(rv_1)
  deallocate(e)
  deallocate(en_1)
  
  stop
 
 
 
 
 10 format ('saved file: ', A) 
 11 format ('time = ', F10.7, '  |  step = ', I7) 
 20 format (6E40.30)
 21 format ('{',E40.30,',',E40.30,',',E40.30,'},')
 30 format ('result.',I0,'.txt')
end program tri   





subroutine ENO(i, j, rnx, rny, la)
	include 'triLRIM.comm.inc'

  real    (kind = 8) :: rnx, rny, la, sina, cosa, lc
	integer (kind = 4) :: i, j, k
	real    (kind = 8) :: RBx,RBy,PBx,PBy,UBx,UBy,VBx,VBy, &
	                      REx,REy,PEx,PEy,UEx,UEy,VEx,VEy

        k = tri_neigh(mod(j+2,3)+1, i)
        xa = x(1, tri_v(mod(j,3)+1, i))
        ya = x(2, tri_v(mod(j,3)+1, i))
        xb = x(1, tri_v(mod(j+1,3)+1, i))
        yb = x(2, tri_v(mod(j+1,3)+1, i))
        xc = (xa+xb)/2.d0
        yc = (ya+yb)/2.d0
        if (k < 0) then
          if (dabs(ya-yb)<=epsilon)  then
            xk=centers(1,i)
            yk=2*ya-centers(2,i)
          else
            xk=2.d0*xb-centers(1,i)
            yk=centers(2,i)
          endif
        else
          xk=centers(1,k)
          yk=centers(2,k)
        endif
        la  = dsqrt((xa-xb)**2+(ya-yb)**2)
        rnx = (xk-centers(1,i))
        rny = (yk-centers(2,i))
        lc  = dsqrt(rnx**2+rny**2)
        rnx = rnx/lc
        rny = rny/lc
        cosa = rnx
        sina = rny
        
        call interpolation(RBx,RBy,PBx,PBy,UBx,UBy,VBx,VBy,i)

        RB = r(i)+RBx*(xc-cm(1,i))+RBy*(yc-cm(2,i))
        PB = p(i)+PBx*(xc-cm(1,i))+PBy*(yc-cm(2,i))
        UB = u(i)+UBx*(xc-cm(1,i))+UBy*(yc-cm(2,i))
        VB = v(i)+VBx*(xc-cm(1,i))+VBy*(yc-cm(2,i))
        WB = 0
        CnB= 0

        if (k < 0) then
          if (dabs(ya-yb)<=epsilon)  then
            if (dabs(ya-ymin)<=epsilon) then
              RE = RO1_
              PE = P1_
              UE = U1_
              VE = V1_
            else
              if (dabs(ya-ymax)<=epsilon) then
                RE =  RB
                UE =  UB
                VE = -VB
                PE =  PB
              else
                RE=  RB
                UE=  UB
                VE= -VB
                PE=  PB
              endif
            endif
          else
            RE=  RB
            UE= -UB
            VE=  VB
            PE=  PB
          endif
        else
        
          call interpolation(REx,REy,PEx,PEy,UEx,UEy,VEx,VEy,k)
          RE = r(k)+REx*(xc-cm(1,k))+REy*(yc-cm(2,k))
          PE = p(k)+PEx*(xc-cm(1,k))+PEy*(yc-cm(2,k))
          UE = u(k)+UEx*(xc-cm(1,k))+UEy*(yc-cm(2,k))
          VE = v(k)+VEx*(xc-cm(1,k))+VEy*(yc-cm(2,k))
        
        endif
        WE= 0
        CnE=0
        
        utmp =  cosa*UB+sina*VB
        vtmp = -sina*UB+cosa*VB
        UB = utmp
        VB = vtmp
        
        utmp =  cosa*UE+sina*VE
        vtmp = -sina*UE+cosa*VE
        UE = utmp
        VE = vtmp
        
        if (RB .NE. RB) goto 333
        if (PB .NE. PB) goto 333
        if (UB .NE. UB) goto 333
        if (VB .NE. VB) goto 333
        if (RE .NE. RE) goto 333
        if (PE .NE. PE) goto 333
        if (UE .NE. UE) goto 333
        if (VE .NE. VE) goto 333
        
        
        return
        
 333    print *, '===== ENO ====='
        print *, 'RB = ', RB
        print *, 'PB = ', PB
        print *, 'UB = ', UB
        print *, 'VB = ', VB
        print *, 'RE = ', RE
        print *, 'PE = ', PE
        print *, 'UE = ', UE
        print *, 'VE = ', VE
        print *
				stop        
end subroutine ENO

subroutine interpolation(Rx,Ry,Px,Py,Ux,Uy,Vx,Vy,ind)
  include 'triLRIM.comm.inc'
	real (kind = 8) :: Rx,Ry,Px,Py,Ux,Uy,Vx,Vy
	real (kind = 8) :: Rx1,Ry1,Px1,Py1,Ux1,Uy1,Vx1,Vy1
	real (kind = 8) :: maxu,maxr,maxp,maxv,cosu,cosv,cosp,cosr,det
	real (kind = 8), dimension(0:2) :: xk1,yk1,rk1,pk1,uk1,vk1
  integer (kind = 4) :: i,j,k,k1,k2,ind 

	do j = 0, 2
	 	k = tri_neigh(mod(j+2,3)+1, ind)
	 
	 	xa = x(1, tri_v(mod(j,3)+1, ind))
	 	ya = x(2, tri_v(mod(j,3)+1, ind))
	 	xb = x(1, tri_v(mod(j+1,3)+1, ind))
	 	yb = x(2, tri_v(mod(j+1,3)+1, ind))
		if (k < 0) then
			if (dabs(ya-yb)<=epsilon)  then
			  xk1(j) = cm(1,ind)
			  yk1(j) = 2*ya-cm(2,ind)
			  if (dabs(ya-ymin)<=epsilon) then
			    rk1(j) = RO1_
			    pk1(j) = P1_
			    uk1(j) = U1_
			    vk1(j) = V1_
			  else
			    rk1(j)=  r(ind) 
			    uk1(j)=  u(ind)
			    vk1(j)= -v(ind)
			    pk1(j)=  p(ind)
			  endif
			else
			  xk1(j) = 2.d0*xa-cm(1,ind)
			  yk1(j) = cm(2,ind)
				rk1(j)=  r(ind)
			  uk1(j)= -u(ind)
			  vk1(j)=  v(ind)
			  pk1(j)=  p(ind)
			endif
	  else
		  xk1(j)=cm(1,k)
		  yk1(j)=cm(2,k)
		  rk1(j)=r(k)
		  uk1(j)=u(k)
		  vk1(j)=v(k)
		  pk1(j)=p(k)
	 	endif
	enddo
  maxr = 0.d0
  maxp = 0.d0
  maxu = 0.d0
  maxv = 0.d0  
  
	do j = 0, 2
 		k1 = mod(j,3)
 		k2 = mod(j+1,3)
 		!U
 		det = ((xk1(k1)-cm(1,ind))*(yk1(k2)-cm(2,ind))-(yk1(k1)-cm(2,ind))*(xk1(k2)-cm(1,ind)))
 		Ux1 = ((uk1(k1)-u(ind))*(yk1(k2)-cm(2,ind))-(yk1(k1)-cm(2,ind))*(uk1(k2)-u(ind)))/det
		Uy1 = ((xk1(k1)-cm(1,ind))*(uk1(k2)-u(ind))-(uk1(k1)-u(ind))*(xk1(k2)-cm(1,ind)))/det
		cosu = 1/dsqrt(1+Ux1*Ux1+Uy1*Uy1)
		!R
 		Rx1 = ((rk1(k1)-r(ind))*(yk1(k2)-cm(2,ind))-(yk1(k1)-cm(2,ind))*(rk1(k2)-r(ind)))/det
		Ry1 = ((xk1(k1)-cm(1,ind))*(rk1(k2)-r(ind))-(rk1(k1)-r(ind))*(xk1(k2)-cm(1,ind)))/det
		cosr = 1/dsqrt(1+Rx1*Rx1+Ry1*Ry1)
		!P
		Px1 = ((pk1(k1)-p(ind))*(yk1(k2)-cm(2,ind))-(yk1(k1)-cm(2,ind))*(pk1(k2)-p(ind)))/det
		Py1 = ((xk1(k1)-cm(1,ind))*(pk1(k2)-p(ind))-(pk1(k1)-p(ind))*(xk1(k2)-cm(1,ind)))/det
		cosp = 1/dsqrt(1+Px1*Px1+Py1*Py1)
		!V
		Vx1 = ((vk1(k1)-v(ind))*(yk1(k2)-cm(2,ind))-(yk1(k1)-cm(2,ind))*(vk1(k2)-v(ind)))/det
		Vy1 = ((xk1(k1)-cm(1,ind))*(vk1(k2)-v(ind))-(vk1(k1)-v(ind))*(xk1(k2)-cm(1,ind)))/det
		cosr = 1/dsqrt(1+Vx1*Vx1+Vy1*Vy1)
		
		if (cosu>=maxu) then
		  maxu = cosu
		  Ux = Ux1
		  Uy = Uy1
 		endif
 		if (cosr>=maxr) then
		  maxr = cosr
		  Rx = Rx1
		  Ry = Ry1
 		endif
 		if (cosp>=maxp) then
		  maxp = cosp
		  Px = Px1
		  Py = Py1
 		endif
 		if (cosv>=maxv) then
		  maxv = cosv
		  Vx = Vx1
		  Vy = Vy1
 		endif
 	enddo
!  	
!  	Rx = 0.d0
!  	Ry = 0.d0
!  	px = 0.d0
!  	py = 0.d0
!  	ux = 0.d0
!  	uy = 0.d0
!  	vx = 0.d0
!  	vy = 0.d0
 	
end subroutine

     
!==========================================================
!    Nikichine
!    module for tube.for /Zmitrenko/
!==========================================================

      SUBROUTINE RIM (RI,EI,PI,UI,VI,WI,CI)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/A/ RB,PB,UB,VB,WB,CnB,RE,PE,UE,VE,WE,CnE
      COMMON/GAMMA/ GAM,AGAM,BGAM,CGAM,DGAM,EGAM,FGAM,GGAM,HGAM,OGAM,&
               PGAM,QGAM,RGAM,SGAM,TGAM

      EPS=1.D-5
      CB= DSQRT(GAM*PB/RB)
      CE= DSQRT(GAM*PE/RE)
      EB= CB**2/SGAM
      EE= CE**2/SGAM
      RCB=RB*CB
      RCE=RE*CE
      DU=UB-UE
      IF (DU.LT.-2.*(CB+CE)/AGAM) THEN
         WRITE(*,'(1X,A)') ' ATTENTION!!!  VACUUM '
         RF=0.
         RS=0.
         EF=0.
         ES=0.
         SBL=UB-CB
         SFL=UB+2.*CB/AGAM
         SSL=UE-2.*CE/AGAM
         SEL=UE+CE
         GOTO 9
      ENDIF
      P=(PB*RCE+PE*RCB+DU*RCB*RCE)/(RCB+RCE)
    5  CONTINUE
           IF (P.LT.EPS) P=EPS

           PPB=P/PB
           IF(PB.GT.P) GO TO 1
           PKB=PGAM*PPB+OGAM
           ZNB=RCB*DSQRT(PKB)
           F1=(P-PB)/ZNB
           FS1=(QGAM*PPB+FGAM)/(RGAM*ZNB*PKB)
           GO TO 2
    1      ZFB=CB*PPB**OGAM
           F1=DGAM*(ZFB-CB)
           FS1=ZFB/(GAM*P)
    2      PPE=P/PE
           IF(PE.GT.P) GO TO 3
           PKE=PGAM*PPE+OGAM
           ZNE=RCE*DSQRT(PKE)
           F2=(P-PE)/ZNE
           FS2=(QGAM*PPE+FGAM)/(RGAM*ZNE*PKE)
           GO TO 4
    3      ZFE=CE*PPE**OGAM
           F2=DGAM*(ZFE-CE)
           FS2=ZFE/(GAM*P)
    4      DP=(DU-F1-F2)/(FS1+FS2)
           P=P+DP
       IF(DABS(DU-F1-F2).GT.EPS) GO TO 5


      PPB=P/PB
      PPE=P/PE

      ZFB=CB*PPB**OGAM
      ZFE=CE*PPE**OGAM
      IF(PB.GT.P) GO TO 6
      D=UB-DSQRT((TGAM*P+HGAM*PB)/RB)
      UBD=UB-D
      RUBD=RB*UBD
      RF=RUBD**2/(PB-P+RUBD*UBD)
      UF=D+RUBD/RF
      EF=P/(AGAM*RF)
      SBL=D
      SFL=D
      GO TO 7
    6 EF=ZFB**2/SGAM
      UF=UB+DGAM*(CB-ZFB)
      RF=P/(AGAM*EF)
      SBL=UB-CB
      SFL=UF-ZFB
    7 IF(PE.GT.P) GO TO 8
      D=UE+DSQRT((TGAM*P+HGAM*PE)/RE)
      UED=UE-D
      RUED=RE*UED
      RS=RUED**2/(PE-P+RUED*UED)
      US=D+RUED/RS
      ES=P/(AGAM*RS)
      SEL=D
      SSL=D
      GO TO 9
    8 ES=ZFE**2/SGAM
      US=UE-DGAM*(CE-ZFE)
      RS=P/(AGAM*ES)
      SSL=US+ZFE
      SEL=UE+CE
    9 CONTINUE

!     compute the interpolation value
      IF (SEL.LE.0.) THEN
         RI= RE
         EI= EE
         UI= UE
         VI= VE
         WI= WE
         CI= CnE
         GO TO 157
      END IF

      IF (SBL.GE.0.) THEN
         RI= RB
         EI= EB
         UI= UB
         VI= VB
         WI= WB
         CI= CnB
         GO TO 157
      END IF

      IF ((SSL.GE.0.).AND.(SFL.LE.0.)) THEN
         IF (US.GE.0.) THEN
            RI= RF
            EI= EF
            UI= UF
            VI= VB
            WI= WB
            CI= CnB
         ELSE
            RI= RS
            EI= ES
            UI= US
            VI= VE
            WI= WE
            CI= CnE
         END IF
         GO TO 157
      END IF

      IF (SFL.GT.0.) THEN
         UI= (UB+DGAM*GGAM*dSQRT(EB))/(1+DGAM)
         VI= VB
         WI= WB
         CI= CnB
         EI= (UI**2)/SGAM
         RI= RB*((EI/EB)**(1/AGAM))
       ELSE
         UI= (UE-DGAM*GGAM*dSQRT(EE))/(1+DGAM)
         VI= VE
         WI= WE
         CI= CnE
         EI= (UI**2)/SGAM
         RI= RE*((EI/EE)**(1/AGAM))
      END IF

  157   CONTINUE
      PI= AGAM*EI*RI

      RETURN
      END
 
subroutine ReadFile(fName)
  include 'triLRIM.comm.inc'
  integer (kind = 4) :: i,j,k
  character*30 fName
  real (kind = 8) la, lb, lc

  open(1, file='tube.1.node')
  read(1, *) nx, dim, na, nb
  if (nx > 0) then
    allocate(x(2, nx))
    allocate(xbound(nx))
    do i = 1, nx
      read(1,*) itmp, x(1,i), x(2,i), xbound(i)
    enddo
  endif
  close(1)
  open(1, file='tube.1.ele')
  read(1, *) ntri, na, nb
  if (ntri > 0) then
    allocate(tri_v(3, ntri))
    do i = 1, ntri
      read(1,*) j, tri_v(1,i), tri_v(2,i), tri_v(3,i)
    enddo
  endif
  close(1)
  open(1, file='tube.1.neigh')
  read(1, *) ntri, itmp
  if (ntri > 0) then
    allocate(tri_neigh(3, ntri))
    do i = 1, ntri
      read(1,*) j, tri_neigh(1,i), tri_neigh(2,i), tri_neigh(3,i)
    enddo
  endif
  close(1)
  allocate(centers(2,ntri))
  allocate(cm(2,ntri))
  allocate(s(ntri))
  do i = 1, ntri
    xa = x(1, tri_v(1, i))
    ya = x(2, tri_v(1, i))
    xb = x(1, tri_v(2, i))
    yb = x(2, tri_v(2, i))
    xc = x(1, tri_v(3, i))
    yc = x(2, tri_v(3, i))
    la = dsqrt((xa-xb)**2+(ya-yb)**2)
    lb = dsqrt((xc-xb)**2+(yc-yb)**2)
    lc = dsqrt((xa-xc)**2+(ya-yc)**2)
    rtmp = (la+lb+lc)/2.d0
    s(i) = dsqrt(rtmp*(rtmp-la)*(rtmp-lb)*(rtmp-lc))
    lc = (xa-xb)*(ya-yc)-(xa-xc)*(ya-yb)
    la = (xa*xa-xb*xb+ya*ya-yb*yb)/2.d0
    lb = (xa*xa-xc*xc+ya*ya-yc*yc)/2.d0
    centers(1, i) = (la*(ya-yc)-lb*(ya-yb))/lc
    centers(2, i) = (lb*(xa-xb)-la*(xa-xc))/lc
    cm(1, i)      = (xa+xb+xc)/3.d0
    cm(2, i)      = (ya+yb+yc)/3.d0
  enddo

  open(1, file='tube_tri.gp')
  do i = 1, ntri
    do j = 1, 3
      write(1,*) x(1, tri_v(j,i)), x(2, tri_v(j,i))
    enddo
    write(1,*) x(1, tri_v(1,i)), x(2, tri_v(1,i))
    write(1,*)
  enddo
  close(1)
  open(1, file='tube_vor.gp')
  do i = 1, ntri
    do j=0, 2
      k = tri_neigh(mod(j+2,3)+1, i)
      if (k>0) then
        write(1,*) centers(1, i), centers(2, i)
        write(1,*) centers(1, k), centers(2, k)
        write(1,*)
      endif
    enddo
  enddo
  close(1)

end subroutine ReadFile



!-----------------------------------------------------------
      REAL*8 FUNCTION BO1(VAR1)
      IMPLICIT REAL*8 (A-H,O-Z),INTEGER (I-N)
          BO1 = 1.D0
      RETURN
      END
!-----------------------------------------------------------
      REAL*8 FUNCTION BO2(VAR1)
        IMPLICIT REAL*8 (A-H,O-Z),INTEGER (I-N)
          IF ((VAR1.GT.3.D0).AND.(VAR1.LT.4.2D0)) THEN
            BO2 = 2.D0-SIN((VAR1-3.D0)*3.14D0/1.2D0)
          ELSE
            BO2 = 2.D0
          ENDIF
        RETURN
      END
!===========================================================
