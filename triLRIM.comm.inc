	implicit none
	real (kind = 8), pointer, dimension(:,:) :: x, xatr, centers, cm
	integer (kind = 4), pointer, dimension(:,:) :: tri_v, tri_neigh
	integer (kind = 4), pointer, dimension(:) :: xbound
	real (kind = 8), pointer, dimension(:) :: u1,s
	real (kind = 8) :: rtmp, t, xa,ya, xb,yb, xc,yc, xk,yk,rk,uk,vk,pk
	integer (kind = 4) :: nx, na, nb, dim, itmp, ntri, step
	logical :: is_bound
	  
	real (kind = 8), pointer, dimension(:) :: r,p,u,v, ro_1, ru_1,rv_1,e, en_1, &
	                                                       ro, ru,rv,en
	real (kind = 8) :: utmp,vtmp,etmp,FR,FU,FV,FE,GR,GU,GV,GE,RB,PB,UB,VB,WB,CnB,   &
	                   EB,RE,PE,UE,VE,WE,CnE,EE,RI,EI,PI,UI,VI,WI,UIn,VIn,WIn,CI
	real (kind = 8) :: GAM,AGAM,BGAM,CGAM,DGAM,EGAM,FGAM,GGAM,HGAM,OGAM,PGAM,QGAM,RGAM,SGAM,TGAM
	
	
	
	common /A/ RB,PB,UB,VB,WB,CnB,RE,PE,UE,VE,WE,CnE
	common /GAMMA/ GAM,AGAM,BGAM,CGAM,DGAM,EGAM,FGAM,GGAM,HGAM,OGAM,PGAM,QGAM,RGAM,SGAM,TGAM
	common /GAS_PARAMS/ r,p,u,v
	common /TRI_DATA/ ntri, nx, x, xatr, centers, cm, tri_v, tri_neigh, xbound, u1, s 
	



  real (kind = 8), parameter :: epsilon = 1.d-5,             &
                                tmax    = 2500.d-4,             &
                                tau     = 1.d-5,            &
                                m_pi    = 3.14159265358979d0
  
  integer, parameter :: PRINT_STEP = 10,    &
                        SAVE_STEP  = 500
 
!   real (kind = 8),parameter  :: RO1_ = 12.09D0, U1_ = 0.D0, P1_ = 2.152D+5, V1_ = 97.76D0
!   real (kind = 8),parameter  :: RO2_ = 6.037D0, U2_ = 0.d0, P2_ = 1.D+5,    V2_ = 0.d0
!   real (kind = 8),parameter  :: RO3_ = 1.198D0, U3_ = 0.d0, P3_ = 1.d+5,    V3_ = 0.d0

  real (kind = 8),parameter :: RO1_=17.6d0, U1_=0.d0, P1_=1.8d+5, V1_=155.556d0
  real (kind = 8),parameter :: RO2_=1.161d0, U2_=0.d0, P2_=1.d+5, V2_=0.d0	
  real (kind = 8),parameter :: RO3_=0.16d0, U3_=0.d0, P3_=1.d+5, V3_=0.d0


  real (kind = 8), parameter :: xmin = 0.d0,   xmax = 7.2D0
  real (kind = 8), parameter :: ymin = 0.D0,  ymax = 10.D0
  real (kind = 8), parameter :: y_ =-2.d-2 

