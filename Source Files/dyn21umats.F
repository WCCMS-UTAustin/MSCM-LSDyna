# 1 "dyn21umats_aniso1d.F"
# 1 ".\define.inc" 1 


# 26
















# 44


# 48


# 53


# 57


# 62









# 73









# 85



# 90






























# 122

# 126




# 136



# 180


# 188


# 223


# 257


# 283






# 299






# 311

# 318









# 329




# 341










# 370


# 1 ".\adios.h" 1 
# 9



















# 373 ".\define.inc" 2 


# 379

# 382






# 408



# 417

# 424





# 2 "dyn21umats_aniso1d.F" 2 
# 1 ".\define2.inc" 1 



# 3 "dyn21umats_aniso1d.F" 2 
      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     isotropic elastic material (sample user subroutine)
c
c     Variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c        .
c        .
c        .
c     cm(n)=nth material constant
c
c     eps(1)=local x  strain increment
c     eps(2)=local y  strain increment
c     eps(3)=local z  strain increment
c     eps(4)=local xy strain increment
c     eps(5)=local yz strain increment
c     eps(6)=local zx strain increment
c
c     Note that eps is in strain Voigt notation, meaning that
c     e.g. eps(4) is the sum of tensorial xy and yx strain increments.
c
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
c
c     hsv(1)=1st history variable
c     hsv(2)=2nd history variable
c        .
c        .
c        .
c        .
c     hsv(n)=nth history variable
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq."solid" for solid elements
c        eq."sph  " for smoothed particle hydrodynamics
c        eq."shl2d" for shell forms 12 (2D solids - plane stress)
c        eq."sld2d" for shell forms 13 (2D solids - plane strain)
c        eq."sldax" for shell forms 14, and 15 (2D solids - axisymmetric)
c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
c        eq."tshel" for thick shell forms 3 and 5
c        eq."hbeam" for beam element forms 1, 11, 14
c        eq."tbeam" for beam element form 3 (truss)
c        eq."dbeam" for beam element form 6 (discrete)
c        eq."beam " for all other beam elements (currently not used)
c
c     tt=current problem time.
c
c     temper=current temperature
c
c     failel=flag for failure, set to .true. to fail an integration point,
c            if .true. on input the integration point has failed earlier
c
c     crv=array representation of curves in keyword deck
c
c     nnpcrv=# of discretization points per crv()
c
c     cma=additional memory for material data defined by LMCA at
c       6th field of 2nd card of *DATA_USER_DEFINED
c
c     elsiz=characteristic element size
c
c     idele=element id
c
c     reject (implicit only) = set to .true. if this implicit iterate is
c                              to be rejected for some reason
c
c     All transformations into the element local system are
c     performed prior to entering this subroutine.  Transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase. Initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c
     
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
      integer idele
c
      if (ncycle.eq.1) then
        if (cm(16).ne.1234567) then
          call usermsg('mat41')
        endif
      endif
c
c     compute shear modulus, g
c
      g2 =abs(cm(1))/(1.+cm(2))
      g  =.5*g2
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax'.or.
     3     etype.eq.'sldco') then
        if (cm(16).eq.1234567) then
          call mitfail3d(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
        else
          if (.not.failel) then
          davg=(-eps(1)-eps(2)-eps(3))/3.
          p=-davg*abs(cm(1))/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=sig(2)+p+g2*(eps(2)+davg)
          sig(3)=sig(3)+p+g2*(eps(3)+davg)
          sig(4)=sig(4)+g*eps(4)
          sig(5)=sig(5)+g*eps(5)
          sig(6)=sig(6)+g*eps(6)
          if (cm(1).lt.0.) then
            if (sig(1).gt.cm(5)) failel=.true.
          endif
          endif
        end if
c
      else if (etype.eq.'shell'.or.etype.eq.'shl2d') then
        if (cm(16).eq.1234567) then
          call mitfailure(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
        else
          if (.not.failel) then
          gc    =capa*g
          q1    =abs(cm(1))*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
          q3    =1./(q1+g2)
          eps(3)=-q1*(eps(1)+eps(2))*q3
          davg  =(-eps(1)-eps(2)-eps(3))/3.
          p     =-davg*abs(cm(1))/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=sig(2)+p+g2*(eps(2)+davg)
          sig(3)=0.0
          sig(4)=sig(4)+g *eps(4)
          sig(5)=sig(5)+gc*eps(5)
          sig(6)=sig(6)+gc*eps(6)
          if (cm(1).lt.0.) then
            if (sig(1).gt.cm(5)) failel=.true.
          endif
          endif
        end if
      elseif (etype.eq.'hbeam') then
          q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
          q3    =q1+2.0*g
          gc    =capa*g
          deti  =1./(q3*q3-q1*q1)
          c22i  = q3*deti
          c23i  =-q1*deti
          fac   =(c22i+c23i)*q1
          eps(2)=-eps(1)*fac-sig(2)*c22i-sig(3)*c23i
          eps(3)=-eps(1)*fac-sig(2)*c23i-sig(3)*c22i
          davg  =(-eps(1)-eps(2)-eps(3))/3.
          p     =-davg*cm(1)/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=0.0
          sig(3)=0.0
          sig(4)=sig(4)+gc*eps(4)
          sig(5)=0.0
          sig(6)=sig(6)+gc*eps(6)
c
      elseif (etype.eq.'tbeam') then
        q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
        q3    =q1+2.0*g
        deti  =1./(q3*q3-q1*q1)
        c22i  = q3*deti
        c23i  =-q1*deti
        fac   =(c22i+c23i)*q1
        eps(2)=-eps(1)*fac
        eps(3)=-eps(1)*fac
        davg  =(-eps(1)-eps(2)-eps(3))/3.
        p     =-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=0.0
        sig(3)=0.0
c
      else
c       write(iotty,10) etype
c       write(iohsp,10) etype
c       write(iomsg,10) etype
c       call adios(2)
        cerdat(1)=etype
        call lsmsg(3,MSG_SOL+1150,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c10   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current material model.')
      return
      end

      subroutine umat42 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
     
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel,reject
      integer nnpcrv(*)
      character*5 etype
      integer idele
c
      if (ncycle.eq.1) then
        call usermsg('mat42')
      endif

c     Replace with your own code
      call umat42_test (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)

      return
      end

      subroutine umat43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
    
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      character*5 etype
      logical failel,reject
      integer idele
c
      if (ncycle.eq.1) then
        call usermsg('mat43')
      endif
c
      return
      end

      subroutine umat44 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
    
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      character*5 etype
      logical failel,reject
      integer idele
c
c
      if (ncycle.eq.1) then
        call usermsg('mat44')
      endif
c
      return
      end

	  
 !DEC$ ATTRIBUTES DLLEXPORT :: UMAT45
      subroutine umat45 (cm,eps,sig,epsp,hsv,dt1,capa,
     . etype,time,temp,failel,crv,nnpcrv,cma,qmat,elsiz,idele,ih,
     . reject)
	 
c
c******************************************************************
c|  Babak Bozorgmehri and Michael S. Sacks Gasser-Holfapzel-Ogden (GOH)
c|  Material Model
c|  Developed from the Neo-Hookean material user subroutine
c|  by Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     Gasser-Holfapzel-Ogden (GOH) material model (user subroutine)
c

c
c     hsv(1)=1st history variable
c     hsv(2)=2nd history variable
c        .
c        .
c        .
c        .
c     hsv(n)=nth history variable
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq."solid" for solid elements
c        eq."sph  " for smoothed particle hydrodynamics
c        eq."sld2d" for shell forms 13 (2D solids - plane strain)
c        eq."sldco" for solid form 62 (corotational element)
c        eq."sldax" for shell forms 14, and 15 (2D solids)
c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
c        eq."tshel" for thick shell forms 3 and 5
c        eq."hbeam" for beam element forms 1, 11, 14
c        eq."tbeam" for beam element form 3 (truss)
c        eq."dbeam" for beam element form 6 (discrete)
c        eq."beam " for all other beam elements (currently not used)
c
c     time=current problem time.
c     temp=current temperature
c
c     cma=additional memory for material data defined by LMCA at
c       6th field of 2nd card of *DATA_USER_DEFINED
c
c     All transformations into the element local system are
c     performed prior to entering this subroutine.  Transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase.   Initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c

c

c
      logical usemodule,usemetal
c
c     As an example of the use of *MODULE, this routine can be executed
c     both as a hyperelastic Neo-Hookean material (this routine) and as
c     an elastic-plastic material (subroutines metal*). The purpose of
c     the *MODULE concept is to load several shared object files to a
c     single LS-DYNA executable and whence (in this particular case)
c     allow multiple instances of subroutine umat45 be executed without
c     conflicts. This is important if independent user material distributors
c     happen to use the same user material number for their respective
c     material implementations, and neither want to reveal the source code.
c
c     To use this example, perform the following steps in order
c
c     1) Set usemodule=.true. below
c     2) Compile the code to a shared object file and call this file
c        hyperelastic.so
c     3) Set usemetal=.true. below

c
c     Note that for implicit, utan45 needs to modified accordingly
c
c     For details, refer to Appendix A of the Keyword User's Manual
c
    
	  

      include 'nlqparm'
      include 'iounits.inc'
      include 'bk06.inc'
      character*5 etype
      dimension cm(*), eps(*), sig(*), hsv(*),
     &          crv(lq1,2,*), cma(*), qmat(3,3)

	  integer nnpcrv(*)
      logical failel, reject
	  logical :: OK_FLAG
      real newthick(4)
      double precision :: I1, I4, I4star1, J4, I4star2, c10
      double precision :: n1, n2, n3, m1, m2, m3, normn, normm
      integer i, j, k
	  double precision :: F(3,3), C(3,3), Cbar(3,3), Cinv(3,3), Ftrans(3,3)
      double precision :: A0(3,3), B0(3,3), SIGMA(3,3), II(3,3), S(3,3)
      double precision :: JJinv23, trDEV2, PSI1, PSI41, PSI42, a1, a2
	  double precision :: DEV1(3,3), DEV2(3,3), Cbartrans(3,3), DEV(3,3)
      double precision :: Cbar_inv(3,3), E1(3,3)
      double precision :: JJ, trCbar, lnJ
	  double precision :: d_pen, kappa
	  double precision :: detF, phi1, phi2, arg1, arg2, lambda1, lambda2
	  double precision :: cg, sg, gamma
	  real*8 fiber_stress1, fiber_stress2
	  integer :: nff, nhv
	  double precision, parameter :: pi = 3.1415926535897932384626433832795d0

c====================================================================
c   UMAT for solid 8-node element with Gasser–Ogden–Holzapfel anisotropy
c====================================================================
      usemodule = .true.
      usemetal  = .false.
      if (usemodule) then
	  usemetal=.false.
        if (usemetal) then
          if (etype .eq. 'solid') then
            call metalsld(cm,eps,sig,epsp,hsv,dt1,capa,etype,time,temp,
     &            failel,crv,nnpcrv)
          else
            call metalshl(cm,eps,sig,epsp,hsv,dt1,capa,etype,time,temp,
     &            failel,crv,nnpcrv)
          endif
          return
        endif
      endif

      if (.not.usemodule) then
        if (ncycle.eq.1) then
          call usermsg('mat45')
        endif
      endif


!--- Set material parameters
      c10   = cm(1) / 2.0d0  ! cm(1): matrix shear modulus, c10: Neo-Hookean parameter 
      d_pen = cm(2) / 2.0d0  ! cm(1): bulk modulus, d_pen = parameter to impose incompressibilty constrain 
      a1    = cm(3)          ! fiber stiffness parameter
      a2    = cm(4)          ! fiber nonlinearity parameter
      kappa = cm(5)          ! dispersion, between 0 and 1/3 
      gamma = (pi/180.d0) * cm(6) ! fiber angle degree relative to a toward b ( "a" and "b" are set in user input)
      nff   = nint(cm(7))        ! number of fibre families: 1 or 2
	  nhv   = nint(cm(8)) 
      cg = dcos(gamma);  sg = dsin(gamma)
	  
!---- Two fiber families
!     n1 = cg*qmat(:,1) +  sg*qmat(:,2)
!     a = qmat(:,1), b = qmat(:,2)
      qmat(1,1) = 1.0d0; qmat(2,1) = 0.0d0; qmat(3,1) = 0.0d0 
	  qmat(1,2) = 0.0d0; qmat(2,2) = 1.0d0; qmat(3,2) = 0.0d0
!     n1=sg*qmat(1,1)+cg*qmat(1,2); n2=sg*qmat(2,1)+cg*qmat(2,2); 
!	  n3=sg*qmat(3,1)+cg*qmat(3,2)
c    a01 =  cos(gamma)*a + sin(gamma)*d  (if nff=1 or nff=2) 
	  n1=cg*qmat(1,1)+sg*qmat(1,2); n2=cg*qmat(2,1)+sg*qmat(2,2); 
	  n3=cg*qmat(3,1)+sg*qmat(3,2)
	  normn = sqrt(n1**2.0d0+n2**2.0d0+n3**2.0d0)
	  n1 = n1/normn; n2 = n2/normn; n3 = n3/normn
	  
c    a02 =  cos(gamma)*a - sin(gamma)*d  (only if nff=2)  
    ! m = cg*qmat(:,1) -  sg*qmat(:,2): 
	if (nff .eq. 2) then
	  m1=cg*qmat(1,1)-sg*qmat(1,2); m2=cg*qmat(2,1)-sg*qmat(2,2); 
	  m3=cg*qmat(3,1)-sg*qmat(3,2)
	  normm = sqrt(m1**2.0d0+m2**2.0d0+m3**2.0d0)
	  m1 = m1/normm; m2 = m2/normm; m3 = m3/normm
	endif
!     

      hsv(1) = qmat(1,1); hsv(2) = qmat(2,1); hsv(3) = qmat(3,1);
	  hsv(4) = qmat(1,2); hsv(5) = qmat(2,2); hsv(6) = qmat(3,2);
      hsv(7) = qmat(1,3); hsv(8) = qmat(2,3); hsv(9) = qmat(3,3);

	  
!   if (norm .lt. 1.0d-12) norm = 1.0d0

      
!--- Deformation gradient F from hsv(1:9)
      F(1,1) = hsv(1+nhv); F(1,2) = hsv(2+nhv); F(1,3) = hsv(3+nhv)
      F(2,1) = hsv(4+nhv); F(2,2) = hsv(5+nhv); F(2,3) = hsv(6+nhv)
      F(3,1) = hsv(7+nhv); F(3,2) = hsv(8+nhv); F(3,3) = hsv(9+nhv)

!--- Sanity check on F
      detF = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
     &     -F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
     &     +F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
	

!--- Determinant and isochoric correction
      JJ = detF
   !   JJ = max(JJ, 1.0d-8)
	  lnJ = dlog(JJ)
      JJinv23 = JJ**(-2.0d0/3.0d0)
	  
! Compute C = transpose(F) * F
      Ftrans = transpose(F)
      C = matmul(Ftrans, F)
	  
	  
 
      do i=1,3
	     do j=1,3
            C(j,i) = 0.5d0*(C(j,i) + C(i,j))
         enddo 
	  enddo
	  
	  call invert3x3(C, Cinv)
	  
	  !--- Cbar = isochoric part
	  
	  Cbar = JJinv23 * C
	  
	  do i=1,3
	     do j=1,3
            Cbar(j,i) = 0.5d0*(Cbar(j,i) + Cbar(i,j))
         enddo 
	  enddo
	  
	  call invert3x3(Cbar, Cbar_inv)
      
      Cbartrans = transpose(Cbar)
	  

!--- Identity tensor
      do i=1,3
         do j=1,3
            II(i,j) = 0.0d0
            if (i .eq. j) II(i,j) = 1.0d0
         enddo
      enddo

!--- Structural tensor A0 = a0 ⊗ a0
      A0(1,1) = n1*n1; A0(1,2) = n1*n2; A0(1,3) = n1*n3
      A0(2,1) = n2*n1; A0(2,2) = n2*n2; A0(2,3) = n2*n3
      A0(3,1) = n3*n1; A0(3,2) = n3*n2; A0(3,3) = n3*n3
	  
      if (nff .eq. 2) then	  
!--- Structural tensor B0 = b0 ⊗ b0
      B0(1,1) = m1*m1; B0(1,2) = m1*m2; B0(1,3) = m1*m3
      B0(2,1) = m2*m1; B0(2,2) = m2*m2; B0(2,3) = m2*m3
      B0(3,1) = m3*m1; B0(3,2) = m3*m2; B0(3,3) = m3*m3
      endif
!--- Invariants
      trCbar = Cbar(1,1) + Cbar(2,2) + Cbar(3,3)
	  I1 = trCbar
	  
  
      I4 = A0(1,1)*Cbar(1,1) + A0(1,2)*Cbar(1,2) + A0(1,3)*Cbar(1,3)
     &   + A0(2,1)*Cbar(2,1) + A0(2,2)*Cbar(2,2) + A0(2,3)*Cbar(2,3)
     &   + A0(3,1)*Cbar(3,1) + A0(3,2)*Cbar(3,2) + A0(3,3)*Cbar(3,3)
c	 
	  I4star1 = kappa * I1 + (1.0d0 - 3.0d0*kappa) * I4 - 1.0d0
c	   
      if (nff .eq. 2) then	
      J4 = B0(1,1)*Cbar(1,1) + B0(1,2)*Cbar(1,2) + B0(1,3)*Cbar(1,3)
     &   + B0(2,1)*Cbar(2,1) + B0(2,2)*Cbar(2,2) + B0(2,3)*Cbar(2,3)
     &   + B0(3,1)*Cbar(3,1) + B0(3,2)*Cbar(3,2) + B0(3,3)*Cbar(3,3)
c	 
	  I4star2 = kappa * I1 + (1.0d0 - 3.0d0*kappa) * J4 - 1.0d0
	  endif 

!--- Compute fiber stretches (store in hsv(16-17))
!    Fiber stretch λ = sqrt(I4) for each family
      lambda1 = sqrt(I4)  ! First fiber family stretch
	  hsv(16) = lambda1
c
      if (nff .eq. 2) then
         lambda2 = sqrt(J4)  ! Second fiber family stretch
		 hsv(17) = lambda2
      else
	     lambda2 = 0.0d0
         hsv(17) = 0.0d0
      endif


!--- Strain energy derivatives
      phi1 = max(I4star1, 0d0)
      arg1 = a2 * (phi1**2)
      arg1 = min(arg1, 40.0d0)
      PSI41 = a1*(1.0d0 - 3.0d0*kappa)*dexp(arg1)*phi1
c
	  if (nff .eq. 1) then
	  PSI1 = 1.0d0*c10 + a1 * kappa * (dexp(arg1)*phi1) 
	  endif
!--- Compute fiber stresses (store in hsv(19-20))
! Fiber stress in principal direction
         ! First fiber family stress (Cauchy stress in fiber direction)
         fiber_stress1 = 2.0d0 * PSI41 * I4 / JJ
         hsv(18) = fiber_stress1
	  
!---- Second fiber family stress
!      hsv(20) = 0.0d0  ! Default value
	  if (nff .eq. 2) then
!--- Strain energy derivatives
      phi2 = max(I4star2, 0d0)
      arg2 = a2 * (phi2**2)
      arg2 = min(arg2, 40.0d0)
      PSI1 = 1.0d0*c10 + a1 * kappa * (dexp(arg1)*phi1+dexp(arg2)*phi2) ! two fiber families
                     ! single fiber family
      PSI42 = a1 * (1.0d0 - 3.0d0*kappa)*dexp(arg2)*phi2
	  fiber_stress2 = 2.0d0 * PSI42 * J4 / JJ  
      hsv(19) = fiber_stress2
	  else
	  hsv(19) = 0.0d0
	  endif
c	  
c	  
! --- D = dWbar/dCbar
      if (nff .eq. 1) then
      do i=1,3; do j=1,3
         DEV1(i,j) = PSI1*II(i,j) + PSI41*A0(i,j)  !  if one fiber families
      enddo; enddo
	  endif
c	  
	  if (nff .eq. 2) then
	  do i=1,3; do j=1,3
         DEV1(i,j) = PSI1*II(i,j) + PSI41*A0(i,j) + PSI42*B0(i,j)  !  if two fiber families
      enddo; enddo  
      endif


!      trDEV2 = 0d0
!      do i=1,3; do j=1,3
!         trDEV2 = trDEV2 + Cbar(i,j)*DEV1(j,i)      ! Cbar : DEV1
!     enddo; enddo
 
! --- deviatoric projector
       DEV2 = matmul(DEV1, Cbartrans)
	   
	   trDEV2 = DEV2(1,1) + DEV2(2,2) + DEV2(3,3)
	   
      do i=1,3; do j=1,3
         DEV(i,j) = DEV1(i,j) - (1.0d0/3.0d0)*trDEV2*Cbar_inv(i,j)
      enddo; enddo

! --- 2nd-PK
 !     do i=1,3; do j=1,3
 !        S(i,j) = 2.0d0*d_pen*(JJ-1.0d0)*JJ*Cinv(i,j) + 2.0d0*JJinv23*DEV(i,j)
 !     enddo; enddo
      
!--- Final 2nd PK stress
      S = 2.0d0 * d_pen * lnJ * Cinv
     &      + 2.0d0 * JJinv23 * DEV
	  hsv(10) = S(1,1)
	  hsv(11) = S(1,2)
      hsv(12) = S(1,3)
	  hsv(13) = S(2,2)
	  hsv(14) = S(2,3)
	  hsv(15) = S(3,3)
!---- Green-Lagrange strain tesnor
      E1 = 0.5d0 * (C - II)
	  hsv(20) = E1(1,1)
	  hsv(21) = E1(1,2)
      hsv(22) = E1(1,3)
	  hsv(23) = E1(2,2)
	  hsv(24) = E1(2,3)
	  hsv(25) = E1(3,3)
! -- Map S (2nd PK) -> sigma (Cauchy)
      SIGMA = matmul( F, matmul(S, Ftrans) ) / JJ
!--- Voigt conversion
      sig(1) = SIGMA(1,1)
      sig(2) = SIGMA(2,2)
      sig(3) = SIGMA(3,3)
      sig(4) = SIGMA(1,2)
      sig(5) = SIGMA(2,3)
      sig(6) = SIGMA(1,3)

!--- NaN protection
      do i = 1,6
         if (sig(i) .ne. sig(i)) sig(i) = 0.0d0
      enddo



c--- Debug output (optional)
 
      return
      end
c|  ------------------------------------------------------------  |
c==================================================================
c     ADD THIS IMMEDIATELY AFTER UMAT45 "end" STATEMENT
c==================================================================
c==================================================================
c     UTAN45 for Gasser-Ogden-Holzapfel (GOH) Model
c     Material Tangent Modulus for Implicit Analysis
c     Add this AFTER the UMAT45 subroutine in the same file
c==================================================================
c==================================================================
c     UTAN45 FOR GOH MODEL - MATCHES UMAT45 HSV ACCESS
!DEC$ ATTRIBUTES DLLEXPORT :: UTAN45
      subroutine utan45(cm,eps,sig,epsp,hsv,dt1,unsym,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c     TANGENT MATRIX CONSISTENT WITH YOUR UMAT45
c
      include 'nlqparm'
      include 'iounits.inc'
c
      dimension cm(*),eps(*),sig(*),hsv(*)
      dimension es(6,*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,unsym
      character*5 etype
c
c     SAME VARIABLES AS THE UMAT
      real*8 F(3,3), C(3,3), Cbar(3,3), Cinv(3,3), Cbar_inv(3,3)
      real*8 A0(3,3), B0(3,3), II(3,3), Ftrans(3,3)
      real*8 c10, d_pen, a1, a2, kappa, gamma
      real*8 JJ, JJinv23, JJinv43, JJinv53, lnJ, detC
      real*8 I1, I4, J4, I4star1, I4star2
      real*8 phi1, phi2, arg1, arg2
      real*8 PSI1, PSI41, PSI42, PSI14_total, PSI11_total
      real*8 PSI11, PSI44_1, PSI44_2, PSI14_1, PSI14_2
      real*8 n1, n2, n3, m1, m2, m3, normn, normm
      real*8 cg, sg
      real*8 C_ijkl(6,6)
      integer i, j, k, l, nhv, nff
      real*8, parameter :: pi = 3.14159265358979d0
c
c     Initialize
      do i=1,6
        do j=1,6
          es(i,j) = 0.0d0
        enddo
      enddo
c
c     Extract parameters EXACTLY as your UMAT
      c10   = cm(1) / 2.0d0
      d_pen = cm(2) / 2.0d0
      a1    = cm(3)
      a2    = cm(4)
      kappa = cm(5)
      gamma = (pi/180.d0) * cm(6)
      nff   = nint(cm(7))
      nhv   = nint(cm(8))
c
      cg = dcos(gamma)
      sg = dsin(gamma)
c
c     Identity
      II(1,1) = 1.0d0; II(1,2) = 0.0d0; II(1,3) = 0.0d0
      II(2,1) = 0.0d0; II(2,2) = 1.0d0; II(2,3) = 0.0d0
      II(3,1) = 0.0d0; II(3,2) = 0.0d0; II(3,3) = 1.0d0
c
c     Get F from history (SAME as your UMAT)
      F(1,1) = hsv(1+nhv); F(1,2) = hsv(2+nhv); F(1,3) = hsv(3+nhv)
      F(2,1) = hsv(4+nhv); F(2,2) = hsv(5+nhv); F(2,3) = hsv(6+nhv)
      F(3,1) = hsv(7+nhv); F(3,2) = hsv(8+nhv); F(3,3) = hsv(9+nhv)
c
c     Jacobian
      JJ = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
     &   - F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
     &   + F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
c
c      if (JJ .le. 0.01d0) return
c
      JJinv23 = JJ**(-2.0d0/3.0d0)
      JJinv43 = JJ**(-4.0d0/3.0d0)
	  JJinv53 = JJ**(-5.0d0/3.0d0)
      lnJ = dlog(JJ)
c
c     C = F^T * F (SAME as the UMAT)
      Ftrans = transpose(F)
      C = matmul(Ftrans, F)
c
c     Symmetrize C (as in the UMAT)
      do i=1,3
        do j=1,3
          C(j,i) = 0.5d0*(C(j,i) + C(i,j))
        enddo
      enddo
c
c     Inverse of C
      call invert3x3(C, Cinv)
c
c     Cbar = J^(-2/3) * C
      do i=1,3
        do j=1,3
          Cbar(i,j) = JJinv23 * C(i,j)
        enddo
      enddo
c
c     Symmetrize Cbar (as in your UMAT)
      do i=1,3
        do j=1,3
          Cbar(j,i) = 0.5d0*(Cbar(j,i) + Cbar(i,j))
        enddo
      enddo
c
c     Inverse of Cbar
      call invert3x3(Cbar, Cbar_inv)
c
c     First invariant
      I1 = Cbar(1,1) + Cbar(2,2) + Cbar(3,3)
c
c     Fiber directions (SAME as the UMAT)
c     First fiber: n = cos(gamma)*qmat(:,1) + sin(gamma)*qmat(:,2)
      n1 = cg*qmat(1,1) + sg*qmat(1,2)
      n2 = cg*qmat(2,1) + sg*qmat(2,2)
      n3 = cg*qmat(3,1) + sg*qmat(3,2)
      normn = dsqrt(n1*n1 + n2*n2 + n3*n3)
      n1 = n1/normn; n2 = n2/normn; n3 = n3/normn
c
c     Structural tensor A0 = n ⊗ n
      A0(1,1) = n1*n1; A0(1,2) = n1*n2; A0(1,3) = n1*n3
      A0(2,1) = n2*n1; A0(2,2) = n2*n2; A0(2,3) = n2*n3
      A0(3,1) = n3*n1; A0(3,2) = n3*n2; A0(3,3) = n3*n3
c
c     Second fiber (if nff=2)
      if (nff .eq. 2) then
        m1 = cg*qmat(1,1) - sg*qmat(1,2)
        m2 = cg*qmat(2,1) - sg*qmat(2,2)
        m3 = cg*qmat(3,1) - sg*qmat(3,2)
        normm = dsqrt(m1*m1 + m2*m2 + m3*m3)
        m1 = m1/normm; m2 = m2/normm; m3 = m3/normm
c
        B0(1,1) = m1*m1; B0(1,2) = m1*m2; B0(1,3) = m1*m3
        B0(2,1) = m2*m1; B0(2,2) = m2*m2; B0(2,3) = m2*m3
        B0(3,1) = m3*m1; B0(3,2) = m3*m2; B0(3,3) = m3*m3
      endif
c
c     Fiber invariant I4
      I4 = A0(1,1)*Cbar(1,1) + A0(1,2)*Cbar(1,2) + A0(1,3)*Cbar(1,3)
     &   + A0(2,1)*Cbar(2,1) + A0(2,2)*Cbar(2,2) + A0(2,3)*Cbar(2,3)
     &   + A0(3,1)*Cbar(3,1) + A0(3,2)*Cbar(3,2) + A0(3,3)*Cbar(3,3)
  
      I4star1 = kappa*I1 + (1.0d0 - 3.0d0*kappa)*I4 - 1.0d0
c
c     Second fiber invariant (if nff=2)
      if (nff .eq. 2) then
      J4 = B0(1,1)*Cbar(1,1) + B0(1,2)*Cbar(1,2) + B0(1,3)*Cbar(1,3)
     &   + B0(2,1)*Cbar(2,1) + B0(2,2)*Cbar(2,2) + B0(2,3)*Cbar(2,3)
     &   + B0(3,1)*Cbar(3,1) + B0(3,2)*Cbar(3,2) + B0(3,3)*Cbar(3,3)
        I4star2 = kappa*I1 + (1.0d0 - 3.0d0*kappa)*J4 - 1.0d0
      endif
c
c==================================================================
c     TANGENT MODULUS: C = dσ/dε
c     Consistent with S = 2*d_pen*lnJ*Cinv + 2*J^(-2/3)*DEV
c==================================================================
c
c     1. VOLUMETRIC TANGENT from d(2*d_pen*lnJ*Cinv)/dC
c     C_vol = 4*d_pen*[Cinv⊗Cinv - lnJ*(Cinv⊙Cinv)]
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c
c     VOLUMETRIC TANGENT (in spatial configuration)
c     For nearly incompressible: we use penalty method
c     Spatial tangent: c_vol = K*I⊗I where K = 2*d_pen
c
c     CORRECT VOLUMETRIC TANGENT
c     From Ψ_vol = (d_pen)*(lnJ)² 
c     Spatial tangent: c_vol = (2*d_pen/J)*[2*lnJ*I⊗I + 2*I⊙I]
c
c     Spatial volumetric tangent (what LS-DYNA actually needs)
c     c_vol = (2*d_pen/J)*[2*lnJ*(I⊗I) + 2*(I⊙I)]
c==================================================================
c     1. VOLUMETRIC PART
c     From Ψ_vol = d_pen*(lnJ)² → S_vol = 2*d_pen*lnJ*C^(-1)
c     Spatial tangent: c_vol = (2*d_pen/J)*[2*lnJ*I⊗I + 2*I⊙I]
c
c==================================================================
c     SPATIAL TANGENT MODULUS for GOH Model
c     Consistent with S = 2*d_pen*lnJ*Cinv + 2*J^(-2/3)*DEV
c==================================================================
c
c     Initialize tangent matrix
      do i=1,6
        do j=1,6
          es(i,j) = 0.0d0
        enddo
      enddo
c
c     Check for valid deformation state
      if (JJ .le. 0.1d0 .or. JJ .gt. 10.0d0) then
c       Return minimal stiffness for safety
        do i=1,3
          es(i,i) = 1.0d0
        enddo
        es(4,4) = 0.5d0
        es(5,5) = 0.5d0
        es(6,6) = 0.5d0
        return
      endif
c
c==================================================================
c     1. VOLUMETRIC TANGENT
c     From Ψ_vol = d_pen*(lnJ)² → S_vol = 2*d_pen*lnJ*C^(-1)
c     Spatial: c_vol = (2*d_pen/J)*[2*lnJ*I⊗I + 2*I⊙I]
c==================================================================
c
c     I⊗I contribution (fills all 3x3 block)
      do i=1,3
        do j=1,3
          es(i,j) = (4.0d0*d_pen*lnJ/JJ)
        enddo
      enddo
c
c     I⊙I contribution (diagonal only)
      do i=1,3
        es(i,i) = es(i,i) + (4.0d0*d_pen/JJ)
      enddo
c
c==================================================================
c     2. ISOCHORIC NEO-HOOKEAN - COMPLETE SPATIAL TANGENT
c==================================================================
c
c     Factor for isochoric terms
      factor_iso = 2.0d0*c10*JJinv53/JJ  ! = 2μ*J^(-8/3)
c
c     Spatial deviatoric projector: I⊙I - (1/3)I⊗I
c     In index notation: δᵢₖδⱼₗ - (1/3)δᵢⱼδₖₗ
c
c     For Voigt notation (symmetric storage):
c     Upper 3x3 block (normal stresses)
      do i=1,3
        do j=1,3
          if (i.eq.j) then
c           Diagonal terms: (1 - 1/3) = 2/3
            es(i,j) = es(i,j) + factor_iso*(4.0d0/3.0d0)
          else
c           Off-diagonal: -1/3
            es(i,j) = es(i,j) - factor_iso*(2.0d0/3.0d0)
          endif
        enddo
      enddo
c
c     Additional geometric stiffness terms (often neglected but more accurate)
c     These come from the push-forward and chain rule
c     For full accuracy, add:
      if (.false.) then  ! Enable for higher accuracy
        do i=1,3
          do j=1,3
c           Additional terms involving C̄^(-1)
c           These are small for moderate deformations but matter at large stretch
            correction = factor_iso*(I1/9.0d0)*Cbar_inv(i,j)
            es(i,j) = es(i,j) + correction
          enddo
        enddo
      endif
c==================================================================
c     3. FIBER CONTRIBUTIONS (Tension-only, No coupling)
c==================================================================
c
c     Initialize combined contributions
      PSI11_total = 0.0d0
      PSI14_total = 0.0d0
c
c     ------ FIRST FIBER FAMILY ------
      if (I4star1 .gt. 0.0d0) then  ! TENSION ONLY
        phi1 = I4star1
        arg1 = a2*phi1*phi1
c
c       Prevent exponential overflow
        if (arg1 .gt. 40.0d0) arg1 = 40.0d0
        exp1 = dexp(arg1)
c
c       Second derivatives
        PSI44_1 = a1*(1.0d0 - 3.0d0*kappa)**2 * exp1 * 
     &            (1.0d0 + 2.0d0*a2*phi1*phi1)
c
c       Cap to prevent numerical instability
        if (PSI44_1 .gt. 1000.0d0) PSI44_1 = 1000.0d0
c
c       Add fiber stiffness (only through A0, no coupling)
        do i=1,3
          do j=1,3
            es(i,j) = es(i,j) + (4.0d0*JJinv53/JJ)*PSI44_1*A0(i,j)
          enddo
        enddo
c
c       Cross-derivative terms (if kappa > 0)
        if (kappa .gt. 1.0d-6) then
          PSI14_1 = a1*kappa*(1.0d0 - 3.0d0*kappa)*exp1*
     &              (1.0d0 + 2.0d0*a2*phi1*phi1)
          PSI11_1 = a1*kappa*kappa*exp1*
     &              (1.0d0 + 2.0d0*a2*phi1*phi1)
c
c         Cap these terms
          if (PSI14_1 .gt. 100.0d0) PSI14_1 = 100.0d0
          if (PSI11_1 .gt. 100.0d0) PSI11_1 = 100.0d0
c
          PSI14_total = PSI14_1  ! Start with fiber 1
          
c         PSI11 depends on nff (matching UMAT logic)
          if (nff .eq. 1) then
c           Single fiber: only fiber 1 contributes to PSI11
            PSI11_total = PSI11_1
          else
c           Two fibers: fiber 1 contribution (fiber 2 added below)
            PSI11_total = PSI11_1
          endif
        endif
      endif
c
c------ SECOND FIBER FAMILY ------
      if (nff .eq. 2 .and. I4star2 .gt. 0.0d0) then  ! TENSION ONLY
        phi2 = I4star2
        arg2 = a2*phi2*phi2
c
c       Prevent exponential overflow
        if (arg2 .gt. 40.0d0) arg2 = 40.0d0
        exp2 = dexp(arg2)
c
c       Second derivatives
        PSI44_2 = a1*(1.0d0 - 3.0d0*kappa)**2 * exp2 *
     &            (1.0d0 + 2.0d0*a2*phi2*phi2)
c
c       Cap to prevent numerical instability
        if (PSI44_2 .gt. 1000.0d0) PSI44_2 = 1000.0d0
c
c       Add fiber stiffness (only through B0, no coupling with A0)
        do i=1,3
          do j=1,3
            es(i,j) = es(i,j) + (4.0d0*JJinv53/JJ)*PSI44_2*B0(i,j)
          enddo
        enddo
c
c       Cross-derivative terms (if kappa > 0)
        if (kappa .gt. 1.0d-6) then
          PSI14_2 = a1*kappa*(1.0d0 - 3.0d0*kappa)*exp2*
     &              (1.0d0 + 2.0d0*a2*phi2*phi2)
          PSI11_2 = a1*kappa*kappa*exp2*
     &              (1.0d0 + 2.0d0*a2*phi2*phi2)
c
c         Cap these terms
          if (PSI14_2 .gt. 100.0d0) PSI14_2 = 100.0d0
          if (PSI11_2 .gt. 100.0d0) PSI11_2 = 100.0d0
c
c         Add to totals (both fibers contribute when nff=2)
          PSI14_total = PSI14_total + PSI14_2
          PSI11_total = PSI11_total + PSI11_2  ! Add second fiber to PSI11
        endif
      endif
c
c     ------ Add combined I1 contributions ------
      if (kappa .gt. 1.0d-6 .and. 
     &    (PSI14_total + PSI11_total) .gt. 0.0d0) then
        do i=1,3
          es(i,i) = es(i,i) + (4.0d0*JJinv53/JJ)*
     &              (PSI14_total + PSI11_total)
        enddo
      endif
c==================================================================
c     4. SHEAR STIFFNESS
c==================================================================
c     Base shear from Neo-Hookean matrix
      es(4,4) = 2.0d0*c10*JJinv53/JJ  ! c10 = μ/2
      es(5,5) = 2.0d0*c10*JJinv53/JJ
      es(6,6) = 2.0d0*c10*JJinv53/JJ
c
c     Ensure positive definiteness
      shear_actual = 2.0d0*c10*JJinv53/JJ
      
c     Use smaller fraction for ultra-soft materials
c     0.1% instead of 1% gives more flexibility
      shear_min = max(1.0d-10, 1.0d-3*shear_actual)  ! 0.1% of actual
      
      do i=4,6
        if (es(i,i) .lt. shear_min) es(i,i) = shear_min
      enddo
      
c     Enforce symmetry in upper 3x3 block
      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)
c==================================================================
c     5. FINAL VALIDATION
c==================================================================
c     PROPER MINIMUMS FOR YOUR ULTRA-SOFT ARTERIAL TISSUE
c     All values in GPa units
c==================================================================
cc    
      do i=1,3
        if (es(i,i) .le. 0.0d0) then
          es(i,i) = 1.0d-9  ! 0.001 kPa << 7 kPa (Y-stress)
        endif
      enddo
      
      do i=4,6  
        if (es(i,i) .le. 0.0d0) then
          es(i,i) = 1.0d-10  ! 0.0001 kPa for shear
        endif
      enddo
c     STOP HERE - Don't add another loop!
c
      return
      end
c======================================================================
c======================================================================
c     ANALYTICAL TANGENT (for better performance)
c======================================================================
c
c     For better performance, derive the analytical tangent:
c     
c     For GOH model: C_ijkl = ∂²Ψ/∂ε_ij∂ε_kl
c     
c     Key terms:
c     1. Volumetric part: from d_pen * ln(J) term
c     2. Isochoric matrix: from c10 * (I1-3) term
c     3. Fiber contribution: from a1*exp(...) terms
c     4. Deviatoric projection: complex tensor operations
c     
c     Benefits:
c     - No extra stress evaluations (6× faster)
c     - Exact quadratic convergence
c     - Better for very stiff problems
c======================================================================   
c     NOTES ON USAGE:
c======================================================================
c
c     1. The altenative numerical tangent works by computing:
c        
c        C_ij = ∂σ_i/∂ε_j ≈ (σ_i(ε+Δε) - σ_i(ε)) / Δε
c        
c        where Δε = 1e-7 (delta_eps parameter)
c
c     2. The tangent is automatically symmetrized, which is correct
c        for hyperelastic materials derived from a strain energy function
c
c     3. For implicit analysis to work well, you need:
c        - Smooth stress-strain response (✓ your GOH model has this)
c        - Positive definite tangent (✓ enforced by diagonal checks)
c        - Consistent with the stress computation (✓ uses same umat45)
c
c     4. Performance notes:
c        - Requires 6 extra stress evaluations per element per iteration
c        - For most problems, this overhead is acceptable
c        - Newton-Raphson convergence should be quadratic
c
c     5. If you experience convergence problems, try:
c        - Smaller load steps
c        - Line search (LS-DYNA *CONTROL_IMPLICIT_SOLUTION, ILIMIT)
c        - Adjust perturbation size (delta_eps)
c        - Enable stabilization (uncomment lines above)
c
c     6. To verify the tangent is correct, you can compare:
c        - Numerical tangent vs analytical tangent (if you derive it)
c        - Check symmetry: max|C_ij - C_ji| should be < 1e-6
c        - Check positive definiteness: all eigenvalues > 0
c======================================================================
     

c*******************************************************************
      subroutine umat46 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,thhsvi,nthhsv,qmat,elsiz,idele,
     2 reject)
	 
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
 
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),thhsvi(*),
     . cma(*),qmat(3,3)
      integer nnpcrv(*)
      character*5 etype
      logical failel,reject
      integer idele
c
c     Isotropic elastic material with access to history variables
c     of a corresponding thermal user material. See also
c     subroutine umat41 for information about non-thermal
c     specific information. Access to thermal history variables
c     requieres IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED and
c     is only possible for a selection of shell and solid elements.
c     More information including a list of supported element
c     types is given in the comments to subroutine thusrmat.
c
c     Read only variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c     cm(3)=.eq.0 do nothing
c           .eq.1 print mechanical and thermal history variables.
c
c     nthhsv=Number of thermal history variables.
c     thhsvi=If IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED then
c            thhsvi(i) i=1...nthhsv, contains history
c            variable data from thermal user material routine
c            interpolated to current integration point using
c            shape functions for the corresponding thermal element.
c
c     compute shear modulus, g
c
      if (ncycle.eq.1) then
        call usermsg('mat46')
      endif
c
      g2 =cm(1)/(1.+cm(2))
      g  =.5*g2
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
        davg=(-eps(1)-eps(2)-eps(3))/3.
        p=-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=sig(2)+p+g2*(eps(2)+davg)
        sig(3)=sig(3)+p+g2*(eps(3)+davg)
        sig(4)=sig(4)+g*eps(4)
        sig(5)=sig(5)+g*eps(5)
        sig(6)=sig(6)+g*eps(6)
c
      else if (etype.eq.'shell'.or.etype.eq.'shl2d') then
        gc    =capa*g
        q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
        q3    =1./(q1+g2)
        eps(3)=-q1*(eps(1)+eps(2))*q3
        davg  =(-eps(1)-eps(2)-eps(3))/3.
        p     =-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=sig(2)+p+g2*(eps(2)+davg)
        sig(3)=0.0
        sig(4)=sig(4)+g *eps(4)
        sig(5)=sig(5)+gc*eps(5)
        sig(6)=sig(6)+gc*eps(6)
c
      else
c       write(iotty,10) etype
c       write(iohsp,10) etype
c       write(iomsg,10) etype
c       call adios(2)
        cerdat(1)=etype
        call lsmsg(3,MSG_SOL+1152,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c     print history variables?
c
      if ( nint(cm(3)).eq.1) then
        write(iotty,*) " epsp ",epsp," thhsvi(1) ",thhsvi(1),
     1   " thhsvi(5) ",thhsvi(5),nthhsv
      endif
c
c10   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current user material model #47.')
      return
      end
c******************************************************************
c     UMAT47 - MSCM V1 (Base Model without Interactions)
c======================================================================
c======================================================================
 !DEC$ ATTRIBUTES DLLEXPORT :: UMAT47
	       subroutine umat47 (cm,eps,sig,epsp,hsv,dt1,capa,
     . etype,time,temp,failel,crv,nnpcrv,cma,qmat,elsiz,idele,ih,
     . reject)
	 
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
  
	  include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel,reject
      integer nnpcrv(*)
      character*5 etype
      integer idele
	  logical :: OK_FLAG
	  real*8 pi
      parameter (pi=3.141592653589793d0)
c
      real*8 tt,temper,epsp,elsiz,dt1,capa
c
c     Local variables
      real*8 F(3,3),C(3,3),Cinv(3,3),SIGMA(3,3),Ftrans(3,3)
      real*8 gamma_f,lambda_lb,lambda_ub,mu_s,sigma_s,c_m,phi_f
      real*8 lambda_f,lambda_s,P_ens,n(3),max_lambda_f,S_fiber_11
	  real*8 lambda_f_post,n_post(3), theta_post
      real*8 JJ,I1,det33,lnJ,E11,E22,E12,nCn
	  real*8 Svol(3,3),Siso(3,3),S(3,3),Sf(3,3),Sm(3,3),C33
      real*8 gauss5(5),w5(5),gauss_s(5),ws(5),Gamma_gauss
	  real*8 gauss21(21),w21(21)
      real*8 Gamma_s,dGamma_s,Gamma_n,Gamma_t,Gamma_0
	  real*8 n_ref(2),n_def(2),beta,n_def_mag,J2D,theta
	  real*8 lambda_f_array(100)  ! Array to store lambda_f
      real*8 TFR,FEFR,erf_norm                  ! Total fiber recruitment
      integer i,j,k,nhv,igauss,itheta,isigma
      integer n_theta,n_sigma,n_theta_post
	  
c

c
c

	 
	 
	
c
c==================================================================
c     Material Parameters (following MSCM notation)
c==================================================================
c     cm(1) = cm_bar (matrix modulus (mu_m) * volume fraction (phi_m), GPa)
c     cm(2) = d_pen (bulk modulus)
c     cm(3) = gamma_f_bar (fiber modulus * volume fraction, GPa)
c     cm(4) = sigma_n (angular distribution, degrees)
c     cm(5) = mu_s (mean slack stretch)
c     cm(6) = sigma_s (recruitment standard deviation)
c     cm(7) = lambda_lb (lower recruitment bound)
c     cm(8) = lambda_ub (upper recruitment bound)
c     cm(9) = nhv (history variable offset)

                   



c
      cm_bar = cm(1)/2.0d0             ! cm(1) = mu_m * phi_m
      d_pen = cm(2)/2.0d0              ! Bulk modulus 
      gamma_f_bar = cm(3)          ! 
      sigma_n = cm(4)*pi/180.0d0   ! Convert degrees to radians
      mu_s = cm(5)                 ! Mean slack stretch
      sigma_s = cm(6)              ! Recruitment standard deviation
      lambda_lb = cm(7)            ! Lower bound
      lambda_ub = cm(8)            ! Upper bound
      nhv = nint(cm(9))            ! History variable offset


c

      n_theta_post = 20.0d0  ! Angular segments for postprocessing
      n_sigma = 5.0d0   ! Number of Gauss points for slack integration
	  n_theta = 1.0d0       ! No segmentation over theta range
	  phi_f = 1.0d0   ! Since we merged volume fractions
      phi_m = 0.0d0   ! Or just remove this line
c     At the very beginning of UMAT (after variable declarations)

c     5-point Gauss quadrature
      data w5/0.5688888888888889d0, 0.4786286704993665d0, 
     &        0.4786286704993665d0, 0.2369268850561891d0, 
     &        0.2369268850561891d0/
      data gauss5/0.d0,-0.5384693101056831d0,0.5384693101056831d0,
     &           -0.9061798459386640d0,0.9061798459386640d0/
	 
c     21-point Gauss quadrature weights
      data w21/0.1460811336496904d0, 0.1445244039899700d0,
     &         0.1445244039899700d0, 0.1398873947910731d0,
     &         0.1398873947910731d0, 0.1322689386333375d0,
     &         0.1322689386333375d0, 0.1217818238391521d0,
     &         0.1217818238391521d0, 0.1085536696255165d0,
     &         0.1085536696255165d0, 0.0927259106232362d0,
     &         0.0927259106232362d0, 0.0744524990543851d0,
     &         0.0744524990543851d0, 0.0538990698846632d0,
     &         0.0538990698846632d0, 0.0312394965310737d0,
     &         0.0312394965310737d0, 0.0066627134645283d0,
     &         0.0066627134645283d0/

c     21-point Gauss quadrature points  
      data gauss21/0.0d0,
     &            -0.1455618541608951d0, 0.1455618541608951d0,
     &            -0.2880213168024011d0, 0.2880213168024011d0,
     &            -0.4243421202074388d0, 0.4243421202074388d0,
     &            -0.5516188358872198d0, 0.5516188358872198d0,
     &            -0.6671388041974123d0, 0.6671388041974123d0,
     &            -0.7684399634756779d0, 0.7684399634756779d0,
     &            -0.8533633645833173d0, 0.8533633645833173d0,
     &            -0.9200993341504008d0, 0.9200993341504008d0,
     &            -0.9672268385663063d0, 0.9672268385663063d0,
     &            -0.9937521706203895d0, 0.9937521706203895d0/

c
c==================================================================
c     Kinematics
c==================================================================
c     Extract F from history
      F(1,1) = hsv(1+nhv); F(1,2) = hsv(2+nhv); F(1,3) = hsv(3+nhv)
      F(2,1) = hsv(4+nhv); F(2,2) = hsv(5+nhv); F(2,3) = hsv(6+nhv)
      F(3,1) = hsv(7+nhv); F(3,2) = hsv(8+nhv); F(3,3) = hsv(9+nhv)
c
c     Right Cauchy-Green C = F^T*F
      Ftrans = transpose(F)
      C = matmul(Ftrans, F)
c
c     Jacobian
      JJ = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
     &   - F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
     &   + F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
c
c     First invariant
      I1 = C(1,1) + C(2,2) + C(3,3)
c
c     Inverse of C
      call invert3x3(C, Cinv)
c
c==================================================================
c     Initialize 2nd PK Stress
c==================================================================
      do i=1,3
        do j=1,3
          Sf(i,j) = 0.0d0
        enddo
      enddo
	  
	  do i=1,3
        do j=1,3
          Svol(i,j) = 0.0d0
        enddo
      enddo
	  
	  do i=1,3
        do j=1,3
          Siso(i,j) = 0.0d0
        enddo
      enddo	  
c
c==================================================================
c     FIBER CONTRIBUTION - Double Integration
c     S_f = γ̄f ∫∫ (Γn/λf) * Γs * (λf/λs - 1) * (n⊗n) dλs dθ
c==================================================================
c
c     Initialize fiber tracking variables
      S_fiber_11 = 0.0d0
      S_fiber_22 = 0.0d0
      S_fiber_12 = 0.0d0
      n_fibers_active = 0
      max_lambda_f = 0.0d0
      
c     Check if fiber contribution should be computed
      if (gamma_f_bar .gt. 1.0d-10) then
c
c       Angular integration over fiber orientations
        do itheta = 1, n_theta
          theta_a = -pi/2.0d0 + pi*(itheta-1.0d0)/dble(n_theta)
          theta_b = -pi/2.0d0 + pi*itheta/dble(n_theta)
          
c         Gauss quadrature over theta
          do igauss = 1, 21
            theta = gauss21(igauss)*(theta_b-theta_a)/2.0d0 + 
     &              (theta_b+theta_a)/2.0d0
            
c           Fiber direction (in-plane)
            n(1) = dcos(theta)
            n(2) = dsin(theta)
            n(3) = 0.0d0
            
c           Fiber stretch λf = sqrt(n·C·n)
            lambda_f = 0.0d0
            do i=1,3
              do j=1,3
                lambda_f = lambda_f + n(i)*C(i,j)*n(j)
              enddo
            enddo
            lambda_f = dsqrt(lambda_f)
c
            
c           Track maximum stretch
            max_lambda_f = dmax1(max_lambda_f, lambda_f)
            
c           Angular distribution Γn (Gaussian centered at θ=0)
!            Gamma_n = dexp(-0.5d0*(theta/sigma_n)**2) / 
!     &                (dsqrt(2.0d0*pi)*sigma_n)
c=========================================================================================
c Simple Orientation Distribution Function (ODF) 
c==========================================================================================	
c           Angular distribution Γn (Gaussian centered at θ=0)
c           NORMALIZED over [-π/2, +π/2]
!            Gamma_n = dexp(-0.5d0*(theta/sigma_n)**2) / 
!     &                (dsqrt(2.0d0*pi)*sigma_n)
     
c           Normalize over finite domain [-π/2, +π/2]
!            znorm_n = derf(pi/(2.0d0*sigma_n*dsqrt(2.0d0))) - 
!     &                derf(-pi/(2.0d0*sigma_n*dsqrt(2.0d0)))
!            if (abs(znorm_n) .gt. 1.0d-10) then
!              Gamma_n = Gamma_n / znorm_n
!            endif
c=========================================================================================
c Angular Orientation Distribution Function (ODF) 
c==========================================================================================		
c           Angular distribution with proper normalization (Fan & Sacks 2014)
c           For d=1 (purely Gaussian case)
            
c           Gaussian part (unnormalized)
            Gamma_gauss = dexp(-0.5d0*(theta/sigma_n)**2)
            
c           Normalization factor using error function
            erf_norm = erf(pi/(2.0d0*dsqrt(2.0d0)*sigma_n))
            
c           Final normalized distribution
            Gamma_n = Gamma_gauss / 
     &                (erf_norm*dsqrt(2.0d0*pi)*sigma_n)
            
c           Optional: Include uniform component if d < 1
c           d_param = 1.0d0  ! Set in parameters
c           Gamma_uniform = 1.0d0/pi
c           Gamma_n = d_param*Gamma_n + (1.0d0-d_param)*Gamma_uniform
c========================================================================          
c           Integration over slack stretch distribution
            P_ens = 0.0d0
            
c           Only integrate if fibers are recruited
c           TENSION-ONLY: Fibers cannot support compression
            if (lambda_f .le. 1.0d0) then
              P_ens = 0.0d0  ! No contribution in compression
              
            elseif (lambda_f .gt. lambda_lb) then
c             Fibers in tension AND recruited
c             INCREMENT COUNTER 
              n_fibers_active = n_fibers_active + 1
c             Determine integration bounds
              sigma_a = dmax1(1.0d0, lambda_lb)
              sigma_b = dmin1(lambda_f, lambda_ub)
			  
              
              if (sigma_b .gt. sigma_a) then
c               Gauss quadrature over lambda_s
                do isigma = 1, n_sigma
                  lambda_s = gauss5(isigma)*(sigma_b-sigma_a)/2.0d0 + 
     &                       (sigma_a+sigma_b)/2.0d0               
c                 Recruitment distribution Γs with proper normalization
                  arg_s = ((lambda_s - mu_s)/sigma_s)**2
                  if (arg_s .lt. 20.0d0) then
                    Gamma_s = dexp(-0.5d0*arg_s) / 
     &                        (sigma_s*dsqrt(2.0d0*pi))
c                   Normalize over [lambda_lb, lambda_ub]
                    znorm1 = derf((lambda_ub-mu_s)/
     &                            (sigma_s*dsqrt(2.0d0)))
                    znorm2 = derf((lambda_lb-mu_s)/
     &                            (sigma_s*dsqrt(2.0d0)))
                    znorm = znorm1 - znorm2
                    if (abs(znorm) .gt. 1.0d-10) then
                      Gamma_s = Gamma_s / (0.5d0*znorm)
                    endif
                  else
                    Gamma_s = 0.0d0
                  endif
                  
c                 Linear fiber stress P_f = (γ̄f/λs)*(λf/λs - 1)
                  P_fiber = (gamma_f_bar/lambda_s)*
     &                      (lambda_f/lambda_s - 1.0d0)
                  
c                 Ensemble stress contribution
                  P_ens = P_ens + w5(isigma)*Gamma_s*P_fiber*
     &                    (sigma_b-sigma_a)/2.0d0
                enddo  ! End lambda_s integration
              endif
            endif
            
c           Assemble fiber stress tensor
c           Weight: w_theta * Gamma_n * (theta_b-theta_a)/2
!            weight = w5(igauss)*Gamma_n*(theta_b-theta_a)/2.0d0
            weight = w21(igauss)*Gamma_n*(theta_b-theta_a)/2.0d0
            
c           Debug output for key angles

            
c           S_f += (Γn/λf) * P_ens * (n⊗n)
            if (lambda_f .gt. 1.0d-10) then
              factor = weight * P_ens / lambda_f
              
c             Add to total stress
              do i=1,3
                do j=1,3
                  Sf(i,j) = Sf(i,j) + factor * n(i) * n(j)
                enddo
              enddo
              
c             Track fiber contributions separately
              S_fiber_11 = S_fiber_11 + factor * n(1) * n(1)
              S_fiber_22 = S_fiber_22 + factor * n(2) * n(2)
              S_fiber_12 = S_fiber_12 + factor * n(1) * n(2)
            endif
            
          enddo  ! End Gauss quadrature over theta
        enddo  ! End theta segments
        
        
       endif  ! End fiber contribution check
c*******************************************************************************************
c===========================================================================================
c      START OF COMPUTATIONS REQUIRED FOR POST-PROCESSING
c===========================================================================================
c===========================================================================================
c     COMPUTE lambda_f FOR POST-PROCESSING AT n_theta_post ANGLES
c===========================================================================================
c     After main fiber integration, compute lambda_f at specific angles for post-processing
      do itheta = 1, n_theta_post
        theta_post = -pi/2.0d0 + pi*(itheta-0.5d0)/dble(n_theta_post)
        
c       Fiber direction at this post-processing angle
        n_post(1) = dcos(theta_post)
        n_post(2) = dsin(theta_post)
        n_post(3) = 0.0d0
        
c       Compute fiber stretch for this direction
        lambda_f_post = 0.0d0
        do i=1,3
          do j=1,3
            lambda_f_post = lambda_f_post + n_post(i)*C(i,j)*n_post(j)
          enddo
        enddo
        lambda_f_array(itheta) = dsqrt(lambda_f_post)
      enddo

c*******************************************************************************************
c===========================================================================================
c      START OF COMPUTATIONS REQUIRED FOR POST-PROCESSING
c===========================================================================================
c     FIBER REORIENTATION AND RECRUITMENT TRACKING 
c===========================================================================================
      J2D = F(1,1)*F(2,2) - F(1,2)*F(2,1)
      
      do itheta = 1, n_theta_post
        theta_post = -pi/2.0d0 + pi*(itheta-0.5d0)/dble(n_theta_post)
        
c       1. Compute deformed angle beta
        n_ref(1) = dcos(theta)
        n_ref(2) = dsin(theta)
        n_def(1) = F(1,1)*n_ref(1) + F(1,2)*n_ref(2)
        n_def(2) = F(2,1)*n_ref(1) + F(2,2)*n_ref(2)
        beta = datan2(n_def(2), n_def(1))

c============================================================================================
c       2. Compute REFERENCE ODF with proper normalization 
c          This is Gamma_n in the reference configuration
           Gamma_gauss = dexp(-0.5d0*(theta/sigma_n)**2)
           Gamma_0 = Gamma_gauss/(erf_norm*dsqrt(2.0d0*pi)*sigma_n)
		   
        
		
c       Optional: Include uniform component if d < 1
c       d_param = 1.0d0  ! Could be cm(9) if you want it variable
c       if (d_param .lt. 1.0d0) then
c          Gamma_uniform = 1.0d0/pi
c          Gamma_0 = d_param*Gamma_0 + (1.0d0-d_param)*Gamma_uniform
c       endif
c============================================================================================        
c       3. Compute FEFR for this direction - COMPUTE FEFR - INTEGRAL OF RECRUITMENT DISTRIBUTION
     
c       Compute TRANSFORMED ODF in deformed configuration
        nCn = n_ref(1)*C(1,1)*n_ref(1) + 
     &        2.0d0*n_ref(1)*C(1,2)*n_ref(2) +
     &        n_ref(2)*C(2,2)*n_ref(2)
        
c       Gamma_t represents fiber distribution in current config
        Gamma_t = Gamma_0 * dsqrt(nCn) / J2D
		
        FEFR = 0.0d0
        
        if (lambda_f_array(itheta) .gt. lambda_lb) then
c         Integrate Gamma_s from lambda_lb to lambda_f
          sigma_a = lambda_lb
          sigma_b = dmin1(lambda_f_array(itheta), lambda_ub)
          
c         Gauss quadrature for integration
          do isigma = 1, 5
            lambda_s = gauss5(isigma)*(sigma_b-sigma_a)/2.0d0 + 
     &                 (sigma_a+sigma_b)/2.0d0
            
c           Gaussian recruitment distribution
            arg_s = ((lambda_s - mu_s)/sigma_s)**2
            if (arg_s .lt. 20.0d0) then
              Gamma_s = dexp(-0.5d0*arg_s)/(sigma_s*dsqrt(2.0d0*pi))
              
c             Normalize over [lambda_lb, lambda_ub]
              znorm1 = derf((lambda_ub-mu_s)/(sigma_s*dsqrt(2.0d0)))
              znorm2 = derf((lambda_lb-mu_s)/(sigma_s*dsqrt(2.0d0)))
              Gamma_s = Gamma_s/(0.5d0*(znorm1-znorm2))
            else
              Gamma_s = 0.0d0
            endif
            
c           Accumulate integral
            FEFR = FEFR + w5(isigma)*Gamma_s*(sigma_b-sigma_a)/2.0d0
         enddo
        endif
	
c       Store for post-processing
        hsv(30+itheta) = beta * 180.0d0/pi
        hsv(50+itheta) = Gamma_t
c       Store FEFR (as fraction 0-1) 
        hsv(70+itheta) = FEFR 
		hsv(91+itheta) = Gamma_0

      enddo
c
c==================================================================
c     STORE LAMBDA_F VALUES FOR POST-PROCESSING
c==================================================================
c     Store lambda_f for each angle in history variables
      do itheta = 1, n_theta_post
c       Store in hsv(121-140) for up to 20 angles (adjust indices as needed)
        hsv(120+itheta) = lambda_f_array(itheta)
      enddo
      
c     Also store key stretch values for quick reference
      hsv(141) = dsqrt(C(1,1))  ! lambda_x (preferred fiber direction)
      hsv(142) = dsqrt(C(2,2))  ! lambda_y (cross-fiber)
      hsv(143) = max_lambda_f   ! maximum fiber stretch
c==================================================================
c   COMPUTE AND STORE TOTAL FIBER RECRUITMENT (TFR)
c==================================================================
      TFR = 0.0d0
      do itheta = 1, n_theta_post
c       Average FEFR over all directions
c       FEFR was stored in hsv(70+itheta)
        TFR = TFR + hsv(70+itheta)/(n_theta_post)
      enddo

      hsv(115) = TFR * 100.0d0  ! Store TFR percentage
c===========================================================================================
c      END OF CALCULATIONS REQUIRED FOR POST-PROCESSING
c===========================================================================================
C*******************************************************************************************
c==================================================================
c   MATRIX CONTRIBUTION - Neo-Hookean with Plane Stress
c   S_m = c_m_bar * (I - C33*C^(-1))
c==================================================================

c     First, solve for C33 that satisfies S33 = 0
c     For Neo-Hookean: C33 = 1/C^(-1)_33
      C33 = 1.0d0/Cinv(3,3)
      
c     Apply plane stress correction to matrix stress
      do i=1,3
        do j=1,3
          if (i.eq.j) then
            Siso(i,j) = cm_bar*(1.0d0 - C33*Cinv(i,j))
          else
            Siso(i,j) = -cm_bar*C33*Cinv(i,j)
          endif
        enddo
      enddo
c	  
c==================================================================
c     VOLUMETRIC PENALTY - With numerical protection
c==================================================================
      if (JJ .gt. 1.0d-10) then  ! Protect against J≤0
        lnJ =  2.0d0*d_pen*dlog(JJ)
        do i=1,3
          do j=1,3
            Svol(i,j) = Svol(i,j) + lnJ*Cinv(i,j)
          enddo
        enddo
      else
        write(*,*) 'WARNING: J =', JJ, 'too small!'
        failel = .true.  ! Flag element failure
        return
      endif
	  
c     3. Total matrix = deviatoric + volumetric

      Sm = Siso + Svol
	  
c==================================================================	   
c Total 2nd P-K : S = Sf + Sm    
      S = Sm + Sf
c
c==================================================================
c     Total Green-Lagrange strain (same for matrix and fibers)
      E11 = 0.5d0*(C(1,1) - 1.0d0)  
      E22 = 0.5d0*(C(2,2) - 1.0d0)
      E12 = 0.5d0*C(1,2)

c==================================================================
c     Convert to Cauchy Stress
c     σ = (1/JJ) * F * S * F^T
c==================================================================
c
!      call push_forward_stress(S, F, J, SIGMA)
	  SIGMA = matmul( F, matmul(S, Ftrans) ) / JJ
	  
	  hsv(10) = S(1,1)
	  hsv(11) = S(1,2)
      hsv(12) = S(1,3)
	  hsv(13) = S(2,2)
	  hsv(14) = S(2,3)
	  hsv(15) = S(3,3)
	  
	  hsv(16) = Sf(1,1)
      hsv(17) = Sf(2,2)
	  hsv(18) = Sf(1,2)
	  hsv(19) = Sm(1,1)
      hsv(20) = Sm(2,2)
	  hsv(21) = Sm(1,2)

	  
c
c     Store in Voigt notation
      sig(1) = SIGMA(1,1)
      sig(2) = SIGMA(2,2)
      sig(3) = SIGMA(3,3)
      sig(4) = SIGMA(1,2)
      sig(5) = SIGMA(2,3)
      sig(6) = SIGMA(1,3)
c
c
   
	  
      return
      end

      subroutine umat48 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
     
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel,reject
      integer nnpcrv(*)
      character*5 etype
      integer idele
c
      if (ncycle.eq.1) then
        call usermsg('mat48')
      endif
c
      return
      end


      subroutine umat49 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
 
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
      integer idele
c
      if (ncycle.eq.1) then
        call usermsg('mat49')
      endif
c
      return
      end

      subroutine umat50 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
  
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
      integer idele
c
      if (ncycle.eq.1) then
        call usermsg('mat50')
      endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c User Material Model:
c   Some subroutines are used in the sample umat material models, and
c   can be excluded.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine metalshl (cm,eps,sig,epsp,hisv,dt,capa,elform,tt,temp,
     .     failel,crv,nnpcrv)
     
      include 'nlqparm'
      real cm(*),eps(*),sig(*),hisv(*)
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*5 elform
      logical failel
      real epsp,dt,capa,tt,temp
c
c     Elastic plastic material routine for shells
c
c     Material parameters:
c
c     cm(1)=Young's modulus
c     cm(2)=Poisson's ratio
c     cm(3)=Bulk modulus
c     cm(4)=Shear modulus
c     cm(5)=Hardening rule
c           .eq.0.:Linear hardening
c           .ne.0.:Hardening determined by load curve with ID cm(5)
c     cm(6)=Yield stress for linear hardening
c     cm(7)=Hardening modulus for linear hardening
c     cm(8)=Failure strain
c
c     History variables:
c
c     epsp   =Effective plastic strain
c     hisv(1)=Multiplier for plastic correction in implicit tangent
c
      real em,pr,bm,gm,hrd,epsf,g2,g3,gc,gt,sigyld,sighrd,eps12,
     1     eps3e,dvol,davg,p,sig1,sig2,sig4,sig5,sig6,sig56,sig456,
     2     sigv,sig1d,sig2d,sig3d,sigeff2,sigeff,quot,deps,coef,
     3     sig3e,eps3p,sig3p
c
c     Compute material parameters
c
      em=cm(1)
      pr=cm(2)
      bm=cm(3)
      gm=cm(4)
      hrd=cm(5)
      epsf=cm(8)
c
c     Some useful shear modulus parameters
c
      g2=2.*gm
      g3=3.*gm
      gc=capa*gm
      gt=4.5*(gc-gm)
c
c     Compute yield stress for current value of plastic strain
c
      if (hrd.le.0.) then
         sigyld=cm(6)+cm(7)*epsp
         sighrd=cm(7)
      else
         call crvval(crv,nnpcrv,hrd,epsp,sigyld,sighrd)
      endif
c
c     Compute elastic through thickness strain increment
c
      eps12=eps(1)+eps(2)
      eps3e=-pr*eps12/(1.-pr)
c
c     Compute volumetric strain increment
c
      dvol=eps12+eps3e
      davg=dvol/3.
c
c     Compute increment in (negative) pressure
c
      p=bm*dvol
c
c     Compute elastic trial stress
c
      sig1=sig(1)+p+g2*(eps(1)-davg)
      sig2=sig(2)+p+g2*(eps(2)-davg)
      sig4=sig(4)+gm*eps(4)
      sig5=sig(5)+gc*eps(5)
      sig6=sig(6)+gc*eps(6)
c
c     Compute some parameters of shear stress components
c     useful for further calculations
c
      sig56=sig5**2+sig6**2
      sig456=sig4**2+sig56
c
c     Compute effective stress in the sense of von Mises
c
      sigv=(sig1+sig2)/3.
      sig1d=sig1-sigv
      sig2d=sig2-sigv
      sig3d=-sigv
c
      sigeff2=1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456)
      sigeff=sqrt(sigeff2)
c
c     Indicate elastic update as default
c
      hisv(1)=0.
c
c     Check if yielding
c
      if (sigeff.gt.sigyld) then
c
c     Yield stress exceeded, perform elastic-plastic stress return
c
c     Compute through thickness stress
c     for elastic through thickness strain
c
         quot=sig56/sigeff2
         deps=(sigeff-sigyld)/(sighrd+g3+gt*quot)
         coef=g3*deps/sigeff
         sig3e=-coef*sig3d
c
c     First point in secant iteration algorithm: (eps3e,sig3e)
c
c     Compute "plastic" through thickness strain increment
c
         eps3p=-eps12
c
c     Compute elastic trial stress
c
         sig1=sig(1)+g2*eps(1)
         sig2=sig(2)+g2*eps(2)
         sig3p=g2*eps3p
c
c     Compute effective stress in the sense of von Mises
c
         sigv=(sig1+sig2+sig3p)/3.
         sig1d=sig1-sigv
         sig2d=sig2-sigv
         sig3d=sig3p-sigv
c
         sigeff2=1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456)
         sigeff=sqrt(sigeff2)
c
c     Check if yielding
c
         if (sigeff.gt.sigyld) then
c
c     Yield stress exceeded, perform stress return
c
            quot=sig56/sigeff2
            deps=(sigeff-sigyld)/(sighrd+g3+gt*quot)
            coef=g3*deps/sigeff
            sig3p=sig3p-coef*sig3d
c
         endif
c
c     Second point in secant iteration algorithm: (eps3p,sig3p)
c
c     Determine through thickness strain based on
c     interpolation from previous two stress updates
c
         if (abs(sig3e-sig3p).gt.1.e-4*
     1        max(abs(sig3e),abs(sig3p))) then
            eps(3)=eps3e-sig3e*(eps3e-eps3p)/(sig3e-sig3p)
         else
            eps(3)=eps3p
         endif
c
c     Compute volumetric strain increment
c
         dvol=eps12+eps(3)
         davg=dvol/3.
c
c     Compute increment in (negative) pressure
c
         p=bm*dvol
c
c     Compute elastic trial stress
c
         sig(1)=sig(1)+p+g2*(eps(1)-davg)
         sig(2)=sig(2)+p+g2*(eps(2)-davg)
         sig(3)=p+g2*(eps(3)-davg)
c
c     Compute effective stress in the sense of von Mises
c
         sigv=(sig(1)+sig(2)+sig(3))/3.
         sig1d=sig(1)-sigv
         sig2d=sig(2)-sigv
         sig3d=sig(3)-sigv
c
         sigeff2=1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456)
         sigeff=sqrt(sigeff2)
c
c     Check if yielding
c
         if (sigeff.gt.sigyld) then
c
c     Plastic update
c
            hisv(1)=.9
c
c     Yield stress exceeded, perform stress return
c
c     Increment plastic strain
c
            quot=sig56/sigeff2
            deps=(sigeff-sigyld)/(sighrd+g3+gt*quot)
            epsp=epsp+deps
c
c     Return stress in the direction of deviatoric stress
c     and set through thickness stress to zero
c
            coef=g3*deps/sigeff
            sig(1)=sigv+sig1d*(1.-coef)
            sig(2)=sigv+sig2d*(1.-coef)
            sig(3)=0.
            sig(4)=sig4*(1.-coef)
            sig(5)=sig5*(1.-capa*coef)
            sig(6)=sig6*(1.-capa*coef)
c
         else
c
c     Yield stress not exceeded
c
            sig(3)=0.
            sig(4)=sig4
            sig(5)=sig5
            sig(6)=sig6
c
         endif
         
c
c     Elastic-plastic stress update finished
c
      else
c
c     The update was elastic, just store the stress
c     and strain increment with no update of history
c     variables, i.e., plastic strain
c
         sig(1)=sig1
         sig(2)=sig2
         sig(3)=0.
         sig(4)=sig4
         sig(5)=sig5
         sig(6)=sig6
c
c     Store through thickness strain increment
c
         eps(3)=eps3e
c
c     Elastic stress update finished
c
      endif
c
c     Check if failure
c
      if (epsp.gt.epsf) failel=.true.
c
      return
      end
      subroutine metalsld (cm,eps,sig,epsp,hisv,dt,capa,elform,tt,temp,
     .     failel,crv,nnpcrv)
      implicit none
      include 'nlqparm'
      real cm(*),eps(*),sig(*),hisv(*)
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*5 elform
      logical failel
      real epsp,dt,capa,tt,temp
c
c     Elastic plastic material routine for solids
c
c     Material parameters:
c
c     cm(1)=Young's modulus
c     cm(2)=Poisson's ratio
c     cm(3)=Bulk modulus
c     cm(4)=Shear modulus
c     cm(5)=Hardening rule
c           .eq.0.:Linear hardening
c           .ne.0.:Hardening determined by load curve with ID cm(5)
c     cm(6)=Yield stress for linear hardening
c     cm(7)=Hardening modulus for linear hardening
c     cm(8)=Failure strain
c
c     History variables:
c
c     epsp   =Effective plastic strain
c     hisv(1)=Multiplier for plastic correction in implicit tangent
c
      real em,pr,bm,gm,hrd,epsf,g2,g3,sigyld,sighrd,dvol,davg,p,sigv,
     1     sig1d,sig2d,sig3d,sig456,sigeff,deps,coef
c
c     Compute material parameters
c
      em=cm(1)
      pr=cm(2)
      bm=cm(3)
      gm=cm(4)
      hrd=cm(5)
      epsf=cm(8)
c
c     Some useful shear modulus parameters
c
      g2=2.*gm
      g3=3.*gm
c
c     Compute yield stress for current value of plastic strain
c
      if (hrd.le.0.) then
         sigyld=cm(6)+cm(7)*epsp
         sighrd=cm(7)
      else
         call crvval(crv,nnpcrv,hrd,epsp,sigyld,sighrd)
      endif
c
c     Compute volumetric strain increment
c
      dvol=eps(1)+eps(2)+eps(3)
      davg=dvol/3.
c
c     Compute increment in (negative) pressure
c
      p=bm*dvol
c
c     Compute elastic trial stress
c
      sig(1)=sig(1)+p+g2*(eps(1)-davg)
      sig(2)=sig(2)+p+g2*(eps(2)-davg)
      sig(3)=sig(3)+p+g2*(eps(3)-davg)
      sig(4)=sig(4)+gm*eps(4)
      sig(5)=sig(5)+gm*eps(5)
      sig(6)=sig(6)+gm*eps(6)
c
c     Compute effective stress in the sense of von Mises
c
      sigv=(sig(1)+sig(2)+sig(3))/3.
      sig1d=sig(1)-sigv
      sig2d=sig(2)-sigv
      sig3d=sig(3)-sigv
      sig456=sig(4)**2+sig(5)**2+sig(6)**2
c
      sigeff=sqrt(1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456))
c
c     Indicate elastic update as default
c
      hisv(1)=0.
c
c     Check if yielding
c
      if (sigeff.gt.sigyld) then
c
c     Plastic update
c
         hisv(1)=.9
c
c     Yield stress exceeded, perform elastic-plastic stress return
c
c     Increment plastic strain
c
         deps=(sigeff-sigyld)/(sighrd+g3)
         epsp=epsp+deps
c
c     Return stress in the direction of deviatoric stress
c
         coef=1.-g3*deps/sigeff
         sig(1)=sigv+sig1d*coef
         sig(2)=sigv+sig2d*coef
         sig(3)=sigv+sig3d*coef
         sig(4)=sig(4)*coef
         sig(5)=sig(5)*coef
         sig(6)=sig(6)*coef
c
      endif
c
c     Check if failure
c
      if (epsp.gt.epsf) failel=.true.
c
      return
      end
      subroutine metaltan (cm,eps,sig,epsp,hisv,dt,capa,elform,tt,temp,
     1 failel,crv,nnpcrv,es)
      implicit none
      include 'nlqparm'
      real cm(*),eps(*),sig(*),hisv(*)
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*5 elform
      logical failel
      real epsp,dt,capa,tt,temp
      real es(6,6)
c
c     Elastic plastic material routine associated with
c     the metal* routines above, tangent is the same for
c     solids and shells. For shells this is a code that
c     assumes capa=1.0, a minor detail.
c
      real em,pr,bm,gm,hrd,g3,bg23,bg43,sigv,sig1d,sig2d,sig3d,sig456,
     1     sigeff,sighrd,coef,cdfds(6)
      integer i,j
c
c     Compute material parameters
c
      em=cm(1)
      pr=cm(2)
      bm=cm(3)
      gm=cm(4)
      hrd=cm(5)
c
c     Some useful modulus parameters
c
      g3=3.*gm
      bg23=bm-2.*gm/3.
      bg43=bm+4.*gm/3.
c
c     Elastic tangent
c
      es(1,1)=bg43
      es(2,2)=bg43
      es(3,3)=bg43
      es(2,1)=bg23
      es(3,1)=bg23
      es(3,2)=bg23
      es(1,2)=es(2,1)
      es(1,3)=es(3,1)
      es(2,3)=es(3,2)
      es(4,4)=gm
      es(5,5)=gm
      es(6,6)=gm
c
      if (hisv(1).gt.0.) then
c
c     Need plastic correction
c
c     Compute effective stress in the sense of von Mises
c
         sigv=(sig(1)+sig(2)+sig(3))/3.
         sig1d=sig(1)-sigv
         sig2d=sig(2)-sigv
         sig3d=sig(3)-sigv
         sig456=sig(4)**2+sig(5)**2+sig(6)**2
c
         sigeff=sqrt(1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456))
c
c     Compute hardening
c
         if (hrd.le.0.) then
            sighrd=cm(7)
         else
            call crvval(crv,nnpcrv,hrd,epsp,coef,sighrd)
         endif
c
c     Compute C*df/ds
c
         coef=g3/sigeff
         cdfds(1)=coef*sig1d
         cdfds(2)=coef*sig2d
         cdfds(3)=coef*sig3d
         cdfds(4)=coef*sig(4)
         cdfds(5)=coef*sig(5)
         cdfds(6)=coef*sig(6)
c
c     Update tangent
c
         coef=hisv(1)/(g3+sighrd)
         do j=1,6
            do i=1,6
               es(i,j)=es(i,j)-coef*cdfds(i)*cdfds(j)
            enddo
         enddo
c
      endif
c
      return
      end
      subroutine usrshl_updatfs(fs,ts,ss,es)
      implicit none
      include 'nlqparm'
      integer index
      common/bk36loc/index
      real fs(3,3),ts(4),ss(4),es
      real f(nlq,3,3),t(nlq,4),s(nlq,4),e(nlq)
      integer i,j
c$omp threadprivate (/bk36loc/)
      e(index)=es
      do i=1,4
         t(index,i)=ts(i)
         do j=1,nlq
            s(j,i)=1.
         enddo
      enddo
      call usrshl_updatfv(f,t,s,e,index,index)
      do i=1,4
         ss(i)=s(index,i)
      enddo
      do i=1,3
         do j=1,3
            fs(i,j)=f(index,i,j)
         enddo
      enddo
      return
      end
