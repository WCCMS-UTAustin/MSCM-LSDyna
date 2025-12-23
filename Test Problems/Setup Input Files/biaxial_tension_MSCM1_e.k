$===============================================================================
$   Solid Three-point Bending Test
$===============================================================================
$ ===   Units: mm ms kg kN GPa                                          ====
$===============================================================================
*KEYWORD
$===============================================================================
*TITLE
Equi-Biaxial Sutured Stretch 
$===============================================================================
$ INCLUDE cards
$===============================================================================
*INCLUDE
full_biaxial_solid_sutured_t1.k
$===============================================================================
$ Parameters
$===============================================================================
*PARAMETER_EXPRESSION
$ Simulation end time
R te      2.0
$ d3plot output time interval
R dtplot  &te/50.
$ binary output time interval
R dtbin   &te/5.
$===============================================================================
$ CONTROL cards
$===============================================================================10.0E-7
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
       &te  
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
     1e-4     0.667         0       0.0      0.00         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     emscl      ihdo
       0.0         0         0                           0.0       0.0         0
$#       -     igado     dtusr    dtdynv
                   1       0.0         0
*CONTROL_SOLID
$# ESORT  FMATRX NIPTETS   SWLOCL  PSFAIL  T10JTOL   ICOH   TET13K
       0       2       5        0       0       0       0        0    
*CONTROL_BULK_VISCOSITY
$#      Q1      Q2    TYPE
       1.5    0.06       1
*CONTROL_HOURGLASS
$#    IHQ      QH      QM
        9     0.1 
$*CONTROL_DYNAMIC_RELAXATION
$    NRCYCK     DRTOL    DRFCTR   DRTERM   TSSFDR   IRELAL   EDTTL   IDRFLG
$     100        0.001     0.995      0        0        0      0.04      3
$    DRPSET
$      1 
$===============================================================================
$ PART cards
$===============================================================================
*SET_PART_LIST
$ sid
   1
$ pid(s)
   1
$===============================================================================
$ ISOPARAMETRIC (IPA) SECTION and SOLID cards
$===============================================================================
*SECTION_SOLID
$#    sid    elform      irl    imass
        1       -18        0        0
$===============================================================================
$ MATERIAL cards
$===============================================================================
*MAT_USER_DEFINED_MATERIAL_MODELS 	   
$#     mid        ro        mt       lmc       nhv    iortho      ibulk         ig
         1   1.41E-6        47        10       145         1         2         1
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0        -1         0         0
$#    aopt      mafc        xp        yp        zp        a1        a2        a3
         2       1.0       0.0       0.0       0.0       1.0       0.0       0.0
$#      v1        v2        v3        d1        d2        d3      beta     ievts
       0.0       0.0       0.0       0.0       1.0       0.0       0.0         0                       
$#     mu         K gamma_fbar   sigma_n     mu_s   Sigma_s       λ_lb      λ_ub      
    1.0E-4      0.10 22.966E-3      33.0    1.118      0.04       1.01      1.24
$#	   nhv     
	   145        
$==============================================================================
$ IMLICIT SOLVER SET UP
$==============================================================================
*CONTROL_IMPLICIT_GENERAL
$   IMFLAG       DT0    IMFORM      NSBS       IGS    CNSTN      FORM    ZERO_V
         0   -0.0001         2         0         1
*CONTROL_IMPLICIT_SOLUTION
$# nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12         1        10       0.1      0.05    1.0E10      0.9   1.0E-08
$#  dnorm   divflag    istif   nlprint    nlnorm    d3itctl
         1         1         1         2         2         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp
$        0         0      0.5         0       0.0
$#   LSMTD     LSDIR      IRAD      SRAD      AWGT      SRED
         5         3         0          0         1        1
*CONTROL_IMPLICIT_SOLVER
$# lsolvr    lprint    negev  ordering       drcm    drcprm   autospc   autotol
         6         1         1         2          4         0        0    1.0E-8
$# lcpack    mtxdmp
         2         0
$
*CONTROL_IMPLICIT_AUTO
$#  iauto    iteopt    itewin     dtmin     dtmax
         1       31              1.E-6     0.001  		 
$
*CONTROL_IMPLICIT_DYNAMICS
$    IMASS     GAMMA      BETA    TDYBIR   TDYDTH   TDYBUR   IRATE    ALPHA
         0       0.5      0.25       0.0     1E28     1E28       0
$===============================================================================
$ BOUNDARY cards
$===============================================================================
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         1         0         0        0         0         0         0         0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         2         0         0        0         0         0         0         0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         3         0         0        0         0         0         0         0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         4         0         0        0         0         0         0         0
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf      vid      death     birth
         2         2         0         1       -1         0    1.0E28       0.0
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf      vid      death     birth
         1         2         0         1        1         0    1.0E28       0.0
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf      vid      death     birth
         3         1         0         1       -1         0    1.0E28       0.0
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf      vid      death     birth
         4         1         0         1        1         0    1.0E28       0.0    
$*LOAD_NODE_SET
$      nsid      dof      lcid        sf       cid       m1         m2        m3
$         4         1         1      1e-7       
$*LOAD_NODE_SET
$      nsid      dof      lcid        sf       cid       m1         m2        m3
$         3         1         1     -1e-7 
$*LOAD_NODE_SET
$      nsid      dof      lcid        sf       cid       m1         m2        m3
$         2         2         1     -1e-3            
$*LOAD_NODE_SET
$      nsid      dof      lcid        sf       cid       m1         m2        m3
$         1         2         1      1e-3    
$==============================================================================
*INCLUDE
four_suture_pin_loc_sets_t1.k
$===============================================================================
$ SET cards
$===============================================================================
$ DEFINE CURVE cards
$===============================================================================
*INCLUDE
load_curve_2s_750stps.k
$===============================================================================
$ DATABASE cards
$===============================================================================
$===============================================================================
*DATABASE_CROSS_SECTION_PLANE
$     psid       xct       yct       zct       xch       ych       zch
       0        0        0.5       0.0       1.0       0.0      0.0
$     xhev      yhev      zhev      lenl      lenm        id     itype
     0.0      0.0        1.0       2.0       1.0         0        0      		
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
     145         0        0         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
       0         0        0         1         1         1         1         0
$# NINTSLD   PKP_SEN     SCLP      HYDRO     MSSCL     THERM    INTOUT    NODOUT
       8         0        1          1         0         0       STRESS   STRESS_GL
$#   DTDT     RESPLT    NEIPB    QUADSLD    CUBSLD   DELERES
       0         0         0         0         0         0
*DATABASE_BINARY_D3PLOT
$#     dt      lcdt      beam     npltc    psetid
   1.0e-2        0         0        0         0
*DATABASE_SPCFORC
$#      dt    binary      lcur     ioopt
   1.0e-1        1         1         0
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2    option3
   1.0e-1        1         1         0         0         0      1         0         0         0     
*END