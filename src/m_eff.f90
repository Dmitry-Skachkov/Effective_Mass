


      Program m_eff
       implicit none
       real(8), parameter     :: pi = 3.141592653589793238462643383279502884197169d0
       real(8), parameter     :: BohrA = 0.52917721067121d0     ! convert Bohr to A
       real(8), parameter     :: hc = 0.1973269804d4            ! Plank constant (h/2pi) * c (eV*A)
       real(8), parameter     :: mec2 = 0.51099895000d6         ! mass of free electron * c^2 (eV)
       real(8)                :: E2k                            ! second derivative d2E/dk2
       real(8)                :: dk(3)                          ! point k close to studied high-symmetry point (crystal representation in units (2*pi/alat)) 
       real(8)                :: dk_a(3)                        ! point k (1/A^2)
       real(8)                :: dk_ad                          ! the absolute value of dk_a
       real(8)                :: E1                             ! energy on the valence (conduction) band in the extreme (eV)
       real(8)                :: E2                             ! energy on the valence (conduction) band close to extreme point by dk (eV)
       real(8)                :: b1(3),b2(3),b3(3)              ! reciprocal vectors 
       real(8)                :: alat                           ! lattice parameter (A)
       real(8)                :: ms                             ! effective mass
       real(8)                :: ms_h                           ! effective mass of hole
       real(8)                :: ms_e                           ! effective mass of electron
       real(8)                :: ms_ex                          ! exciton reduced mass
       real(8)                :: dkx(3,3)                       ! three directions
       integer                :: k
       integer                :: kb
       real(8)                :: E1x(2)
       real(8)                :: E2x(3,2)        

       alat = 3.743d0   ! GaSe                                          ! alat in A 

       dkx(1:3,1) = (/ 0.02d0,   0.00d0,   0.00d0 /)                    ! dk_x
       dkx(1:3,2) = (/ 0.00d0,   0.02d0,   0.00d0 /)                    ! dk_y
       dkx(1:3,3) = (/ 0.00d0,   0.00d0,   0.02d0 /)                    ! dk_z

       E1x(1) =    1.0288d0                                             ! VB at extreme
       E1x(2) =    3.0783d0                                             ! CB at extreme
       
       E2x(1,1) =  1.0176d0                                             ! VB at dk_x
       E2x(2,1) =  1.0115d0                                             ! VB at dk_y
       E2x(3,1) =  1.0288d0                                             ! VB at dk_z

       E2x(1,2) =  3.1069d0                                             ! CB at dk_x
       E2x(2,2) =  3.1069d0                                             ! CB at dk_y
       E2x(3,2) =  3.0783d0                                             ! CB at dk_z

                                                                        ! reciprocal axes: (cart. coord. in units 2 pi/alat)
       b1(1:3) = (/  1.000000d0,  0.577350d0,  0.000000d0 /)            ! GaSe ML 
       b2(1:3) = (/  0.000000d0,  1.154701d0,  0.000000d0 /)  
       b3(1:3) = (/  0.000000d0,  0.000000d0,  0.187150d0 /)  

       b1 = b1*(2.d0*pi/alat)                                           ! convert to (1/A)
       b2 = b2*(2.d0*pi/alat)
       b3 = b3*(2.d0*pi/alat)
       
       do k=1,3                                                          ! 3 directions in k-space
       
        dk(1:3) = dkx(1:3,k)
        print *
        print 1,dk(1:3)
       
        do kb=1,2                                                       ! 1 - VB and 2 - CB
       
         E1 = E1x(kb)                                                   ! energy at extreme
         E2 = E2x(k,kb)                                                 ! energy at nearest point
       
         dk_a(1:3) = dk(1)*b1(1:3)+dk(2)*b2(1:3)+dk(3)*b3(1:3)          ! convert dk to (1/A)
         dk_ad = dsqrt(dk_a(1)**2 + dk_a(2)**2 + dk_a(3)**2)            ! the length of dk

         E2k = 2.d0*(E2 - E1)/(dk_ad**2)                                ! second derivative of energy  (eV*A^2)

         ms = dabs(((hc**2)/E2k)/mec2)                                  ! effective mass in units of free electron mass
         if(kb==1) then
          ms_h = ms                                                     ! hole (in valence band)
         elseif(kb==2) then
          ms_e = ms                                                     ! electron (in conduction band)
         endif        
        enddo

        ms_ex = 1.d0/(1.d0/ms_h+1.d0/ms_e)                               ! reduced mass of exciton
        print 2,ms_e,ms_h,ms_ex
       enddo
1      format('dk=',3F15.4)
2      format(' m_e =',F15.4,'  m_h =',F15.4,'  m_ex =',F15.4)
      end program m_eff





