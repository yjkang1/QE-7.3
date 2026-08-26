! WIY
!
!-----------------------------------------------------------------------
SUBROUTINE comp_dspinldav() 
  !---------------------------------------------------------------------
  !! Initializes the spin rotation matrix for each symmetry operation.
  !! It is needed when symmetrizing the occupation matrix with
  !! noncollinear case.
  !! Based on 'comp_dspinldau' in 'plus_u_full.f90'  
  !
  USE kinds,     ONLY : DP
  USE ldaV,      ONLY : d_spin_ldav
  USE symm_base, ONLY : nsym, sr, t_rev, sname
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  COMPLEX(DP) :: a, b
  INTEGER     :: isym
  INTEGER     :: is1, is2
  REAL(DP), PARAMETER    :: eps=1e-6_DP
  !
  d_spin_ldav = 0.0_DP
  !  
  DO isym = 1, nsym  
     !        
     CALL find_u( sr(1,1,isym), d_spin_ldav(1,1,isym) )  
     !       
     ! ... if time-reversal:  d_spin_ldau --> i sigma_y d_spin_ldau^*         
     !
     DO is1 = 1, 2
     DO is2 = 1, 2
       IF( ABS(d_spin_ldav(is1,is2,isym)**2) < eps) d_spin_ldav(is1,is2,isym) = 0.0_DP
     ENDDO
     ENDDO
     !
     IF (t_rev(isym)==1) THEN        
       a = CONJG( d_spin_ldav(1,1,isym) )        
       b = CONJG( d_spin_ldav(1,2,isym) )        
       d_spin_ldav(1,1,isym) = CONJG( d_spin_ldav(2,1,isym) )        
       d_spin_ldav(1,2,isym) = CONJG( d_spin_ldav(2,2,isym) )        
       d_spin_ldav(2,1,isym) = -a        
       d_spin_ldav(2,2,isym) = -b        
     ENDIF         
     !  
  ENDDO 
  !  
  RETURN
  !
END SUBROUTINE comp_dspinldav   
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_v_nc( ik, wfcatom )
  !---------------------------------------------------------------------
  !! For noncollinear case: builds up the superposition (for a k-point 
  !! \(\text{ik}\)) of pure spin-up or spin-down atomic wavefunctions.
  !! Refer to atomic_wfc.f90
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi, fpi, pi
  USE cell_base,        ONLY : omega, tpiba
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,            ONLY : natomwfc
  USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,            ONLY : xk, igk_k, ngk
  USE wvfct,            ONLY : npwx
  USE uspp_param,       ONLY : upf, nwfcm
  USE noncollin_module, ONLY : noncolin, domag, npol, angle1, angle2, &
                               starting_spin_angle
  USE upf_spinorb,      ONLY : rot_ylm, lmaxx
  !
  implicit none
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom( npwx, npol, natomwfc )
  !! Superposition of atomic wavefunctions
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, npw
  REAL(DP),    ALLOCATABLE :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  COMPLEX(DP), ALLOCATABLE :: sk (:), aux(:)
  COMPLEX(DP) :: kphase, lphase
  REAL(DP)    :: arg

  REAL(DP) :: j
  REAL(DP), ALLOCATABLE :: chiaux(:)
  INTEGER :: nc, ib


  CALL start_clock( 'atomic_wfc' )

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, ntyp
     lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  END DO
  !
  npw = ngk(ik)
  !
  ALLOCATE( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
       gk(3,npw), qg(npw) )
  !
  DO ig = 1, npw
     iig = igk_k (ig,ik)
     gk (1,ig) = xk(1, ik) + g(1,iig)
     gk (2,ig) = xk(2, ik) + g(2,iig)
     gk (3,ig) = xk(3, ik) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  END DO
  !
  !  ylm = spherical harmonics
  !
  CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )
  !
  ! set now q=|k+G| in atomic units
  !
  DO ig = 1, npw
     qg(ig) = SQRT( qg(ig) )*tpiba
  END DO
  !
  CALL interp_atwfc ( npw, qg, nwfcm, chiq )
  !
  DEALLOCATE( qg, gk )
  ALLOCATE( aux(npw), sk(npw) )
  ALLOCATE(chiaux(npw))
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  n_starting_wfc = 0
  !
  DO na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
     !
     !     sk is the structure factor
     !
     DO ig = 1, npw
        iig = igk_k (ig,ik)
        sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                           eigts2 (mill (2,iig), na) * &
                           eigts3 (mill (3,iig), na)
     END DO
     !
     nt = ityp (na)
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           l = upf(nt)%lchi(nb)
           lphase = (0.d0,1.d0)**l
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           j = upf(nt)%jchi(nb)
           !
           !  This routine creates two functions only in the case j=l+1/2 or exit in the
           !  other case 
           !    
           IF (ABS(j-l+0.5_DP)<1.d-4) CYCLE
           !
           !  Find the functions j=l-1/2
           !
           IF (l == 0)  THEN
              chiaux(:)=chiq(:,nb,nt)
           ELSE
              DO ib=1, upf(nt)%nwfc
                 IF ((upf(nt)%lchi(ib) == l).AND. &
                              (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
                    nc=ib
                    EXIT
                 ENDIF
              ENDDO
              !
              !  Average the two functions
              !
              chiaux(:)=(chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
              !
           ENDIF 
           !
           !  and construct the starting wavefunctions as in the noncollinear case.
           !
           DO m = 1, 2 * l + 1
              lm = l**2 + m
              n_starting_wfc = n_starting_wfc + 1
              IF ( n_starting_wfc + 2*l+1 > natomwfc ) CALL errore &
                    ('atomic_wfc_so_mag', 'internal error: too many wfcs', 1)
              DO ig = 1, npw
                 aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
              END DO
              !
              DO ig = 1, npw
                 !
                 ! Now, build the orthogonal wfc
                 !
                 wfcatom(ig,1,n_starting_wfc) = lphase * aux(ig)
                 wfcatom(ig,2,n_starting_wfc) = (0.0_DP, 0.0_DP)
                 !
                 !
                 wfcatom(ig,1,n_starting_wfc+2*l+1) = (0.0_DP, 0.0_DP)
                 !wfcatom(ig,2,n_starting_wfc+2*l+1) = (0.0_DP, 1.0_DP) * lphase * aux(ig)
                 wfcatom(ig,2,n_starting_wfc+2*l+1) = lphase * aux(ig)
                 !
              END DO
           END DO
           n_starting_wfc = n_starting_wfc + 2*l+1
           !
        END IF
        !
     END DO
     !
  END DO

  IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
       'internal error: some wfcs were lost ', 1 )

  DEALLOCATE( aux, sk, chiq, ylm )
  DEALLOCATE( chiaux )

  CALL stop_clock( 'atomic_wfc' )
  
  RETURN

END SUBROUTINE atomic_wfc_v_nc
!
! WIY
SUBROUTINE atomic_wfc_nc_updown_v( ik, wfcatom )
  !-----------------------------------------------------------------------
  !! For noncollinear case: builds up the superposition (for a k-point 
  !! \(\text{ik}\)) of pure spin-up or spin-down atomic wavefunctions.
  !
  !! Based on 'atomic_wfc.f90'
  !! Refer to SUBROUTINE 'atomic_wfc_so_mag' and 'atomic_wfc_nc'
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : tpi, fpi, pi
  USE cell_base,         ONLY : tpiba
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,             ONLY : natomwfc
  USE gvect,             ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,             ONLY : xk, ngk, igk_k
  USE wvfct,             ONLY : npwx, nbnd
  USE uspp_param,        ONLY : upf, nwfcm
  USE noncollin_module,  ONLY : noncolin, npol, angle1, angle2, domag, lspinorb
  USE upf_spinorb,       ONLY : rot_ylm
  USE constants,         ONLY : sqrt2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! the superposition of atomic wavefunctions (up or down)
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, npw
  REAL(DP), ALLOCATABLE :: qg(:), ylm(:,:), chiq(:,:,:), gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk(:), aux(:)
  COMPLEX(DP) :: kphase, lphase
  REAL(DP) :: arg, px, ux, vx, wx
  !
  CALL start_clock( 'atomic_wfc' )
  !
  ! ... calculate max angular momentum required in wavefunctions
  !
  lmax_wfc = 0
  do nt = 1, ntyp
     lmax_wfc = MAX( lmax_wfc, MAXVAL(upf(nt)%lchi(1:upf(nt)%nwfc)) )
  ENDDO
  !
  npw = ngk(ik)
  ALLOCATE ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
       gk(3,npw), qg(npw) )
  !
  DO ig = 1, npw
     iig = igk_k(ig,ik)
     gk(1,ig) = xk(1,ik) + g(1,iig)
     gk(2,ig) = xk(2,ik) + g(2,iig)
     gk(3,ig) = xk(3,ik) + g(3,iig)
     qg(ig) = gk(1,ig)**2 +  gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
  ! ... ylm = spherical harmonics
  !
  CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )
  !
  ! ... set now q=|k+G| in atomic units
  !
  DO ig = 1, npw
     qg(ig) = SQRT(qg(ig))*tpiba
  ENDDO
  !
  CALL interp_atwfc ( npw, qg, nwfcm, chiq )
  !
  DEALLOCATE( qg, gk )
  ALLOCATE( aux(npw), sk(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  n_starting_wfc = 0
  !
  DO na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX(COS(arg), - SIN(arg), KIND=DP)
     !
     ! ... sk is the structure factor
     !
     DO ig = 1, npw
        iig = igk_k(ig,ik)
        sk(ig) = kphase * eigts1 (mill(1,iig), na) * &
                          eigts2 (mill(2,iig), na) * &
                          eigts3 (mill(3,iig), na)
     ENDDO
     !
     nt = ityp (na)
     DO nb = 1, upf(nt)%nwfc
        IF (upf(nt)%oc(nb) >= 0.d0) THEN
           l = upf(nt)%lchi(nb)
           lphase = (0.d0,1.d0)**l
           !
           ! the factor i^l MUST BE PRESENT in order to produce
           ! wavefunctions for k = 0 that are real in real space
           !
           IF ( upf(nt)%has_so ) THEN
              !
              CALL wfc_atom( .TRUE. )
              !
           ELSE
              !
              CALL wfc_atom( .FALSE. )
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  IF (n_starting_wfc /= natomwfc) CALL errore( 'atomic_wfc_nc_updown', &
                             'internal error: some wfcs were lost ', 1 )
  !
  DEALLOCATE( aux, sk, chiq, ylm )
  !
  CALL stop_clock ( 'atomic_wfc' )
  !
  RETURN
  !
CONTAINS
   !
   !--------------------------
   SUBROUTINE wfc_atom( soc )
      !---------------------------
      !
      LOGICAL :: soc
      !! .TRUE. if the fully-relativistic pseudo
      !
      ! ... local variables
      !
      REAL(DP) :: j, gamman
      REAL(DP), ALLOCATABLE :: chiaux(:)
      INTEGER :: nc, ib 
      !
      ! ... If SOC go on only if j=l+1/2
      IF (soc) j = upf(nt)%jchi(nb)
      IF (soc .AND. ABS(j-l+0.5_DP)<1.d-4 ) RETURN
      !
      ALLOCATE( chiaux(npw) )
      !
      IF (soc) THEN 
        !
        ! ... Find the index for j=l-1/2
        !
        IF (l == 0)  THEN
           chiaux(:)=chiq(:,nb,nt)
        ELSE
           DO ib=1, upf(nt)%nwfc
              IF ((upf(nt)%lchi(ib) == l) .AND. &
                           (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
                 nc=ib
                 exit
              ENDIF
           ENDDO
           !
           ! ... Average the two radial functions 
           !
           chiaux(:) = (chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
        ENDIF
        !
      ELSE
        !
        chiaux(:) = chiq(:,nb,nt)
        !
      ENDIF
      !
      DO m = 1, 2*l+1
         lm = l**2 + m
         n_starting_wfc = n_starting_wfc + 1
         !
         IF (n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
               ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
         !
         aux=(0.d0,0.d0)
         !
         ! temp
         !IF ( soc == .true. ) THEN
         !    !
         !    IF ( l == 1 ) THEN
         !       IF ( m == 2 ) THEN
         !           aux(:) = CMPLX((-1.d0)/sqrt2, 0.d0, KIND=DP) * ylm(:,lm) + & 
         !                    CMPLX(0.d0, (-1.d0)/sqrt2, KIND=DP) * ylm(:,lm+1) 
         !       ELSE IF ( m == 3 ) THEN
         !           aux(:) = CMPLX(1.d0/sqrt2, 0.d0, KIND=DP) * ylm(:,lm-1) + & 
         !                    CMPLX(0.d0, (-1.d0)/sqrt2, KIND=DP) * ylm(:,lm)
         !       ELSE
         !           aux(:) = ylm(:,lm)
         !       ENDIF
         !       !
         !    ELSE IF ( l == 2 ) THEN
         !       IF ( m == 2 ) THEN
         !           aux(:) = CMPLX((-1.d0)/sqrt2, 0.d0, KIND=DP) * ylm(:,lm) + & 
         !                    CMPLX(0.d0, (-1.d0)/sqrt2, KIND=DP) * ylm(:,lm+1)    
         !       ELSE IF ( m == 3 ) THEN
         !           aux(:) = CMPLX(1.d0/sqrt2, 0.d0, KIND=DP) * ylm(:,lm-1) + &
         !                    CMPLX(0.d0, (-1.d0)/sqrt2, KIND=DP) * ylm(:,lm)
         !       ELSE IF ( m == 4 ) THEN
         !           aux(:) = CMPLX(1.d0/sqrt2, 0.d0, KIND=DP) * ylm(:,lm) + & 
         !                    CMPLX(0.d0, (1.d0)/sqrt2, KIND=DP) * ylm(:,lm+1)
         !       ELSE IF ( m == 5 ) THEN
         !           aux(:) = CMPLX(1.d0/sqrt2, 0.d0, KIND=DP) * ylm(:,lm-1) + &
         !                    CMPLX(0.d0, (-1.d0)/sqrt2, KIND=DP) * ylm(:,lm)
         !       ELSE
         !           aux(:) = ylm(:,lm)
         !       ENDIF
         !    ELSE
         !        aux(:) = ylm(:,lm)
         !    ENDIF
         !ELSE
         !   aux(:) = ylm(:,lm)
         !ENDIF
         ! temp
         !
         DO ig = 1, npw
            aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
         ENDDO
         ! 
         DO ig = 1, npw
            !
            wfcatom(ig,1,n_starting_wfc) = lphase*aux(ig)
            !wfcatom(ig,1,n_starting_wfc) = lphase*sk(ig)*aux(ig)*chiaux(ig)
            wfcatom(ig,2,n_starting_wfc) = (0.0_dp, 0.0_dp)
            !
            wfcatom(ig,1,n_starting_wfc+2*l+1) = (0.0_dp, 0.0_dp)
            !wfcatom(ig,2,n_starting_wfc+2*l+1) = lphase*sk(ig)*CONJG(aux(ig))*chiaux(ig)
            wfcatom(ig,2,n_starting_wfc+2*l+1) = lphase * aux(ig)
            !
         ENDDO
      ENDDO
      !
      n_starting_wfc = n_starting_wfc + 2*l+1
      !
      DEALLOCATE( chiaux )
      !
   END SUBROUTINE wfc_atom
   !
   !
END SUBROUTINE atomic_wfc_nc_updown_v
!

