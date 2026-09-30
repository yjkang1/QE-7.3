MODULE ldaV

  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  USE basis,      ONLY : natomwfc
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE uspp_param, ONLY : upf
  USE io_global,  ONLY : stdout

  SAVE

  COMPLEX(DP), ALLOCATABLE ::   &
    wfcv(:,:),                  & ! atomic wavefunctions having V correction
    ! WIY
    d_spin_ldav(:,:,:)            ! the rotations in spin space for all symmetries
    ! WIY

  REAL(DP) ::                   &
    ehubene,                    & ! Hubbard contribution to total energy
    ehubene_old,                &
    dehubene,                   &
    ehub_nn_distance,           & ! considered as the nearest neighbors
                                  ! if the distance is less than "ehub_nn_distance"
    ehub_conv_thr,              &
    ehub_mixing

  REAL(DP), ALLOCATABLE ::      &
    ehub_v(:,:,:,:,:),          & ! Hubbard inter-site interaction energy V
                                  ! ehub_v(na1, na2, m1, m2, ncells)
                                  ! na1, na2: atom indices
                                  ! m1, m2: orbital indices
                                  ! ncells: nearest-neighbor cell indices
    ehub_u(:,:,:),              & ! ehub_u(na, m1, m2)
    ehub_v_diff(:,:,:,:,:),     &
    ehub_u_diff(:,:,:),         &
    ehub_u_raw(:,:,:),          & ! raw (unmixed) one-shot U of the current SCF step [rawprint]
    ehub_v_raw(:,:,:,:,:),      & ! raw (unmixed) one-shot V of the current SCF step [rawprint]
    ehub_nn_leng(:,:,:,:,:),    &
    ehub_nn_leng_all(:,:,:,:,:) 

  INTEGER ::                    &
    nwfcv,                      & ! total # of atomic wavefunctions having V correction
    ehub_lm(ntypx),             & ! sum of angular momenta of V-corrected orbitals of each types
    ehub_lm_max,                & ! maximum value of "Hubbard_lm"
    ehub_l_choice(ntypx,6),     &
    ehub_l_choice_inp(ntypx,6), &
    remove_ehub_u(ntypx,6),     &
    remove_ehub_v(ntypx,6),     &
    acbn0_type,                 & ! 1: Hubbard, 2: extended Hubbard, 3: orbital-dependent extended Hubbard
    ehub_nn_dim(3),             &
    ehub_nnnc

  ! [hubdat] labels of the per-call blocks in outdir/eACBN0_hub_stock_*.dat
  INTEGER ::                    &
    ehub_scf_iter = 0,          & ! current SCF iteration, set in electrons_scf
    ehub_dat_iter = -1,         & ! SCF iteration of the last block written
    ehub_dat_step = -1,         & ! ionic step (istep) of the last block written
    ehub_dat_pass = 0,          & ! kernel calls so far inside that iteration
    ehub_dat_nblk = 0             ! blocks written in this run (0 = replace files)

  INTEGER, ALLOCATABLE ::       &
    ofst(:), ofstv(:),          & ! offset of atomic wavfunctions used for projections
    ehub_nn_vect(:,:,:,:,:),    & ! nearest neighbor(1) or not(0)
                                  ! ex) ehub_nn_vect(a, b, c, na1, na2)
                                  ! a, b, c: cell indices
                                  ! na1: atom located at home unit cell
                                  ! na2: atom located at (a,b,c) cell
    ehub_l(:,:),                &
    counter_lm(:,:)               ! counter for calculating occupation number in noncollinear calculations

  INTEGER, PARAMETER ::         &
    ehub_csize = 10

  LOGICAL ::                    &
    lda_plus_v,                 & ! .TRUE. if extended Hubbard model calculation is performed
    lacbn0,                     & ! .TRUE. if ACBN0 functional is used
    pot_fix,                    & ! .TRUE. for fixed Hubbard potential calculation
    read_ehub_ns_file = .FALSE.,&
    read_ehub_uv_file = .FALSE.,&
    stop_acbn0        = .FALSE.

CONTAINS

  SUBROUTINE init_lda_plus_v( noncolin, nspin )
 
    IMPLICIT NONE
 
    LOGICAL, INTENT(IN) :: noncolin

    INTEGER, INTENT(IN) :: nspin
 
    LOGICAL             :: neq

    INTEGER             :: na1, na2, nt1, nt2, m1, m2
    INTEGER             :: nw, na, nt, m, l, j, n, counter, counter_
    INTEGER             :: n1, n2, n3, nnc, ics
    INTEGER             :: nmax1, nmax2, nmax3
    INTEGER             :: l_choice(ntypx,6)
    INTEGER, PARAMETER  :: l_wfc = 6 ! 1 = s, 2 = p, 3 = d, 4 = f, 5 = d for semicore

    REAL(DP), PARAMETER :: eps = 1.0D-2
    REAL(DP), PARAMETER :: eps2 = 1.0D-6

    CHARACTER(LEN=2)    :: orbital
    CHARACTER(LEN=1)    :: spdf(0:3) = ['S','P','D','F']

    IF ( .NOT. lda_plus_v ) THEN
       ehub_lm_max = 0
       RETURN
    ENDIF
    !
    IF ( nspin == 1 ) &
       CALL errore('init_lda_plus_v', 'nspin == 1 is not supported yet', 1)
    ! WIY
    ! Following Do loops for checking pp is scalar or
    ! full relativistic.
    ! Full relativistic pp is not working on calculations
    ! without lspinorb = .true. in case of eACBN0 calculation.
    !
    DO nt = 1, ntyp
        !
        DO nw = upf(nt)%nwfc, 1, -1
            !
            orbital = trim(upf(nt)%els(nw))
            !
            l = upf(nt)%lchi(nw)
            !
            IF(spdf(l) /= orbital(2:)) THEN
                CALL errore('init_lda_plus_v', 'Full relativistic pp can be used with lspinorb = .true. in case of eACBN0 calculation.', 1)
            ENDIF
            !
        ENDDO
    ENDDO
    !
    IF ( noncolin ) THEN
       IF ( .NOT. ALLOCATED (d_spin_ldav) ) ALLOCATE( d_spin_ldav(2,2,48))
       CALL comp_dspinldav()
    ENDIF
    ! WIY
    ! ... reset ehub_l_choice
    !
    l_choice(:,:) = 0
    !
    ! WIY
    DO nt = 1, ntyp
      !
      DO nw = 1, l_wfc
         !
         n = 0
         !
         IF ( nw == 3 ) THEN ! 3 is for d orbital
            !
            DO na = upf(nt)%nwfc, 1, -1
               !
               IF ( upf(nt)%lchi(na) /= nw-1 ) CYCLE
               !
               orbital = trim(upf(nt)%els(na)) ! label in PP_CHI
               !
               IF ( n < ICHAR(orbital(1:1)) ) THEN
                  !
                  n = ICHAR(orbital(1:1))
                  !
                  counter_ = na
                  !
               ENDIF
               !
            ENDDO
            !
            IF ( ehub_l_choice(nt,nw) == 1) THEN
            !    
                IF ( noncolin ) THEN
                    l_choice(nt,nw) = counter_
                ELSE
                    l_choice(nt,counter_) = 1
                ENDIF
            !
            ENDIF
            !
         !ENDIF
         !!
         ELSEIF ( ehub_l_choice(nt,nw) == 1 .AND. nw <= 4 .AND. nw /= 3 ) THEN
            !
            DO na = upf(nt)%nwfc, 1, -1
               !
               IF ( upf(nt)%lchi(na) /= nw-1 ) CYCLE !.OR. nw == 3 ) CYCLE
               !
               orbital = trim(upf(nt)%els(na)) ! label in PP_CHI
               !
               IF ( n < ICHAR(orbital(1:1)) ) THEN
                  !
                  n = ICHAR(orbital(1:1))
                  !
                  counter = na
                  !
               ENDIF
               !
            ENDDO
            !
            !IF ( nw == 3 ) counter = counter_
            !
            IF ( noncolin ) THEN
               l_choice(nt,nw) = counter
            ELSE
               l_choice(nt,counter) = 1
            ENDIF
            !
            ! ... if valence d & semicore d are coexisted in pseudofile, for
            ! ... semicore d orbital -> ehub_l_choice(type,5)
            !
         ELSEIF ( ehub_l_choice(nt,nw) == 1 .AND. nw == 5 ) THEN
            !
            DO na = upf(nt)%nwfc, 1, -1
               !
               IF ( upf(nt)%lchi(na) /= nw-3 ) CYCLE
               !
               orbital = trim(upf(nt)%els(na))
               !
               IF ( n < ICHAR(orbital(1:1)) ) THEN
                  !
                  n = ICHAR(orbital(1:1))
                  !
                  IF ( na /= counter_ ) counter = na
                  !
               ENDIF
               !
            ENDDO
            !
            IF ( noncolin ) THEN
               l_choice(nt,nw) = counter
            ELSE
               l_choice(nt,counter) = 1
            ENDIF
            !
         ENDIF
         !
      ENDDO
      !
    ENDDO
    !
    ehub_l_choice_inp(:,:) = ehub_l_choice(:,:)
    ehub_l_choice(:,:) = l_choice(:,:)
    ! ... ehub_lm max
    !
    ehub_lm_max = -1
    !
    ehub_lm(:) = 0
    !
    DO nt = 1, ntyp
       !
       DO nw = 1, upf(nt)%nwfc
          !
          l = upf(nt)%lchi(nw)
          !
          IF ( noncolin ) THEN
            IF ( upf(nt)%oc(nw) >= 0.0_DP .AND. ehub_l_choice(nt, l+1) == nw ) THEN
               ! Below part is for the noncollinear ACBN0 functionals.
               !IF ( upf(nt)%has_so ) THEN
               !   ! j = l-1/2, degeneracy 2*l
               !   ehub_lm(nt) = ehub_lm(nt) + 2*l
               !   ! j = l+1/2, degeneracy 2*l+2
               !   IF ( ABS( upf(nt)%jchi(nw)-l-0.5_DP) < eps2 ) &
               !   & ehub_lm(nt) = ehub_lm(nt) + 2
               !   !
               !   IF ( (nw-1 > 0) .AND. ( TRIM(upf(nt)%els(nw)) == TRIM(upf(nt)%els(nw-1)))) THEN
               !      ! 
               !      ehub_lm(nt) = ehub_lm(nt) + 2*l
               !      !
               !      IF ( ABS( upf(nt)%jchi(nw-1)-l-0.5_DP) < eps2 ) &
               !      & ehub_lm(nt) = ehub_lm(nt) + 2
               !      !
               !   ENDIF                     
               !ELSE
               !   ehub_lm(nt) = ehub_lm(nt) + 2*(2*l+1)
               !ENDIF 
               !
               ehub_lm(nt) = ehub_lm(nt) + 2*l + 1
            ENDIF            
          ELSE 
            IF ( upf(nt)%oc(nw) >= 0.0_DP .AND. ehub_l_choice(nt,nw) == 1 ) THEN
               !
               ehub_lm(nt) = ehub_lm(nt) + 2*l + 1
               !
            ENDIF
          ENDIF
          !
       ENDDO
       !
       ehub_lm_max = MAX(ehub_lm_max, ehub_lm(nt))
       !
    ENDDO
    !
    ! ... ehub_l
    !
    IF ( .NOT. ALLOCATED(ehub_l) ) ALLOCATE ( ehub_l(ehub_lm_max, ntyp) ) ; ehub_l(:,:) = 0
    !
    DO nt = 1, ntyp
       !
       counter = 1
       !
       DO nw = 1, upf(nt)%nwfc
          !
          l = upf(nt)%lchi(nw)
          !
          IF ( noncolin ) THEN
            IF ( upf(nt)%oc(nw) >= 0.0_DP .AND. ehub_l_choice(nt, l+1) == nw ) THEN
               ! Below part is for the noncollinear ACBN0 functionals.
               !IF ( upf(nt)%has_so ) THEN
               !   !
               !   IF ( ABS( upf(nt)%jchi(nw)-l-0.5_DP) < eps2 ) THEN
               !      ! j = l+1/2, degeneracy 2*l+2
               !      DO m = counter, counter + 2*l + 2
               !         ehub_l(m,nt) = l+1
               !      ENDDO
               !      !
               !      counter = counter + 2*l + 2
               !      !
               !   ELSE
               !      ! j = l-1/2, degeneracy 2*l
               !      DO m = counter, counter + 2*l 
               !         ehub_l(m,nt) = l+1
               !      ENDDO
               !      !
               !      counter = counter 2*l
               !      !
               !   ENDIF
               !   !
               !   IF ( (nw-1 > 0) .AND. ( TRIM(upf(nt)%els(nw)) == TRIM(upf(nt)%els(nw-1)))) THEN
               !      !
               !      IF ( ABS( upf(nt)%jchi(nw-1)-l-0.5_DP) < eps2 ) THEN
               !         ! j = l+1/2, degeneracy 2*l+2
               !         DO m = counter, counter + 2*l + 2
               !            ehub_l(m,nt) = l+1
               !         ENDDO
               !         !
               !         counter = counter + 2*l + 2
               !         !
               !      ELSE
               !         ! j = l-1/2, degeneracy 2*l
               !         DO m = counter, counter + 2*l 
               !            ehub_l(m,nt) = l+1
               !         ENDDO
               !         !
               !         counter = counter + 2*l 
               !         !
               !      ENDIF
               !      !
               !   ENDIF
               !ELSE 
               !   !
               !   DO m = counter, counter + 2*(2*l + 1)
               !      ehub_l(m,nt) = l + 1
               !   ENDDO
               !   !
               !   counter = counter + 2*(2*l + 1)
               !   !
               !ENDIF 
               !
               DO m = counter, counter + 2*l
                  ehub_l(m,nt) = l + 1
               ENDDO 
               !
               counter = counter + 2*l + 1
               !
            ENDIF
            !
          ELSE
            !
            IF ( upf(nt)%oc(nw) >= 0.0_DP .AND. ehub_l_choice(nt,nw) == 1 ) THEN
               !
               DO m = counter, counter + 2*l
                  ehub_l(m,nt) = l + 1
               ENDDO
               !
               counter = counter + 2*l + 1
               !
            ENDIF
            !
          ENDIF
       ENDDO
       !
    ENDDO
    !
    ! WIY
    IF (noncolin) THEN
      !
      IF ( .NOT. ALLOCATED(counter_lm) ) ALLOCATE ( counter_lm(ehub_lm_max,ntyp) ); counter_lm(:,:) = 0
      !
      DO nt = 1, ntyp
         !
         DO m = 2, ehub_lm(nt)
            !
            counter_lm(m,nt) = counter_lm(m-1,nt)
            !
            IF (ehub_l(m,nt) /= ehub_l(m-1,nt)) &
            counter_lm(m,nt) = counter_lm(m,nt) + 2*(ehub_l(m-1,nt)-1)+1
            !
         ENDDO
         !
      ENDDO
      !
    ENDIF
    ! WIY
    IF ( ehub_lm_max == 0 ) CALL errore('init_lda_plus_v', 'set ehub_l_choice', 1)
    !
    ! ... ofstv
    !
    IF ( .NOT. ALLOCATED(ofst) ) ALLOCATE ( ofst(nat) )
    !
    CALL offset_atom_wfc_v(.false., ofst, nwfcv)
    !
    IF ( nwfcv .NE. natomwfc ) &
       CALL errore('offset_atom_wfc_v', 'wrong number of wavefunctions', 1)
    !
    IF ( .NOT. ALLOCATED(ofstv) ) ALLOCATE ( ofstv(nat) )
    !
    CALL offset_atom_wfc_v(.true., ofstv, nwfcv)
    !
    ics = ehub_csize + 1
    !
    IF ( .NOT. lacbn0 .OR. .NOT. acbn0_type == 1 ) THEN
       !
       ! ... ehub_nn_vect, ehub_nn_leng
       !
       IF ( .NOT. ALLOCATED(ehub_nn_vect) ) &
          ALLOCATE ( ehub_nn_vect(2*ehub_csize+1,2*ehub_csize+1,2*ehub_csize+1,nat,nat) )
       IF ( .NOT. ALLOCATED(ehub_nn_leng) ) &
          ALLOCATE ( ehub_nn_leng(2*ehub_csize+1,2*ehub_csize+1,2*ehub_csize+1,nat,nat) )
       IF ( .NOT. ALLOCATED(ehub_nn_leng_all) ) &
          ALLOCATE ( ehub_nn_leng_all(2*ehub_csize+1,2*ehub_csize+1,2*ehub_csize+1,nat,nat) )
       !
       CALL hubbard_nns()
       !
       IF (MAXVAL(ehub_nn_dim(:)) > 10) &
          CALL errore('hubbard_nns', 'reduce the ehub_nn_distance (ehub_nn_dim > 10)', 1)
       !
       ! Print the nearest neighbors vectors
       !
       WRITE ( stdout,'(a)' ) '     -----------------------------------------------------------'
       WRITE ( stdout,'(a)' ) '     List of the nearest neighbors for inter-site interactions V'
       WRITE ( stdout,'(a)' ) '     -----------------------------------------------------------'
       !
       DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
       DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
       DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
          !
          DO na1 = 1, nat
          DO na2 = 1, nat
             !
             IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
             !
             WRITE ( stdout,'(a,i2,a,i2,a,3i6)' ) &
                   '     atom 1: ', na1, ',  atom 2: ', na2, ',',  n1, n2, n3
          ENDDO ! na2
          ENDDO ! na1
          !
       ENDDO ! n3
       ENDDO ! n2
       ENDDO ! n1
       !
       WRITE( stdout,'(a)' ) ' '
       !
    ENDIF
    !
    IF ( lacbn0 .AND. stop_acbn0 ) THEN
       IF ( acbn0_type == 1 ) THEN
          ALLOCATE ( ehub_u_diff(nat,ehub_lm_max,ehub_lm_max) )
          ehub_u_diff = 1000.0_DP
       ELSE
          ALLOCATE ( ehub_v_diff(nat,nat,ehub_lm_max,ehub_lm_max,ehub_nnnc) )
          ehub_v_diff = 1000.0_DP
       ENDIF
    ENDIF
    !
    IF ( lda_plus_v .AND. .NOT.lacbn0 .AND. .NOT.read_ehub_ns_file .AND. .NOT.read_ehub_uv_file ) THEN
       !
       ! ... initialize U & V values
       !
       !
       ! Set U for same type of atoms
       !
       na1 = 1
       neq = .false.
       !
       DO nt1 = 1, ntyp
          DO na2 = na1, nat
             !
             IF ( nt1 == ityp(na2) ) THEN
                !
                DO m1 = 1, ehub_lm(nt1)
                DO m2 = 1, ehub_lm(nt1)
                   ehub_v(na2,na2,m1,m2,1) = ehub_v(na1,na1,m1,m2,1)
                ENDDO
                ENDDO
                !
             ELSE
                !
                IF (.not. neq) THEN
                   na = na2
                   neq = .true.
                ENDIF
                !
             ENDIF
             !
          ENDDO
          !
          na1 = na
          neq = .false.
          !
       ENDDO
       !
       ! Hermiticity of Hubbard interaction parameters U & V (in the unit cell)
       !
       DO na1 = 1, nat
          nt1 = ityp(na1)
          !
       DO na2 = na1, nat
          nt2 = ityp(na2)
          !
          DO m1 = 1, ehub_lm(nt1)
          DO m2 = 1, ehub_lm(nt2)
             ehub_v(na2,na1,m2,m1,1) = ehub_v(na1,na2,m1,m2,1)
          ENDDO
          ENDDO
          !
       ENDDO
       ENDDO
       !
       ! Find the first nearest neighbors within adjacent negibor cells
       ! and set the inter-site interaction parameters V
       !
       nnc = 2
       !
       DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
       DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
       DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
          !
          IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
          !
          DO na1 = 1, nat
             nt1 = ityp(na1)
             !
          DO na2 = 1, nat
             nt2 = ityp(na2)
             !
             IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
             !
             IF ( (ehub_nn_leng_all(n1+ics,n2+ics,n3+ics,na1,na2) - MINVAL(ehub_nn_leng(:,:,:,na1,na2))) > eps .OR. na1 == na2 ) CYCLE
             !
             DO m1 = 1, ehub_lm(nt1)
             DO m2 = 1, ehub_lm(nt2)
                ehub_v(na1,na2,m1,m2,nnc) = ehub_v(na1,na2,m1,m2,1)
             ENDDO ! m2
             ENDDO ! m1
             !
          ENDDO ! na2
          ENDDO ! na1
          !
          IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
          !
       ENDDO ! n3
       ENDDO ! n2
       ENDDO ! n1
       !
       DO na1 = 1, nat
       DO na2 = 1, nat
          IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) &
             ehub_v(na1,na2,:,:,1) = 0.0_DP
       ENDDO
       ENDDO
       !
    ENDIF
    !
  END SUBROUTINE init_lda_plus_v
  !
  SUBROUTINE deallocate_ldaV(flag)

  LOGICAL, INTENT (IN) :: flag

  IF ( flag ) THEN
     IF ( ALLOCATED(ofst) )             DEALLOCATE ( ofst )
     IF ( ALLOCATED(ofstv) )            DEALLOCATE ( ofstv )
     IF ( ALLOCATED(ehub_nn_vect) )     DEALLOCATE ( ehub_nn_vect )
     IF ( ALLOCATED(ehub_nn_leng) )     DEALLOCATE ( ehub_nn_leng )
     IF ( ALLOCATED(ehub_nn_leng_all) ) DEALLOCATE ( ehub_nn_leng_all ) 
     IF ( ALLOCATED(ehub_v) )           DEALLOCATE ( ehub_v )
     IF ( ALLOCATED(ehub_u) )           DEALLOCATE ( ehub_u )
     IF ( ALLOCATED(ehub_u_raw) )       DEALLOCATE ( ehub_u_raw )
     IF ( ALLOCATED(ehub_v_raw) )       DEALLOCATE ( ehub_v_raw )
     IF ( ALLOCATED(ehub_l) )           DEALLOCATE ( ehub_l )
     IF ( ALLOCATED(ehub_v_diff) )      DEALLOCATE ( ehub_v_diff )
     IF ( ALLOCATED(ehub_u_diff) )      DEALLOCATE ( ehub_u_diff )
     IF ( ALLOCATED(d_spin_ldav) )      DEALLOCATE ( d_spin_ldav )
  ENDIF
  IF (ALLOCATED(wfcv)) DEALLOCATE (wfcv)
 

  END SUBROUTINE deallocate_ldaV


  SUBROUTINE copy_v_wfc( swfcatom, noncolin )
  !
  !  Copy (orthogonalized) atomic wavefunctions "swfcatom"
  !  having a extended Hubbard correction to array "wfcv"
  !
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)           :: swfcatom(:,:)
  LOGICAL    , INTENT(IN), OPTIONAL :: noncolin
  INTEGER                           :: na, nt, m1, m2

  DO na = 1, nat
     nt = ityp(na)
     !
     m1 = 1
     m2 = ehub_lm(nt)
     !
     IF ( PRESENT(noncolin) ) THEN
        IF ( noncolin ) m2 = 2*m2
     ENDIF
     !
     wfcv(:,ofstv(na)+m1:ofstv(na)+m2) = swfcatom(:,ofst(na)+m1:ofst(na)+m2)
     !
  ENDDO

  END SUBROUTINE copy_v_wfc


  SUBROUTINE write_ehub_uv

  IMPLICIT NONE

  INTEGER :: na, na1, na2
  INTEGER :: nt, nt1, nt2
  INTEGER :: m1, m2
  INTEGER :: n1, n2, n3, nnc, ics

  IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
     !
     ! Write U matrix
     !
     WRITE ( stdout,'(a)' ) '                          '
     WRITE ( stdout,'(a)' ) '     ---------------------'
     WRITE ( stdout,'(a)' ) '     U[atom i, #orb, #orb]'
     WRITE ( stdout,'(a)' ) '     ---------------------'
 
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        WRITE ( stdout,'(3(a,i2),a)' ) &
              '     [ ', na, ', ', 2*MAXVAL(ehub_l(:,nt))-1, ', ', 2*MAXVAL(ehub_l(:,nt))-1, '  ]'
        !
        DO m1 = 1, ehub_lm(nt)
           WRITE ( stdout,'(4x,16f9.4)' ) &
                 (ehub_u(na,m1,m2) * rytoev, m2 = 1, ehub_lm(nt))
        ENDDO
        !
        WRITE ( stdout,* ) ''
        !
     ENDDO ! na
     !
  ELSE
     !
     ! Write V matrices
     !
     WRITE ( stdout,'(a)' ) '                                        '
     WRITE ( stdout,'(a)' ) '     -----------------------------------'
     WRITE ( stdout,'(a)' ) '     V[atom i, atom j, #orb(i), #orb(j)]'
     WRITE ( stdout,'(a)' ) '     -----------------------------------'
     !
     ics = ehub_csize + 1
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = 1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        WRITE ( stdout,'(4(a,i2),a)' ) &
              '     [ ', na1, ', ', na2, ', ', ehub_lm(nt1), ', ', ehub_lm(nt2), '  ]'
        !
        WRITE ( stdout,'(a,i2,2i6,a)' ) '     [ ', 0, 0, 0, '  ]'
        !
        DO m1 = 1, ehub_lm(nt1)
           WRITE ( stdout,'(4x,16f9.4)' ) &
                 (ehub_v(na1,na2,m1,m2,1) * rytoev, m2 = 1, ehub_lm(nt2))
        ENDDO
        !
        WRITE ( stdout,* ) ''
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     nnc = 2
     !
     DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
     DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
     DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
        !
        IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
        !
        DO na1 = 1, nat
           nt1 = ityp(na1)
           !
           IF ( ehub_lm(nt1) == 0 ) CYCLE
           !
        DO na2 = 1, nat
           nt2 = ityp(na2)
           !
           IF ( ehub_lm(nt2) == 0 ) CYCLE
           !
           IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           !
           WRITE ( stdout,'(4(a,i2),a)' ) &
                 '     [ ', na1, ', ', na2, ', ', ehub_lm(nt1), ', ', ehub_lm(nt2), '  ]'
           !
           WRITE ( stdout,'(a,i2,2i6,a)' ) '     [ ', n1, n2, n3, '  ]'
           !
           DO m1 = 1, ehub_lm(nt1)
              WRITE ( stdout,'(4x,16f9.4)' ) &
                    (ehub_v(na1,na2,m1,m2,nnc) * rytoev, m2 = 1, ehub_lm(nt2))
           ENDDO
           !
           WRITE ( stdout,* ) ''
           !
        ENDDO ! na2
        ENDDO ! na1
        !
        IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
        !
     ENDDO ! n3
     ENDDO ! n2
     ENDDO ! n1
  ENDIF


  ! [rawprint] reference print: mixed (module) value next to the raw one-shot
  ! value of this step, full precision, eV.  Only meaningful with lacbn0.
  !
  IF ( lacbn0 ) THEN
     IF ( acbn0_type == 1 ) THEN
        IF ( ALLOCATED(ehub_u_raw) ) THEN
           WRITE ( stdout,'(/5x,a,f6.3,a)' ) &
              'U (eV) MIXED   : U_n = (1-m)*U_{n-1} + m*U_raw,  m = ehub_mixing = ', ehub_mixing
           CALL rawprint_u ( ehub_u )
           WRITE ( stdout,'(/5x,a)' ) &
              'U (eV) RAW     : one-shot U_raw from the current wavefunctions, no mixing'
           CALL rawprint_u ( ehub_u_raw )
        ENDIF
     ELSE
        IF ( ALLOCATED(ehub_v_raw) ) THEN
           WRITE ( stdout,'(/5x,a,f6.3,a)' ) &
              'V (eV) MIXED   : V_n = (1-m)*V_{n-1} + m*V_raw,  m = ehub_mixing = ', ehub_mixing
           CALL rawprint_v ( ehub_v )
           WRITE ( stdout,'(/5x,a)' ) &
              'V (eV) RAW     : one-shot V_raw from the current wavefunctions, no mixing'
           CALL rawprint_v ( ehub_v_raw )
        ENDIF
     ENDIF
  ENDIF

  CONTAINS

  SUBROUTINE rawprint_u ( u )
  REAL(DP), INTENT(IN) :: u(:,:,:)
  INTEGER :: na, nt, m1, m2
  DO na = 1, nat
     nt = ityp(na)
     IF ( ehub_lm(nt) == 0 ) CYCLE
     WRITE ( stdout,'(3(a,i2),a)' ) &
           '     [ ', na, ', ', 2*MAXVAL(ehub_l(:,nt))-1, ', ', 2*MAXVAL(ehub_l(:,nt))-1, '  ]'
     DO m1 = 1, ehub_lm(nt)
        WRITE ( stdout,'(4x,16f16.10)' ) (u(na,m1,m2) * rytoev, m2 = 1, ehub_lm(nt))
     ENDDO
     WRITE ( stdout,* ) ''
  ENDDO
  END SUBROUTINE rawprint_u

  SUBROUTINE rawprint_v ( v )
  REAL(DP), INTENT(IN) :: v(:,:,:,:,:)
  INTEGER :: na1, na2, nt1, nt2, m1, m2, n1, n2, n3, nnc, ics
  ics = ehub_csize + 1
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( ehub_lm(nt1) == 0 ) CYCLE
  DO na2 = 1, nat
     nt2 = ityp(na2)
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     WRITE ( stdout,'(4(a,i2),a)' ) &
           '     [ ', na1, ', ', na2, ', ', ehub_lm(nt1), ', ', ehub_lm(nt2), '  ]'
     WRITE ( stdout,'(a,i2,2i6,a)' ) '     [ ', 0, 0, 0, '  ]'
     DO m1 = 1, ehub_lm(nt1)
        WRITE ( stdout,'(4x,16f16.10)' ) (v(na1,na2,m1,m2,1) * rytoev, m2 = 1, ehub_lm(nt2))
     ENDDO
     WRITE ( stdout,* ) ''
  ENDDO ! na2
  ENDDO ! na1
  nnc = 2
  DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
  DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
  DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     DO na1 = 1, nat
        nt1 = ityp(na1)
        IF ( ehub_lm(nt1) == 0 ) CYCLE
     DO na2 = 1, nat
        nt2 = ityp(na2)
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
        WRITE ( stdout,'(4(a,i2),a)' ) &
              '     [ ', na1, ', ', na2, ', ', ehub_lm(nt1), ', ', ehub_lm(nt2), '  ]'
        WRITE ( stdout,'(a,i2,2i6,a)' ) '     [ ', n1, n2, n3, '  ]'
        DO m1 = 1, ehub_lm(nt1)
           WRITE ( stdout,'(4x,16f16.10)' ) (v(na1,na2,m1,m2,nnc) * rytoev, m2 = 1, ehub_lm(nt2))
        ENDDO
        WRITE ( stdout,* ) ''
     ENDDO ! na2
     ENDDO ! na1
     IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
  ENDDO ! n3
  ENDDO ! n2
  ENDDO ! n1
  END SUBROUTINE rawprint_v

  END SUBROUTINE write_ehub_uv

  !-----------------------------------------------------------------------
  ! [hubdat]  Stock (Lee-Son) counterpart of the dev hp.x file eACBN0_hub_dev.dat
  !-----------------------------------------------------------------------
  !
  ! Every call of the collinear ACBN0 kernel (new_ns_a0 for acbn0_type = 1,
  ! new_ns_e0_1 for acbn0_type = 2) appends one block to two files in outdir:
  !
  !   eACBN0_hub_stock_w_mixing.dat   U, V AFTER mixing (ehub_u / ehub_v),
  !                                   i.e. the values that enter the Hamiltonian
  !   eACBN0_hub_stock_wo_mixing.dat  one-shot U, V of the same call, no mixing
  !                                   (ehub_u_raw / ehub_v_raw)
  !
  ! The block body copies hp_acbn0_write / hp_acbn0_write_v of the dev hp.x
  ! (qe-7.3_thiswork/HP/src/eACBN0.f90), so every NUMERIC line of the last
  ! block of the wo_mixing file matches eACBN0_hub_dev.dat written by hp.x on
  ! the same outdir; the header and the energy prose differ on purpose.  The
  ! energies use the loops of v_hubbard_a0 / v_hubbard_v with the U, V of the
  ! block and the output occupations ns of the same call.  Blocks are written
  ! only for lacbn0 runs, only from inside electrons_scf (ehub_scf_iter > 0),
  ! and only when write_scf writes prefix.save (io_level > -2), so the last
  ! block always matches ehub_uv.txt / ehub_uv_raw.txt.  The first block of a
  ! run replaces the files, later blocks are appended; only ionode writes.
  ! acbn0_type = 3 and noncollinear runs are not covered (dev hp.x scope).
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_ehub_dat_u( ns )
    !-----------------------------------------------------------------------
    USE control_flags, ONLY : io_level
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ns(:,:,:,:)      ! (ehub_lm_max,ehub_lm_max,nspin,nat)
    LOGICAL              :: replace
    !
    IF ( .NOT. lacbn0 .OR. ehub_scf_iter <= 0 .OR. io_level < -1 ) RETURN
    IF ( .NOT. ALLOCATED(ehub_u) .OR. .NOT. ALLOCATED(ehub_u_raw) ) RETURN
    CALL ehub_dat_next( replace )
    CALL ehub_dat_block_u( 'eACBN0_hub_stock_w_mixing.dat',  .TRUE.,  ehub_u,     ns, replace )
    CALL ehub_dat_block_u( 'eACBN0_hub_stock_wo_mixing.dat', .FALSE., ehub_u_raw, ns, replace )
    !
  END SUBROUTINE write_ehub_dat_u
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_ehub_dat_v( ns, updated )
    !-----------------------------------------------------------------------
    USE control_flags, ONLY : io_level
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ns(:,:,:,:,:,:)  ! (ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
    LOGICAL,  INTENT(IN) :: updated          ! .FALSE. when stop_acbn0 skipped the update
    LOGICAL              :: replace
    !
    IF ( .NOT. lacbn0 .OR. ehub_scf_iter <= 0 .OR. io_level < -1 ) RETURN
    IF ( .NOT. ALLOCATED(ehub_v) .OR. .NOT. ALLOCATED(ehub_v_raw) ) RETURN
    CALL ehub_dat_next( replace )
    CALL ehub_dat_block_v( 'eACBN0_hub_stock_w_mixing.dat',  .TRUE.,  ehub_v,     ns, replace, updated )
    CALL ehub_dat_block_v( 'eACBN0_hub_stock_wo_mixing.dat', .FALSE., ehub_v_raw, ns, replace, updated )
    !
  END SUBROUTINE write_ehub_dat_v
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ehub_dat_next( replace )
    !-----------------------------------------------------------------------
    ! Advance the (ionic step, SCF iteration, pass) label of the next block.
    ! Runs on every rank so the bookkeeping never depends on ionode.
    USE control_flags, ONLY : istep
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: replace
    !
    replace = ( ehub_dat_nblk == 0 )
    IF ( ehub_scf_iter == ehub_dat_iter .AND. istep == ehub_dat_step ) THEN
       ehub_dat_pass = ehub_dat_pass + 1
    ELSE
       ehub_dat_pass = 1
    ENDIF
    ehub_dat_iter = ehub_scf_iter
    ehub_dat_step = istep
    ehub_dat_nblk = ehub_dat_nblk + 1
    !
  END SUBROUTINE ehub_dat_next
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ehub_dat_open( fname, replace, lmixed, updated, iunit )
    !-----------------------------------------------------------------------
    ! Open fname in outdir (replace on the first block of the run, append
    ! otherwise), write the file header once, then the block separator.
    ! iunit = -1 if the file cannot be opened: the block is skipped, pw.x runs on.
    USE io_files,      ONLY : tmp_dir, prefix
    USE control_flags, ONLY : lmd, lbfgs
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: fname
    LOGICAL,          INTENT(IN)  :: replace, lmixed, updated
    INTEGER,          INTENT(OUT) :: iunit
    INTEGER                       :: ios
    INTEGER, EXTERNAL             :: find_free_unit
    !
    ! find_free_unit returns a positive unit, so iunit = -1 below is an
    ! unambiguous failure flag (NEWUNIT= would return a negative unit).
    iunit = find_free_unit()
    IF ( replace ) THEN
       OPEN( UNIT = iunit, FILE = TRIM(tmp_dir)//fname, FORM = 'formatted', &
             STATUS = 'replace', IOSTAT = ios )
    ELSE
       OPEN( UNIT = iunit, FILE = TRIM(tmp_dir)//fname, FORM = 'formatted', &
             STATUS = 'unknown', POSITION = 'append', IOSTAT = ios )
    ENDIF
    IF ( ios /= 0 ) THEN
       CALL infomsg( 'ehub_dat_open', 'cannot open '//TRIM(tmp_dir)//fname//', block skipped' )
       iunit = -1
       RETURN
    ENDIF
    !
    IF ( replace ) THEN
       WRITE( iunit, '(/2x,79("="))')
       IF ( acbn0_type == 1 ) THEN
          WRITE( iunit, '(/25x,"Hubbard parameters from ACBN0")')
       ELSE
          WRITE( iunit, '(/24x,"Hubbard parameters from eACBN0")')
       ENDIF
       WRITE( iunit, '(/13x,"computed by pw.x, stock QE (Lee-Son)  (acbn0_type = ",i1,")")') acbn0_type
       WRITE( iunit, '(/2x,79("="),/)')
       WRITE( iunit, '(2x,"prefix     : ",a)') TRIM(prefix)
       WRITE( iunit, '(2x,"projectors : ortho-atomic")')
       WRITE( iunit, '(2x,"bookkeeping: MODULE ldaV (Lee-Son patch)")')
       IF ( lmixed ) THEN
          WRITE( iunit, '(2x,"mixing     : ",f5.3,"  U_n = (1-m)*U_{n-1} + m*U_raw;  U, V AFTER mixing")') &
                ehub_mixing
       ELSE
          WRITE( iunit, '(2x,"mixing     : none   one-shot U_raw, V_raw;  ehub_mixing = ",f5.3," NOT applied")') &
                ehub_mixing
       ENDIF
       IF ( acbn0_type == 2 ) &
          WRITE( iunit, '(2x,"neighbours : cutoff d_IJ = ",f6.3," A,  nnnc = ",i4)') ehub_nn_distance, ehub_nnnc
       IF ( acbn0_type == 1 ) THEN
          WRITE( iunit, '(2x,"estimator  : pseudo-hybrid ACBN0 ratio, NOT linear response")')
       ELSE
          WRITE( iunit, '(2x,"estimator  : pseudo-hybrid eACBN0 ratio, NOT linear response")')
       ENDIF
       WRITE( iunit, '(2x,"blocks     : one per call of the ACBN0 kernel, i.e. one per SCF iteration.")')
       WRITE( iunit, '(2x,"             A 2nd pass inside one iteration is the re-diagonalisation with")')
       WRITE( iunit, '(2x,"             lowered ethr (iteration 1) or the ehub_pot_fix re-evaluation")')
       WRITE( iunit, '(2x,"             after convergence.  Occupations = output ns of that call.")')
       IF ( lmd .OR. lbfgs ) &
          WRITE( iunit, '(2x,"             ionic step = istep of pw.x (0 = starting geometry).")')
       IF ( lmixed ) THEN
          WRITE( iunit, '(2x,"             The LAST block holds the values of prefix.save/ehub_uv.txt.")')
       ELSE
          WRITE( iunit, '(2x,"             The LAST block holds the values of prefix.save/ehub_uv_raw.txt.")')
       ENDIF
    ENDIF
    !
    WRITE( iunit, '(/2x,79("#"))')
    IF ( lmd .OR. lbfgs .OR. ehub_dat_step > 0 ) THEN
       WRITE( iunit, '(2x,"#  ionic step ",i4,"   SCF iteration ",i4,"   pass ",i2,"   (block ",i5,")")') &
             ehub_dat_step, ehub_dat_iter, ehub_dat_pass, ehub_dat_nblk
    ELSE
       WRITE( iunit, '(2x,"#  SCF iteration ",i4,"   pass ",i2,"   (block ",i5,")")') &
             ehub_dat_iter, ehub_dat_pass, ehub_dat_nblk
    ENDIF
    IF ( .NOT. updated ) &
       WRITE( iunit, '(2x,"#  stop_acbn0: U, V NOT recomputed in this call; values repeat the last update")')
    WRITE( iunit, '(2x,79("#"))')
    !
  END SUBROUTINE ehub_dat_open
  !
  !-----------------------------------------------------------------------
  FUNCTION ehub_dat_chanlab( nt, c ) RESULT( lab )
    !-----------------------------------------------------------------------
    ! Label ('4s', '3d', ...) of channel c (= l+1) of type nt.  The principal
    ! quantum number is the first character of the PP_CHI label of the
    ! wavefunction init_lda_plus_v selected for that l (collinear case:
    ! ehub_l_choice(nt,nw) = 1 marks the selected wavefunction nw).
    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: nt, c
    CHARACTER(LEN=2)            :: lab
    CHARACTER(LEN=1), PARAMETER :: spdf(4) = ['s','p','d','f']
    INTEGER                     :: nw
    !
    lab = '??'
    IF ( c < 1 .OR. c > 4 ) RETURN
    lab = ' ' // spdf(c)
    DO nw = 1, MIN( upf(nt)%nwfc, SIZE(ehub_l_choice,2) )
       IF ( ehub_l_choice(nt,nw) /= 1 )   CYCLE
       IF ( upf(nt)%lchi(nw) /= c-1 )     CYCLE
       IF ( upf(nt)%oc(nw) < 0.0_DP )     CYCLE
       IF ( INDEX('123456789', upf(nt)%els(nw)(1:1)) > 0 ) lab = upf(nt)%els(nw)(1:1) // spdf(c)
       RETURN
    ENDDO
    !
  END FUNCTION ehub_dat_chanlab
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ehub_dat_block_u( fname, lmixed, hub_u, ns, replace )
    !-----------------------------------------------------------------------
    ! One acbn0_type = 1 block.  Body = hp_acbn0_write of the dev hp.x.
    USE io_global, ONLY : ionode
    USE ions_base, ONLY : atm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL,          INTENT(IN) :: lmixed, replace
    REAL(DP),         INTENT(IN) :: hub_u(:,:,:), ns(:,:,:,:)
    INTEGER                      :: na, nt, m1, m2, is, iunit, nspin
    REAL(DP)                     :: occ, e_on
    !
    IF ( .NOT. ionode ) RETURN
    nspin = SIZE(ns,3)
    CALL ehub_dat_open( fname, replace, lmixed, .TRUE., iunit )
    IF ( iunit < 0 ) RETURN
    !
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"Summary")')
    WRITE( iunit, '(2x,79("-"))')
    WRITE( iunit, '(/4x,"site n.  type  label  channel   occupation   U_ACBN0 (eV)")')
    WRITE( iunit, '( 4x,"(one row per channel; occupation and U of that channel only)")')
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       DO m1 = 1, ehub_lm(nt)
          IF ( ANY( ehub_l(1:m1-1,nt) == ehub_l(m1,nt) ) ) CYCLE     ! first m of this channel
          occ = 0.0_DP
          DO is = 1, nspin
             DO m2 = 1, ehub_lm(nt)
                IF ( ehub_l(m2,nt) == ehub_l(m1,nt) ) occ = occ + ns(m2,m2,is,na)
             ENDDO
          ENDDO
          WRITE( iunit, '(6x,i3,4x,i3,3x,a6,4x,a2,6x,f9.4,4x,f10.4)') &
                na, nt, TRIM(atm(nt)), ehub_dat_chanlab(nt,ehub_l(m1,nt)), occ, hub_u(na,m1,m1)*rytoev
       ENDDO
    ENDDO
    !
    ! ... energy: the v_hubbard_a0 loop with this block's U and this call's ns
    e_on = 0.0_DP
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       DO is = 1, nspin
          DO m1 = 1, ehub_lm(nt)
             e_on = e_on + 0.5_DP*hub_u(na,m1,m1)*ns(m1,m1,is,na)
             DO m2 = 1, ehub_lm(nt)
                e_on = e_on - 0.5_DP*hub_u(na,m1,m2)*ns(m1,m2,is,na)*ns(m2,m1,is,na)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    IF ( nspin == 1 ) e_on = 2.0_DP*e_on
    !
    WRITE( iunit, '(/2x,79("-"))')
    IF ( lmixed ) THEN
       WRITE( iunit, '(2x,"Hubbard energy with the MIXED U")')
    ELSE
       WRITE( iunit, '(2x,"One-shot Hubbard energy  (NOT self-consistent)")')
    ENDIF
    WRITE( iunit, '(2x,79("-"))')
    WRITE( iunit, '(/4x,"E(+U) : ",f16.9," Ry  = ",f12.6," eV")') e_on, e_on*rytoev
    IF ( lmixed ) THEN
       WRITE( iunit, '( 4x,"Evaluated with the mixed U above and the output occupations of this call")')
       WRITE( iunit, '( 4x,"(v_hubbard_a0 in v_hubbard_v.f90).  In the converged iteration this is the")')
       WRITE( iunit, '( 4x,"E(+U) pw.x prints (not with ehub_pot_fix); earlier ones use the mixed ns.")')
    ELSE
       WRITE( iunit, '( 4x,"Evaluated with the one-shot U above and the occupations it was computed")')
       WRITE( iunit, '( 4x,"from (v_hubbard_a0 in v_hubbard_v.f90); not a total-energy contribution.")')
    ENDIF
    !
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"Full U matrix  U[na, m1, m2]  (eV)")')
    WRITE( iunit, '(2x,79("-"))')
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"   (",i2," x",i2,")")') &
             na, TRIM(atm(nt)), ehub_lm(nt), ehub_lm(nt)
       DO m1 = 1, ehub_lm(nt)
          WRITE( iunit, '(6x,16f10.4)') (hub_u(na,m1,m2)*rytoev, m2 = 1, ehub_lm(nt))
       ENDDO
    ENDDO
    !
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"Occupation matrix  ns[na, m1, m2, spin]")')
    WRITE( iunit, '(2x,79("-"))')
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       DO is = 1, nspin
          WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"  spin ",i1)') na, TRIM(atm(nt)), is
          DO m1 = 1, ehub_lm(nt)
             WRITE( iunit, '(6x,16f10.4)') (ns(m1,m2,is,na), m2 = 1, ehub_lm(nt))
          ENDDO
       ENDDO
    ENDDO
    !
    WRITE( iunit, '(/2x,79("="),/)')
    CLOSE( iunit )
    !
  END SUBROUTINE ehub_dat_block_u
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ehub_dat_block_v( fname, lmixed, hub_v, ns, replace, updated )
    !-----------------------------------------------------------------------
    ! One acbn0_type = 2 block.  Body = hp_acbn0_write_v of the dev hp.x; the
    ! cell -> nnc map is the counting used by write_ehub_uv and v_hubbard_v.
    USE constants, ONLY : bohr_radius_angs
    USE cell_base, ONLY : alat, at
    USE io_global, ONLY : ionode
    USE ions_base, ONLY : atm, tau
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL,          INTENT(IN) :: lmixed, replace, updated
    REAL(DP),         INTENT(IN) :: hub_v(:,:,:,:,:), ns(:,:,:,:,:,:)
    INTEGER                      :: na, na1, na2, nt, nt1, nt2, m1, m2, is, iunit
    INTEGER                      :: npair, ics, n1, n2, n3, nnc, thisnnc, ix, nspin
    REAL(DP)                     :: occ, e_on, e_int, dtau(3), r(3), r2
    REAL(DP), PARAMETER          :: veps = 1.0D-8   ! Ry; below this a pair carries no V
    !
    IF ( .NOT. ionode ) RETURN
    CALL ehub_dat_open( fname, replace, lmixed, updated, iunit )
    IF ( iunit < 0 ) RETURN
    ics   = ehub_csize + 1
    nspin = SIZE(ns,3)
    !
    ! ---- on-site U ------------------------------------------------------
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"On-site U   = hub_v[na, na, m1, m2, cell 0]")')
    WRITE( iunit, '(2x,79("-"))')
    WRITE( iunit, '(/4x,"site n.  type  label  channel   occupation    U (eV)")')
    WRITE( iunit, '( 4x,"(one row per channel; occupation and U of that channel only)")')
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       DO m1 = 1, ehub_lm(nt)
          IF ( ANY( ehub_l(1:m1-1,nt) == ehub_l(m1,nt) ) ) CYCLE
          occ = 0.0_DP
          DO is = 1, nspin
             DO m2 = 1, ehub_lm(nt)
                IF ( ehub_l(m2,nt) == ehub_l(m1,nt) ) occ = occ + ns(m2,m2,is,na,na,1)
             ENDDO
          ENDDO
          WRITE( iunit, '(6x,i3,4x,i3,3x,a6,4x,a2,6x,f9.4,4x,f10.4)') &
                na, nt, TRIM(atm(nt)), ehub_dat_chanlab(nt,ehub_l(m1,nt)), occ, hub_v(na,na,m1,m1,1)*rytoev
       ENDDO
    ENDDO
    !
    ! ---- inter-site V, with the cell vector and distance -----------------
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"Inter-site V")')
    WRITE( iunit, '(2x,79("-"))')
    WRITE( iunit, '(/4x,"na1  type1     na2  type2      cell [n1 n2 n3]", &
          &"    d (A)    block      V (eV)")')
    WRITE( iunit, '( 4x,"(one row per (channel of na1, channel of na2) block; ", &
          &"V is constant inside a block)")')
    npair = 0
    nnc   = 2
    DO n1 = -ehub_nn_dim(1), ehub_nn_dim(1)
    DO n2 = -ehub_nn_dim(2), ehub_nn_dim(2)
    DO n3 = -ehub_nn_dim(3), ehub_nn_dim(3)
       IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) THEN
          thisnnc = 1
       ELSE
          thisnnc = nnc
       ENDIF
       DO na1 = 1, nat
          nt1 = ityp(na1)
          IF ( ehub_lm(nt1) == 0 ) CYCLE
          DO na2 = 1, nat
             nt2 = ityp(na2)
             IF ( ehub_lm(nt2) == 0 ) CYCLE
             IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
             IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 .AND. na1 == na2 ) CYCLE
             IF ( MAXVAL(ABS(hub_v(na1,na2,1:ehub_lm(nt1),1:ehub_lm(nt2),thisnnc))) < veps ) CYCLE
             npair = npair + 1
             ! distance at the CURRENT geometry, same expression as hubbard_nns.f90
             ! (ehub_nn_leng keeps the setup geometry during relax / md)
             dtau(:) = tau(:,na1) - tau(:,na2)
             r2 = 0.0_DP
             DO ix = 1, 3
                r(ix) = n1 * at(ix,1) + n2 * at(ix,2) + n3 * at(ix,3) - dtau(ix)
                r2 = r2 + r(ix) * r(ix)
             ENDDO
             DO m1 = 1, ehub_lm(nt1)
                IF ( ANY( ehub_l(1:m1-1,nt1) == ehub_l(m1,nt1) ) ) CYCLE
                DO m2 = 1, ehub_lm(nt2)
                   IF ( ANY( ehub_l(1:m2-1,nt2) == ehub_l(m2,nt2) ) ) CYCLE
                   WRITE( iunit, '(4x,i3,2x,a6,4x,i3,2x,a6,4x,"[",3i3," ]",4x,f7.4,3x,a2,"-",a2,4x,f10.4)') &
                         na1, TRIM(atm(nt1)), na2, TRIM(atm(nt2)), n1, n2, n3, &
                         SQRT(r2)*alat*bohr_radius_angs, &
                         ehub_dat_chanlab(nt1,ehub_l(m1,nt1)), ehub_dat_chanlab(nt2,ehub_l(m2,nt2)), &
                         hub_v(na1,na2,m1,m2,thisnnc)*rytoev
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF ( .NOT.( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) ) THEN
          IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    IF ( npair == 0 ) WRITE( iunit, '(/4x,"(none - check ehub_nn_distance)")')
    WRITE( iunit, '(/4x,i5," pair(s) with |V| > ",es8.1," Ry")') npair, veps
    !
    ! ---- energy: the v_hubbard_v loops with this block's U, V and this call's ns
    e_on = 0.0_DP ; e_int = 0.0_DP
    DO na1 = 1, nat
       nt1 = ityp(na1)
       IF ( ehub_lm(nt1) == 0 ) CYCLE
       DO is = 1, nspin
          DO m1 = 1, ehub_lm(nt1)
             e_on = e_on + 0.5_DP*hub_v(na1,na1,m1,m1,1)*ns(m1,m1,is,na1,na1,1)
             DO m2 = 1, ehub_lm(nt1)
                e_on = e_on - 0.5_DP*hub_v(na1,na1,m1,m2,1) &
                            * ns(m1,m2,is,na1,na1,1)*ns(m2,m1,is,na1,na1,1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DO na1 = 1, nat
       nt1 = ityp(na1)
       IF ( ehub_lm(nt1) == 0 ) CYCLE
       DO na2 = 1, nat
          nt2 = ityp(na2)
          IF ( ehub_lm(nt2) == 0 ) CYCLE
          IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
          DO is = 1, nspin
             DO m1 = 1, ehub_lm(nt1)
                DO m2 = 1, ehub_lm(nt2)
                   e_int = e_int - 0.5_DP*hub_v(na1,na2,m1,m2,1) &
                                 * ns(m1,m2,is,na1,na2,1)*ns(m2,m1,is,na2,na1,1)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    nnc = 2
    DO n1 = -ehub_nn_dim(1), ehub_nn_dim(1)
    DO n2 = -ehub_nn_dim(2), ehub_nn_dim(2)
    DO n3 = -ehub_nn_dim(3), ehub_nn_dim(3)
       IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
       DO na1 = 1, nat
          nt1 = ityp(na1)
          IF ( ehub_lm(nt1) == 0 ) CYCLE
          DO na2 = 1, nat
             nt2 = ityp(na2)
             IF ( ehub_lm(nt2) == 0 ) CYCLE
             IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
             DO is = 1, nspin
                DO m1 = 1, ehub_lm(nt1)
                   DO m2 = 1, ehub_lm(nt2)
                      ! NOTE: ns* = ns  (v_hubbard_v.f90)
                      e_int = e_int - 0.5_DP*hub_v(na1,na2,m1,m2,nnc) &
                                    * ns(m1,m2,is,na1,na2,nnc)*ns(m1,m2,is,na1,na2,nnc)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
    ENDDO
    ENDDO
    ENDDO
    IF ( nspin == 1 ) THEN
       e_on = 2.0_DP*e_on ; e_int = 2.0_DP*e_int
    ENDIF
    !
    WRITE( iunit, '(/2x,79("-"))')
    IF ( lmixed ) THEN
       WRITE( iunit, '(2x,"Extended Hubbard energy with the MIXED U, V")')
    ELSE
       WRITE( iunit, '(2x,"One-shot extended Hubbard energy  (NOT self-consistent, see below)")')
    ENDIF
    WRITE( iunit, '(2x,79("-"))')
    WRITE( iunit, '(/4x,"on-site   (U)  : ",f16.9," Ry  = ",f12.6," eV")') e_on,        e_on*rytoev
    WRITE( iunit, '( 4x,"inter-site (V) : ",f16.9," Ry  = ",f12.6," eV")') e_int,       e_int*rytoev
    WRITE( iunit, '( 4x,"E(+U+V)        : ",f16.9," Ry  = ",f12.6," eV")') e_on+e_int, (e_on+e_int)*rytoev
    IF ( lmixed ) THEN
       WRITE( iunit, '(/4x,"Evaluated with the mixed U, V above and the output occupations of this")')
       WRITE( iunit, '( 4x,"call (v_hubbard_v loops).  In the converged iteration this is the E(+U+V)")')
       WRITE( iunit, '( 4x,"pw.x prints (not with ehub_pot_fix); earlier iterations use the mixed ns.")')
    ELSEIF ( updated ) THEN
       WRITE( iunit, '(/4x,"Evaluated with the one-shot U, V above and the occupations they were")')
       WRITE( iunit, '( 4x,"computed from (v_hubbard_v loops); U, V and the density are not")')
       WRITE( iunit, '( 4x,"mutually self-consistent, so this is not a total-energy contribution.")')
    ELSE
       WRITE( iunit, '(/4x,"stop_acbn0 skip: U, V are those of the last update, the occupations")')
       WRITE( iunit, '( 4x,"those of this call (v_hubbard_v loops); not a consistent pair and not")')
       WRITE( iunit, '( 4x,"a total-energy contribution.")')
    ENDIF
    !
    ! ---- full on-site matrices ------------------------------------------
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"Full on-site U matrix  (eV)")')
    WRITE( iunit, '(2x,79("-"))')
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"   (",i2," x",i2,")")') &
             na, TRIM(atm(nt)), ehub_lm(nt), ehub_lm(nt)
       DO m1 = 1, ehub_lm(nt)
          WRITE( iunit, '(6x,16f10.4)') (hub_v(na,na,m1,m2,1)*rytoev, m2 = 1, ehub_lm(nt))
       ENDDO
    ENDDO
    !
    ! ---- on-site occupation matrices ------------------------------------
    WRITE( iunit, '(/2x,79("-"))')
    WRITE( iunit, '(2x,"On-site occupation matrix  ns[na, m1, m2, spin]")')
    WRITE( iunit, '(2x,79("-"))')
    DO na = 1, nat
       nt = ityp(na)
       IF ( ehub_lm(nt) == 0 ) CYCLE
       DO is = 1, nspin
          WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"  spin ",i1)') na, TRIM(atm(nt)), is
          DO m1 = 1, ehub_lm(nt)
             WRITE( iunit, '(6x,16f10.4)') (ns(m1,m2,is,na,na,1), m2 = 1, ehub_lm(nt))
          ENDDO
       ENDDO
    ENDDO
    !
    WRITE( iunit, '(/2x,79("="),/)')
    CLOSE( iunit )
    !
  END SUBROUTINE ehub_dat_block_v
  
END MODULE ldaV
