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

  END SUBROUTINE write_ehub_uv
  
END MODULE ldaV
