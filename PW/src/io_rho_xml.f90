!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE io_rho_xml
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE io_files,    ONLY : create_directory
  USE io_base,     ONLY : write_rhog, read_rhog
  !
  PRIVATE
  ! by LSH WIY
  PUBLIC :: write_scf, read_scf, read_ehub_uv, read_ehub_ns
  ! by LSH WIY
  !
  ! {read|write}_rho: read or write the charge density
  ! {read|write}_scf: as above, plus ldaU ns, PAW becsum, meta-GGA
  !
  CONTAINS

    SUBROUTINE write_scf ( rho, nspin )
      !
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, hub_back, lda_plus_u_kind, nsg
      USE two_chem,         ONLY : twochem
      ! by LSH WIY
      USE ldaV,             ONLY : lda_plus_v, lacbn0, acbn0_type, &
                                   ehub_u, ehub_v, ehub_lm, ehub_l, ehub_nn_vect, &
                                   ehub_nn_dim, ehub_csize, ehub_lm_max, ehub_nnnc
      ! by LSH WIY
      USE xc_lib,           ONLY : xclib_dft_is
      USE noncollin_module, ONLY : noncolin, domag
      USE scf,              ONLY : scf_type
      !
      USE cell_base,        ONLY : bg, tpiba
      USE gvect,            ONLY : ig_l2g, mill
      USE control_flags,    ONLY : gamma_only
      USE io_files,         ONLY : restart_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_pools,         ONLY : my_pool_id
      USE mp_bands,         ONLY : my_bgrp_id, root_bgrp_id, &
                                   root_bgrp, intra_bgrp_comm
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast
      ! by KSJ
      USE constants,  ONLY : rytoev
      USE kinds,      ONLY : DP
      USE parameters, ONLY : ntypx
      USE ions_base,  ONLY : nat, ntyp => nsp, ityp
      USE uspp_param,  ONLY : upf
      ! by KSJ

      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(IN)           :: rho
      INTEGER,          INTENT(IN)           :: nspin
      !
      CHARACTER (LEN=256) :: dirname
      INTEGER :: nspin_, iunocc, iunpaw, iuntwochem, ierr
      ! by KSJ LSH WIY
      INTEGER :: na, na1, na2
      INTEGER :: nt, nt1, nt2
      INTEGER :: m1, m2
      INTEGER :: n1, n2, n3, nnc, ics
      INTEGER :: iunehubu, iunehubv, iunoccv, iunoccuv
      INTEGER :: lm1, lm2, spin_idx
      ! by KSJ LSH WIY
      INTEGER, EXTERNAL :: find_free_unit

      dirname = restart_dir ( )
      CALL create_directory( dirname )
      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         nspin_ = 1
      ELSE
         nspin_ = nspin
      ENDIF
      ! Write G-space density
      IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
           CALL write_rhog( TRIM(dirname) // "charge-density", &
           root_bgrp, intra_bgrp_comm, &
           bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
           gamma_only, mill, ig_l2g, rho%of_g(:,1:nspin_) )
      !
      ! Write kinetic energy density density
      IF ( xclib_dft_is('meta') ) THEN 
         IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
              CALL write_rhog( TRIM(dirname) // "ekin-density", &
              root_bgrp, intra_bgrp_comm, &
              bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
              gamma_only, mill, ig_l2g, rho%kin_g(:,1:nspin_) )
         WRITE(stdout,'(5x,"Writing meta-gga kinetic term")')
      ENDIF

      ! Then write the other terms to separate files

      IF ( lda_plus_u ) THEN
         !
         IF ( ionode ) THEN
            OPEN ( NEWUNIT=iunocc, FILE = TRIM(dirname) // 'occup.txt', &
                 FORM='formatted', STATUS='unknown' )
            IF (lda_plus_u_kind.EQ.0) THEN
               IF (noncolin) THEN
                   WRITE( iunocc, * , iostat = ierr) rho%ns_nc
               ELSE   
                  WRITE( iunocc, * , iostat = ierr) rho%ns
               ENDIF
               IF (hub_back) WRITE( iunocc, * , iostat = ierr) rho%nsb
            ELSEIF (lda_plus_u_kind.EQ.1) THEN
               IF (noncolin) THEN
                  WRITE( iunocc, * , iostat = ierr) rho%ns_nc
               ELSE
                  WRITE( iunocc, * , iostat = ierr) rho%ns
               ENDIF
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               WRITE( iunocc, * , iostat = ierr) nsg
            ENDIF
         ENDIF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing ldaU ns', 1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
         ENDIF
         !
      END IF
      !
      ! by LSH
      IF ( lda_plus_v ) THEN
         ! by KSJ
         iunoccv = find_free_unit ()
         iunoccuv = iunoccv+1
!         WRITE( stdout, *) iunoccv, iunoccuv
         ! by KSJ
         IF ( ionode ) THEN
            !
            OPEN ( NEWUNIT = iunoccv, FILE = TRIM(dirname)//'occupation_uv.txt', FORM = 'FORMATTED', STATUS = 'unknown')
            ! by KSJ
            OPEN ( NEWUNIT = iunoccuv, FILE = TRIM(dirname)//'occupation_uv0.txt', FORM = 'FORMATTED', STATUS = 'unknown')
            ! by KSJ
            !
            IF ( noncolin ) THEN
               !
               IF ( lacbn0 ) THEN
                  IF ( acbn0_type == 1 ) THEN
                     WRITE ( iunoccv, *, IOSTAT = ierr ) rho%na0_nc
                  ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                      ! by KSJ
!                     WRITE ( stdout,*) '=================================='
!                     WRITE ( stdout,*) '     occupation_uv information    '
!                     WRITE ( stdout,*) SIZE(rho%ne0_nc), SHAPE(rho%ne0_nc)
!                     WRITE ( stdout,*) '=================================='             
                      ! by KSJ
                      WRITE ( iunoccv, *, IOSTAT = ierr ) rho%ne0_nc
                      ! by KSJ
                      WRITE ( iunoccuv,'(a)' ) '# ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc'
                      WRITE ( iunoccuv, '(A, I3, 5I4)') '#', ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc  
                      DO lm1 = 1, ehub_lm_max
                         DO lm2 = 1, ehub_lm_max
                            DO spin_idx = 1, nspin
                               DO na1 = 1, nat
                                  DO na2 = 1, nat
                                     DO nnc = 1, ehub_nnnc
                                     !
                                     IF (rho%ne0_nc(lm1, lm2, spin_idx, na1, na2, nnc) /= 0.) WRITE ( iunoccuv, '(6i4, 2f25.15)') lm1, lm2, spin_idx, na1, na2, nnc, rho%ne0_nc(lm1, lm2, spin_idx, na1, na2, nnc)
                                     !
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO 
                      ! by KSJ
                  ENDIF
               ELSE
                  WRITE ( iunoccv, *, IOSTAT = ierr ) rho%ne0_nc
                  ! by KSJ
                  WRITE ( iunoccuv,'(a)' ) '# ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc'
                  WRITE ( iunoccuv, '(A, I3, 5I4)') '#', ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc  
                  DO lm1 = 1, ehub_lm_max
                     DO lm2 = 1, ehub_lm_max
                        DO spin_idx = 1, nspin
                           DO na1 = 1, nat
                              DO na2 = 1, nat
                                 DO nnc = 1, ehub_nnnc
                                 !
                                 IF (rho%ne0_nc(lm1, lm2, spin_idx, na1, na2, nnc) /= 0.) WRITE ( iunoccuv, '(6i4, 2f25.15)') lm1, lm2, spin_idx, na1, na2, nnc, rho%ne0_nc(lm1, lm2, spin_idx, na1, na2, nnc)
                                 !
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO 
                  ! by KSJ
               ENDIF
               !
            ELSE
               !
               IF ( lacbn0 ) THEN
                  IF ( acbn0_type == 1 ) THEN
                     WRITE ( iunoccv, *, IOSTAT = ierr ) rho%na0
                  ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                     ! by KSJ
!                     WRITE ( stdout,*) '=================================='
!                     WRITE ( stdout,*) '     occupation_uv information    '
!                     WRITE ( stdout,*) SIZE(rho%ne0), SHAPE(rho%ne0)
!                     WRITE ( stdout,*) '=================================='
                     ! by KSJ
                     !
                     WRITE ( iunoccv, *, IOSTAT = ierr ) rho%ne0
                     !
                     ! by KSJ
                     WRITE ( iunoccuv,'(a)' ) '# ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc'
                     WRITE ( iunoccuv, '(A, I3, 5I4)') '#', ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc  
                     DO lm1 = 1, ehub_lm_max
                        DO lm2 = 1, ehub_lm_max
                           DO spin_idx = 1, nspin
                              DO na1 = 1, nat
                                 DO na2 = 1, nat
                                    DO nnc = 1, ehub_nnnc
                                    !
                                    IF (rho%ne0(lm1, lm2, spin_idx, na1, na2, nnc) /= 0.) WRITE ( iunoccuv, '(6i4, 1f25.15)') lm1, lm2, spin_idx, na1, na2, nnc, rho%ne0(lm1, lm2, spin_idx, na1, na2, nnc)
                                    !
                                    ENDDO
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO 
                     !
                     ! by KSJ
                     !
                  ENDIF
               ELSE
                  WRITE ( iunoccv, *, IOSTAT = ierr ) rho%ne0
                  !
                  ! by KSJ
                  ! rho%ne0(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
!                  WRITE ( stdout,*) '=================================='
!                  WRITE ( stdout,*) '     occupation_uv information    '
!                  WRITE ( stdout,*) SIZE(rho%ne0), SHAPE(rho%ne0)
!                  WRITE ( stdout,*) '=================================='
                  !
                  WRITE ( iunoccuv,'(a)' ) '# ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc'
                  WRITE ( iunoccuv, '(A, I3, 5I4)') '#', ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc
                  DO lm1 = 1, ehub_lm_max
                     DO lm2 = 1, ehub_lm_max
                        DO spin_idx = 1, nspin
                           DO na1 = 1, nat
                              DO na2 = 1, nat
                                 DO nnc = 1, ehub_nnnc
                                 !
                                 IF (rho%ne0(lm1, lm2, spin_idx, na1, na2, nnc) /= 0.) WRITE ( iunoccuv, '(6i4, 1f25.15)') lm1, lm2, spin_idx, na1, na2, nnc, rho%ne0(lm1, lm2, spin_idx, na1, na2, nnc)
                                 !
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  ! by KSJ
               ENDIF
            ENDIF
            !
         ENDIF
         !
         CALL mp_bcast(ierr, ionode_id, intra_image_comm)
         !
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing occupation matrices ns for DFT+U+V calculation', 1)
         !
         IF ( ionode ) THEN
            CLOSE ( UNIT = iunoccv, STATUS = 'KEEP' )
            ! by KSJ
            CLOSE ( UNIT = iunoccuv, STATUS = 'KEEP' )
            ! by KSJ
         ENDIF
         ! by KSJ
         iunoccv = find_free_unit ()
         ! by KSJ
         !
         IF ( ionode ) THEN
            !
            OPEN ( NEWUNIT = iunoccv, FILE = TRIM(dirname)//'ehub_uv.txt', FORM = 'FORMATTED', STATUS = 'unknown')
            !
            ! by KSJ
!            WRITE ( stdout,*) '=================================='
!            WRITE ( stdout,*) '        ehub_uv information       '
!            WRITE ( stdout,*) '=================================='
!            WRITE ( stdout,*) SIZE(ehub_v), SHAPE(ehub_v)
!            WRITE ( stdout,*) '=================================='
            ! by KSJ
            IF ( lacbn0 ) THEN
               IF ( acbn0_type == 1 ) THEN
                  WRITE ( iunoccv, *, IOSTAT = ierr ) ehub_u
               ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                  WRITE ( iunoccv, *, IOSTAT = ierr ) ehub_v
               ENDIF
            ELSE
               WRITE ( iunoccv, *, IOSTAT = ierr ) ehub_v
            ENDIF
            !
         ENDIF
         !
         CALL mp_bcast(ierr, ionode_id, intra_image_comm)
         !
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing Hubbard parameters U or V', 1)
         !
         IF ( ionode ) THEN
            CLOSE( UNIT = iunoccv, STATUS = 'KEEP' )
         ENDIF
         !
         ! by KSJ
         IF ( lacbn0 .AND. acbn0_type == 1 ) THEN

            iunehubu = find_free_unit ()

            IF ( ionode ) THEN
               OPEN (UNIT= iunehubu, FILE = TRIM(dirname)//'ehub_u0.txt', FORM='FORMATTED', STATUS = 'unknown')

               !
               ! Write U matrix
               !
               WRITE ( iunehubu,'(a)' ) '                          '
               WRITE ( iunehubu,'(a)' ) '     ---------------------'
               WRITE ( iunehubu,'(a)' ) '     U[atom i, #orb, #orb]'
               WRITE ( iunehubu,'(a)' ) '     ---------------------'
 
               DO na = 1, nat
                  nt = ityp(na)
                  !
                  IF ( ehub_lm(nt) == 0 ) CYCLE
                  !
                  WRITE ( iunehubu,'(3(a,i2),a)' ) &
                        '     [ ', na, ', ', 2*MAXVAL(ehub_l(:,nt))-1, ', ', 2*MAXVAL(ehub_l(:,nt))-1, '  ]'
                  !
                  DO m1 = 1, ehub_lm(nt)
                     WRITE ( iunehubu,'(4x,16f9.4)' ) &
                           (ehub_u(na,m1,m2) * rytoev, m2 = 1, ehub_lm(nt))
                  ENDDO
                  !
                  WRITE ( iunehubu,* ) ''
                 !
               ENDDO ! na
               !
            ENDIF

            CALL mp_bcast( ierr, ionode_id, intra_image_comm)

            IF (ierr/=0) CALL errore('write_ehub_uv_file', 'Writing Hubbard parameters U', 1)

            IF ( ionode ) THEN
               CLOSE( UNIT=iunehubu, STATUS = 'KEEP')
            ENDIF

         ELSE

            iunehubv = find_free_unit ()

            IF ( ionode ) THEN

               OPEN (UNIT= iunehubv, FILE = TRIM(dirname)//'ehub_uv0.txt', FORM='FORMATTED', STATUS = 'unknown')
               !
               ! Write V matrices
               !
               !WRITE ( iunehubv,'(a)' ) '                                        '
               !WRITE ( iunehubv,'(a)' ) '     -----------------------------------'
               WRITE ( iunehubv,'(a)' ) '     # V[atom i, atom j, #orb(i), #orb(j), l, m, n, nnc]'
               !WRITE ( iunehubv,'(a)' ) '     -----------------------------------'
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
                  !WRITE ( iunehubv,'(4(a,i2),a)' ) &
                  !      '     [ ', na1, ', ', na2, ', ', ehub_lm(nt1), ', ', ehub_lm(nt2),'  ]'
                  !
                  !WRITE ( iunehubv,'(7i4)' ) na1, na2, ehub_lm(nt1), ehub_lm(nt2), 0, 0, 0
                  !
                  !
                  !WRITE ( iunehubv,'(a,i2,2i6,a)' ) '     [ ', 0, 0, 0, '  ]'
                  !
                  DO m1 = 1, ehub_lm(nt1)
                     DO m2 = 1, ehub_lm(nt2)
                        WRITE ( iunehubv,'(8i4, 1f25.15)' ) na1, na2, m1, m2, 0, 0, 0, 1, ehub_v(na1,na2,m1,m2,1)
                        !WRITE ( iunehubv,'(4x,16f9.4)' ) &
                        !      (ehub_v(na1,na2,m1,m2,1) * rytoev, m2 = 1, ehub_lm(nt2))
                     ENDDO
                  ENDDO
                  !
                  !WRITE ( iunehubv,* ) ''
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
                  !WRITE ( iunehubv, '(4(a,i2),a)' ) &
                  !      '     [ ', na1, ', ', na2, ', ', ehub_lm(nt1), ', ', ehub_lm(nt2), '  ]'
                  !
                  !WRITE ( iunehubv, '(a,i2,2i6,a)' ) '     [ ', n1, n2, n3, '  ]'
                  !
                  DO m1 = 1, ehub_lm(nt1)
                     DO m2 = 1, ehub_lm(nt2)
                        WRITE ( iunehubv,'(8i4, 1f25.15)' ) na1, na2, m1, m2, n1, n2, n3, nnc, ehub_v(na1,na2,m1,m2,nnc)
                        !WRITE ( iunehubv,'(4x,16f9.4)' ) &
                        !      (ehub_v(na1,na2,m1,m2,nnc) * rytoev, m2 = 1, ehub_lm(nt2))
                     ENDDO
                  ENDDO
                  !
                  !WRITE ( iunehubv,* ) ''
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

         CALL mp_bcast(ierr, ionode_id, intra_image_comm)

         IF ( ierr/=0 ) CALL errore('write_ehub_uv_file', 'Writing Hubbard parameters U or V', 1)

         IF ( ionode ) THEN 
            CLOSE( UNIT = iunehubv, STATUS = 'KEEP' )
         ENDIF
      ENDIF
         ! by KSJ
      ENDIF
      ! by LSH
      !
      IF ( okpaw ) THEN
         !
         IF ( ionode ) THEN
            OPEN ( NEWUNIT=iunpaw, FILE = TRIM(dirname) // 'paw.txt', &
                 FORM='formatted', STATUS='unknown' )
            WRITE( iunpaw, * , iostat = ierr) rho%bec
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing PAW becsum',1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunpaw, STATUS = 'KEEP' )
         ENDIF
         !
      END IF
      !
      RETURN
    END SUBROUTINE write_scf

    SUBROUTINE read_scf ( rho, nspin, gamma_only )
      !
      USE scf,              ONLY : scf_type
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, starting_ns, hub_back, &
                                   lda_plus_u_kind, nsg
      ! by LSH
      USE ldaV,             ONLY : lda_plus_v, lacbn0, acbn0_type, &
                                   ehub_u, ehub_v
      ! by LSH
      USE noncollin_module, ONLY : noncolin, domag
      USE gvect,            ONLY : ig_l2g
      USE xc_lib,           ONLY : xclib_dft_is
      USE io_files,         ONLY : restart_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_bands,         ONLY : root_bgrp, intra_bgrp_comm
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(INOUT)        :: rho
      INTEGER,          INTENT(IN)           :: nspin
      LOGICAL, OPTIONAL,INTENT(IN)           :: gamma_only
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL :: lexist
      INTEGER :: nspin_, iunocc, iunpaw, ierr
      ! by LSH
      INTEGER :: iunoccv
      ! by LSH
      INTEGER, EXTERNAL :: find_free_unit

      dirname = restart_dir ( )
      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         nspin_=1
      ELSE
         nspin_=nspin
      ENDIF
      ! read charge density
      CALL read_rhog( TRIM(dirname) // "charge-density", &
           root_bgrp, intra_bgrp_comm, &
           ig_l2g, nspin_, rho%of_g, gamma_only )
      IF ( nspin > nspin_) rho%of_g(:,nspin_+1:nspin) = (0.0_dp, 0.0_dp)
      !
      ! read kinetic energy density
      IF ( xclib_dft_is('meta') ) THEN
         CALL read_rhog( TRIM(dirname) // "ekin-density", &
              root_bgrp, intra_bgrp_comm, &
              ig_l2g, nspin_, rho%kin_g, gamma_only, ierr )
         IF ( ierr == 0 ) THEN
            WRITE(stdout,'(5x,"Reading meta-gga kinetic term")')
         ELSE
            rho%kin_g(:,:) = (0.0_dp, 0.0_dp)
            WRITE(stdout,'(5x,"BEWARE: kinetic-energy density file not found,",&
                    & " Kinetic-energy density set to 0")')
         ENDIF
      END IF

      IF ( lda_plus_u ) THEN
         !
         ! The occupations ns also need to be read in order to build up
         ! the potential
         !
         IF ( ionode ) THEN
            OPEN ( NEWUNIT=iunocc, FILE = TRIM(dirname) // 'occup.txt', &
                 FORM='formatted', STATUS='old', IOSTAT=ierr )
            IF (lda_plus_u_kind.EQ.0) THEN
               IF (noncolin) THEN
                  READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns_nc
               ELSE        
                  READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
               ENDIF   
               IF (hub_back) READ( UNIT = iunocc, FMT = * , iostat = ierr) rho%nsb
            ELSEIF (lda_plus_u_kind.EQ.1) THEN
               IF (noncolin) THEN
                  READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns_nc
               ELSE
                  READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
               ENDIF
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               READ( UNIT = iunocc, FMT = * , iostat = ierr) nsg 
            ENDIF
         ENDIF
         !
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('read_scf', 'Reading ldaU ns', 1)
         !
         IF ( ionode ) THEN
            CLOSE( UNIT = iunocc, STATUS = 'KEEP')
         ELSE
            IF (lda_plus_u_kind.EQ.0) THEN
               IF (noncolin) THEN
                  rho%ns_nc(:,:,:,:) = 0.D0      
               ELSE        
                  rho%ns(:,:,:,:) = 0.D0
               ENDIF
               IF (hub_back) rho%nsb(:,:,:,:) = 0.D0
            ELSEIF (lda_plus_u_kind.EQ.1) THEN
               IF (noncolin) THEN
                  rho%ns_nc(:,:,:,:) = 0.D0
               ELSE
                  rho%ns(:,:,:,:) = 0.D0
               ENDIF 
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               nsg(:,:,:,:,:) = (0.d0, 0.d0) 
            ENDIF
         ENDIF
         !
         IF (lda_plus_u_kind.EQ.0) THEN
            IF (noncolin) THEN
                CALL mp_sum(rho%ns_nc, intra_image_comm)
            ELSE    
                CALL mp_sum(rho%ns, intra_image_comm)
            ENDIF    
            IF (hub_back) CALL mp_sum(rho%nsb, intra_image_comm)   
         ELSEIF (lda_plus_u_kind.EQ.1) THEN
            IF (noncolin) THEN
               CALL mp_sum(rho%ns_nc, intra_image_comm)
            ELSE
               CALL mp_sum(rho%ns, intra_image_comm)
            ENDIF
         ELSEIF (lda_plus_u_kind.EQ.2) THEN
            CALL mp_sum(nsg, intra_image_comm)
         ENDIF
         !
         ! If projections on Hubbard manifold are read from file, there is no
         ! need to set starting values: reset them to prevent any problem
         starting_ns = -1.0_dp
         !
      ENDIF
      !
      IF ( okpaw ) THEN
         !
         ! Also the PAW coefficients are needed:
         !
         IF ( ionode ) THEN
            OPEN ( NEWUNIT=iunpaw, FILE = TRIM(dirname) // 'paw.txt', &
                 FORM='formatted', STATUS='old', IOSTAT=ierr )
            READ( UNIT = iunpaw, FMT = *, iostat=ierr ) rho%bec
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('read_scf', 'Reading PAW becsum',1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunpaw, STATUS = 'KEEP')
         ELSE
            rho%bec(:,:,:) = 0.D0
         END IF
         CALL mp_sum(rho%bec, intra_image_comm)
         !
      END IF
      !
      ! by LSH
      IF ( lda_plus_v ) THEN
         !
         ! The occupations ns also need to be read in order to build up
         ! the potential
         !
         CALL read_ehub_ns( rho )
         CALL read_ehub_uv
         !
      ENDIF
      ! by LSH
      !
      RETURN
    END SUBROUTINE read_scf
    !
    ! by LSH
    SUBROUTINE read_ehub_ns( rho )
      !
      USE scf,              ONLY : scf_type
      USE ldaV,             ONLY : lda_plus_v, lacbn0, acbn0_type
      USE noncollin_module, ONLY : noncolin
      USE io_files,         ONLY : restart_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      !
      IMPLICIT NONE

      TYPE(scf_type), INTENT(INOUT) :: rho

      CHARACTER(LEN=256)            :: dirname
      LOGICAL                       :: lexist
      INTEGER                       :: iunoccv, ierr
      INTEGER, EXTERNAL             :: find_free_unit

      dirname = restart_dir ( )
      !
      iunoccv = find_free_unit ()
      !
      IF ( ionode ) THEN
         !
         OPEN (UNIT = iunoccv, FILE = TRIM(dirname)//'occupation_uv.txt', FORM = 'FORMATTED', STATUS = 'old', IOSTAT = ierr )
         !
         IF ( noncolin ) THEN
            !
            IF ( lacbn0 ) THEN
               IF ( acbn0_type == 1 ) THEN
                  READ ( UNIT = iunoccv, FMT = *, IOSTAT = ierr ) rho%na0_nc
               ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                  READ ( UNIT = iunoccv, FMT = *, IOSTAT = ierr ) rho%ne0_nc
               ENDIF
            ELSE
               READ ( UNIT = iunoccv, FMT = *, IOSTAT = ierr ) rho%ne0_nc
            ENDIF
            !
         ELSE
            !
            IF ( lacbn0 ) THEN
               IF ( acbn0_type == 1 ) THEN
                  READ ( UNIT = iunoccv, FMT = *, IOSTAT = ierr ) rho%na0
               ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                  READ ( UNIT = iunoccv, FMT = *, IOSTAT = ierr ) rho%ne0
               ENDIF
            ELSE
               READ( UNIT = iunoccv, FMT = *, IOSTAT = ierr ) rho%ne0
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(ierr, ionode_id, intra_image_comm)
      !
      IF ( ierr/=0 ) CALL errore('read_scf', 'Reading occupation matrices ns', 1)
      !
      IF ( ionode ) THEN
         CLOSE ( UNIT = iunoccv, STATUS = 'KEEP' )
      ELSE
         !
         IF ( noncolin ) THEN
            !
            IF ( lacbn0 ) THEN
               IF ( acbn0_type == 1 ) THEN
                  rho%na0_nc = 0.D0
               ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                  rho%ne0_nc = 0.D0
               ENDIF
            ELSE
               rho%ne0_nc = 0.D0
            ENDIF
            !
         ELSE
            !
            IF ( lacbn0 ) THEN
               IF ( acbn0_type == 1 ) THEN
                  rho%na0 = 0.D0
               ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
                  rho%ne0 = 0.D0
               ENDIF
            ELSE
               rho%ne0 = 0.D0
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
      IF ( noncolin ) THEN
         !
         IF ( lacbn0 ) THEN
            IF ( acbn0_type == 1 ) THEN
               CALL mp_sum(rho%na0_nc, intra_image_comm)
            ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
               CALL mp_sum(rho%ne0_nc, intra_image_comm)
            ENDIF
         ELSE
            CALL mp_sum(rho%ne0_nc, intra_image_comm)
         ENDIF
         !
      ELSE
         !
         IF ( lacbn0 ) THEN
            IF ( acbn0_type == 1 ) THEN
               CALL mp_sum(rho%na0, intra_image_comm)
            ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
               CALL mp_sum(rho%ne0, intra_image_comm)
            ENDIF
         ELSE
            CALL mp_sum(rho%ne0, intra_image_comm)
         ENDIF
         !
      ENDIF
      !
      RETURN
      !
    END SUBROUTINE read_ehub_ns
      !
    SUBROUTINE read_ehub_uv
      !
      USE ldaV,             ONLY : lda_plus_v, lacbn0, acbn0_type, &
                                   ehub_u, ehub_v
      USE io_files,         ONLY : restart_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      !
      IMPLICIT NONE

      CHARACTER(LEN=256)            :: dirname
      LOGICAL                       :: lexist
      INTEGER                       :: iunoccv, ierr
      INTEGER, EXTERNAL             :: find_free_unit

      dirname = restart_dir ( )
      !
      iunoccv = find_free_unit ()
      !
      IF ( ionode ) THEN
         !
         OPEN (UNIT = iunoccv, FILE = TRIM(dirname)//'ehub_uv.txt', FORM = 'FORMATTED', STATUS = 'old', IOSTAT = ierr )
         !
         IF ( lacbn0 ) THEN
            IF ( acbn0_type == 1 ) THEN
               READ ( UNIT = iunoccv, FMT = *, iostat = ierr ) ehub_u
            ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
               READ ( UNIT = iunoccv, FMT = *, iostat = ierr ) ehub_v
            ENDIF
         ELSE
            READ( UNIT = iunoccv, FMT = *, iostat = ierr ) ehub_v
         ENDIF
         !
      ENDIF
      !
      CALL mp_bcast(ierr, ionode_id, intra_image_comm)
      !
      IF ( ierr/=0 ) CALL errore('read_scf', 'Reading Hubbard U or V', 1)
      !
      IF ( ionode ) THEN
          CLOSE ( UNIT = iunoccv, STATUS = 'KEEP' )
      ELSE
         IF ( lacbn0 ) THEN
            IF ( acbn0_type == 1 ) THEN
               ehub_u = 0.D0
            ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
               ehub_v = 0.D0
            ENDIF
         ELSE
            ehub_v = 0.D0
         ENDIF
      ENDIF
      !
      IF ( lacbn0 ) THEN
         IF ( acbn0_type == 1 ) THEN
            CALL mp_sum(ehub_u, intra_image_comm)
         ELSEIF ( acbn0_type == 2 .OR. acbn0_type == 3 ) THEN
            CALL mp_sum(ehub_v, intra_image_comm)
         ENDIF
      ELSE
         CALL mp_sum(ehub_v, intra_image_comm)
      ENDIF
      !
      RETURN
      !
    END SUBROUTINE read_ehub_uv
      !
      ! by LSH
END MODULE io_rho_xml
