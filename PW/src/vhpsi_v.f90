SUBROUTINE vhpsi_v( ldap, np, mps, psip, hpsi )
  !! This routine computes the Hubbard potential applied to the electronic 
  !! structure of the current k-point. The result is added to hpsi.
  !
  USE kinds,         ONLY : DP
  USE becmod,        ONLY : bec_type, calbec,  &
                            allocate_bec_type, &
                            deallocate_bec_type
  USE ldaV,          ONLY : ehub_lm_max, ehub_lm,      &
                            ehub_nn_vect, ehub_nn_dim, &
                            ehub_csize, nwfcv, wfcv, ofstv
  USE constants,     ONLY : tpi
  USE cell_base,     ONLY : at
  USE klist,         ONLY : xk
  USE wvfct,         ONLY : current_k
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ityp
  USE control_flags, ONLY : gamma_only

  IMPLICIT NONE

  INTEGER,     INTENT(IN)    :: ldap, np, mps
  COMPLEX(DP), INTENT(IN)    :: psip(ldap,mps)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ldap,mps)

  TYPE(bec_type)             :: proj

  INTEGER                    :: na1, na2, nt1, nt2, ldim1, ldim2
  INTEGER                    :: n1, n2, n3, nnc, ics

  REAL(DP),    ALLOCATABLE   :: rtemp(:,:)
  REAL(DP)                   :: arg

  COMPLEX(DP), ALLOCATABLE   :: ctemp(:,:), vtemp(:,:)
  COMPLEX(DP)                :: kphase

  CALL start_clock('vhpsi_v')
  !
  CALL allocate_bec_type(nwfcv, mps, proj)
  !
  CALL calbec(np, wfcv, psip, proj)
  !
  ics = ehub_csize + 1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     ldim1 = ehub_lm(nt1)
     IF ( ldim1 == 0 ) CYCLE
     !
     IF ( gamma_only ) THEN
        ALLOCATE ( rtemp(ldim1,mps) ) ; rtemp = 0.0_DP
     ELSE
        ALLOCATE ( ctemp(ldim1,mps) ) ; ctemp = 0.0_DP
     ENDIF
     !
     DO na2 = 1, nat
        nt2 = ityp(na2)
        !
        ldim2 = ehub_lm(nt2)
        IF ( ldim2 == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        IF ( gamma_only ) THEN
           !
           CALL DGEMM('n', 'n', ldim1, mps, ldim2,                      &
                1.0_dp, v%ne0(1,1,current_spin,na1,na2,1), ehub_lm_max, &
                proj%r(ofstv(na2)+1,1), nwfcv, 0.0_dp, rtemp, ldim1)
           CALL DGEMM('n','n', 2*np, mps, ldim1,      &
                1.0_dp, wfcv(1,ofstv(na1)+1), 2*ldap, &
                rtemp, ldim1, 1.0_dp, hpsi, 2*ldap)
           !
        ELSE
           !
           ALLOCATE ( vtemp(ldim1,ldim2) ) ; vtemp = 0.0_DP
           !
           vtemp(:,:) = v%ne0(1:ldim1,1:ldim2,current_spin,na1,na2,1)
           !
           CALL ZGEMM('n', 'n', ldim1, mps, ldim2, &
                (1.0_dp,0.0_dp), vtemp, ldim1,     &
                proj%k(ofstv(na2)+1,1), nwfcv,     &
                (0.0_dp,0.0_dp), ctemp, ldim1)
           CALL ZGEMM('n', 'n', np, mps, ldim1,              &
                (1.0_dp,0.0_dp), wfcv(1,ofstv(na1)+1), ldap, &
                ctemp, ldim1, (1.0_dp,0.0_dp), hpsi, ldap)
           !
           DEALLOCATE ( vtemp )
           !
        ENDIF
        !
     ENDDO
     !
     IF ( gamma_only ) THEN
        DEALLOCATE ( rtemp )
     ELSE
        DEALLOCATE ( ctemp )
     ENDIF
     !
  ENDDO

  nnc = 2
  !
  DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
  DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
  DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
     !
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     !
     arg = tpi * (                                                     &
      xk(1,current_k) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
      xk(2,current_k) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
      xk(3,current_k) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
     !
     kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        ldim1 = ehub_lm(nt1)
        IF ( ldim1 == 0 ) CYCLE
        !
        IF ( gamma_only ) THEN
           ALLOCATE ( rtemp(ldim1,mps) ) ; rtemp = 0.0_DP
        ELSE
           ALLOCATE ( ctemp(ldim1,mps) ) ; ctemp = 0.0_DP
        ENDIF
        !
        DO na2 = 1, nat
           nt2 = ityp(na2)
           !
           ldim2 = ehub_lm(nt2)
           IF ( ldim2 == 0 ) CYCLE
           !
           IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           !
           IF ( gamma_only ) THEN
              !
              CALL DGEMM('n', 'n', ldim1, mps, ldim2,                        &
                   1.0_dp, v%ne0(1,1,current_spin,na1,na2,nnc), ehub_lm_max, &
                   proj%r(ofstv(na2)+1,1), nwfcv, 0.0_dp, rtemp, ldim1)
              CALL DGEMM('n','n', 2*np, mps, ldim1,      &
                   1.0_dp, wfcv(1,ofstv(na1)+1), 2*ldap, &
                   rtemp, ldim1, 1.0_dp, hpsi, 2*ldap)
              !
           ELSE
              !
              ALLOCATE ( vtemp(ldim1,ldim2) ) ; vtemp = 0.0_DP
              !
              vtemp(:,:) = v%ne0(1:ldim1,1:ldim2,current_spin,na1,na2,nnc)
              !
              CALL ZGEMM('n', 'n', ldim1, mps, ldim2, &
                   kphase, vtemp, ldim1,              &
                   proj%k(ofstv(na2)+1,1), nwfcv,     &
                   (0.0_dp,0.0_dp), ctemp, ldim1)
              CALL ZGEMM('n', 'n', np, mps, ldim1,              &
                   (1.0_dp,0.0_dp), wfcv(1,ofstv(na1)+1), ldap, &
                   ctemp, ldim1, (1.0_dp,0.0_dp), hpsi, ldap)
              !
              DEALLOCATE ( vtemp )
              !
           ENDIF
           !
        ENDDO ! na2
        !
        IF ( gamma_only ) THEN
           DEALLOCATE ( rtemp )
        ELSE
           DEALLOCATE ( ctemp )
        ENDIF
        !
     ENDDO ! na1
     !
     IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
     !
  ENDDO ! n3
  ENDDO ! n2
  ENDDO ! n1
  !
  CALL deallocate_bec_type(proj)
  !
  CALL stop_clock('vhpsi_v')
  !
  RETURN

END SUBROUTINE vhpsi_v
!
SUBROUTINE vhpsi_a0( ldap, np, mps, psip, hpsi )

  USE kinds,         ONLY : DP
  USE becmod,        ONLY : bec_type, calbec,  &
                            allocate_bec_type, &
                            deallocate_bec_type
  USE ldaV,          ONLY : ehub_lm_max, ehub_lm, &
                            nwfcv, wfcv, ofstv
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ityp
  USE control_flags, ONLY : gamma_only

  IMPLICIT NONE

  INTEGER,     INTENT(IN)    :: ldap, np, mps
  COMPLEX(DP), INTENT(IN)    :: psip(ldap,mps)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ldap,mps)

  TYPE(bec_type)             :: proj

  INTEGER                    :: na, nt, ldim

  REAL(DP),    ALLOCATABLE   :: rtemp(:,:)

  COMPLEX(DP), ALLOCATABLE   :: ctemp(:,:), vtemp(:,:)

  CALL start_clock('vhpsi_a0')
  !
  CALL allocate_bec_type(nwfcv, mps, proj)
  !
  CALL calbec(np, wfcv, psip, proj)
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     ldim = ehub_lm(nt)
     IF ( ldim == 0 ) CYCLE
     !
     IF ( gamma_only ) THEN
        ALLOCATE ( rtemp(ldim,mps) ) ; rtemp = 0.0_DP
     ELSE
        ALLOCATE ( ctemp(ldim,mps) ) ; ctemp = 0.0_DP
        ALLOCATE ( vtemp(ldim,ldim) ) ; vtemp = 0.0_DP
     ENDIF
     !
     IF ( gamma_only ) THEN
        !
        CALL DGEMM('n', 'n', ldim, mps, ldim,                 &
             1.0_dp, v%na0(1,1,current_spin,na), ehub_lm_max, &
             proj%r(ofstv(na)+1,1), nwfcv, 0.0_dp, rtemp, ldim)
        CALL DGEMM('n','n', 2*np, mps, ldim,      &
             1.0_dp, wfcv(1,ofstv(na)+1), 2*ldap, &
             rtemp, ldim, 1.0_dp, hpsi, 2*ldap)
        !
     ELSE
        !
        vtemp(:,:) =  v%na0(1:ldim,1:ldim,current_spin,na)
        !
        CALL ZGEMM('n', 'n', ldim, mps, ldim, &
             (1.0_dp,0.0_dp), vtemp, ldim,    &
             proj%k(ofstv(na)+1,1), nwfcv,    &
             (0.0_dp,0.0_dp), ctemp, ldim)
        CALL ZGEMM('n', 'n', np, mps, ldim,              &
             (1.0_dp,0.0_dp), wfcv(1,ofstv(na)+1), ldap, &
             ctemp, ldim, (1.0_dp,0.0_dp), hpsi, ldap)
        !
     ENDIF
     !
     IF (gamma_only) THEN
        DEALLOCATE (rtemp)
     ELSE
        DEALLOCATE (ctemp, vtemp)
     ENDIF
     !
  ENDDO
  !
  CALL deallocate_bec_type(proj)
  !
  CALL stop_clock('vhpsi_a0')
  !
  RETURN

END SUBROUTINE vhpsi_a0
!
! by WIY
SUBROUTINE vhpsi_v_nc( ldap, np, mps, psip, hpsi )
   !! Noncollinear version of vhpsi_v routine following vhpsi_nc routine.
   !
   USE kinds,              ONLY : DP
   USE becmod,             ONLY : bec_type, calbec,  &
                                  allocate_bec_type, &
                                  deallocate_bec_type
   USE ldaV,               ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                                  ehub_nn_vect, ehub_nn_dim,    &
                                  ehub_csize, nwfcv, wfcv,      &
                                  ofstv, counter_lm
   USE constants,          ONLY : tpi
   USE cell_base,          ONLY : at
   USE klist,              ONLY : xk
   USE wvfct,              ONLY : current_k
   USE lsda_mod,           ONLY : current_spin
   USE scf,                ONLY : v
   USE ions_base,          ONLY : nat, ityp
   USE noncollin_module,   ONLY : npol, noncolin
   USE mp_pools,           ONLY : inter_pool_comm
   USE mp_bands,           ONLY : intra_bgrp_comm
   USE mp,                 ONLY : mp_sum
   USE io_global,          ONLY : stdout
   !
   IMPLICIT NONE
   !
   INTEGER,     INTENT(IN)    :: ldap
   !! leading dimension of arrays psip, hpsi
   INTEGER,     INTENT(IN)    :: np
   !! true dimension of psip, hpsi
   INTEGER,     INTENT(IN)    :: mps
   !! number of states psip
   COMPLEX(DP), INTENT(IN)    :: psip(ldap*npol,mps)
   !! the wavefunction
   COMPLEX(DP), INTENT(INOUT) :: hpsi(ldap*npol,mps)
   !! Hamiltonian dot psi
   !
   COMPLEX(DP), ALLOCATABLE   :: proj(:,:)
   !
   INTEGER                    :: na, na1, na2, nt1, nt2, ldim1, ldim2
   INTEGER                    :: n1, n2, n3, nnc, ics
   INTEGER                    :: is1, is2
   INTEGER                    :: ibnd, m1, m2
   !
   REAL(DP)                   :: arg
   !
   COMPLEX(DP)                :: temp
   COMPLEX(DP)                :: kphase
   !
   CALL start_clock('vhpsi_v_nc')
   !
   ALLOCATE( proj(nwfcv, mps))
   !
   !! proj = <wfcU|psip>
   !
   Do ibnd = 1, mps
      Do na = 1, nwfcv
         proj(na, ibnd) = DOT_PRODUCT( wfcv(1:ldap*npol, na), psip(1:ldap*npol, ibnd))
      ENDDO
   ENDDO
   !
   CALL mp_sum( proj, intra_bgrp_comm)
   !
   ics = ehub_csize + 1
   !
   ! Compute teh action of the Hubbard potential on the KS wave functons:
   ! V_Hub |psip > = \sum v%ne0 |wfcv><wfcv|psip>
   ! where v%ne0 is computed in v_hubbard_v.f90
   !
   DO ibnd = 1, mps
      DO na1 = 1, nat
         nt1 = ityp(na1)
         !
         ldim1 = ehub_lm(nt1)
         IF ( ldim1 == 0 ) CYCLE
      !
      DO na2 = 1, nat
         nt2 = ityp(na2)
         !
         ldim2 = ehub_lm(nt2)
         IF ( ldim2 == 0 ) CYCLE
         !
         IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
         !
         !! Since the order of proj is interwined with quantum lm number and spin index,
         !! ZGEMM is hard to be used.
         !
         DO m1 = 1, ldim1
         DO is1 = 1, npol
            !
            temp = 0.0_DP
            !
            DO m2 = 1, ldim2
            DO is2 = 1, npol
               temp = temp + v%ne0_nc(m1, m2, npol*(is1-1)+is2, na1, na2, 1) * &
                              proj(ofstv(na2)+m2+counter_lm(m2,nt2) +          &
                                  (2*(ehub_l(m2,nt2)-1)+1)*(is2-1), ibnd)
            ENDDO ! is2
            ENDDO ! m2
               CALL ZAXPY(ldap*npol, temp, wfcv(1, ofstv(na1)+m1+counter_lm(m1,nt1) + &
                                               (2*(ehub_l(m1,nt1)-1)+1)*(is1-1)),     &
                                               1, hpsi(1, ibnd), 1)
         !
         ENDDO ! is1
         ENDDO ! m1
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
         arg = tpi * (                                                       &
            xk(1,current_k) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
            xk(2,current_k) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
            xk(3,current_k) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
         !
         kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
         !
         DO na1 = 1, nat
            nt1 = ityp(na1)
            !
            ldim1 = ehub_lm(nt1)
            IF( ldim1 == 0 ) CYCLE 
            !
            DO na2 = 1, nat
               nt2 = ityp(na2)
               !
               ldim2 = ehub_lm(nt2)
               IF( ldim2 == 0 ) CYCLE
               !
               IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
               !
               DO m1 = 1, ldim1
               DO is1 = 1, npol
                  !
                  temp = 0.0_DP
                  !
                  DO m2 = 1, ldim2
                  DO is2 = 1, npol
                     temp = temp + kphase * v%ne0_nc(m1, m2, npol*(is1-1)+is2, na1, na2, nnc) * &
                                             proj(ofstv(na2)+m2+counter_lm(m2,nt2) +            &
                                                 (2*(ehub_l(m2,nt2)-1)+1)*(is2-1), ibnd)
                  ENDDO ! is2
                  ENDDO ! m2
                     CALL ZAXPY(ldap*npol, temp, wfcv(1, ofstv(na1)+m1+counter_lm(m1,nt1) + &
                                                     (2*(ehub_l(m1,nt1)-1)+1)*(is1-1)),     &
                                                     1, hpsi(1, ibnd), 1)
               !
               ENDDO ! is1
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
   ENDDO ! mps
   !
   DEALLOCATE(proj)
   CALL stop_clock('vhpsi_v_nc')
   !
   RETURN
 
 END SUBROUTINE vhpsi_v_nc
 ! by WIY
