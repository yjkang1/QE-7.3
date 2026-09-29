!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_ns_v(ns)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at
  USE klist,                ONLY : nks, ngk, xk
  USE ldaV,                 ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                                   ehub_nn_vect, ehub_nn_dim,    &
                                   ehub_nnnc, ehub_csize, wfcv,  &
                                   nwfcv, ofstv
  USE ldaU,                 ONLY : Hubbard_projectors
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcv, iunhubv
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf

  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)

  TYPE occmat_5
     REAL(DP), ALLOCATABLE :: O(:,:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:,:)
  END TYPE occmat_5

  TYPE(bec_type)        :: proj

  COMPlEX(DP)           :: kphase

  REAL(DP)              :: arg

  TYPE(occmat_5)        :: n(ehub_nnnc)

  INTEGER               :: ik, ibnd, is, i, npw
  INTEGER               :: na1, na2, nt1, nt2
  INTEGER               :: m1, m2, m3, m4
  INTEGER               :: n1, n2, n3, nnc, ics

  CALL start_clock('new_ns_v')
  !
  CALL allocate_bec_type(nwfcv, nbnd, proj)
  !
  DO i = 1, ehub_nnnc
     ALLOCATE ( n(i)%O(ehub_lm_max,ehub_lm_max,nspin,nat,nat) ) ; n(i)%O = 0.0_DP
  ENDDO
  !
  ics = ehub_csize + 1
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcv | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)     
     IF ( Hubbard_projectors == 'pseudo' ) THEN
        CALL errore('new_ns_v', 'Hubbard_projectors = pseudo is not supported', 1)
     ELSE
        IF ( nks > 1 ) CALL get_buffer(wfcv, nwordwfcv, iunhubv, ik)
        CALL calbec(npw, wfcv, evc, proj)
     ENDIF
     !
     ! calculate the occupation matrices n
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                         proj%r(ofstv(na2)+m2,ibnd) * &
                         proj%r(ofstv(na1)+m1,ibnd)
              ENDDO
           ELSE
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                    DBLE(proj%k(ofstv(na2)+m2,ibnd) * &
                   CONJG(proj%k(ofstv(na1)+m1,ibnd)))
              ENDDO
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           n(1)%O(m2,m1,current_spin,na2,na1) = n(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
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
        arg = tpi * (                                              &
         xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
         xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
         xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
        !
        kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
           DO m1 = 1, ehub_lm(nt1)
           DO m2 = 1, ehub_lm(nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                 DBLE(kphase*wg(ibnd,ik) * &
                              proj%r(ofstv(na2)+m2,ibnd) * &
                              proj%r(ofstv(na1)+m1,ibnd))
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                 DBLE(kphase*wg(ibnd,ik) * &
                              proj%k(ofstv(na2)+m2,ibnd) * &
                        CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
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
  ENDDO ! ik
  !
  CALL deallocate_bec_type(proj)
  !
  DO i = 1, ehub_nnnc
     CALL mp_sum(n(i)%O, inter_pool_comm)
  ENDDO
  !
  DO i = 1, ehub_nnnc
     IF ( nspin == 1 ) n(i)%O = 0.5_DP * n(i)%O
  ENDDO
  !
  CALL symmetrize_n5(n)
  !
  ns = 0.0_DP
  !
  DO i = 1, ehub_nnnc
     ns(:,:,:,:,:,i) = n(i)%O(:,:,:,:,:)
     DEALLOCATE ( n(i)%O )
  ENDDO
  !
  CALL stop_clock('new_ns_v')
  !
  RETURN
  !
END SUBROUTINE new_ns_v


SUBROUTINE new_ns_a0(ns)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE cell_base,            ONLY : at, alat
  USE klist,                ONLY : nks, ngk, xk
  USE ldaV,                 ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                                   ehub_mixing, ehub_u, wfcv,    &
                                   nwfcv, ofstv, remove_ehub_u,  &
                                   write_ehub_uv, ehub_u_raw
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, wg
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcv, iunhubv
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf

  IMPLICIT NONE

  REAL(DP),    INTENT(OUT) :: ns(ehub_lm_max,ehub_lm_max,nspin,nat)

  TYPE hubbard_ene_3
     REAL(DP), ALLOCATABLE :: U(:,:,:)
     REAL(DP), ALLOCATABLE :: J(:,:,:)
  END TYPE hubbard_ene_3

  TYPE occmat_4
     REAL(DP), ALLOCATABLE :: O(:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:)
  END TYPE occmat_4

  TYPE(bec_type)           :: proj

  TYPE(hubbard_ene_3)      :: A0, HF

  TYPE(occmat_4)           :: num(nbnd), n, p

  REAL(DP),    ALLOCATABLE :: ehub_u_(:,:,:)

  REAL(DP), EXTERNAL       :: compute_eri
  REAL(DP)                 :: c1(6), c2(6), c3(6), c4(6)

  INTEGER                  :: i, ik, ibnd, is, npw
  INTEGER                  :: na, nt, l, m
  INTEGER                  :: na1, na2
  INTEGER                  :: nt1, nt2
  INTEGER                  :: m1, m2, m3, m4
  INTEGER                  :: l1, l2, l3, l4

  CALL start_clock('new_ns_a0')
  !
  CALL allocate_bec_type(nwfcv, nbnd, proj)
  !
  ALLOCATE ( n%O(ehub_lm_max,ehub_lm_max,nspin,nat) ) ; n%O = 0.0_DP
  ALLOCATE ( p%O(ehub_lm_max,ehub_lm_max,nspin,nat) ) ; p%O = 0.0_DP
  ALLOCATE ( A0%U(4,4,nat), A0%J(4,4,nat) )
  ALLOCATE ( HF%U(4,4,nat), HF%J(4,4,nat) )
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcv | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( nks > 1 ) CALL get_buffer(wfcv, nwordwfcv, iunhubv, ik)
     CALL calbec(npw, wfcv, evc, proj)
     !
     ! calculate the occupation matrices n
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt)
        DO m2 = 1, ehub_lm(nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n%O(m1,m2,current_spin,na) = &
                 n%O(m1,m2,current_spin,na) + &
                                wg(ibnd,ik) * &
                  proj%r(ofstv(na)+m2,ibnd) * &
                  proj%r(ofstv(na)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 n%O(m1,m2,current_spin,na) = &
                 n%O(m1,m2,current_spin,na) + &
                                wg(ibnd,ik) * &
             DBLE(proj%k(ofstv(na)+m2,ibnd) * &
            CONJG(proj%k(ofstv(na)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na
     !
     ! calculate the renormalized occupations N(ACBN0)
     !
     DO i = 1, nbnd
        ALLOCATE ( num(i)%O(4,4,nspin,ntyp) ) ; num(i)%O = 0.0_DP
        ALLOCATE ( num(i)%T(4,4,nspin,ntyp) ) ; num(i)%T = 0.0_DP
     ENDDO
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m = 1, ehub_lm(nt)
           !
           l = ehub_l(m,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt) + &
                         proj%r(ofstv(na)+m,ibnd) * &
                         proj%r(ofstv(na)+m,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt) + &
                    DBLE(proj%k(ofstv(na)+m,ibnd) * &
                   CONJG(proj%k(ofstv(na)+m,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m
        !
     ENDDO ! na
     !
     DO i = 1, nbnd
        !
        DO nt1 = 1, ntyp
        DO nt2 = 1, ntyp
           !
           IF ( upf(nt1)%psd /= upf(nt2)%psd ) CYCLE
           !
           num(i)%T(:,:,:,nt1) = num(i)%T(:,:,:,nt1) &
                               + num(i)%O(:,:,:,nt2)
           !
        ENDDO ! na2
        ENDDO ! na1
        !
     ENDDO ! i
     !
     DO i = 1, nbnd
        num(i)%O = num(i)%T
        DEALLOCATE ( num(i)%T )
     ENDDO
     !
     ! calculate the renormalized occupation matrices P(ACBN0)
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt)
           !
           l1 = ehub_l(m1,nt)
           !
        DO m2 = 1, ehub_lm(nt)
           !
           l2 = ehub_l(m2,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 p%O(m1,m2,current_spin,na) = &
                 p%O(m1,m2,current_spin,na) + &
                         num(ibnd)%O(l1,l2,   &
                           current_spin,nt) * &
                                wg(ibnd,ik) * &
                  proj%r(ofstv(na)+m2,ibnd) * &
                  proj%r(ofstv(na)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 p%O(m1,m2,current_spin,na) = &
                 p%O(m1,m2,current_spin,na) + &
                         num(ibnd)%O(l1,l2,   &
                           current_spin,nt) * &
                                wg(ibnd,ik) * &
             DBLE(proj%k(ofstv(na)+m2,ibnd) * &
            CONJG(proj%k(ofstv(na)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na
     !
     DO i = 1, nbnd
        DEALLOCATE ( num(i)%O )
     ENDDO

  ENDDO ! ik

  CALL deallocate_bec_type(proj)

  CALL mp_sum(n%O, inter_pool_comm)
  CALL mp_sum(p%O, inter_pool_comm)

  IF ( nspin == 1 ) n%O = 0.5_DP * n%O
  !
  CALL symmetrize_n4(n)
  !
  ns = 0.0_DP
  !
  ns = n%O
  !
  DEALLOCATE (n%O)
  !
  IF ( nspin == 1 ) p%O = 0.5_DP * p%O
  !
  CALL symmetrize_n4(p)
  !
  ! calculate ACBN0 U & J
  !
  A0%U = 0.0_DP ; A0%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, ehub_lm(nt)
        CALL orberi(ehub_l(m1,nt), l1)
        CALL pao_3g(upf(nt)%psd, ehub_l(m1,nt), c1)
        !
     l2 = 0
     DO m2 = 1, ehub_lm(nt)
        CALL orberi(ehub_l(m2,nt), l2)
        CALL pao_3g(upf(nt)%psd, ehub_l(m2,nt), c2)
        !
     l3 = 0
     DO m3 = 1, ehub_lm(nt)
        IF ( ehub_l(m3,nt) /= ehub_l(m1,nt) ) CYCLE
        CALL orberi(ehub_l(m3,nt), l3)
        CALL pao_3g(upf(nt)%psd, ehub_l(m3,nt), c3)
        !
     l4 = 0
     DO m4 = 1, ehub_lm(nt)
        IF ( ehub_l(m4,nt) /= ehub_l(m2,nt) ) CYCLE
        CALL orberi(ehub_l(m4,nt), l4)
        CALL pao_3g(upf(nt)%psd, ehub_l(m4,nt), c4)
        !
        IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
              A0%U(ehub_l(m1,nt),ehub_l(m2,nt),na)  = &
              A0%U(ehub_l(m1,nt),ehub_l(m2,nt),na)  + &
                  (p%O(m3,m1,1,na)*p%O(m4,m2,2,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,1,na)) * &
                          compute_eri(l3,l1,l4,l2,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                                      c3,c1,c4,c2)
        ELSE
              A0%U(ehub_l(m1,nt),ehub_l(m2,nt),na)  = &
              A0%U(ehub_l(m1,nt),ehub_l(m2,nt),na)  + &
                  (p%O(m3,m1,1,na)*p%O(m4,m2,1,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,2,na)  + &
                   p%O(m3,m1,1,na)*p%O(m4,m2,2,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,1,na)) * &
                          compute_eri(l3,l1,l4,l2,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                                      c3,c1,c4,c2)
              A0%J(ehub_l(m1,nt),ehub_l(m2,nt),na)  = &
              A0%J(ehub_l(m1,nt),ehub_l(m2,nt),na)  + &
                  (p%O(m3,m1,1,na)*p%O(m4,m2,1,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,2,na)) * &
                          compute_eri(l3,l2,l4,l1,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                                      c3,c2,c4,c1)
        ENDIF
        !
     ENDDO ! m4
     ENDDO ! m3
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  DEALLOCATE (p%O)
  !
  HF%U = 0.0_DP ; HF%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt)
     DO m2 = 1, ehub_lm(nt)
        !
        IF ( m1 == m2 ) THEN
             HF%U(ehub_l(m1,nt),ehub_l(m2,nt),na) = &
             HF%U(ehub_l(m1,nt),ehub_l(m2,nt),na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,2,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,1,na)
        ELSE
             HF%U(ehub_l(m1,nt),ehub_l(m2,nt),na) = &
             HF%U(ehub_l(m1,nt),ehub_l(m2,nt),na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,1,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,2,na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,2,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,1,na)
 
             HF%J(ehub_l(m1,nt),ehub_l(m2,nt),na) = &
             HF%J(ehub_l(m1,nt),ehub_l(m2,nt),na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,1,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,2,na)
        ENDIF
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  ALLOCATE (ehub_u_(nat,ehub_lm_max,ehub_lm_max)) ; ehub_u_ = 0.0_DP
  IF ( .NOT.ALLOCATED(ehub_u_raw) ) THEN   ! [rawprint]
     ALLOCATE ( ehub_u_raw(SIZE(ehub_u,1),SIZE(ehub_u,2),SIZE(ehub_u,3)) )
     ehub_u_raw = 0.0_DP
  ENDIF
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt)
     DO m2 = 1, ehub_lm(nt)
        !
        IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
           ehub_u_(na,m1,m2) =                          &
                 A0%U(ehub_l(m1,nt),ehub_l(m2,nt),na) / &
                 HF%U(ehub_l(m1,nt),ehub_l(m2,nt),na) * &
                      remove_ehub_u(nt,ehub_l(m1,nt)) * &
                      remove_ehub_u(nt,ehub_l(m2,nt))
        ELSE
           ehub_u_(na,m1,m2) =                          &
                (A0%U(ehub_l(m1,nt),ehub_l(m2,nt),na) / &
                 HF%U(ehub_l(m1,nt),ehub_l(m2,nt),na) - &
                 A0%J(ehub_l(m1,nt),ehub_l(m2,nt),na) / &
                 HF%J(ehub_l(m1,nt),ehub_l(m2,nt),na))* &
                      remove_ehub_u(nt,ehub_l(m1,nt)) * &
                      remove_ehub_u(nt,ehub_l(m2,nt))
        ENDIF
        !
        ehub_u_raw(na,m1,m2) = ehub_u_(na,m1,m2)   ! [rawprint]
        !
        ehub_u(na,m1,m2) =                               &
              (1.0_DP - ehub_mixing)*ehub_u (na,m1,m2) + &
              (         ehub_mixing)*ehub_u_(na,m1,m2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = m1, ehub_lm(nt1)
        !
        ehub_u(na1,m1,m2) = 0.5D0 * (ehub_u(na1,m1,m2) + &
                                     ehub_u(na1,m2,m1))
        !
        ehub_u(na1,m2,m1) = ehub_u(na1,m1,m2)
        ! [rawprint]
        ehub_u_raw(na1,m1,m2) = 0.5D0 * (ehub_u_raw(na1,m1,m2) + ehub_u_raw(na1,m2,m1))
        ehub_u_raw(na1,m2,m1) = ehub_u_raw(na1,m1,m2)
        !
     ENDDO
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE ( ehub_u_ )
  DEALLOCATE ( A0%U, A0%J )
  DEALLOCATE ( HF%U, HF%J )

  CALL write_ehub_uv

  CALL stop_clock('new_ns_a0')

  RETURN

END SUBROUTINE new_ns_a0


SUBROUTINE new_ns_e0_1(ns)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi, rytoev
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE cell_base,            ONLY : at, alat
  USE klist,                ONLY : nks, ngk, xk
  USE ldaV,                 ONLY : ehub_lm_max, ehub_lm, ehub_l,                      &
                                   ehub_nn_vect, ehub_nn_leng, ehub_nn_dim,           &
                                   ehub_csize, ehub_nnnc, ehub_mixing, ehub_conv_thr, &
                                   ehub_v, ehub_v_diff, remove_ehub_u, write_ehub_uv, &
                                   wfcv, nwfcv, ofstv, stop_acbn0, ehub_v_raw
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcv, iunhubv
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf

  IMPLICIT NONE

  REAL(DP), INTENT(OUT)    :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)

  TYPE hubbard_ene_4
     REAL(DP), ALLOCATABLE :: U(:,:,:,:)
     REAL(DP), ALLOCATABLE :: V(:,:,:,:)
     REAL(DP), ALLOCATABLE :: J(:,:,:,:)
  END TYPE hubbard_ene_4

  TYPE occmat_4
     REAL(DP), ALLOCATABLE :: O(:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:)
  END TYPE occmat_4

  TYPE occmat_5
     REAL(DP), ALLOCATABLE :: O(:,:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:,:)
  END TYPE occmat_5

  TYPE(bec_type)           :: proj

  TYPE(hubbard_ene_4)      :: E0(ehub_nnnc)
  TYPE(hubbard_ene_4)      :: HF(ehub_nnnc)

  TYPE(occmat_5)           :: n(ehub_nnnc)
  TYPE(occmat_5)           :: p(ehub_nnnc)
  TYPE(occmat_5)           :: num(nbnd)

  TYPE(occmat_4)           :: dm0
  TYPE(occmat_4)           :: dm1(ntyp,4)

  COMPLEX(DP)              :: kphase

  REAL(DP), ALLOCATABLE    :: ehub_v_(:,:,:,:)

  REAL(DP)                 :: arg
  REAL(DP)                 :: ntau(3,nat)
  REAL(DP)                 :: compute_eri
  REAL(DP)                 :: c1(6), c2(6), c3(6), c4(6)

  INTEGER                  :: npw, nw
  INTEGER                  :: i, ik, ibnd, is, ix
  INTEGER                  :: na, nt, m, l
  INTEGER                  :: na1, na2
  INTEGER                  :: nt1, nt2
  INTEGER                  :: m1, m2, m3, m4
  INTEGER                  :: l1, l2, l3, l4
  INTEGER                  :: n1, n2, n3
  INTEGER                  :: nnc, ics

  CALL start_clock('new_ns_e0_#1')

  IF ( stop_acbn0 ) THEN
     !
     IF ( MAXVAL(ehub_v_diff(:,:,:,:,:)) < SQRT(ehub_conv_thr) ) THEN
        CALL new_ns_v(ns)
        GOTO 1
     ENDIF
     !
     ehub_v_diff = 0.0_DP
     !
  ENDIF

  CALL allocate_bec_type(nwfcv, nbnd, proj)
  !
  ALLOCATE ( dm0%O(ehub_lm_max,ehub_lm_max,nspin,nat) ) ; dm0%O = 0.0_DP
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        ALLOCATE ( dm1(nt,l)%O(ehub_lm_max,ehub_lm_max,nspin,nat) ) ; dm1(nt,l)%O = 0.0_DP
     ENDDO
  ENDDO
  !
  DO i = 1, ehub_nnnc
     ALLOCATE ( n(i)%O(ehub_lm_max,ehub_lm_max,nspin,nat,nat) ) ; n(i)%O = 0.0_DP
     ALLOCATE ( p(i)%O(ehub_lm_max,ehub_lm_max,nspin,nat,nat) ) ; p(i)%O = 0.0_DP
     ALLOCATE ( E0(i)%V(4,4,nat,nat), E0(i)%J(4,4,nat,nat) )
     ALLOCATE ( HF(i)%V(4,4,nat,nat), HF(i)%J(4,4,nat,nat) )
  ENDDO

  ics = ehub_csize + 1
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcv | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( nks > 1 ) CALL get_buffer(wfcv, nwordwfcv, iunhubv, ik)
     CALL calbec(npw, wfcv, evc, proj)
     !
     ! calculate the occupation matrices n
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2 ) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                         proj%r(ofstv(na2)+m2,ibnd) * &
                         proj%r(ofstv(na1)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                    DBLE(proj%k(ofstv(na2)+m2,ibnd) * &
                   CONJG(proj%k(ofstv(na1)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           n(1)%O(m2,m1,current_spin,na2,na1) = n(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
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
        arg = tpi * (                                              &
         xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
         xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
         xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
        !
        kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
           DO m1 = 1, ehub_lm(nt1)
           DO m2 = 1, ehub_lm(nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                             wg(ibnd,ik) * &
                              proj%r(ofstv(na2)+m2,ibnd) * &
                              proj%r(ofstv(na1)+m1,ibnd)
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                             wg(ibnd,ik) * &
                                           DBLE((kphase) * &
                              proj%k(ofstv(na2)+m2,ibnd) * &
                        CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
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
     ! calculate the renormalized occupations N(ACBN0)
     !
     DO i = 1, nbnd
        ALLOCATE ( num(i)%O(4,4,nspin,ntyp,ntyp) ) ; num(i)%O = 0.0_DP
        ALLOCATE ( num(i)%T(4,4,nspin,ntyp,ntyp) ) ; num(i)%T = 0.0_DP
     ENDDO
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m = 1, ehub_lm(nt)
           !
           l = ehub_l(m,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt,nt) + &
                            proj%r(ofstv(na)+m,ibnd) * &
                            proj%r(ofstv(na)+m,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt,nt) + &
                       DBLE(proj%k(ofstv(na)+m,ibnd) * &
                      CONJG(proj%k(ofstv(na)+m,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m
        !
     ENDDO ! na
     !
     DO ibnd = 1, nbnd
        !
        DO nt1 = 1, ntyp
        DO nt2 = 1, ntyp
           !
           IF ( upf(nt1)%psd /= upf(nt2)%psd ) CYCLE
           !
           num(ibnd)%T(:,:,:,nt1,nt1) = num(ibnd)%T(:,:,:,nt1,nt1) + &
                                        num(ibnd)%O(:,:,:,nt2,nt2)
           !
        ENDDO ! nt2
        ENDDO ! nt1
        !
     ENDDO ! ibnd
     !
     DO ibnd = 1, nbnd
        !
        DO nt1 = 1, ntyp
        DO nt2 = 1, ntyp
           !
           DO l1 = 1, 4
           DO l2 = 1, 4
              !
              ! Our definition
              num(ibnd)%T(l1,l2,:,nt1,nt2) = num(ibnd)%T(l1,l1,:,nt1,nt1) + &
                                             num(ibnd)%T(l2,l2,:,nt2,nt2)
              !
              IF ( upf(nt1)%psd == upf(nt2)%psd )                           &
              num(ibnd)%T(l1,l2,:,nt1,nt2) = num(ibnd)%T(l1,l2,:,nt1,nt2) / &
                                             2.0_DP
              !
              ! The definition of Nicolas Tancogne-Dejean & Angel Rubio
!             num(ibnd)%T(l1,l2,:,nt1,nt2) = SQRT(num(ibnd)%T(l1,l1,:,nt1,nt1) * &
!                                                 num(ibnd)%T(l2,l2,:,nt2,nt2))
              !
           ENDDO ! l2
           ENDDO ! l1
           !
        ENDDO ! nt2
        ENDDO ! nt1
        !
     ENDDO ! ibnd
     !
     DO ibnd = 1, nbnd
        num(ibnd)%O = num(ibnd)%T
        DEALLOCATE ( num(ibnd)%T )
     ENDDO
     !
     ! calculate the renormalized occupation matrices P(ACBN0)
     !
     ! 1) P^II (dm0, dm1)
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt)
           !
           l1 = ehub_l(m1,nt)
           !
        DO m2 = 1, ehub_lm(nt)
           !
           l2 = ehub_l(m2,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 dm0%O(m1,m2,current_spin,na) = &
                 dm0%O(m1,m2,current_spin,na) + &
                           num(ibnd)%O(l1,l2,   &
                          current_spin,nt,nt) * &
                                  wg(ibnd,ik) * &
                    proj%r(ofstv(na)+m2,ibnd) * &
                    proj%r(ofstv(na)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 dm0%O(m1,m2,current_spin,na) = &
                 dm0%O(m1,m2,current_spin,na) + &
                           num(ibnd)%O(l1,l2,   &
                          current_spin,nt,nt) * &
                                  wg(ibnd,ik) * &
               DBLE(proj%k(ofstv(na)+m2,ibnd) * &
              CONJG(proj%k(ofstv(na)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na

     DO nt = 1, ntyp
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO l = 1, 4
           !
           DO na1 = 1, nat
              nt1 = ityp(na1)
              !
              IF ( ehub_lm(nt1) == 0 ) CYCLE
              !
              DO m1 = 1, ehub_lm(nt1)
              DO m2 = 1, ehub_lm(nt1)
                 !
                 IF ( gamma_only ) THEN
                    DO ibnd = 1, nbnd
                       dm1(nt,l)%O(m1,m2,current_spin,na1) = &
                       dm1(nt,l)%O(m1,m2,current_spin,na1) + &
                             num(ibnd)%O(ehub_l(m1,nt1),l,   &
                                      current_spin,nt1,nt) * &
                                               wg(ibnd,ik) * &
                                proj%r(ofstv(na1)+m2,ibnd) * &
                                proj%r(ofstv(na1)+m1,ibnd)
                    ENDDO ! ibnd
                 ELSE
                    DO ibnd = 1, nbnd
                       dm1(nt,l)%O(m1,m2,current_spin,na1) = &
                       dm1(nt,l)%O(m1,m2,current_spin,na1) + &
                             num(ibnd)%O(ehub_l(m1,nt1),l,   &
                                      current_spin,nt1,nt) * &
                                               wg(ibnd,ik) * &
                           DBLE(proj%k(ofstv(na1)+m2,ibnd) * &
                          CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                    ENDDO ! ibnd
                 ENDIF
                 !
              ENDDO ! m2
              ENDDO ! m1
              !
           ENDDO ! na1
           !
        ENDDO ! l
        !
     ENDDO ! nt
     !
     ! 2) P^IJ
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
           !
           l1 = ehub_l(m1,nt1)
           !
        DO m2 = 1, ehub_lm(nt2)
           !
           l2 = ehub_l(m2,nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                   p(1)%O(m1,m2,current_spin,na1,na2) = &
                   p(1)%O(m1,m2,current_spin,na1,na2) + &
                                   num(ibnd)%O(l1,l2,   &
                                current_spin,nt1,nt2) * &
                                          wg(ibnd,ik) * &
                           proj%r(ofstv(na2)+m2,ibnd) * &
                           proj%r(ofstv(na1)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                   p(1)%O(m1,m2,current_spin,na1,na2) = &
                   p(1)%O(m1,m2,current_spin,na1,na2) + &
                                   num(ibnd)%O(l1,l2,   &
                                current_spin,nt1,nt2) * &
                                          wg(ibnd,ik) * &
                      DBLE(proj%k(ofstv(na2)+m2,ibnd) * &
                     CONJG(proj%k(ofstv(na1)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           p(1)%O(m2,m1,current_spin,na2,na1) = p(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
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
        arg = tpi * (                                              &
         xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
         xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
         xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
        !
        kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
           DO m1 = 1, ehub_lm(nt1)
              !
              l1 = ehub_l(m1,nt1)
              !
           DO m2 = 1, ehub_lm(nt2)
              !
              l2 = ehub_l(m2,nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    p(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    p(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                      num(ibnd)%O(l1,l2,   &
                                   current_spin,nt1,nt2) * &
                                             wg(ibnd,ik) * &
                              proj%r(ofstv(na2)+m2,ibnd) * &
                              proj%r(ofstv(na1)+m1,ibnd)
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    p(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    p(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                      num(ibnd)%O(l1,l2,   &
                                   current_spin,nt1,nt2) * &
                                             wg(ibnd,ik) * &
                                           DBLE((kphase) * &
                              proj%k(ofstv(na2)+m2,ibnd) * &
                        CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
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
     DO i = 1, nbnd
        DEALLOCATE ( num(i)%O )
     ENDDO

  ENDDO ! ik

  CALL deallocate_bec_type(proj)

  CALL mp_sum(dm0%O, inter_pool_comm)
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        CALL mp_sum(dm1(nt,l)%O, inter_pool_comm)
     ENDDO
  ENDDO
  !
  DO i = 1, ehub_nnnc
     CALL mp_sum(n(i)%O, inter_pool_comm)
     CALL mp_sum(p(i)%O, inter_pool_comm)
  ENDDO

  DO i = 1, ehub_nnnc
     IF ( nspin == 1 ) n(i)%O = 0.5_DP * n(i)%O
  ENDDO
  !
  CALL symmetrize_n5(n)
  !
  ns = 0.0_DP
  !
  DO i = 1, ehub_nnnc
     ns(:,:,:,:,:,i) = n(i)%O(:,:,:,:,:)
     DEALLOCATE ( n(i)%O )
  ENDDO
  !
  IF ( nspin == 1 ) dm0%O = 0.5_DP * dm0%O
  !
  CALL symmetrize_n4(dm0)
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        IF ( nspin == 1 ) dm1(nt,l)%O = 0.5_DP * dm1(nt,l)%O
     ENDDO
  ENDDO
  !
  CALL symmetrize_n4m(dm1)
  !
  DO i = 1, ehub_nnnc
     IF ( nspin == 1 ) p(i)%O = 0.5_DP * p(i)%O
  ENDDO
  !
  CALL symmetrize_n5(p)
  !
  ! on-site interaction: U_ij^II
  !
  E0(1)%V = 0.0_DP ; E0(1)%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, ehub_lm(nt)
        CALL orberi(ehub_l(m1,nt), l1)
        CALL pao_3g(upf(nt)%psd, ehub_l(m1,nt), c1)
        !
     l2 = 0
     DO m2 = 1, ehub_lm(nt)
        CALL orberi(ehub_l(m2,nt), l2)
        CALL pao_3g(upf(nt)%psd, ehub_l(m2,nt), c2)
        !
     l3 = 0
     DO m3 = 1, ehub_lm(nt)
        IF ( ehub_l(m3,nt) /= ehub_l(m1,nt) ) CYCLE
        CALL orberi(ehub_l(m3,nt), l3)
        CALL pao_3g(upf(nt)%psd, ehub_l(m3,nt), c3)
        !
     l4 = 0
     DO m4 = 1, ehub_lm(nt)
        IF ( ehub_l(m4,nt) /= ehub_l(m2,nt) ) CYCLE
        CALL orberi(ehub_l(m4,nt), l4)
        CALL pao_3g(upf(nt)%psd, ehub_l(m4,nt), c4)
        !
        IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,1,na)) * &
                             compute_eri(l3,l1,l4,l2,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c1,c4,c2)
        ELSE
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,1,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,1,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,1,na)) * &
                             compute_eri(l3,l1,l4,l2,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c1,c4,c2)

           E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
           E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,1,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,2,na)) * &
                             compute_eri(l3,l2,l4,l1,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c2,c4,c1)
        ENDIF
        !
     ENDDO ! m4
     ENDDO ! m3
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  HF(1)%V = 0.0_DP ; HF(1)%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt)
     DO m2 = 1, ehub_lm(nt)
        !
        IF ( m1 == m2 ) THEN
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,1,na,na,1)
        ELSE
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,1,na,na,1)

            HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
            HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,2,na,na,1)
        ENDIF
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  ALLOCATE ( ehub_v_(nat,nat,ehub_lm_max,ehub_lm_max) ) ; ehub_v_ = 0.0_DP
  IF ( .NOT.ALLOCATED(ehub_v_raw) ) THEN   ! [rawprint]
     ALLOCATE ( ehub_v_raw(SIZE(ehub_v,1),SIZE(ehub_v,2),SIZE(ehub_v,3),SIZE(ehub_v,4),SIZE(ehub_v,5)) )
     ehub_v_raw = 0.0_DP
  ENDIF
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt)
     DO m2 = 1, ehub_lm(nt)
        !
        IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
           ehub_v_(na,na,m1,m2) =                           &
               E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) / &
               HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) * &
                          remove_ehub_u(nt,ehub_l(m1,nt)) * &
                          remove_ehub_u(nt,ehub_l(m2,nt))
        ELSE
           ehub_v_(na,na,m1,m2) =                            &
              (E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  / &
               HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  - &
               E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  / &
               HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)) * &
                          remove_ehub_u(nt,ehub_l(m1,nt))  * &
                          remove_ehub_u(nt,ehub_l(m2,nt))
        ENDIF
        !
        IF ( stop_acbn0 ) ehub_v_diff(na,na,m1,m2,1) = ABS( ehub_v (na,na,m1,m2,1) &
                                                     -      ehub_v_(na,na,m1,m2  ) )
        !
        ehub_v_raw(na,na,m1,m2,1) = ehub_v_(na,na,m1,m2)   ! [rawprint]
        !
        ehub_v(na,na,m1,m2,1)                                 &
            = (1.0_DP - ehub_mixing)*ehub_v (na,na,m1,m2,1) + &
              (         ehub_mixing)*ehub_v_(na,na,m1,m2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na

  !
  ! inter-site interaction: U_ij^IJ (I /= J)
  !
  ! ... in the unit cell
  !
  E0(1)%V = 0.0_DP
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, ehub_lm(nt1)
        CALL orberi(ehub_l(m1,nt1), l1)
        CALL pao_3g(upf(nt1)%psd, ehub_l(m1,nt1), c1)
        !
     l2 = 0
     DO m2 = 1, ehub_lm(nt2)
        CALL orberi(ehub_l(m2,nt2), l2)
        CALL pao_3g(upf(nt2)%psd, ehub_l(m2,nt2), c2)
        !
     l3 = 0
     DO m3 = 1, ehub_lm(nt1)
        IF ( ehub_l(m3,nt1) /= ehub_l(m1,nt1) ) CYCLE
        CALL orberi(ehub_l(m3,nt1), l3)
        CALL pao_3g(upf(nt1)%psd, ehub_l(m3,nt1), c3)
        !
     l4 = 0
     DO m4 = 1, ehub_lm(nt2)
        IF ( ehub_l(m4,nt2) /= ehub_l(m2,nt2) ) CYCLE
        CALL orberi(ehub_l(m4,nt2), l4)
        CALL pao_3g(upf(nt2)%psd, ehub_l(m4,nt2), c4)
        !
        E0(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  = &
        E0(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  + &
               (dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  + &
                dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  - &
                               p(1)%O(m3,m2,1,na1,na2)  * &
                               p(1)%O(m4,m1,1,na2,na1)  - &
                               p(1)%O(m3,m2,2,na1,na2)  * &
                               p(1)%O(m4,m1,2,na2,na1)) * &
                              compute_eri(l3,l1,l4,l2,    &
                      tau(:,na1)*alat,tau(:,na1)*alat,    &
                      tau(:,na2)*alat,tau(:,na2)*alat,    &
                                          c3,c1,c4,c2)
        !
     ENDDO ! m4
     ENDDO ! m3
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        !
        E0(1)%V(ehub_l(m2,nt2),ehub_l(m1,nt1),na2,na1) = E0(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  HF(1)%V = 0.0_DP
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        !
        HF(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2) = &
        HF(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2) + &
         ns(m1,m1,1,na1,na1,1) * ns(m2,m2,1,na2,na2,1) + &
         ns(m1,m1,2,na1,na1,1) * ns(m2,m2,2,na2,na2,1) + &
         ns(m1,m1,1,na1,na1,1) * ns(m2,m2,2,na2,na2,1) + &
         ns(m1,m1,2,na1,na1,1) * ns(m2,m2,1,na2,na2,1) - &
         ns(m1,m2,1,na1,na2,1) * ns(m2,m1,1,na2,na1,1) - &
         ns(m1,m2,2,na1,na2,1) * ns(m2,m1,2,na2,na1,1)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        HF(1)%V(ehub_l(m2,nt2),ehub_l(m1,nt1),na2,na1) = HF(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
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
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
        !
        l1 = ehub_l(m1,nt1)
        !
     DO m2 = 1, ehub_lm(nt2)
        !
        l2 = ehub_l(m2,nt2)
        !
        ehub_v_(na1,na2,m1,m2) = &
               (E0(1)%V(l1,l2,na1,na2) / HF(1)%V(l1,l2,na1,na2)) / 2.0_DP
        !
        IF ( stop_acbn0 ) ehub_v_diff(na1,na2,m1,m2,1) = ABS( ehub_v (na1,na2,m1,m2,1) &
                                                       -      ehub_v_(na1,na2,m1,m2  ) )
        !
        ehub_v_raw(na1,na2,m1,m2,1) = ehub_v_(na1,na2,m1,m2)   ! [rawprint]
        !
        ehub_v(na1,na2,m1,m2,1)                                 &
            = (1.0_DP - ehub_mixing)*ehub_v (na1,na2,m1,m2,1) + &
              (         ehub_mixing)*ehub_v_(na1,na2,m1,m2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  ! ... out of the unit cell
  !
  DO i = 2, ehub_nnnc
     E0(i)%V = 0.0_DP
  ENDDO
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
        DO ix = 1, 3
           ntau(ix,na2) = n1 * at(ix,1) + n2 * at(ix,2) + n3 * at(ix,3) + tau(ix,na2)
        ENDDO
        !
        l1 = 0
        DO m1 = 1, ehub_lm(nt1)
           CALL orberi(ehub_l(m1,nt1), l1)
           CALL pao_3g(upf(nt1)%psd, ehub_l(m1,nt1), c1)
           !
        l2 = 0
        DO m2 = 1, ehub_lm(nt2)
           CALL orberi(ehub_l(m2,nt2), l2)
           CALL pao_3g(upf(nt2)%psd, ehub_l(m2,nt2), c2)
           !
        l3 = 0
        DO m3 = 1, ehub_lm(nt1)
           IF ( ehub_l(m3,nt1) /= ehub_l(m1,nt1) ) CYCLE
           CALL orberi(ehub_l(m3,nt1), l3)
           CALL pao_3g(upf(nt1)%psd, ehub_l(m3,nt1), c3)
           !
        l4 = 0
        DO m4 = 1, ehub_lm(nt2)
           IF ( ehub_l(m4,nt2) /= ehub_l(m2,nt2) ) CYCLE
           CALL orberi(ehub_l(m4,nt2), l4)
           CALL pao_3g(upf(nt2)%psd, ehub_l(m4,nt2), c4)
           !
           E0(nnc)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  = &
           E0(nnc)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  + &
                    (dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                     dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  + &
                     dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                     dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                     dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                     dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                     dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                     dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  - &
                                  p(nnc)%O(m3,m2,1,na1,na2)  * &
                                  p(nnc)%O(m1,m4,1,na1,na2)  - &
                                  p(nnc)%O(m3,m2,2,na1,na2)  * &
                                  p(nnc)%O(m1,m4,2,na1,na2)) * &
                                   compute_eri(l3,l1,l4,l2,    &
                          tau(:,na1)*alat, tau(:,na1)*alat,    &
                         ntau(:,na2)*alat,ntau(:,na2)*alat,    &
                                               c3,c1,c4,c2)
           !
        ENDDO ! m4
        ENDDO ! m3
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
  DEALLOCATE ( dm0%O )
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        DEALLOCATE ( dm1(nt,l)%O )
     ENDDO
  ENDDO
  !
  DO i = 2, ehub_nnnc
     HF(i)%V = 0.0_DP
  ENDDO
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
        DO m1 = 1, ehub_lm(nt1)
           !
           l1 = ehub_l(m1,nt1)
           !
        DO m2 = 1, ehub_lm(nt2)
           !
           l2 = ehub_l(m2,nt2)
           !
           HF(nnc)%V(l1,l2,na1,na2) = HF(nnc)%V(l1,l2,na1,na2) + &
               ns(m1,m1,1,na1,na1,  1)*ns(m2,m2,1,na2,na2,  1) + &
               ns(m1,m1,2,na1,na1,  1)*ns(m2,m2,2,na2,na2,  1) + &
               ns(m1,m1,1,na1,na1,  1)*ns(m2,m2,2,na2,na2,  1) + &
               ns(m1,m1,2,na1,na1,  1)*ns(m2,m2,1,na2,na2,  1) - &
               ns(m1,m2,1,na1,na2,nnc)*ns(m1,m2,1,na1,na2,nnc) - &
               ns(m1,m2,2,na1,na2,nnc)*ns(m1,m2,2,na1,na2,nnc)
           !
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
        DO m1 = 1, ehub_lm(nt1)
           !
           l1 = ehub_l(m1,nt1)
           !
        DO m2 = 1, ehub_lm(nt2)
           !
           l2 = ehub_l(m2,nt2)
           !
           ehub_v_(na1,na2,m1,m2) = &
                  (E0(nnc)%V(l1,l2,na1,na2) / HF(nnc)%V(l1,l2,na1,na2)) / 2.0_DP
           !
           IF ( stop_acbn0 ) ehub_v_diff(na1,na2,m1,m2,nnc) = ABS( ehub_v (na1,na2,m1,m2,nnc) &
                                                            -      ehub_v_(na1,na2,m1,m2    ) )
           !
           ehub_v_raw(na1,na2,m1,m2,nnc) = ehub_v_(na1,na2,m1,m2)   ! [rawprint]
           !
           ehub_v(na1,na2,m1,m2,nnc) =                               &
                 (1.0_DP - ehub_mixing)*ehub_v (na1,na2,m1,m2,nnc) + &
                 (         ehub_mixing)*ehub_v_(na1,na2,m1,m2    )
           !
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
  DEALLOCATE ( ehub_v_ )
  !
  DO i = 1, ehub_nnnc
     DEALLOCATE ( E0(i)%V, E0(i)%J )
     DEALLOCATE ( HF(i)%V, HF(i)%J )
  ENDDO
  !
  ! Hermiticity of Hubbard parameters U & V
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        ehub_v(na1,na2,m1,m2,1) = 0.5D0 * (ehub_v(na1,na2,m1,m2,1) + &
                                           ehub_v(na2,na1,m2,m1,1))
        !
        ehub_v(na2,na1,m2,m1,1) = ehub_v(na1,na2,m1,m2,1)
        ! [rawprint]
        ehub_v_raw(na1,na2,m1,m2,1) = 0.5D0 * (ehub_v_raw(na1,na2,m1,m2,1) + ehub_v_raw(na2,na1,m2,m1,1))
        ehub_v_raw(na2,na1,m2,m1,1) = ehub_v_raw(na1,na2,m1,m2,1)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
  DO na2 = 1, nat
     IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) THEN
        ehub_v(na1,na2,:,:,1) = 0.0_DP
        ehub_v_raw(na1,na2,:,:,1) = 0.0_DP   ! [rawprint]
     ENDIF
  ENDDO ! na2
  ENDDO ! na1

1 CALL write_ehub_uv

  CALL stop_clock('new_ns_e0_#1')

  RETURN

END SUBROUTINE new_ns_e0_1


SUBROUTINE new_ns_e0_2(ns)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi, rytoev
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE cell_base,            ONLY : at, alat
  USE klist,                ONLY : nks, ngk, xk
  USE ldaV,                 ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                                   ehub_nn_vect, ehub_nn_leng,   &
                                   ehub_nn_dim, ehub_csize,      &
                                   ehub_nnnc, ehub_mixing,       &
                                   ehub_v, remove_ehub_u,        &
                                   wfcv, nwfcv, ofstv,           &
                                   write_ehub_uv, ehub_v_raw
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcv, iunhubv
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf

  IMPLICIT NONE

  REAL(DP), INTENT(OUT)    :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)

  TYPE hubbard_ene_4
     REAL(DP), ALLOCATABLE :: U(:,:,:,:)
     REAL(DP), ALLOCATABLE :: V(:,:,:,:)
     REAL(DP), ALLOCATABLE :: J(:,:,:,:)
  END TYPE hubbard_ene_4

  TYPE occmat_4
     REAL(DP), ALLOCATABLE :: O(:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:)
  END TYPE occmat_4

  TYPE occmat_5
     REAL(DP), ALLOCATABLE :: O(:,:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:,:)
  END TYPE occmat_5

  TYPE(bec_type)           :: proj

  TYPE(hubbard_ene_4)      :: E0(ehub_nnnc)
  TYPE(hubbard_ene_4)      :: HF(ehub_nnnc)

  TYPE(occmat_5)           :: n(ehub_nnnc)
  TYPE(occmat_5)           :: p(ehub_nnnc)
  TYPE(occmat_5)           :: num(nbnd)

  TYPE(occmat_4)           :: dm0
  TYPE(occmat_4)           :: dm1(ntyp,4)

  COMPLEX(DP)              :: kphase

  REAL(DP), ALLOCATABLE    :: ehub_v_(:,:,:,:)

  REAL(DP)                 :: arg
  REAL(DP)                 :: ntau(3,nat)
  REAL(DP)                 :: compute_eri
  REAL(DP)                 :: c1(6), c2(6), c3(6), c4(6)

  INTEGER                  :: npw, nw
  INTEGER                  :: i, ik, ibnd, is, ix
  INTEGER                  :: na, nt, m, l
  INTEGER                  :: na1, na2
  INTEGER                  :: nt1, nt2
  INTEGER                  :: m1, m2, m3, m4
  INTEGER                  :: l1, l2, l3, l4
  INTEGER                  :: n1, n2, n3
  INTEGER                  :: nnc, ics

  CALL start_clock('new_ns_e0_#2')

  CALL allocate_bec_type(nwfcv, nbnd, proj)
  !
  ALLOCATE ( dm0%O(ehub_lm_max,ehub_lm_max,nspin,nat) ) ; dm0%O = 0.0_DP
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        ALLOCATE ( dm1(nt,l)%O(ehub_lm_max,ehub_lm_max,nspin,nat) ) ; dm1(nt,l)%O = 0.0_DP
     ENDDO
  ENDDO
  !
  DO i = 1, ehub_nnnc
     ALLOCATE ( n(i)%O(ehub_lm_max,ehub_lm_max,nspin,nat,nat) ) ; n(i)%O = 0.0_DP
     ALLOCATE ( p(i)%O(ehub_lm_max,ehub_lm_max,nspin,nat,nat) ) ; p(i)%O = 0.0_DP
     ALLOCATE ( E0(i)%V(ehub_lm_max,ehub_lm_max,nat,nat) )
     ALLOCATE ( E0(i)%J(ehub_lm_max,ehub_lm_max,nat,nat) )
     ALLOCATE ( HF(i)%V(ehub_lm_max,ehub_lm_max,nat,nat) )
     ALLOCATE ( HF(i)%J(ehub_lm_max,ehub_lm_max,nat,nat) )
  ENDDO

  ics = ehub_csize + 1
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcv | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( nks > 1 ) CALL get_buffer(wfcv, nwordwfcv, iunhubv, ik)
     CALL calbec(npw, wfcv, evc, proj)
     !
     ! calculate the occupation matrices n
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2 ) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                         proj%r(ofstv(na2)+m2,ibnd) * &
                         proj%r(ofstv(na1)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                    DBLE(proj%k(ofstv(na2)+m2,ibnd) * &
                   CONJG(proj%k(ofstv(na1)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           !
           n(1)%O(m2,m1,current_spin,na2,na1) = n(1)%O(m1,m2,current_spin,na1,na2)
           !
        ENDDO ! m2
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
        arg = tpi * (                                              &
         xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
         xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
         xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
        !
        kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
           DO m1 = 1, ehub_lm(nt1)
           DO m2 = 1, ehub_lm(nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                             wg(ibnd,ik) * &
                              proj%r(ofstv(na2)+m2,ibnd) * &
                              proj%r(ofstv(na1)+m1,ibnd)
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                             wg(ibnd,ik) * &
                                           DBLE((kphase) * &
                              proj%k(ofstv(na2)+m2,ibnd) * &
                        CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
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
     ! calculate the renormalized occupations N(ACBN0)
     !
     DO i = 1, nbnd
        ALLOCATE ( num(i)%O(4,4,nspin,ntyp,ntyp) ) ; num(i)%O = 0.0_DP
        ALLOCATE ( num(i)%T(4,4,nspin,ntyp,ntyp) ) ; num(i)%T = 0.0_DP
     ENDDO
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m = 1, ehub_lm(nt)
           !
           l = ehub_l(m,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt,nt) + &
                            proj%r(ofstv(na)+m,ibnd) * &
                            proj%r(ofstv(na)+m,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt,nt) + &
                       DBLE(proj%k(ofstv(na)+m,ibnd) * &
                      CONJG(proj%k(ofstv(na)+m,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m
        !
     ENDDO ! na
     !
     DO ibnd = 1, nbnd
        !
        DO nt1 = 1, ntyp
        DO nt2 = 1, ntyp
           !
           IF ( upf(nt1)%psd /= upf(nt2)%psd ) CYCLE
           !
           num(ibnd)%T(:,:,current_spin,nt1,nt1) = &
           num(ibnd)%T(:,:,current_spin,nt1,nt1) + &
           num(ibnd)%O(:,:,current_spin,nt2,nt2)
           !
        ENDDO ! nt2
        ENDDO ! nt1
        !
     ENDDO ! ibnd
     !
     DO ibnd = 1, nbnd
        !
        DO nt1 = 1, ntyp
        DO nt2 = 1, ntyp
           !
           DO l1 = 1, 4
           DO l2 = 1, 4
              !
              num(ibnd)%T(l1,l2,current_spin,nt1,nt2) = &
              num(ibnd)%T(l1,l1,current_spin,nt1,nt1) + &
              num(ibnd)%T(l2,l2,current_spin,nt2,nt2)
              !
              IF ( upf(nt1)%psd == upf(nt2)%psd )       &
              num(ibnd)%T(l1,l2,current_spin,nt1,nt2) = &
              num(ibnd)%T(l1,l2,current_spin,nt1,nt2) / 2.0_DP
              !
           ENDDO ! l2
           ENDDO ! l1
           !
        ENDDO ! nt2
        ENDDO ! nt1
        !
     ENDDO ! ibnd
     !
     DO ibnd = 1, nbnd
        !
        num(ibnd)%O = num(ibnd)%T
        !
        DEALLOCATE ( num(ibnd)%T )
        !
     ENDDO
     !
     ! calculate the renormalized occupation matrices P(ACBN0)
     !
     ! 1) P^II (dm0, dm1)
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt)
           !
           l1 = ehub_l(m1,nt)
           !
        DO m2 = 1, ehub_lm(nt)
           !
           l2 = ehub_l(m2,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 dm0%O(m1,m2,current_spin,na) = &
                 dm0%O(m1,m2,current_spin,na) + &
                           num(ibnd)%O(l1,l2,   &
                          current_spin,nt,nt) * &
                                  wg(ibnd,ik) * &
                    proj%r(ofstv(na)+m2,ibnd) * &
                    proj%r(ofstv(na)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 dm0%O(m1,m2,current_spin,na) = &
                 dm0%O(m1,m2,current_spin,na) + &
                           num(ibnd)%O(l1,l2,   &
                          current_spin,nt,nt) * &
                                  wg(ibnd,ik) * &
               DBLE(proj%k(ofstv(na)+m2,ibnd) * &
              CONJG(proj%k(ofstv(na)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na

     DO nt = 1, ntyp
        !
        IF ( ehub_lm(nt) == 0 ) CYCLE
        !
        DO l = 1, 4
           !
           DO na1 = 1, nat
              nt1 = ityp(na1)
              !
              IF ( ehub_lm(nt1) == 0 ) CYCLE
              !
              DO m1 = 1, ehub_lm(nt1)
              DO m2 = 1, ehub_lm(nt1)
                 !
                 IF ( gamma_only ) THEN
                    DO ibnd = 1, nbnd
                       dm1(nt,l)%O(m1,m2,current_spin,na1) = &
                       dm1(nt,l)%O(m1,m2,current_spin,na1) + &
                             num(ibnd)%O(ehub_l(m1,nt1),l,   &
                                      current_spin,nt1,nt) * &
                                               wg(ibnd,ik) * &
                                proj%r(ofstv(na1)+m2,ibnd) * &
                                proj%r(ofstv(na1)+m1,ibnd)
                    ENDDO ! ibnd
                 ELSE
                    DO ibnd = 1, nbnd
                       dm1(nt,l)%O(m1,m2,current_spin,na1) = &
                       dm1(nt,l)%O(m1,m2,current_spin,na1) + &
                             num(ibnd)%O(ehub_l(m1,nt1),l,   &
                                      current_spin,nt1,nt) * &
                                               wg(ibnd,ik) * &
                           DBLE(proj%k(ofstv(na1)+m2,ibnd) * &
                          CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                    ENDDO ! ibnd
                 ENDIF
                 !
              ENDDO ! m2
              ENDDO ! m1
              !
           ENDDO ! na1
           !
        ENDDO ! l
        !
     ENDDO ! nt
     !
     ! 2) P^IJ
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
           !
           l1 = ehub_l(m1,nt1)
           !
        DO m2 = 1, ehub_lm(nt2)
           !
           l2 = ehub_l(m2,nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                   p(1)%O(m1,m2,current_spin,na1,na2) = &
                   p(1)%O(m1,m2,current_spin,na1,na2) + &
                                   num(ibnd)%O(l1,l2,   &
                                current_spin,nt1,nt2) * &
                                          wg(ibnd,ik) * &
                           proj%r(ofstv(na2)+m2,ibnd) * &
                           proj%r(ofstv(na1)+m1,ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                   p(1)%O(m1,m2,current_spin,na1,na2) = &
                   p(1)%O(m1,m2,current_spin,na1,na2) + &
                                   num(ibnd)%O(l1,l2,   &
                                current_spin,nt1,nt2) * &
                                          wg(ibnd,ik) * &
                      DBLE(proj%k(ofstv(na2)+m2,ibnd) * &
                     CONJG(proj%k(ofstv(na1)+m1,ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( ehub_lm(nt2) == 0 ) CYCLE
        !
        IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           p(1)%O(m2,m1,current_spin,na2,na1) = p(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
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
        arg = tpi * (                                              &
         xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
         xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
         xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
        !
        kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
           DO m1 = 1, ehub_lm(nt1)
              !
              l1 = ehub_l(m1,nt1)
              !
           DO m2 = 1, ehub_lm(nt2)
              !
              l2 = ehub_l(m2,nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    p(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    p(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                      num(ibnd)%O(l1,l2,   &
                                   current_spin,nt1,nt2) * &
                                             wg(ibnd,ik) * &
                              proj%r(ofstv(na2)+m2,ibnd) * &
                              proj%r(ofstv(na1)+m1,ibnd)
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    p(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    p(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                      num(ibnd)%O(l1,l2,   &
                                   current_spin,nt1,nt2) * &
                                             wg(ibnd,ik) * &
                                           DBLE((kphase) * &
                              proj%k(ofstv(na2)+m2,ibnd) * &
                        CONJG(proj%k(ofstv(na1)+m1,ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
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
     DO i = 1, nbnd
        DEALLOCATE ( num(i)%O )
     ENDDO

  ENDDO ! ik

  CALL deallocate_bec_type(proj)

  CALL mp_sum(dm0%O, inter_pool_comm)
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        CALL mp_sum(dm1(nt,l)%O, inter_pool_comm)
     ENDDO
  ENDDO
  !
  DO i = 1, ehub_nnnc
     CALL mp_sum(n(i)%O, inter_pool_comm)
     CALL mp_sum(p(i)%O, inter_pool_comm)
  ENDDO

  DO i = 1, ehub_nnnc
     IF ( nspin == 1 ) n(i)%O = 0.5_DP * n(i)%O
  ENDDO
  !
  CALL symmetrize_n5(n)
  !
  ns = 0.0_DP
  !
  DO i = 1, ehub_nnnc
     ns(:,:,:,:,:,i) = n(i)%O(:,:,:,:,:)
     DEALLOCATE ( n(i)%O )
  ENDDO
  !
  IF ( nspin == 1 ) dm0%O = 0.5_DP * dm0%O
  !
  CALL symmetrize_n4(dm0)
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        IF ( nspin == 1 ) dm1(nt,l)%O = 0.5_DP * dm1(nt,l)%O
     ENDDO
  ENDDO
  !
  CALL symmetrize_n4m(dm1)
  !
  DO i = 1, ehub_nnnc
     IF ( nspin == 1 ) p(i)%O = 0.5_DP * p(i)%O
  ENDDO
  !
  CALL symmetrize_n5(p)

  !
  ! on-site interaction: U_ij^II
  !
  E0(1)%V = 0.0_DP ; E0(1)%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, ehub_lm(nt)
        CALL orberi(ehub_l(m1,nt), l1)
        CALL pao_3g(upf(nt)%psd, ehub_l(m1,nt), c1)
        !
     l2 = 0
     DO m2 = 1, ehub_lm(nt)
        CALL orberi(ehub_l(m2,nt), l2)
        CALL pao_3g(upf(nt)%psd, ehub_l(m2,nt), c2)
        !
     l3 = 0
     DO m3 = 1, ehub_lm(nt)
        IF ( ehub_l(m3,nt) /= ehub_l(m1,nt) ) CYCLE
        CALL orberi(ehub_l(m3,nt), l3)
        CALL pao_3g(upf(nt)%psd, ehub_l(m3,nt), c3)
        !
     l4 = 0
     DO m4 = 1, ehub_lm(nt)
        IF ( ehub_l(m4,nt) /= ehub_l(m2,nt) ) CYCLE
        CALL orberi(ehub_l(m4,nt), l4)
        CALL pao_3g(upf(nt)%psd, ehub_l(m4,nt), c4)
        !
        IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,1,na)) * &
                             compute_eri(l3,l1,l4,l2,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c1,c4,c2)
        ELSE
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
           E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,1,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,1,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,1,na)) * &
                             compute_eri(l3,l1,l4,l2,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c1,c4,c2)

           E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
           E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,1,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,2,na)) * &
                             compute_eri(l3,l2,l4,l1,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c2,c4,c1)
        ENDIF
        !
     ENDDO ! m4
     ENDDO ! m3
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  HF(1)%V = 0.0_DP ; HF(1)%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt)
     DO m2 = 1, ehub_lm(nt)
        !
        IF ( m1 == m2 ) THEN
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,1,na,na,1)
        ELSE
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
            HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,1,na,na,1)

            HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
            HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,2,na,na,1)
        ENDIF
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  ALLOCATE ( ehub_v_(nat,nat,ehub_lm_max,ehub_lm_max) ) ; ehub_v_ = 0.0_DP
  IF ( .NOT.ALLOCATED(ehub_v_raw) ) THEN   ! [rawprint]
     ALLOCATE ( ehub_v_raw(SIZE(ehub_v,1),SIZE(ehub_v,2),SIZE(ehub_v,3),SIZE(ehub_v,4),SIZE(ehub_v,5)) )
     ehub_v_raw = 0.0_DP
  ENDIF
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt)
     DO m2 = 1, ehub_lm(nt)
        !
        IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
           ehub_v_(na,na,m1,m2) =                           &
               E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) / &
               HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) * &
                          remove_ehub_u(nt,ehub_l(m1,nt)) * &
                          remove_ehub_u(nt,ehub_l(m2,nt))
        ELSE
           ehub_v_(na,na,m1,m2) =                            &
              (E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  / &
               HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  - &
               E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  / &
               HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)) * &
                          remove_ehub_u(nt,ehub_l(m1,nt))  * &
                          remove_ehub_u(nt,ehub_l(m2,nt))
        ENDIF
        !
        ehub_v_raw(na,na,m1,m2,1) = ehub_v_(na,na,m1,m2)   ! [rawprint]
        !
        ehub_v(na,na,m1,m2,1)                                 &
            = (1.0_DP - ehub_mixing)*ehub_v (na,na,m1,m2,1) + &
              (         ehub_mixing)*ehub_v_(na,na,m1,m2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na

  !
  ! inter-site interaction: U_ij^IJ (I /= J)
  !
  ! ... in the unit cell
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, ehub_lm(nt1)
        CALL orberi(ehub_l(m1,nt1), l1)
        CALL pao_3g(upf(nt1)%psd, ehub_l(m1,nt1), c1)
        !
     l2 = 0
     DO m2 = 1, ehub_lm(nt2)
        CALL orberi(ehub_l(m2,nt2), l2)
        CALL pao_3g(upf(nt2)%psd, ehub_l(m2,nt2), c2)
        !
     l3 = 0
     DO m3 = 1, ehub_lm(nt1)
        IF ( ehub_l(m3,nt1) /= ehub_l(m1,nt1) ) CYCLE
        CALL orberi(ehub_l(m3,nt1), l3)
        CALL pao_3g(upf(nt1)%psd, ehub_l(m3,nt1), c3)
        !
     l4 = 0
     DO m4 = 1, ehub_lm(nt2)
        IF ( ehub_l(m4,nt2) /= ehub_l(m2,nt2) ) CYCLE
        CALL orberi(ehub_l(m4,nt2), l4)
        CALL pao_3g(upf(nt2)%psd, ehub_l(m4,nt2), c4)
        !
        E0(1)%V(m1,m2,na1,na2) = E0(1)%V(m1,m2,na1,na2)  + &
                (dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                 dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  + &
                 dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                 dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                 dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                 dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                 dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                 dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  - &
                                p(1)%O(m3,m2,1,na1,na2)  * &
                                p(1)%O(m4,m1,1,na2,na1)  - &
                                p(1)%O(m3,m2,2,na1,na2)  * &
                                p(1)%O(m4,m1,2,na2,na1)) * &
                               compute_eri(l3,l1,l4,l2,    &
                       tau(:,na1)*alat,tau(:,na1)*alat,    &
                       tau(:,na2)*alat,tau(:,na2)*alat,    &
                                           c3,c1,c4,c2)
        !
     ENDDO ! m4
     ENDDO ! m3
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        !
        E0(1)%V(m2,m1,na2,na1) = E0(1)%V(m1,m2,na1,na2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  HF(1)%V = 0.0_DP
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        !
        HF(1)%V(m1,m2,na1,na2) = HF(1)%V(m1,m2,na1,na2) + &
         ns(m1,m1,1,na1,na1,1) * ns(m2,m2,1,na2,na2,1)  + &
         ns(m1,m1,2,na1,na1,1) * ns(m2,m2,2,na2,na2,1)  + &
         ns(m1,m1,1,na1,na1,1) * ns(m2,m2,2,na2,na2,1)  + &
         ns(m1,m1,2,na1,na1,1) * ns(m2,m2,1,na2,na2,1)  - &
         ns(m1,m2,1,na1,na2,1) * ns(m2,m1,1,na2,na1,1)  - &
         ns(m1,m2,2,na1,na2,1) * ns(m2,m1,2,na2,na1,1)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        HF(1)%V(m2,m1,na2,na1) = HF(1)%V(m1,m2,na1,na2)
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
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
     IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
        !
        l1 = ehub_l(m1,nt1)
        !
     DO m2 = 1, ehub_lm(nt2)
        !
        l2 = ehub_l(m2,nt2)
        !
        ehub_v_(na1,na2,m1,m2) = &
               (E0(1)%V(m1,m2,na1,na2) / HF(1)%V(m1,m2,na1,na2)) / 2.0_DP
        !
        ehub_v_raw(na1,na2,m1,m2,1) = ehub_v_(na1,na2,m1,m2)   ! [rawprint]
        !
        ehub_v(na1,na2,m1,m2,1)                                 &
            = (1.0_DP - ehub_mixing)*ehub_v (na1,na2,m1,m2,1) + &
              (         ehub_mixing)*ehub_v_(na1,na2,m1,m2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  ! ... out of the unit cell
  !
  DO i = 2, ehub_nnnc
     E0(i)%V = 0.0_DP
  ENDDO
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
        DO ix = 1, 3
           ntau(ix,na2) = n1 * at(ix,1) + n2 * at(ix,2) + n3 * at(ix,3) + tau(ix,na2)
        ENDDO
        !
        l1 = 0
        DO m1 = 1, ehub_lm(nt1)
           CALL orberi(ehub_l(m1,nt1), l1)
           CALL pao_3g(upf(nt1)%psd, ehub_l(m1,nt1), c1)
           !
        l2 = 0
        DO m2 = 1, ehub_lm(nt2)
           CALL orberi(ehub_l(m2,nt2), l2)
           CALL pao_3g(upf(nt2)%psd, ehub_l(m2,nt2), c2)
           !
        l3 = 0
        DO m3 = 1, ehub_lm(nt1)
           IF ( ehub_l(m3,nt1) /= ehub_l(m1,nt1) ) CYCLE
           CALL orberi(ehub_l(m3,nt1), l3)
           CALL pao_3g(upf(nt1)%psd, ehub_l(m3,nt1), c3)
           !
        l4 = 0
        DO m4 = 1, ehub_lm(nt2)
           IF ( ehub_l(m4,nt2) /= ehub_l(m2,nt2) ) CYCLE
           CALL orberi(ehub_l(m4,nt2), l4)
           CALL pao_3g(upf(nt2)%psd, ehub_l(m4,nt2), c4)
           !
           E0(nnc)%V(m1,m2,na1,na2) = E0(nnc)%V(m1,m2,na1,na2)  + &
                       (dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                        dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  + &
                        dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                        dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                        dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,na1)  * &
                        dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,na2)  + &
                        dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,na1)  * &
                        dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,na2)  - &
                                     p(nnc)%O(m3,m2,1,na1,na2)  * &
                                     p(nnc)%O(m1,m4,1,na1,na2)  - &
                                     p(nnc)%O(m3,m2,2,na1,na2)  * &
                                     p(nnc)%O(m1,m4,2,na1,na2)) * &
                                      compute_eri(l3,l1,l4,l2,    &
                             tau(:,na1)*alat, tau(:,na1)*alat,    &
                            ntau(:,na2)*alat,ntau(:,na2)*alat,    &
                                                  c3,c1,c4,c2)
           !
        ENDDO ! m4
        ENDDO ! m3
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
  DEALLOCATE ( dm0%O )
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        DEALLOCATE ( dm1(nt,l)%O )
     ENDDO
  ENDDO
  !
  DO i = 2, ehub_nnnc
     HF(i)%V = 0.0_DP
  ENDDO
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
        DO m1 = 1, ehub_lm(nt1)
           !
           l1 = ehub_l(m1,nt1)
           !
        DO m2 = 1, ehub_lm(nt2)
           !
           l2 = ehub_l(m2,nt2)
           !
           HF(nnc)%V(m1,m2,na1,na2) = HF(nnc)%V(m1,m2,na1,na2) + &
               ns(m1,m1,1,na1,na1,  1)*ns(m2,m2,1,na2,na2,  1) + &
               ns(m1,m1,2,na1,na1,  1)*ns(m2,m2,2,na2,na2,  1) + &
               ns(m1,m1,1,na1,na1,  1)*ns(m2,m2,2,na2,na2,  1) + &
               ns(m1,m1,2,na1,na1,  1)*ns(m2,m2,1,na2,na2,  1) - &
               ns(m1,m2,1,na1,na2,nnc)*ns(m1,m2,1,na1,na2,nnc) - &
               ns(m1,m2,2,na1,na2,nnc)*ns(m1,m2,2,na1,na2,nnc)
           !
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
        DO m1 = 1, ehub_lm(nt1)
           !
           l1 = ehub_l(m1,nt1)
           !
        DO m2 = 1, ehub_lm(nt2)
           !
           l2 = ehub_l(m2,nt2)
           !
           ehub_v_(na1,na2,m1,m2) = &
                  (E0(nnc)%V(m1,m2,na1,na2) / HF(nnc)%V(m1,m2,na1,na2)) / 2.0_DP
           !
           ehub_v_raw(na1,na2,m1,m2,nnc) = ehub_v_(na1,na2,m1,m2)   ! [rawprint]
           !
           ehub_v(na1,na2,m1,m2,nnc) =                               &
                 (1.0_DP - ehub_mixing)*ehub_v (na1,na2,m1,m2,nnc) + &
                 (         ehub_mixing)*ehub_v_(na1,na2,m1,m2    )
           !
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
  DEALLOCATE ( ehub_v_ )
  !
  DO i = 1, ehub_nnnc
     DEALLOCATE ( E0(i)%V, E0(i)%J )
     DEALLOCATE ( HF(i)%V, HF(i)%J )
  ENDDO
  !
  ! Hermiticity of Hubbard parameters U & V
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        ehub_v(na1,na2,m1,m2,1) = 0.5D0 * (ehub_v(na1,na2,m1,m2,1) + &
                                           ehub_v(na2,na1,m2,m1,1))
        !
        ehub_v(na2,na1,m2,m1,1) = ehub_v(na1,na2,m1,m2,1)
        ! [rawprint]
        ehub_v_raw(na1,na2,m1,m2,1) = 0.5D0 * (ehub_v_raw(na1,na2,m1,m2,1) + ehub_v_raw(na2,na1,m2,m1,1))
        ehub_v_raw(na2,na1,m2,m1,1) = ehub_v_raw(na1,na2,m1,m2,1)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
  DO na2 = 1, nat
     IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) THEN
        ehub_v(na1,na2,:,:,1) = 0.0_DP
        ehub_v_raw(na1,na2,:,:,1) = 0.0_DP   ! [rawprint]
     ENDIF
  ENDDO ! na2
  ENDDO ! na1

  CALL write_ehub_uv

  CALL stop_clock('new_ns_e0_#2')

  RETURN

END SUBROUTINE new_ns_e0_2


SUBROUTINE symmetrize_n4( n )

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE ldaV,                 ONLY : ehub_lm, ehub_l, ehub_lm_max, &
                                   ehub_nn_vect, ehub_nn_dim,    &
                                   ehub_csize, ehub_nnnc
  USE lsda_mod,             ONLY : nspin
  USE symm_base,            ONLY : d1, d2, d3, nsym, irt

  IMPLICIT NONE

  TYPE occmat_4
     REAL(DP), ALLOCATABLE :: O(:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:)
  END TYPE occmat_4

  TYPE(occmat_4), INTENT(INOUT) :: n

  INTEGER                       :: is, isym, l, m, nt
  INTEGER                       :: na1, na2
  INTEGER                       :: nt1, nt2
  INTEGER                       :: m1, m2, m3, m4
  INTEGER                       :: m1_, m2_, m3_, m4_
  INTEGER                       :: temp, l_separation(5,ntyp)

  ALLOCATE (n%T(ehub_lm_max,ehub_lm_max,nspin,nat)) ; n%T = 0.0_DP

  DO nt = 1, ntyp
     !
     temp = ehub_l(1,nt)
     l_separation(ehub_l(1,nt),nt) = 1
     !
     DO m = 2, ehub_lm(nt)
        !
        IF ( ehub_l(m,nt) /= temp ) THEN
           temp = ehub_l(m,nt)
           l_separation(ehub_l(m,nt),nt) = m
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  DO is = 1, nspin
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
           m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
           !
        DO m2 = 1, ehub_lm(nt1)
           m2_ = m2 - l_separation(ehub_l(m2,nt1),nt1) + 1
           !
           DO isym = 1, nsym
              !
              na2 = irt(isym, na1)
              !
              DO m3 = 1, ehub_lm(nt1)
                 IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                 m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                 !
              DO m4 = 1, ehub_lm(nt1)
                 IF ( ehub_l(m2,nt1) /= ehub_l(m4,nt1) ) CYCLE
                 m4_ = m4 - l_separation(ehub_l(m4,nt1),nt1) + 1
                 !
                 IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)
  
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)

                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)

                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE
                    CALL errore('new_ns_v','angular momentum not implemented',1)
                 ENDIF
                 !
              ENDDO ! m4
              ENDDO ! m3
              !
           ENDDO ! isym
           !
        ENDDO ! m2
        ENDDO ! m1
        ! 
     ENDDO ! na1
     !
  ENDDO ! is
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
     DO is = 1, nspin
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = m1, ehub_lm(nt1)
           !
           n%T(m1,m2,is,na1) = 0.5D0 * (n%T(m1,m2,is,na1) + &
                                        n%T(m2,m1,is,na1))
           !
           n%T(m2,m1,is,na1) = n%T(m1,m2,is,na1)
           !
        ENDDO
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  n%O = n%T
  DEALLOCATE (n%T)

END SUBROUTINE symmetrize_n4


SUBROUTINE symmetrize_n4m( n )

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE ldaV,                 ONLY : ehub_lm, ehub_l, ehub_lm_max, &
                                   ehub_nn_vect, ehub_nn_dim,    &
                                   ehub_csize, ehub_nnnc
  USE lsda_mod,             ONLY : nspin
  USE symm_base,            ONLY : d1, d2, d3, nsym, irt

  IMPLICIT NONE

  TYPE occmat_4
     REAL(DP), ALLOCATABLE :: O(:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:)
  END TYPE occmat_4

  TYPE(occmat_4), INTENT(INOUT) :: n(ntyp,4)

  INTEGER                       :: is, isym, i, l, m, nt
  INTEGER                       :: na1, na2
  INTEGER                       :: nt1, nt2
  INTEGER                       :: m1, m2, m3, m4
  INTEGER                       :: m1_, m2_, m3_, m4_
  INTEGER                       :: temp, l_separation(5,ntyp)

  DO i = 1, ntyp
     DO l = 1, 4
        ALLOCATE (n(i,l)%T(ehub_lm_max,ehub_lm_max,nspin,nat)) ; n(i,l)%T = 0.0_DP
     ENDDO
  ENDDO
  !
  DO nt = 1, ntyp
     !
     temp = ehub_l(1,nt)
     l_separation(ehub_l(1,nt),nt) = 1
     !
     DO m = 2, ehub_lm(nt)
        !
        IF ( ehub_l(m,nt) /= temp ) THEN
           temp = ehub_l(m,nt)
           l_separation(ehub_l(m,nt),nt) = m
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  DO nt = 1, ntyp
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO l = 1, 4
        !
        DO is = 1, nspin
           !
           DO na1 = 1, nat
              nt1 = ityp(na1)
              !
              IF ( ehub_lm(nt1) == 0 ) CYCLE
              !
              DO m1 = 1, ehub_lm(nt1)
                 m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
                 !
              DO m2 = 1, ehub_lm(nt1)
                 m2_ = m2 - l_separation(ehub_l(m2,nt1),nt1) + 1
                 !
                 DO isym = 1, nsym
                    !
                    na2 = irt(isym, na1)
                    !
                    DO m3 = 1, ehub_lm(nt1)
                       IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                       m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                       !
                    DO m4 = 1, ehub_lm(nt1)
                       IF ( ehub_l(m2,nt1) /= ehub_l(m4,nt1) ) CYCLE
                       m4_ = m4 - l_separation(ehub_l(m4,nt1),nt1) + 1
                       !
                       IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)
      
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)

                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)

                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE
                          CALL errore('new_ns_v','angular momentum not implemented',1)
                       ENDIF
                       !
                    ENDDO ! m4
                    ENDDO ! m3
                    !
                 ENDDO ! isym
                 !
              ENDDO ! m2
              ENDDO ! m1
              ! 
           ENDDO ! na1
           !
        ENDDO ! is
        !
     ENDDO ! l
     !
  ENDDO ! nt
  !
  DO nt = 1, ntyp
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO l = 1, 4
        !
        DO na1 = 1, nat
           nt1 = ityp(na1)
           !
           IF ( ehub_lm(nt1) == 0 ) CYCLE
           !
           DO is = 1, nspin
              !
              DO m1 = 1, ehub_lm(nt1)
              DO m2 = m1, ehub_lm(nt1)
                 !
                 n(nt,l)%T(m1,m2,is,na1) = 0.5D0 * (n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%T(m2,m1,is,na1))
                 !
                 n(nt,l)%T(m2,m1,is,na1) = n(nt,l)%T(m1,m2,is,na1)
                 !
              ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        n(nt,l)%O = n(nt,l)%T
        DEALLOCATE (n(nt,l)%T)
     ENDDO
  ENDDO

END SUBROUTINE symmetrize_n4m


SUBROUTINE symmetrize_n5( n )

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE ldaV,                 ONLY : ehub_lm, ehub_l, ehub_lm_max, &
                                   ehub_nn_vect, ehub_nn_dim,    &
                                   ehub_csize, ehub_nnnc
  USE lsda_mod,             ONLY : nspin
  USE symm_base,            ONLY : d1, d2, d3, nsym, irt

  IMPLICIT NONE

  TYPE occmat_5
     REAL(DP), ALLOCATABLE :: O(:,:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:,:)
  END TYPE occmat_5

  TYPE(occmat_5), INTENT(INOUT) :: n(ehub_nnnc)

  INTEGER                       :: is, isym, i, m, nt
  INTEGER                       :: na1, na2, na3, na4
  INTEGER                       :: nt1, nt2, nt3, nt4
  INTEGER                       :: m1, m2, m3, m4
  INTEGER                       :: m1_, m2_, m3_, m4_
  INTEGER                       :: n0(3), n1, n2, n3
  INTEGER                       :: nnc0, nnc, ics
  INTEGER                       :: temp, l_separation(5,ntyp)

  DO i = 1, ehub_nnnc
     ALLOCATE (n(i)%T(ehub_lm_max,ehub_lm_max,nspin,nat,nat)) ; n(i)%T = 0.0_DP
  ENDDO
  !
  DO nt = 1, ntyp
     !
     temp = ehub_l(1,nt)
     l_separation(ehub_l(1,nt),nt) = 1
     !
     DO m = 2, ehub_lm(nt)
        !
        IF ( ehub_l(m,nt) /= temp ) THEN
           temp = ehub_l(m,nt)
           l_separation(ehub_l(m,nt),nt) = m
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ics = ehub_csize + 1
  !
  DO is = 1, nspin
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
        DO m1 = 1, ehub_lm(nt1)
           m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
           !
        DO m2 = 1, ehub_lm(nt1)
           m2_ = m2 - l_separation(ehub_l(m2,nt1),nt1) + 1
           !
           DO isym = 1, nsym
              !
              na2 = irt(isym, na1)
              !
              DO m3 = 1, ehub_lm(nt1)
                 IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                 m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                 !
              DO m4 = 1, ehub_lm(nt1)
                 IF ( ehub_l(m2,nt1) /= ehub_l(m4,nt1) ) CYCLE
                 m4_ = m4 - l_separation(ehub_l(m4,nt1),nt1) + 1
                 !
                 IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)

                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)

                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)

                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE
                    CALL errore('new_ns_v','angular momentum not implemented',1)
                 ENDIF
                 !
              ENDDO ! m4
              ENDDO ! m3
              !
           ENDDO ! isym
           !
        ENDDO ! m2
        ENDDO ! m1
        ! 
     ENDDO ! na1
     !
  ENDDO ! is

  DO is = 1, nspin
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( ehub_lm(nt1) == 0 ) CYCLE
        !
        DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
        DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
        DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
           !
           n0(1) = n1 ; n0(2) = n2 ; n0(3) = n3
           !
           CALL hubbard_get_nnc(n0,nnc0)
           !
           DO na2 = 1, nat
              nt2 = ityp(na2)
              !
              IF ( ehub_lm(nt2) == 0 ) CYCLE
              !
              IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
              !
              DO m1 = 1, ehub_lm(nt1)
                 m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
                 !
              DO m2 = 1, ehub_lm(nt2)
                 m2_ = m2 - l_separation(ehub_l(m2,nt2),nt2) + 1
                 !
                 DO isym = 1, nsym
                    !
                    CALL hubbard_symmetry(na1,na2,na3,na4,n0,nnc,isym)
                    !
                    DO m3 = 1, ehub_lm(nt1)
                       IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                       m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                       !
                    DO m4 = 1, ehub_lm(nt2)
                       IF ( ehub_l(m2,nt2) /= ehub_l(m4,nt2) ) CYCLE
                       m4_ = m4 - l_separation(ehub_l(m4,nt2),nt2) + 1
                       !
                       IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)

                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)

                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                            (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)

                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE
                          CALL errore('new_ns_v','angular momentum not implemented',1)
                       ENDIF
                       !
                    ENDDO ! m4
                    ENDDO ! m3
                    !
                 ENDDO ! isym
                 !
              ENDDO ! m2
              ENDDO ! m1
              !
           ENDDO ! na2
           ! 
        ENDDO ! n3
        ENDDO ! n2
        ENDDO ! n1
        !
     ENDDO ! na1
     !
  ENDDO ! is
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( ehub_lm(nt2) == 0 ) CYCLE
     !
     DO m1 = 1, ehub_lm(nt1)
     DO m2 = 1, ehub_lm(nt2)
        !
        DO is = 1, nspin
           n(1)%T(m1,m2,is,na1,na2) = 0.5D0 * (n(1)%T(m1,m2,is,na1,na2) + &
                                               n(1)%T(m2,m1,is,na2,na1))
           !
           n(1)%T(m2,m1,is,na2,na1) = n(1)%T(m1,m2,is,na1,na2)
        ENDDO
        !
     ENDDO
     ENDDO
     !
  ENDDO
  ENDDO
  !
  DO na1 = 1, nat
  DO na2 = 1, nat
     !
     DO is = 1, nspin
        IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) &
           n(1)%T(:,:,is,na1,na2) = 0.0_DP
     ENDDO
     !
  ENDDO
  ENDDO
  !
  DO i = 1, ehub_nnnc
     n(i)%O = n(i)%T
     DEALLOCATE (n(i)%T)
  ENDDO

END SUBROUTINE symmetrize_n5


SUBROUTINE orberi(ehub_l, eri_l)

  ! ehub_l = 1: s orbital [eri_l = 1(s)]
  ! ehub_l = 2: p orbital [eri_l = 2(pz), 3(px), 4(py)]
  ! ehub_l = 3: d orbital [eri_l = 5(dz^2-r^2), 6(dxz), 7(dyz), 8(dx^2-y^2),
  ! 9(dxy)]

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: ehub_l
  INTEGER, INTENT(INOUT) :: eri_l

  IF (ehub_l == 1) THEN
     eri_l = 1
  ELSEIF (ehub_l == 2) THEN
     IF (eri_l == 2) THEN
        eri_l = 3
     ELSEIF (eri_l == 3) THEN
        eri_l = 4
     ELSE
        eri_l = 2
     ENDIF
  ELSEIF (ehub_l == 3) THEN
     IF (eri_l == 5) THEN
        eri_l = 6
     ELSEIF (eri_l == 6) THEN
        eri_l = 7
     ELSEIF (eri_l == 7) THEN
        eri_l = 8
     ELSEIF (eri_l == 8) THEN
        eri_l = 9
     ELSE
        eri_l = 5
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE orberi
!
SUBROUTINE new_ns_v_nc(ns)

   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi
   USE ions_base,            ONLY : nat, ityp, tau, ntyp => nsp
   USE cell_base,            ONLY : at
   USE klist,                ONLY : nks, ngk, xk
   USE ldaV,                 ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                                    ehub_nn_vect, ehub_nn_dim,    &
                                    ehub_nnnc, ehub_csize, wfcv,  &
                                    nwfcv, ofstv, counter_lm
   USE ldaU,                 ONLY : Hubbard_projectors
   USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
   USE noncollin_module,     ONLY : npol, noncolin, domag, lspinorb
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE io_global,            ONLY : stdout
   USE wavefunctions,        ONLY : evc
   USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcv, iunhubv
   USE buffers,              ONLY : get_buffer
   USE mp_pools,             ONLY : inter_pool_comm
   USE mp,                   ONLY : mp_sum
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE becmod,               ONLY : bec_type, calbec,  &
                                    allocate_bec_type, &
                                    deallocate_bec_type
   USE uspp_param,           ONLY : upf
 
   IMPLICIT NONE
 
   COMPLEX(DP), INTENT(OUT) :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
 
   TYPE occmat_5_nc
      COMPLEX(DP), ALLOCATABLE :: O(:,:,:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: T(:,:,:,:,:,:)
   END TYPE occmat_5_nc
  
   COMPlEX(DP)           :: kphase
   COMPLEX(DP), ALLOCATABLE           :: proj(:,:)
 
   REAL(DP)              :: arg, psum
 
   TYPE(occmat_5_nc)        :: n(ehub_nnnc)
 
   INTEGER               :: ik, ibnd, is1, is2, i, j, k, npw, npw_
   INTEGER               :: na1, na2, nt1, nt2
   INTEGER               :: m1, m2, m3, m4, m1_, m2_
   INTEGER               :: n1, n2, n3, nnc, ics

   LOGICAL               :: save_flag
 
   CALL start_clock('new_ns_v_nc')
   !
   ALLOCATE( proj(nwfcv,nbnd))
   !
   proj(:,:) = (0.0_DP,0.0_DP)
   !
   DO i = 1, ehub_nnnc
      ALLOCATE ( n(i)%O(ehub_lm_max,ehub_lm_max,npol,npol,nat,nat) ) ; n(i)%O = 0.0_DP
   ENDDO
   !
   ics = ehub_csize + 1
   !
   DO ik = 1, nks
      !
      npw = ngk(ik)
      npw_ = npwx * npol
      !
      ! projection < wfcv | evc >
      !
      IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)     
      IF ( Hubbard_projectors == 'pseudo' ) THEN
         CALL errore('new_ns_v_nc', 'Hubbard_projectors = pseudo is not supported', 1)
      ELSE
         IF ( nks > 1 ) CALL get_buffer(wfcv, nwordwfcv, iunhubv, ik)
         !
         ! Refer to "projwave" in projwfc.f90
         !
         CALL calbec(npw_, wfcv, evc, proj)
      ENDIF
      !
      ! WIY
      !
      !
      ! calculate the occupation matrices n
      !
      DO na1 = 1, nat
         nt1 = ityp(na1)
         !
         IF ( ehub_lm(nt1) == 0 ) CYCLE
         !
      DO na2 = na1, nat
         nt2 = ityp(na2)
         !
         IF ( ehub_lm(nt2) == 0 ) CYCLE
         !
         IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
         !
         DO m1 = 1, ehub_lm(nt1)
            !
         DO m2 = 1, ehub_lm(nt2)
            !
            DO is1 = 1, npol
            DO is2 = 1, npol
               !
               DO ibnd = 1, nbnd
                n(1)%O(m1,m2,is1,is2,na1,na2) = &
                n(1)%O(m1,m2,is1,is2,na1,na2) + &
                                    wg(ibnd,ik) * &
                        proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd)* &
                  CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))
               ENDDO
               !
            ENDDO
            ENDDO
            !
         ENDDO ! m2
         ENDDO ! m1
         !
      ENDDO ! na2
      ENDDO ! na1
      !
      DO na1 = 1, nat
         nt1 = ityp(na1)
         !
         IF ( ehub_lm(nt1) == 0 ) CYCLE
         !
      DO na2 = na1, nat
         nt2 = ityp(na2)
         !
         IF ( ehub_lm(nt2) == 0 ) CYCLE
         !
         IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
         !
         DO m1 = 1, ehub_lm(nt1)
         DO m2 = 1, ehub_lm(nt2)
            DO is1 = 1, npol
            DO is2 = 1, npol
                n(1)%O(m2,m1,is2,is1,na2,na1) = CONJG(n(1)%O(m1,m2,is1,is2,na1,na2))
            ENDDO !is2
            ENDDO !is1
         ENDDO ! m2
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
         arg = tpi * (                                              &
          xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
          xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
          xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
         !
         kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
            IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
            !
            DO m1 = 1, ehub_lm(nt1)
               !
            DO m2 = 1, ehub_lm(nt2)
               !
               DO is1 = 1, npol
               DO is2 = 1, npol
                  !
                  DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,is1,is2,na1,na2) = &
                    n(nnc)%O(m1,m2,is1,is2,na1,na2) + &
                        kphase*wg(ibnd,ik) * &
                        proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) * &
                  CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))
                  ENDDO ! ibnd
                  !
               ENDDO 
               ENDDO
               !
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
   ENDDO ! ik
   !
   DO i = 1, ehub_nnnc
      CALL mp_sum(n(i)%O, inter_pool_comm)
   ENDDO
   !
   CALL symmetrize_n5_nc(n)
   !
   ns = 0.0_DP
   !
   DO i = 1, ehub_nnnc
    !
    DO is1 = 1, npol
    DO is2 = 1, npol
        j = npol*(is1-1) + is2
        ns(:,:,j,:,:,i) = n(i)%O(:,:,is1,is2,:,:)
    ENDDO
    ENDDO
    !
   ENDDO ! i
   !   
   ! checking ns Hermitean
   !
   !! temp
   !!DO i = 1, ehub_nnnc
   !   DO na1 = 1, nat
   !      nt1 = ityp(na1)
   !      !
   !      IF ( ehub_lm(nt1) == 0 ) CYCLE
   !      !
   !      DO na2 = 1, nat
   !         nt2 = ityp(na2)
   !         !
   !         IF ( ehub_lm(nt2) == 0 ) CYCLE
   !         !
   !         !IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
   !         !
   !         DO is1 = 1, npol
   !         DO is2 = 1, npol
   !            j = npol * (is1 - 1) + is2
   !            k = is1 + npol * (is2 - 1)
   !            DO m2 = 1, ehub_lm(nt2)
   !            DO m1 = 1, ehub_lm(nt1)
   !               !psum = ABS( ns(m1,m2,j,na1,na2,i) - CONJG(ns(m2,m1,k,na2,na1,i)))
   !               !IF (psum > 1.d-10) CALL errore( 'new_ns_v_nc', 'non Hermitean matrix', 1)
   !            ENDDO !m2
   !            ENDDO !m1
   !         ENDDO !is2
   !         ENDDO !is1
   !      ENDDO !na2
   !   ENDDO !na1
   !!ENDDO !i
   !! temp
   DO i = 1, ehub_nnnc
      DEALLOCATE ( n(i)%O )
   ENDDO
   !
   DEALLOCATE(proj)
   !
   CALL stop_clock('new_ns_v_nc')
   !
   RETURN
   !
END SUBROUTINE new_ns_v_nc
!
SUBROUTINE new_ns_e0_nc_1(ns)

   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi, rytoev
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
   USE cell_base,            ONLY : at, alat
   USE klist,                ONLY : nks, ngk, xk
   USE ldaV,                 ONLY : ehub_lm_max, ehub_lm, ehub_l,                      &
                                    ehub_nn_vect, ehub_nn_leng, ehub_nn_dim,           &
                                    ehub_csize, ehub_nnnc, ehub_mixing, ehub_conv_thr, &
                                    ehub_v, ehub_v_diff, remove_ehub_u, write_ehub_uv, &
                                    wfcv, nwfcv, ofstv, stop_acbn0, counter_lm, ehub_v_raw
   USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
   USE noncollin_module,     ONLY : npol, noncolin, domag, lspinorb
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE io_global,            ONLY : stdout
   USE control_flags,        ONLY : gamma_only
   USE wavefunctions,        ONLY : evc
   USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcv, iunhubv
   USE buffers,              ONLY : get_buffer
   USE mp_pools,             ONLY : inter_pool_comm
   USE mp,                   ONLY : mp_sum
   USE becmod,               ONLY : bec_type, calbec,  &
                                    allocate_bec_type, &
                                    deallocate_bec_type
   USE uspp_param,           ONLY : upf
 
   IMPLICIT NONE
 
   COMPLEX(DP), INTENT(OUT)    :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
 
   TYPE hubbard_ene_4
      REAL(DP), ALLOCATABLE :: U(:,:,:,:)
      REAL(DP), ALLOCATABLE :: V(:,:,:,:)
      REAL(DP), ALLOCATABLE :: J(:,:,:,:)
   END TYPE hubbard_ene_4
 
   TYPE occmat_4_nc
      COMPLEX(DP), ALLOCATABLE :: O(:,:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: T(:,:,:,:,:)
   END TYPE occmat_4_nc
 
   TYPE occmat_5_nc
      COMPLEX(DP), ALLOCATABLE :: O(:,:,:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: T(:,:,:,:,:,:)
   END TYPE occmat_5_nc
   
   TYPE occmat_4
      REAL(DP), ALLOCATABLE    :: O(:,:,:,:)
      REAL(DP), ALLOCATABLE    :: T(:,:,:,:)
   END TYPE occmat_4

   COMPLEX(DP), ALLOCATABLE    :: proj(:,:)
 
   TYPE(hubbard_ene_4)      :: E0(ehub_nnnc)
   TYPE(hubbard_ene_4)      :: HF(ehub_nnnc)
 
   TYPE(occmat_5_nc)           :: n(ehub_nnnc)
   TYPE(occmat_5_nc)           :: p(ehub_nnnc)
   ! temp 
   TYPE(occmat_4)           :: num(nbnd)
   ! temp
   TYPE(occmat_4_nc)           :: dm0
   TYPE(occmat_4_nc)           :: dm1(ntyp,4)
 
   COMPLEX(DP)              :: kphase
   COMPLEX(DP)              :: n11, n12, n21, n22
 
   REAL(DP), ALLOCATABLE    :: ehub_v_(:,:,:,:)
 
   REAL(DP)                 :: arg
   REAL(DP)                 :: ntau(3,nat)
   REAL(DP)                 :: compute_eri
   REAL(DP)                 :: c1(6), c2(6), c3(6), c4(6)
 
   INTEGER                  :: npw, nw, npw_
   INTEGER                  :: i, j, ik, ibnd, is, ix
   INTEGER                  :: na, nt, m, l
   INTEGER                  :: na1, na2
   INTEGER                  :: nt1, nt2
   INTEGER                  :: m1, m2, m3, m4
   INTEGER                  :: l1, l2, l3, l4
   INTEGER                  :: n1, n2, n3
   INTEGER                  :: nnc, ics
   INTEGER                  :: is1, is2
 
   CALL start_clock('new_ns_e0_nc')
 
   IF ( stop_acbn0 ) THEN
      !
      IF ( MAXVAL(ehub_v_diff(:,:,:,:,:)) < SQRT(ehub_conv_thr) ) THEN
         CALL new_ns_v_nc(ns)
         GOTO 1
      ENDIF
      !
      ehub_v_diff = 0.0_DP
      !
   ENDIF
 
   ALLOCATE( proj(nwfcv,nbnd))
   !
   ALLOCATE ( dm0%O(ehub_lm_max,ehub_lm_max,npol,npol,nat) ) ; dm0%O = 0.0_DP
   !
   DO nt = 1, ntyp
      DO l = 1, 4
         ALLOCATE ( dm1(nt,l)%O(ehub_lm_max,ehub_lm_max,npol,npol,nat) ) ; dm1(nt,l)%O = 0.0_DP
      ENDDO
   ENDDO
   !
   DO i = 1, ehub_nnnc
      ALLOCATE ( n(i)%O(ehub_lm_max,ehub_lm_max,npol,npol,nat,nat) ) ; n(i)%O = 0.0_DP
      ALLOCATE ( p(i)%O(ehub_lm_max,ehub_lm_max,npol,npol,nat,nat) ) ; p(i)%O = 0.0_DP
      ALLOCATE ( E0(i)%V(4,4,nat,nat), E0(i)%J(4,4,nat,nat) )
      ALLOCATE ( HF(i)%V(4,4,nat,nat), HF(i)%J(4,4,nat,nat) )
   ENDDO
 
   ics = ehub_csize + 1
   !
   DO ik = 1, nks
      !
      npw = ngk(ik)
      npw_ = npwx * npol
      !
      ! projection < wfcv | evc >
      !
      IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      IF ( nks > 1 ) CALL get_buffer(wfcv, nwordwfcv, iunhubv, ik)
      CALL calbec(npw_, wfcv, evc, proj)
      !
      ! calculate the occupation matrices n
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
         IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2 ) == 0 ) CYCLE
         !
         DO m1 = 1, ehub_lm(nt1)
         DO m2 = 1, ehub_lm(nt2)
            !
            DO is1 = 1, npol
            DO is2 = 1, npol
               !
               DO ibnd = 1, nbnd
                  n(1)%O(m1,m2,is1,is2,na1,na2) = &
                  n(1)%O(m1,m2,is1,is2,na1,na2) + &
                                         wg(ibnd,ik) * &
                     proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) * &
               CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))
               ENDDO ! ibnd
               !
            ENDDO
            ENDDO
            !
         ENDDO ! m2
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
         arg = tpi * (                                              &
          xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
          xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
          xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
         !
         kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
            DO m1 = 1, ehub_lm(nt1)
            DO m2 = 1, ehub_lm(nt2)
               !
               DO is1 = 1, npol
               DO is2 = 1, npol
                  !
                  DO ibnd = 1, nbnd
                     n(nnc)%O(m1,m2,is1,is2,na1,na2) = &
                     n(nnc)%O(m1,m2,is1,is2,na1,na2) + &
                           kphase*wg(ibnd,ik) * &
                           proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) * &
                     CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))
                  ENDDO ! ibnd
                  !
               ENDDO ! is2
               ENDDO ! is1
               !
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
      ! calculate the renormalized occupations N(ACBN0)
      !
      DO i = 1, nbnd
         ALLOCATE ( num(i)%O(4,4,ntyp,ntyp) ) ; num(i)%O = 0.0_DP
         ALLOCATE ( num(i)%T(4,4,ntyp,ntyp) ) ; num(i)%T = 0.0_DP
      ENDDO
      !
      DO na = 1, nat
         nt = ityp(na)
         !
         IF ( ehub_lm(nt) == 0 ) CYCLE
         !
         DO m = 1, ehub_lm(nt)
            !
            l = ehub_l(m,nt)
            !
            DO is1 = 1, npol
            !DO is2 = 1, npol
            !   !
            !   DO ibnd = 1, nbnd
            !      ! question: How to consider renormalization nubmer in
            !      ! nonlinear case?
            !      ! The definition of Nicolas Tancogne-Dejean & Angel Rubio
            !      num(ibnd)%O(l,l,nt,nt) = &
            !      num(ibnd)%O(l,l,nt,nt) + &
            !       DBLE(proj(ofstv(na)+m+counter_lm(m,nt)+(2*(ehub_l(m,nt)-1)+1)*(is2-1),ibnd) * &
            !      CONJG(proj(ofstv(na)+m+counter_lm(m,nt)+(2*(ehub_l(m,nt)-1)+1)*(is1-1),ibnd)))
            !   ENDDO ! ibnd
            !   !
            !ENDDO ! is2
            ! Consider total charge in each band.
               !
               DO ibnd = 1, nbnd
                  num(ibnd)%O(l,l,nt,nt) = &
                  num(ibnd)%O(l,l,nt,nt) + &
                   DBLE(proj(ofstv(na)+m+counter_lm(m,nt)+(2*(ehub_l(m,nt)-1)+1)*(is1-1),ibnd) * &
                  CONJG(proj(ofstv(na)+m+counter_lm(m,nt)+(2*(ehub_l(m,nt)-1)+1)*(is1-1),ibnd)))
               ENDDO
               !
            ENDDO ! is1
            !
         ENDDO ! m
         !
      ENDDO ! na
      !
      DO ibnd = 1, nbnd
         !
         DO nt1 = 1, ntyp
         DO nt2 = 1, ntyp
            !
            IF ( upf(nt1)%psd /= upf(nt2)%psd ) CYCLE
            !
            num(ibnd)%T(:,:,nt1,nt1) = num(ibnd)%T(:,:,nt1,nt1) + &
                                           num(ibnd)%O(:,:,nt2,nt2)
            !
         ENDDO ! nt2
         ENDDO ! nt1
         !
      ENDDO ! ibnd
      !
      DO ibnd = 1, nbnd
         !
         DO nt1 = 1, ntyp
         DO nt2 = 1, ntyp
            !
            DO l1 = 1, 4
            DO l2 = 1, 4
               !
               ! Our definition
               num(ibnd)%T(l1,l2,nt1,nt2) = num(ibnd)%T(l1,l1,nt1,nt1) + &
                                                num(ibnd)%T(l2,l2,nt2,nt2)
               !
               IF ( upf(nt1)%psd == upf(nt2)%psd )                       &
               num(ibnd)%T(l1,l2,nt1,nt2) = num(ibnd)%T(l1,l2,nt1,nt2) / &
                                                2.0_DP
               !
               ! The definition of Nicolas Tancogne-Dejean & Angel Rubio
 !             num(ibnd)%T(l1,l2,:,nt1,nt2) = SQRT(num(ibnd)%T(l1,l1,:,nt1,nt1) * &
 !                                                 num(ibnd)%T(l2,l2,:,nt2,nt2))
               !
            ENDDO ! l2
            ENDDO ! l1
            !
         ENDDO ! nt2
         ENDDO ! nt1
         !
      ENDDO ! ibnd
      !
      DO ibnd = 1, nbnd
         num(ibnd)%O = num(ibnd)%T
         DEALLOCATE ( num(ibnd)%T )
      ENDDO
      !
      ! calculate the renormalized occupation matrices P(ACBN0)
      !
      ! 1) P^II (dm0, dm1)
      !
      DO na = 1, nat
         nt = ityp(na)
         !
         IF ( ehub_lm(nt) == 0 ) CYCLE
         !
         DO m1 = 1, ehub_lm(nt)
            !
            l1 = ehub_l(m1,nt)
            !
         DO m2 = 1, ehub_lm(nt)
            !
            l2 = ehub_l(m2,nt)
            !
            DO is1 = 1, npol
            DO is2 = 1, npol
               !
               DO ibnd = 1, nbnd
                  dm0%O(m1,m2,is1,is2,na) = &
                  dm0%O(m1,m2,is1,is2,na) + &
                            num(ibnd)%O(l1,l2,nt,nt) * wg(ibnd,ik) * &
                  proj(ofstv(na)+m2+counter_lm(m2,nt)+(2*(ehub_l(m2,nt)-1)+1)*(is2-1),ibnd) * &
            CONJG(proj(ofstv(na)+m1+counter_lm(m1,nt)+(2*(ehub_l(m1,nt)-1)+1)*(is1-1),ibnd))
               ENDDO ! ibnd
               !
            ENDDO ! is1
            ENDDO ! is2
            !
         ENDDO ! m2
         ENDDO ! m1
         !
      ENDDO ! na
 
      DO nt = 1, ntyp
         !
         IF ( ehub_lm(nt) == 0 ) CYCLE
         !
         DO l = 1, 4
            !
            DO na1 = 1, nat
               nt1 = ityp(na1)
               !
               IF ( ehub_lm(nt1) == 0 ) CYCLE
               !
               DO m1 = 1, ehub_lm(nt1)
               DO m2 = 1, ehub_lm(nt1)
                  !
                  DO is1 = 1, npol
                  DO is2 = 1, npol
                     !
                     DO ibnd = 1, nbnd
                        dm1(nt,l)%O(m1,m2,is1,is2,na1) = &
                        dm1(nt,l)%O(m1,m2,is1,is2,na1) + &
                           num(ibnd)%O(ehub_l(m1,nt1),l,nt1,nt) * wg(ibnd,ik) * &
                           proj(ofstv(na1)+m2+counter_lm(m2,nt1)+(2*(ehub_l(m2,nt1)-1)+1)*(is2-1),ibnd) * &
                     CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))                           
                     ENDDO ! ibnd
                     !
                  ENDDO !is1
                  ENDDO !is2
                  !
               ENDDO ! m2
               ENDDO ! m1
               !
            ENDDO ! na1
            !
         ENDDO ! l
         !
      ENDDO ! nt
      !
      ! 2) P^IJ
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
         IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
         !
         DO m1 = 1, ehub_lm(nt1)
            !
            l1 = ehub_l(m1,nt1)
            !
         DO m2 = 1, ehub_lm(nt2)
            !
            l2 = ehub_l(m2,nt2)
            !
            DO is1 = 1, npol
            DO is2 = 1, npol
               !
               DO ibnd = 1, nbnd
                    p(1)%O(m1,m2,is1,is2,na1,na2) = &
                    p(1)%O(m1,m2,is1,is2,na1,na2) + &
                      num(ibnd)%O(l1,l2,nt1,nt2) * wg(ibnd,ik) * &
                      proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) * &
                CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd)) 
               ENDDO ! ibnd
               !
            ENDDO
            ENDDO
            !
         ENDDO ! m2
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
         arg = tpi * (                                              &
          xk(1,ik) * (n1 * at(1,1) + n2 * at(1,2) + n3 * at(1,3)) + &
          xk(2,ik) * (n1 * at(2,1) + n2 * at(2,2) + n3 * at(2,3)) + &
          xk(3,ik) * (n1 * at(3,1) + n2 * at(3,2) + n3 * at(3,3)))
         !
         kphase = CMPLX(COS(arg), SIN(arg), KIND=DP)
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
            DO m1 = 1, ehub_lm(nt1)
               !
               l1 = ehub_l(m1,nt1)
               !
            DO m2 = 1, ehub_lm(nt2)
               !
               l2 = ehub_l(m2,nt2)
               !
               DO is1 = 1, npol
               DO is2 = 1, npol
                  !
                  DO ibnd = 1, nbnd
                     p(nnc)%O(m1,m2,is1,is2,na1,na2) = &
                     p(nnc)%O(m1,m2,is1,is2,na1,na2) + &
                         kphase * num(ibnd)%O(l1,l2,nt1,nt2) * wg(ibnd,ik) * &
                         proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) * &
                   CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))                          
                  ENDDO ! ibnd
                  !
               ENDDO ! is1
               ENDDO ! is2
               !
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
      DO i = 1, nbnd
         DEALLOCATE ( num(i)%O )
      ENDDO
 
   ENDDO ! ik
 
   CALL mp_sum(dm0%O, inter_pool_comm)
   !
   DO nt = 1, ntyp
      DO l = 1, 4
         CALL mp_sum(dm1(nt,l)%O, inter_pool_comm)
      ENDDO
   ENDDO
   !
   DO i = 1, ehub_nnnc
      CALL mp_sum(n(i)%O, inter_pool_comm)
      CALL mp_sum(p(i)%O, inter_pool_comm)
   ENDDO
   !
   CALL symmetrize_n5_nc(n)
   !
   ns = 0.0_DP
   !
   DO i = 1, ehub_nnnc
      DO is1 = 1, npol
      DO is2 = 1, npol
         j = npol*(is1-1) + is2
         ns(:,:,j,:,:,i) = n(i)%O(:,:,is1,is2,:,:)
      ENDDO
      ENDDO 
      !
      DEALLOCATE ( n(i)%O )
      !
   ENDDO
   !
   DEALLOCATE(proj)
   !
   CALL symmetrize_n4_nc(dm0)
   !
   CALL symmetrize_n4m_nc(dm1)
   !
   CALL symmetrize_n5_nc(p)
   !
   ! on-site interaction: U_ij^II
   !
   E0(1)%V = 0.0_DP ; E0(1)%J = 0.0_DP
   !
   DO na = 1, nat
      nt = ityp(na)
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      l1 = 0
      DO m1 = 1, ehub_lm(nt)
         CALL orberi(ehub_l(m1,nt), l1)
         CALL pao_3g(upf(nt)%psd, ehub_l(m1,nt), c1)
         !
      l2 = 0
      DO m2 = 1, ehub_lm(nt)
         CALL orberi(ehub_l(m2,nt), l2)
         CALL pao_3g(upf(nt)%psd, ehub_l(m2,nt), c2)
         !
      l3 = 0
      DO m3 = 1, ehub_lm(nt)
         IF ( ehub_l(m3,nt) /= ehub_l(m1,nt) ) CYCLE
         CALL orberi(ehub_l(m3,nt), l3)
         CALL pao_3g(upf(nt)%psd, ehub_l(m3,nt), c3)
         !
      l4 = 0
      DO m4 = 1, ehub_lm(nt)
         IF ( ehub_l(m4,nt) /= ehub_l(m2,nt) ) CYCLE
         CALL orberi(ehub_l(m4,nt), l4)
         CALL pao_3g(upf(nt)%psd, ehub_l(m4,nt), c4)
         !
         IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
            E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
            E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
          DBLE(dm0%O(m3,m1,1,1,na)*dm0%O(m4,m2,2,2,na)  + &
               dm0%O(m3,m1,2,2,na)*dm0%O(m4,m2,1,1,na)) * &
                              compute_eri(l3,l1,l4,l2,    &
                        tau(:,na)*alat,tau(:,na)*alat,    &
                        tau(:,na)*alat,tau(:,na)*alat,    &
                                          c3,c1,c4,c2)
         ELSE
            E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
            E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
          DBLE(dm0%O(m3,m1,1,1,na)*dm0%O(m4,m2,1,1,na)  + &
               dm0%O(m3,m1,2,2,na)*dm0%O(m4,m2,2,2,na)  + &
               dm0%O(m3,m1,1,1,na)*dm0%O(m4,m2,2,2,na)  + &
               dm0%O(m3,m1,2,2,na)*dm0%O(m4,m2,1,1,na)) * &
                              compute_eri(l3,l1,l4,l2,    &
                        tau(:,na)*alat,tau(:,na)*alat,    &
                        tau(:,na)*alat,tau(:,na)*alat,    &
                                          c3,c1,c4,c2)
 
            E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  = &
            E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  + &
          DBLE(dm0%O(m3,m1,1,1,na)*dm0%O(m4,m2,1,1,na)  + &
               dm0%O(m3,m1,2,2,na)*dm0%O(m4,m2,2,2,na)  + &
               dm0%O(m3,m1,2,1,na)*dm0%O(m4,m2,1,2,na)  + &
               dm0%O(m3,m1,1,2,na)*dm0%O(m4,m2,2,1,na)) * &
                              compute_eri(l3,l2,l4,l1,    &
                        tau(:,na)*alat,tau(:,na)*alat,    &
                        tau(:,na)*alat,tau(:,na)*alat,    &
                                          c3,c2,c4,c1)
         ENDIF
         !
      ENDDO ! m4
      ENDDO ! m3
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na
   !
   HF(1)%V = 0.0_DP ; HF(1)%J = 0.0_DP
   !
   DO na = 1, nat
      nt = ityp(na)
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt)
      DO m2 = 1, ehub_lm(nt)
         !
         IF ( m1 == m2 ) THEN
             HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
             HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
         DBLE(ns(m1,m1,1,na,na,1) * ns(m2,m2,4,na,na,1) + &
              ns(m1,m1,4,na,na,1) * ns(m2,m2,1,na,na,1) - &
              ns(m1,m1,2,na,na,1) * ns(m2,m2,3,na,na,1) - &
              ns(m1,m1,3,na,na,1) * ns(m2,m2,2,na,na,1))
         ELSE
             HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
             HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
         DBLE(ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
              ns(m1,m1,4,na,na,1) * ns(m2,m2,4,na,na,1) + &
              ns(m1,m1,1,na,na,1) * ns(m2,m2,4,na,na,1) + &
              ns(m1,m1,4,na,na,1) * ns(m2,m2,1,na,na,1))
 
             HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na) = &
             HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na) + &
         DBLE(ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
              ns(m1,m1,4,na,na,1) * ns(m2,m2,4,na,na,1) + &
              ns(m1,m1,2,na,na,1) * ns(m2,m2,3,na,na,1) + &
              ns(m1,m1,3,na,na,1) * ns(m2,m2,2,na,na,1)) 
         ENDIF
         !
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na
   !
   ALLOCATE ( ehub_v_(nat,nat,ehub_lm_max,ehub_lm_max) ) ; ehub_v_ = 0.0_DP
   IF ( .NOT.ALLOCATED(ehub_v_raw) ) THEN   ! [rawprint]
      ALLOCATE ( ehub_v_raw(SIZE(ehub_v,1),SIZE(ehub_v,2),SIZE(ehub_v,3),SIZE(ehub_v,4),SIZE(ehub_v,5)) )
      ehub_v_raw = 0.0_DP
   ENDIF
   !
   DO na = 1, nat
      nt = ityp(na)
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt)
      DO m2 = 1, ehub_lm(nt)
         !
         IF ( ehub_l(m1,nt) == 1 .AND. ehub_l(m2,nt) == 1 ) THEN
            ehub_v_(na,na,m1,m2) =                           &
                E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) / &
                HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na) * &
                           remove_ehub_u(nt,ehub_l(m1,nt)) * &
                           remove_ehub_u(nt,ehub_l(m2,nt))
         ELSE
            ehub_v_(na,na,m1,m2) =                            &
               (E0(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  / &
                HF(1)%V(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  - &
                E0(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)  / &
                HF(1)%J(ehub_l(m1,nt),ehub_l(m2,nt),na,na)) * &
                           remove_ehub_u(nt,ehub_l(m1,nt))  * &
                           remove_ehub_u(nt,ehub_l(m2,nt))
         ENDIF
         !
         IF ( stop_acbn0 ) ehub_v_diff(na,na,m1,m2,1) = ABS( ehub_v (na,na,m1,m2,1) &
                                                      -      ehub_v_(na,na,m1,m2  ) )
         !
         ehub_v_raw(na,na,m1,m2,1) = ehub_v_(na,na,m1,m2)   ! [rawprint]
         !
         ehub_v(na,na,m1,m2,1)                                 &
             = (1.0_DP - ehub_mixing)*ehub_v (na,na,m1,m2,1) + &
               (         ehub_mixing)*ehub_v_(na,na,m1,m2)
         !
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na
 
   !
   ! inter-site interaction: U_ij^IJ (I /= J)
   !
   ! ... in the unit cell
   !
   E0(1)%V = 0.0_DP
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
   DO na2 = na1, nat
      nt2 = ityp(na2)
      !
      IF ( ehub_lm(nt2) == 0 ) CYCLE
      !
      IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
      !
      l1 = 0
      DO m1 = 1, ehub_lm(nt1)
         CALL orberi(ehub_l(m1,nt1), l1)
         CALL pao_3g(upf(nt1)%psd, ehub_l(m1,nt1), c1)
         !
      l2 = 0
      DO m2 = 1, ehub_lm(nt2)
         CALL orberi(ehub_l(m2,nt2), l2)
         CALL pao_3g(upf(nt2)%psd, ehub_l(m2,nt2), c2)
         !
      l3 = 0
      DO m3 = 1, ehub_lm(nt1)
         IF ( ehub_l(m3,nt1) /= ehub_l(m1,nt1) ) CYCLE
         CALL orberi(ehub_l(m3,nt1), l3)
         CALL pao_3g(upf(nt1)%psd, ehub_l(m3,nt1), c3)
         !
      l4 = 0
      DO m4 = 1, ehub_lm(nt2)
         IF ( ehub_l(m4,nt2) /= ehub_l(m2,nt2) ) CYCLE
         CALL orberi(ehub_l(m4,nt2), l4)
         CALL pao_3g(upf(nt2)%psd, ehub_l(m4,nt2), c4)
         !
         E0(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  = &
         E0(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  + &
          DBLE(dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,1,na1)  * &
               dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,1,na2)  + &
               dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,2,na1)  * &
               dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,2,na2)  + &
               dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,1,na1)  * &
               dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,2,na2)  + &
               dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,2,na1)  * &
               dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,1,na2)  - &
                              p(1)%O(m3,m2,1,1,na1,na2)  * &
                        CONJG(p(1)%O(m1,m4,2,2,na1,na2))  - &
                              p(1)%O(m3,m2,2,2,na1,na2)  * &
                        CONJG(p(1)%O(m1,m4,1,1,na1,na2))  - &
                              p(1)%O(m3,m2,1,2,na1,na2)  * &
                        CONJG(p(1)%O(m1,m4,1,2,na1,na2))  - &
                              p(1)%O(m3,m2,2,1,na1,na2)  * &
                        CONJG(p(1)%O(m1,m4,2,1,na1,na2)))* &
                               compute_eri(l3,l1,l4,l2,    &
                       tau(:,na1)*alat,tau(:,na1)*alat,    &
                       tau(:,na2)*alat,tau(:,na2)*alat,    &
                                           c3,c1,c4,c2)
         !
      ENDDO ! m4
      ENDDO ! m3
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
   DO na2 = na1, nat
      nt2 = ityp(na2)
      !
      IF ( ehub_lm(nt2) == 0 ) CYCLE
      !
      IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt2)
         !
         E0(1)%V(ehub_l(m2,nt2),ehub_l(m1,nt1),na2,na1) = E0(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)
         !
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
   !
   HF(1)%V = 0.0_DP
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
   DO na2 = na1, nat
      nt2 = ityp(na2)
      !
      IF ( ehub_lm(nt2) == 0 ) CYCLE
      !
      IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt2)
         !
         HF(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2) = &
         HF(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2) + &
     DBLE(ns(m1,m1,1,na1,na1,1) * ns(m2,m2,1,na2,na2,1) + &
          ns(m1,m1,4,na1,na1,1) * ns(m2,m2,4,na2,na2,1) + &
          ns(m1,m1,1,na1,na1,1) * ns(m2,m2,4,na2,na2,1) + &
          ns(m1,m1,4,na1,na1,1) * ns(m2,m2,1,na2,na2,1) - &
          ns(m1,m2,1,na1,na2,1) *CONJG(ns(m1,m2,4,na1,na2,1)) - &
          ns(m1,m2,4,na1,na2,1) *CONJG(ns(m1,m2,1,na1,na2,1)) - &
          ns(m1,m2,2,na1,na2,1) *CONJG(ns(m1,m2,2,na1,na2,1)) - &
          ns(m1,m2,3,na1,na2,1) *CONJG(ns(m1,m2,3,na1,na2,1)))
         !
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
   DO na2 = na1, nat
      nt2 = ityp(na2)
      !
      IF ( ehub_lm(nt2) == 0 ) CYCLE
      !
      IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt2)
         HF(1)%V(ehub_l(m2,nt2),ehub_l(m1,nt1),na2,na1) = HF(1)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
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
      IF ( ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
         !
         l1 = ehub_l(m1,nt1)
         !
      DO m2 = 1, ehub_lm(nt2)
         !
         l2 = ehub_l(m2,nt2)
         !
         ehub_v_(na1,na2,m1,m2) = &
                (E0(1)%V(l1,l2,na1,na2) / HF(1)%V(l1,l2,na1,na2)) / 2.0_DP
         !
         IF ( stop_acbn0 ) ehub_v_diff(na1,na2,m1,m2,1) = ABS( ehub_v (na1,na2,m1,m2,1) &
                                                        -      ehub_v_(na1,na2,m1,m2  ) )
         !
         ehub_v_raw(na1,na2,m1,m2,1) = ehub_v_(na1,na2,m1,m2)   ! [rawprint]
         !
         ehub_v(na1,na2,m1,m2,1)                                 &
             = (1.0_DP - ehub_mixing)*ehub_v (na1,na2,m1,m2,1) + &
               (         ehub_mixing)*ehub_v_(na1,na2,m1,m2)
         !
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
   !
   ! ... out of the unit cell
   !
   DO i = 2, ehub_nnnc
      E0(i)%V = 0.0_DP
   ENDDO
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
         DO ix = 1, 3
            ntau(ix,na2) = n1 * at(ix,1) + n2 * at(ix,2) + n3 * at(ix,3) + tau(ix,na2)
         ENDDO
         !
         l1 = 0
         DO m1 = 1, ehub_lm(nt1)
            CALL orberi(ehub_l(m1,nt1), l1)
            CALL pao_3g(upf(nt1)%psd, ehub_l(m1,nt1), c1)
            !
         l2 = 0
         DO m2 = 1, ehub_lm(nt2)
            CALL orberi(ehub_l(m2,nt2), l2)
            CALL pao_3g(upf(nt2)%psd, ehub_l(m2,nt2), c2)
            !
         l3 = 0
         DO m3 = 1, ehub_lm(nt1)
            IF ( ehub_l(m3,nt1) /= ehub_l(m1,nt1) ) CYCLE
            CALL orberi(ehub_l(m3,nt1), l3)
            CALL pao_3g(upf(nt1)%psd, ehub_l(m3,nt1), c3)
            !
         l4 = 0
         DO m4 = 1, ehub_lm(nt2)
            IF ( ehub_l(m4,nt2) /= ehub_l(m2,nt2) ) CYCLE
            CALL orberi(ehub_l(m4,nt2), l4)
            CALL pao_3g(upf(nt2)%psd, ehub_l(m4,nt2), c4)
            !
            E0(nnc)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  = &
            E0(nnc)%V(ehub_l(m1,nt1),ehub_l(m2,nt2),na1,na2)  + &
                 DBLE(dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,1,na1)  * &
                      dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,1,na2)  + &
                      dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,2,na1)  * &
                      dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,2,na2)  + &
                      dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,1,1,na1)  * &
                      dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,2,2,na2)  + &
                      dm1(nt2,ehub_l(m2,nt2))%O(m3,m1,2,2,na1)  * &
                      dm1(nt1,ehub_l(m1,nt1))%O(m4,m2,1,1,na2)  - &
                                   p(nnc)%O(m3,m2,1,1,na1,na2)  * &
                             CONJG(p(nnc)%O(m1,m4,2,2,na1,na2)) - &
                                   p(nnc)%O(m3,m2,2,2,na1,na2)  * &
                             CONJG(p(nnc)%O(m1,m4,1,1,na1,na2)) - &
                                   p(nnc)%O(m3,m2,1,2,na1,na2)  * &
                             CONJG(p(nnc)%O(m1,m4,1,2,na1,na2)) - &
                                   p(nnc)%O(m3,m2,2,1,na1,na2)  * &
                             CONJG(p(nnc)%O(m1,m4,2,1,na1,na2)))* &
                                    compute_eri(l3,l1,l4,l2,    &
                           tau(:,na1)*alat, tau(:,na1)*alat,    &
                          ntau(:,na2)*alat,ntau(:,na2)*alat,    &
                                                c3,c1,c4,c2)
            !
         ENDDO ! m4
         ENDDO ! m3
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
   DEALLOCATE ( dm0%O )
   !
   DO nt = 1, ntyp
      DO l = 1, 4
         DEALLOCATE ( dm1(nt,l)%O )
      ENDDO
   ENDDO
   !
   DO i = 2, ehub_nnnc
      HF(i)%V = 0.0_DP
   ENDDO
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
         DO m1 = 1, ehub_lm(nt1)
            !
            l1 = ehub_l(m1,nt1)
            !
         DO m2 = 1, ehub_lm(nt2)
            !
            l2 = ehub_l(m2,nt2)
            !
            HF(nnc)%V(l1,l2,na1,na2) = HF(nnc)%V(l1,l2,na1,na2) + &
           DBLE(ns(m1,m1,1,na1,na1,  1)*ns(m2,m2,1,na2,na2,  1) + &
                ns(m1,m1,4,na1,na1,  1)*ns(m2,m2,4,na2,na2,  1) + &
                ns(m1,m1,1,na1,na1,  1)*ns(m2,m2,4,na2,na2,  1) + &
                ns(m1,m1,4,na1,na1,  1)*ns(m2,m2,1,na2,na2,  1) - &          
                ns(m1,m2,1,na1,na2,nnc)*CONJG(ns(m1,m2,4,na1,na2,nnc)) - &
                ns(m1,m2,4,na1,na2,nnc)*CONJG(ns(m1,m2,1,na1,na2,nnc)) - &
                ns(m1,m2,2,na1,na2,nnc)*CONJG(ns(m1,m2,2,na1,na2,nnc)) - &
                ns(m1,m2,3,na1,na2,nnc)*CONJG(ns(m1,m2,3,na1,na2,nnc)))
            !
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
         DO m1 = 1, ehub_lm(nt1)
            !
            l1 = ehub_l(m1,nt1)
            !
         DO m2 = 1, ehub_lm(nt2)
            !
            l2 = ehub_l(m2,nt2)
            !
            ehub_v_(na1,na2,m1,m2) = &
                   (E0(nnc)%V(l1,l2,na1,na2) / HF(nnc)%V(l1,l2,na1,na2)) / 2.0_DP
            !
            IF ( stop_acbn0 ) ehub_v_diff(na1,na2,m1,m2,nnc) = ABS( ehub_v (na1,na2,m1,m2,nnc) &
                                                             -      ehub_v_(na1,na2,m1,m2    ) )
            !
            ehub_v_raw(na1,na2,m1,m2,nnc) = ehub_v_(na1,na2,m1,m2)   ! [rawprint]
            !
            ehub_v(na1,na2,m1,m2,nnc) =                               &
                  (1.0_DP - ehub_mixing)*ehub_v (na1,na2,m1,m2,nnc) + &
                  (         ehub_mixing)*ehub_v_(na1,na2,m1,m2    )
            !
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
   DEALLOCATE ( ehub_v_ )
   !
   DO i = 1, ehub_nnnc
      DEALLOCATE ( E0(i)%V, E0(i)%J )
      DEALLOCATE ( HF(i)%V, HF(i)%J )
   ENDDO
   !
   ! Hermiticity of Hubbard parameters U & V
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
   DO na2 = na1, nat
      nt2 = ityp(na2)
      !
      IF ( ehub_lm(nt2) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt2)
         ehub_v(na1,na2,m1,m2,1) = 0.5D0 * (ehub_v(na1,na2,m1,m2,1) + &
                                            ehub_v(na2,na1,m2,m1,1))
         !
         ehub_v(na2,na1,m2,m1,1) = ehub_v(na1,na2,m1,m2,1)
         ! [rawprint]
         ehub_v_raw(na1,na2,m1,m2,1) = 0.5D0 * (ehub_v_raw(na1,na2,m1,m2,1) + ehub_v_raw(na2,na1,m2,m1,1))
         ehub_v_raw(na2,na1,m2,m1,1) = ehub_v_raw(na1,na2,m1,m2,1)
         !
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
   !
   DO na1 = 1, nat
   DO na2 = 1, nat
      IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) THEN
         ehub_v(na1,na2,:,:,1) = 0.0_DP
         ehub_v_raw(na1,na2,:,:,1) = 0.0_DP   ! [rawprint]
      ENDIF
   ENDDO ! na2
   ENDDO ! na1
 
 1 CALL write_ehub_uv
 
   CALL stop_clock('new_ns_e0_nc')
 
   RETURN
 
END SUBROUTINE new_ns_e0_nc_1
!
SUBROUTINE symmetrize_n4_nc( n )

   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp
   USE ldaV,                 ONLY : ehub_lm, ehub_l, ehub_lm_max, &
                                    ehub_nn_vect, ehub_nn_dim,    &
                                    ehub_csize, ehub_nnnc, d_spin_ldav
   USE lsda_mod,             ONLY : nspin
   USE noncollin_module,     ONLY : npol
   USE symm_base,            ONLY : d1, d2, d3, nsym, irt, time_reversal, t_rev, sr
 
   IMPLICIT NONE
 
   TYPE occmat_4_nc
      COMPLEX(DP), ALLOCATABLE :: O(:,:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: T(:,:,:,:,:)
   END TYPE occmat_4_nc
 
   TYPE(occmat_4_nc), INTENT(INOUT) :: n
 
   INTEGER                       :: isym, l, m, nt
   INTEGER                       :: na1, na2
   INTEGER                       :: nt1, nt2
   INTEGER                       :: m1, m2, m3, m4
   INTEGER                       :: m1_, m2_, m3_, m4_
   INTEGER                       :: temp, l_separation(5,ntyp)
   integer                       :: is1, is2, is3, is4
   COMPLEX(DP)                   :: s_spin(2,2)
   COMPLEX(DP)                   :: n11, n12, n21, n22
 
   ALLOCATE (n%T(ehub_lm_max,ehub_lm_max,npol,npol,nat)) ; n%T = 0.0_DP
 
   DO nt = 1, ntyp
      !
      temp = ehub_l(1,nt)
      l_separation(ehub_l(1,nt),nt) = 1
      !
      DO m = 2, ehub_lm(nt)
         !
         IF ( ehub_l(m,nt) /= temp ) THEN
            temp = ehub_l(m,nt)
            l_separation(ehub_l(m,nt),nt) = m
         ENDIF
         !
      ENDDO
      !
   ENDDO
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
         m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
         !
      DO m2 = 1, ehub_lm(nt1)
         m2_ = m2 - l_separation(ehub_l(m2,nt1),nt1) + 1
         !
         DO is1 = 1, npol
         DO is2 = 1, npol
            !
            DO isym = 1, nsym
               !
               CALL find_u(sr(:,:,isym), s_spin(:,:))
               !
               na2 = irt(isym, na1)
               !
               DO m3 = 1, ehub_lm(nt1)
                  IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                  m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                  !
               DO m4 = 1, ehub_lm(nt1)
                  IF ( ehub_l(m2,nt1) /= ehub_l(m4,nt1) ) CYCLE
                  m4_ = m4 - l_separation(ehub_l(m4,nt1),nt1) + 1
                  !
                  DO is3 = 1, npol
                  DO is4 = 1, npol
                  !
                  IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))           
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))           
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))       
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))       
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF 
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                      CONJG(n%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n%T(m1,m2,is1,is2,na1) = n%T(m1,m2,is1,is2,na1) + &
                                            d_spin_ldav(is3,is1,isym) * &
                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                            n%O(m3,m4,is3,is4,na2) / (nsym) * &
                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                      CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF 
                  ELSE
                     CALL errore('new_ns_v','angular momentum not implemented',1)
                  ENDIF
                  !
                  ENDDO ! is4
                  ENDDO ! is3
                  !
               ENDDO ! m4
               ENDDO ! m3
               !
            ENDDO ! isym
            !
         ENDDO ! is2
         ENDDO ! is1
         ! 
      ENDDO ! m2
      ENDDO ! m1
      ! 
   ENDDO ! na1
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO is1 = 1, npol
      DO is2 = 1, npol
         !
         DO m1 = 1, ehub_lm(nt1)
         DO m2 = m1, ehub_lm(nt1)
            !
            n%T(m1,m2,is1,is2,na1) = 0.5D0 * (n%T(m1,m2,is1,is2,na1) + &
                                        CONJG(n%T(m2,m1,is2,is1,na1)))
            !
            n%T(m2,m1,is2,is1,na1) = CONJG(n%T(m1,m2,is1,is2,na1))
            !
         ENDDO
         ENDDO
         !
      ENDDO ! is2
      ENDDO ! is1
      !
   ENDDO
   !
   !IF ( time_reversal) THEN
   !!
   !n11 = 0.0_DP
   !n12 = 0.0_DP
   !n21 = 0.0_DP
   !n22 = 0.0_DP
   !!
   !DO na1 = 1, nat
   !   nt1 = ityp(na1)
   !   !
   !   IF (ehub_lm(nt1) == 0 ) CYCLE
   !   !
   !   DO m1 = 1, ehub_lm(nt1)
   !   DO m2 = 1, ehub_lm(nt1)
   !   !
   !   n11 = n%T(m1,m2,1,1,na1)
   !   n11 = n%T(m1,m2,1,2,na1)
   !   n11 = n%T(m1,m2,2,1,na1)
   !   n11 = n%T(m1,m2,2,2,na1)
   !   !
   !   n%T(m1,m2,1,1,na1) = 0.5_DP*(n11+CONJG(n22))
   !   n%T(m1,m2,1,2,na1) = 0.5_DP*(n12-CONJG(n21))
   !   n%T(m1,m2,2,1,na1) = -CONJG(n%T(m1,m2,1,2,na1))
   !   n%T(m1,m2,2,2,na1) = CONJG(n%T(m1,m2,1,1,na1))
   !   !
   !   ENDDO ! m2
   !   ENDDO ! m1
   !   !
   !ENDDO ! na1
   !!
   !ENDIF
   !
   n%O = n%T
   DEALLOCATE (n%T)
 
END SUBROUTINE symmetrize_n4_nc
!
SUBROUTINE symmetrize_n4m_nc( n )

   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp
   USE ldaV,                 ONLY : ehub_lm, ehub_l, ehub_lm_max, &
                                    ehub_nn_vect, ehub_nn_dim,    &
                                    ehub_csize, ehub_nnnc, d_spin_ldav
   USE lsda_mod,             ONLY : nspin
   USE noncollin_module,     ONLY : npol
   USE symm_base,            ONLY : d1, d2, d3, nsym, irt, time_reversal, t_rev, sr
 
   IMPLICIT NONE
 
   TYPE occmat_4_nc
      COMPLEX(DP), ALLOCATABLE :: O(:,:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: T(:,:,:,:,:)
   END TYPE occmat_4_nc
 
   TYPE(occmat_4_nc), INTENT(INOUT) :: n(ntyp,4)
 
   INTEGER                       :: isym, i, l, m, nt
   INTEGER                       :: na1, na2
   INTEGER                       :: nt1, nt2
   INTEGER                       :: m1, m2, m3, m4
   INTEGER                       :: m1_, m2_, m3_, m4_
   INTEGER                       :: temp, l_separation(5,ntyp)
   INTEGER                       :: is1, is2, is3, is4
   COMPLEX(DP)                   :: s_spin(2,2)
   COMPLEX(DP)                   :: n11, n12, n21, n22
 
   DO i = 1, ntyp
      DO l = 1, 4
         ALLOCATE (n(i,l)%T(ehub_lm_max,ehub_lm_max,npol,npol,nat)) ; n(i,l)%T = 0.0_DP
      ENDDO
   ENDDO
   !
   DO nt = 1, ntyp
      !
      temp = ehub_l(1,nt)
      l_separation(ehub_l(1,nt),nt) = 1
      !
      DO m = 2, ehub_lm(nt)
         !
         IF ( ehub_l(m,nt) /= temp ) THEN
            temp = ehub_l(m,nt)
            l_separation(ehub_l(m,nt),nt) = m
         ENDIF
         !
      ENDDO
      !
   ENDDO
   !
   DO nt = 1, ntyp
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      DO l = 1, 4
         !
         DO na1 = 1, nat
            nt1 = ityp(na1)
            !
            IF ( ehub_lm(nt1) == 0 ) CYCLE
            !
            DO m1 = 1, ehub_lm(nt1)
               m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
               !
            DO m2 = 1, ehub_lm(nt1)
               m2_ = m2 - l_separation(ehub_l(m2,nt1),nt1) + 1
               !
               DO is1 = 1, npol
               DO is2 = 1, npol
                  !
                  DO isym = 1, nsym
                     !
                     CALL find_u(sr(:,:,isym), s_spin(:,:))
                     !
                     na2 = irt(isym, na1)
                     !
                     DO m3 = 1, ehub_lm(nt1)
                        IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                        m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                        !
                     DO m4 = 1, ehub_lm(nt1)
                        IF ( ehub_l(m2,nt1) /= ehub_l(m4,nt1) ) CYCLE
                        m4_ = m4 - l_separation(ehub_l(m4,nt1),nt1) + 1
                        !
                        DO is3 = 1, npol
                        DO is4 = 1, npol
                        !
                        IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))      
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))    
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF 
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF 
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF 
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF 
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF 
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF 
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                    CONJG(n(nt,l)%O(m3,m4,is3,is4,na2)) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nt,l)%T(m1,m2,is1,is2,na1) = n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                          d_spin_ldav(is3,is1,isym) * &
                                                    cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &
                                                          n(nt,l)%O(m3,m4,is3,is4,na2) / (nsym) * &
                                                    cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) * &
                                                    CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE
                           CALL errore('new_ns_v','angular momentum not implemented',1)
                        ENDIF
                        !
                        ENDDO ! is4
                        ENDDO ! is3
                        !
                     ENDDO ! m4
                     ENDDO ! m3
                     !
                  ENDDO ! isym
                  !
               ENDDO ! is2
               ENDDO ! is1
               !
            ENDDO ! m2
            ENDDO ! m1
            ! 
         ENDDO ! na1
         !
      ENDDO ! l
      !
   ENDDO ! nt
   !
   DO nt = 1, ntyp
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      DO l = 1, 4
         !
         DO na1 = 1, nat
            nt1 = ityp(na1)
            !
            IF ( ehub_lm(nt1) == 0 ) CYCLE
            !
            DO is1 = 1, npol
            DO is2 = 1, npol
               !
               DO m1 = 1, ehub_lm(nt1)
               DO m2 = m1, ehub_lm(nt1)
                  !
                  n(nt,l)%T(m1,m2,is1,is2,na1) = 0.5D0 * (n(nt,l)%T(m1,m2,is1,is2,na1) + &
                                                    CONJG(n(nt,l)%T(m2,m1,is2,is1,na1)))
                  !
                  n(nt,l)%T(m2,m1,is2,is1,na1) = CONJG(n(nt,l)%T(m1,m2,is1,is2,na1))
                  !
               ENDDO
               ENDDO
               !
            ENDDO
            ENDDO 
            !
         ENDDO
         !
      ENDDO
      !
   ENDDO
   !
   DO nt = 1, ntyp
      DO l = 1, 4
         n(nt,l)%O = n(nt,l)%T
         DEALLOCATE (n(nt,l)%T)
      ENDDO
   ENDDO
 
END SUBROUTINE symmetrize_n4m_nc
!
SUBROUTINE symmetrize_n5_nc( n )

   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp
   USE ldaV,                 ONLY : ehub_lm, ehub_l, ehub_lm_max, &
                                    ehub_nn_vect, ehub_nn_dim,    &
                                    ehub_csize, ehub_nnnc, d_spin_ldav
   USE lsda_mod,             ONLY : nspin
   USE noncollin_module,     ONLY : npol
   USE symm_base,            ONLY : d1, d2, d3, nsym, irt, time_reversal, t_rev, sr
 
   IMPLICIT NONE
 
   TYPE occmat_5_nc
      COMPLEX(DP), ALLOCATABLE :: O(:,:,:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: T(:,:,:,:,:,:)
   END TYPE occmat_5_nc
 
   TYPE(occmat_5_nc), INTENT(INOUT) :: n(ehub_nnnc)
 
   INTEGER                       :: isym, i, m, nt
   INTEGER                       :: na1, na2, na3, na4
   INTEGER                       :: nt1, nt2, nt3, nt4
   INTEGER                       :: m1, m2, m3, m4
   INTEGER                       :: m1_, m2_, m3_, m4_
   INTEGER                       :: n0(3), n1, n2, n3
   INTEGER                       :: nnc0, nnc, ics
   INTEGER                       :: temp, l_separation(5,ntyp)
   INTEGER                       :: is1, is2, is3, is4
   COMPLEX(DP)                   :: s_spin(2,2)
   COMPLEX(DP)                   :: n11, n12, n21, n22

   DO i = 1, ehub_nnnc
      ALLOCATE (n(i)%T(ehub_lm_max,ehub_lm_max,npol,npol,nat,nat)) ; n(i)%T = 0.0_DP
   ENDDO
   !
   DO nt = 1, ntyp
      !
      temp = ehub_l(1,nt)
      l_separation(ehub_l(1,nt),nt) = 1
      !
      DO m = 2, ehub_lm(nt)
         !
         IF ( ehub_l(m,nt) /= temp ) THEN
            temp = ehub_l(m,nt)
            l_separation(ehub_l(m,nt),nt) = m
         ENDIF
         !
      ENDDO
      !
   ENDDO
   !
   ics = ehub_csize + 1
   !
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
         m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
         !
      DO m2 = 1, ehub_lm(nt1)
         m2_ = m2 - l_separation(ehub_l(m2,nt1),nt1) + 1
         !
         DO is1 = 1, npol
         DO is2 = 1, npol
            DO isym = 1, nsym
               !
               CALL find_u(sr(:,:,isym), s_spin(:,:))
               !
               na2 = irt(isym, na1)
               !
               DO m3 = 1, ehub_lm(nt1)
                  IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                  m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                  !
               DO m4 = 1, ehub_lm(nt1)
                  IF ( ehub_l(m2,nt1) /= ehub_l(m4,nt1) ) CYCLE
                  m4_ = m4 - l_separation(ehub_l(m4,nt1),nt1) + 1
                  !
                  DO is3 = 1, npol
                  DO is4 = 1, npol
                  !
                  ! Below routines are from new_ns_nc in new_ns.f90
                  ! Or d_matrix_nc in d_matrix_nc.f90 in PP
                  ! With time reversal symmetry: spin 1 -> 2*, 2 -> -1* 
                  !
                  IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &    
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &    
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &                  
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) * &                  
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG((n(1)%O(m3,m4,is3,is4,na2,na2))) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym)) 
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF 
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1 )THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF 
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 1 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 2 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 3 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt1) == 4 ) THEN
                     IF (t_rev(isym) == 1) THEN
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                             CONJG(n(1)%O(m3,m4,is3,is4,na2,na2)) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ELSE
                     n(1)%T(m1,m2,is1,is2,na1,na1) = n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                   d_spin_ldav(is3,is1,isym) *  &
                                             cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                   n(1)%O(m3,m4,is3,is4,na2,na2) * &
                                             cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                             CONJG(d_spin_ldav(is4,is2,isym))
                     ENDIF
                  ELSE
                     CALL errore('new_ns_v','angular momentum not implemented',1)
                  ENDIF
                  !
                  ENDDO ! is4
                  ENDDO ! is3
                  !
               ENDDO ! m4
               ENDDO ! m3
               !
            ENDDO ! isym
         ENDDO ! is1
         ENDDO ! is2
         !
      ENDDO ! m2
      ENDDO ! m1
      ! 
   ENDDO ! na1
      !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
      DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
      DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
         !
         n0(1) = n1 ; n0(2) = n2 ; n0(3) = n3
         !
         CALL hubbard_get_nnc(n0,nnc0)
         !
         DO na2 = 1, nat
            nt2 = ityp(na2)
            !
            IF ( ehub_lm(nt2) == 0 ) CYCLE
            !
            IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
            !
            DO m1 = 1, ehub_lm(nt1)
               m1_ = m1 - l_separation(ehub_l(m1,nt1),nt1) + 1
               !
            DO m2 = 1, ehub_lm(nt2)
               m2_ = m2 - l_separation(ehub_l(m2,nt2),nt2) + 1
               !
               DO is1 = 1, npol
               DO is2 = 1, npol
                  !
                  DO isym = 1, nsym
                     !
                     CALL find_u(sr(1,1,isym), s_spin)
                     !
                     CALL hubbard_symmetry(na1,na2,na3,na4,n0,nnc,isym)
                     !
                     DO m3 = 1, ehub_lm(nt1)
                        IF ( ehub_l(m1,nt1) /= ehub_l(m3,nt1) ) CYCLE
                        m3_ = m3 - l_separation(ehub_l(m3,nt1),nt1) + 1
                        !
                     DO m4 = 1, ehub_lm(nt2)
                        IF ( ehub_l(m2,nt2) /= ehub_l(m4,nt2) ) CYCLE
                        m4_ = m4 - l_separation(ehub_l(m4,nt2),nt2) + 1
                        !
                        DO is3 = 1, npol
                        DO is4 = 1, npol
                        !
                        IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))

                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &                   
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &                   
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &                   
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &                   
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 1 .AND. ehub_l(m2,nt2) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &                   
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) * &                   
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))  
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))  
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 2 .AND. ehub_l(m2,nt2) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d1(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 3 .AND. ehub_l(m2,nt2) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d2(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 1 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym))
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 2 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d1(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 3 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d2(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE IF ( ehub_l(m1,nt1) == 4 .AND. ehub_l(m2,nt2) == 4 ) THEN
                           IF (t_rev(isym) == 1) THEN
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                      CONJG(n(nnc)%O(m3,m4,is3,is4,na3,na4)) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ELSE
                           n(nnc0)%T(m1,m2,is1,is2,na1,na2) = n(nnc0)%T(m1,m2,is1,is2,na1,na2) + &
                                                            d_spin_ldav(is3,is1,isym) *  &
                                                      cmplx(d3(m3_,m1_,isym), 0.d0, kind=DP) * &                                                                      
                                                            n(nnc)%O(m3,m4,is3,is4,na3,na4) * &
                                                      cmplx(d3(m4_,m2_,isym), 0.d0, kind=DP) / (nsym) * &
                                                      CONJG(d_spin_ldav(is4,is2,isym)) 
                           ENDIF
                        ELSE
                           CALL errore('new_ns_v','angular momentum not implemented',1)
                        ENDIF
                        !
                        ENDDO ! is4
                        ENDDO ! is3
                     !
                     ENDDO ! m4
                     ENDDO ! m3
                  !
                  ENDDO ! isym
                  !
               ENDDO ! is2
               ENDDO ! is1
               !
            ENDDO ! m2
            ENDDO ! m1
            !
         ENDDO ! na2
         ! 
      ENDDO ! n3
      ENDDO ! n2
      ENDDO ! n1
      !
   ENDDO ! na1
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
   !DO na2 = na1, nat
   !   nt2 = ityp(na2)
   !   !
   !   IF ( ehub_lm(nt2) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt1)
         !
         DO is1 = 1, npol
         DO is2 = 1, npol
            n(1)%T(m1,m2,is1,is2,na1,na1) = 0.5_DP * (n(1)%T(m1,m2,is1,is2,na1,na1) + &
                                                CONJG(n(1)%T(m2,m1,is2,is1,na1,na1)))
            !
            n(1)%T(m2,m1,is2,is1,na1,na1) = CONJG(n(1)%T(m1,m2,is1,is2,na1,na1))
         ENDDO
         ENDDO
         !
      ENDDO
      ENDDO
      !
   !ENDDO
   ENDDO
   !
   IF ( time_reversal ) THEN
   ! in a nonmagnetic case, considering the time reversal
   !
   DO i = 1, ehub_nnnc 
      !
      n11 = 0.0_DP
      n12 = 0.0_DP
      n21 = 0.0_DP
      n22 = 0.0_DP
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
         DO m1 = 1, ehub_lm(nt1)
         DO m2 = 1, ehub_lm(nt2)
         !
         n11 = n(i)%T(m1,m2,1,1,na1,na2)
         n12 = n(i)%T(m1,m2,1,2,na1,na2)
         n21 = n(i)%T(m1,m2,2,1,na1,na2)
         n22 = n(i)%T(m1,m2,2,2,na1,na2)
         !
         n(i)%T(m1,m2,1,1,na1,na2) = 0.5_DP*(n11+CONJG(n22))
         n(i)%T(m1,m2,1,2,na1,na2) = 0.5_DP*(n12-CONJG(n21))
         n(i)%T(m1,m2,2,1,na1,na2) = -CONJG(n(i)%T(m1,m2,1,2,na1,na2))
         n(i)%T(m1,m2,2,2,na1,na2) = CONJG(n(i)%T(m1,m2,1,1,na1,na2))
         !
         ENDDO ! m2
         ENDDO ! m1
         !
      ENDDO ! na2
      ENDDO ! na1
      !
   ENDDO ! ehub_nnnc
   !
   ENDIF
   !
   DO na1 = 1, nat
   DO na2 = 1, nat
      !
      DO is1 = 1, npol
      Do is2 = 1, npol
         IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) &
            n(1)%T(:,:,is1,is2,na1,na2) = 0.0_DP
      ENDDO
      ENDDO
      !
   ENDDO
   ENDDO
   !
   DO i = 1, ehub_nnnc
      n(i)%O = n(i)%T
      DEALLOCATE (n(i)%T)
   ENDDO
 
END SUBROUTINE symmetrize_n5_nc
!


