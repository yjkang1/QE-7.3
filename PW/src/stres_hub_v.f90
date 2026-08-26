!
! Copyright (C) 2002-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stres_hub_v (sigmah)
   !----------------------------------------------------------------------
   ! This routines computes the Hubbard contribution to the internal stress
   ! tensor. It gives in output the array sigmah(i,j) which corresponds to
   ! the quantity \( -(1/\omega)dE_h/\epsilon_{i,j} \)
   !
   USE kinds,                ONLY : DP
   USE wavefunctions,        ONLY : evc
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp
   USE cell_base,            ONLY : omega, at, bg
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE klist,                ONLY : nks, xk, ngk, igk_k
   USE basis,                ONLY : natomwfc, wfcatom, swfcatom
   USE io_files,             ONLY : nwordwfc, iunwfc
   USE buffers,              ONLY : get_buffer
   USE scf,                  ONLY : v, rho
   USE lsda_mod,             ONLY : nspin
   USE symme,                ONLY : symmatrix
   USE io_global,            ONLY : stdout
   USE mp_pools,             ONLY : inter_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE control_flags,        ONLY : gamma_only, offload_type
   USE ldaU,                 ONLY : Hubbard_projectors
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, ehub_csize,     &
                                    ehub_nn_vect, ehub_nn_dim, ehub_nnnc, &
                                    wfcv, nwfcv, ofstv, copy_v_wfc,       &
                                    ehub_l_choice, ehub_l, acbn0_type
   USE mp_bands,             ONLY : use_bgrp_in_hpsi, intra_bgrp_comm
   USE noncollin_module,     ONLY : noncolin, npol
   USE force_mod,            ONLY : eigenval, eigenvect, overlap_inv, at_dy, &
                                    at_dj, us_dy, us_dj
   USE uspp_init,            ONLY : init_us_2, gen_us_dj, gen_us_dy 
   IMPLICIT NONE
   !
   ! the Hubbard contribution to stresses
   !
   REAL(DP), INTENT(OUT) :: sigmah(3,3) 
   !
   ! ... local variables
   !
   TYPE(bec_type) :: proj
   COMPLEX(DP), ALLOCATABLE :: proj_nc(:,:)

   INTEGER :: ipol, jpol, is, ik, npw, i
   INTEGER :: m1, m2, na1, na2, nt1, nt2
   INTEGER :: ics, nnc
   INTEGER :: n1, n2, n3
   INTEGER :: nb_s, nb_e, mykey
   INTEGER :: ehub_l_max_
   INTEGER :: ehub_l_choice_(ntyp,6)
   INTEGER :: ehub_lm_(ntyp)

   REAL(DP) :: xyz(3,3)

   COMPLEX(DP), ALLOCATABLE :: dns(:,:,:,:,:,:)
   COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
   LOGICAL :: save_flag
   !
   CALL start_clock('stres_hub_v')
   !
   save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
   !
   IF ( Hubbard_projectors /= "ortho-atomic" ) CALL errore("stres_hub", &
                   "stress for this Hubbard_projectors type not implemented",1)
   IF( .NOT.(acbn0_type == 2 .or. acbn0_type == 3)) CALL errore("stres_hub_v", &
         " acbn0_type should be assigned to 2 or 3 ",1)
   IF( gamma_only ) CALL errore("stres_hub_v", " gamma_only calculation is not implemented ",1)
   !
   ehub_l_max_ = MAXVAL(ehub_l(:,:))
   ehub_l_choice_ = ehub_l_choice
   ehub_lm_ = ehub_lm
   sigmah(:,:) = 0.D0
   !
   ALLOCATE ( wfcatom(npwx*npol,natomwfc) )
   ALLOCATE ( spsi(npwx*npol,nbnd) )
   ALLOCATE ( at_dy(npwx*npol,natomwfc), at_dj(npwx*npol, natomwfc))
   IF (okvan) ALLOCATE (us_dy(npwx,nkb), us_dj(npwx,nkb))
   ALLOCATE ( swfcatom(npwx*npol, natomwfc))
   ALLOCATE ( eigenval(natomwfc))
   ALLOCATE ( eigenvect(natomwfc, natomwfc))
   ALLOCATE ( overlap_inv(natomwfc,natomwfc))
   !
   IF (noncolin) THEN
     ALLOCATE( proj_nc(nwfcv, nbnd))
   ELSE
     CALL allocate_bec_type(nwfcv, nbnd, proj)
   ENDIF
   !
   ! poor-man parallelization over bands
   ! - if nproc_pool=1     : nb_s=1, nb_e=nbnd, mykey=0
   ! - if nproc_pool<=nbnd : each processor calculates band nb_s to nb_e; mykey=0
   ! - if nproc_pool>nbnd  : each processor takes care of band na_s=nb_e;
   !   mykey labels how many times each band appears (mykey=0 first time etc.)
   !
   CALL block_distribute(nbnd, me_pool, nproc_pool, nb_s, nb_e, mykey)
   !
   ALLOCATE ( dns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc) )
   !
   ics = ehub_csize + 1
   !
   ! We start a loop on k points
   !
   DO ik = 1, nks
      !
      IF ( lsda ) current_spin = isk(ik)
      npw = ngk(ik)
      !
      IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE.)
      CALL allocate_bec_type ( nkb, nbnd, becp)
      CALL calbec(offload_type, npw, vkb, evc, becp)
      CALL s_psi(npwx, npw, nbnd, evc, spsi)
      CALL deallocate_bec_type(becp)
      !
      CALL orthowfcv2(ik)
      ! proj=<wfcv|S|evc> - no need to read S*wfcv from buffer
      IF (noncolin) THEN
         CALL ZGEMM ('C', 'N', nwfcv, nbnd, npwx*npol, (1.0_DP, 0.0_DP), wfcv, &
                    npwx*npol, spsi, npwx*npol, (0.0_DP, 0.0_DP),  proj_nc, nwfcv)
         CALL mp_sum( proj_nc( :, 1:nbnd ), intra_bgrp_comm )
      ELSE
         !$acc data create(proj)
         CALL calbec(offload_type, npw, wfcv, spsi, proj)
         !$acc update self(projkd)
         !$acc end data
      ENDIF
      !
      ! Compute derivatives of spherical harmonics and spherical Bessel fnctions
      !
      ! xyz are three unit vectors in the x,y,z, directions
      xyz(:,:) = 0.d0
      DO i =1,3
        xyz(i,i) = 1.d0
      ENDDO
      ! The derivative of spherical Bessel functions (for atomic functions)
      CALL gen_at_dj ( ik, at_dj)
      ! The derivative of spherical Bessel functions (for beta functions)
      IF (okvan) CALL gen_us_dj( ik, us_dj)
      !
      ! NB: both ipol and jpol must run from 1 to 3 because this stress
      !     contribution is not in general symmetric when computed only
      !     from k-points in the irreducible wedge of the BZ.
      !     It is (must be) symmetric after symmetrization but this requires
      !     the full stress tensor not only its upper triangular part.
      !
      DO ipol = 1, 3
      !
      ! The derivative of the spherical harmonic
      !
      CALL gen_at_dy (ik, xyz(1,ipol), at_dy)
      ! The derivative of spherical harmonics (for beta functions)
      IF ( okvan) CALL gen_us_dy (ik, xyz(1,ipol), us_dy)
      !
      DO jpol = 1, 3
         !
         ! Compute the derivative of the occupation matrix w.r.t epsilon
         !
         IF ( noncolin ) THEN
            CALL dndepsilon_k_v(ipol,jpol,proj_nc,spsi,ik,nb_s,nb_e,mykey,dns)
         ELSE
            CALL dndepsilon_k_v(ipol,jpol,proj%k,spsi,ik,nb_s,nb_e,mykey,dns)
         ENDIF
         !
         DO is = 1, nspin
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
                  IF (noncolin) THEN
                     sigmah(ipol,jpol) = sigmah(ipol,jpol) &
                                       - DBLE(v%ne0_nc(m1,m2,is,na1,na2,1) &
                                       * dns(m1,m2,is,na1,na2,1))
                  ELSE
                     sigmah(ipol,jpol) = sigmah(ipol,jpol) &
                                       - DBLE(v%ne0(m1,m2,is,na1,na2,1) &
                                       * dns(m1,m2,is,na1,na2,1))
                  ENDIF
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
                  IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2 ) == 0 ) CYCLE
                  !
                  DO m1 = 1, ehub_lm(nt1)
                  DO m2 = 1, ehub_lm(nt2)
                     !
                     IF (noncolin) THEN
                        sigmah(ipol,jpol) = sigmah(ipol,jpol)  &
                                          - DBLE(v%ne0_nc(m1,m2,is,na1,na2,nnc) &
                                          * dns(m1,m2,is,na1,na2,nnc))
                     ELSE
                        sigmah(ipol,jpol) = sigmah(ipol,jpol)  &
                                          - DBLE(v%ne0(m1,m2,is,na1,na2,nnc) &
                                          * dns(m1,m2,is,na1,na2,nnc))
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
         ENDDO ! is
         ! 
      ENDDO ! jpol
      ENDDO ! ipol
      !
   ENDDO ! ik
   !
   CALL mp_sum(sigmah, inter_pool_comm)
   !
   sigmah(:,:) = sigmah(:,:) / omega
   !
   IF ( nspin == 1 ) sigmah(:,:) = 2.D0 * sigmah(:,:)
   !
   ! Symmetrization
   ! 
   CALL symmatrix(sigmah)
   !
   ! Impose symmetry s(i,j) = s(j,i) to the stress tensor.
   ! It should NOT be needed, but let us do it for safety.
   !
   DO ipol = 1, 3
   DO jpol = ipol, 3
      IF ( ABS(sigmah(ipol,jpol)-sigmah(jpol,ipol)) > 1.D-6 ) THEN
         WRITE ( stdout,'(2i3,2f12.7)' ) ipol,jpol,sigmah(ipol,jpol), &
                                                   sigmah(jpol,ipol)
         CALL errore('stres_hub_v','non-symmetric stress contribution',1)
      ENDIF
      sigmah(ipol,jpol) = 0.5D0 * (sigmah(ipol,jpol) + sigmah(jpol,ipol))
      sigmah(jpol,ipol) = sigmah(ipol,jpol)
   ENDDO
   ENDDO
   !
   IF (noncolin) THEN
     DEALLOCATE(proj_nc)
   ELSE
     CALL deallocate_bec_type(proj)
   ENDIF
   IF ( ALLOCATED(dns) )  DEALLOCATE (dns)
   DEALLOCATE ( wfcatom )
   DEALLOCATE ( spsi )
   DEALLOCATE ( at_dy, at_dj)
   IF (okvan) DEALLOCATE (us_dy, us_dj)
   DEALLOCATE ( swfcatom)
   DEALLOCATE ( eigenval)
   DEALLOCATE ( eigenvect)
   DEALLOCATE ( overlap_inv)
   !
   use_bgrp_in_hpsi = save_flag
   !
   CALL stop_clock('stres_hub_v')
   !
   RETURN
   !
END  SUBROUTINE stres_hub_v
!
!
SUBROUTINE dndepsilon_k_v (ipol,jpol,proj,spsi,ik,nb_s,nb_e,mykey,dns)
   !-----------------------------------------------------------------------
   ! This routine computes the derivative of the ns atomic occupations with
   ! respect to the strain epsilon(ipol,jpol) used to obtain the generalized
   ! Hubbard contribution to the internal stres tensor.
   !
   USE kinds,             ONLY : DP
   USE constants,         ONLY : tpi
   USE ions_base,         ONLY : nat, ityp
   USE cell_base,         ONLY : at
   USE klist,             ONLY : ngk, xk
   USE lsda_mod,          ONLY : nspin, current_spin
   USE wvfct,             ONLY : nbnd, npwx, wg
   USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type
   USE mp_pools,          ONLY : intra_pool_comm
   USE mp,                ONLY : mp_sum
   USE ldaV,              ONLY : ehub_lm, ehub_lm_max, ehub_csize,     &
                                 ehub_nn_vect, ehub_nn_dim, ehub_nnnc, &
                                 nwfcv, ofstv, ehub_l, counter_lm
   USE noncollin_module,  ONLY : noncolin, npol

   IMPLICIT NONE
   !
   ! I/O variables
   !
   INTEGER, INTENT(IN)       :: ipol, jpol, ik
   INTEGER, INTENT(IN)       :: nb_s, nb_e, mykey

   COMPLEX (DP), INTENT(IN)  :: proj(nwfcv,nbnd)
   COMPLEX (DP), INTENT(IN)  :: spsi(npwx*npol,nbnd)

   COMPLEX (DP), INTENT(OUT) :: dns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
   !
   ! ... local variables
   !
   TYPE (bec_type)           :: dproj
   COMPLEX(DP), ALLOCATABLE  :: dproj_nc(:,:)

   COMPLEX(DP)               :: kphase

   REAL(DP)                  :: arg

   INTEGER                   :: ibnd, is, npw
   INTEGER                   :: ics, nnc
   INTEGER                   :: m1, m2
   INTEGER                   :: na1, na2, nt1, nt2
   INTEGER                   :: n1, n2, n3
   INTEGER                   :: is1, is2
   !
   IF (noncolin) THEN
     ALLOCATE (dproj_nc(nwfcv, nbnd))
   ELSE
     CALL allocate_bec_type (nwfcv, nbnd, dproj)
   ENDIF
   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
   !
   dns(:,:,:,:,:,:) = (0.0_dp, 0.0_dp)
   !
   npw = ngk(ik)
   !
   ! Calculate the first derivative of proj with respect to epsilon(ipol,jpol)
   !
   IF (noncolin) THEN
     CALL dprojdepsilon_k_v (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj_nc)
   ELSE
     CALL dprojdepsilon_k_v (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj%k)
   ENDIF
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !

   ics = ehub_csize +1

   IF ( mykey /= 0 ) GO TO 10
   !
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
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt2)
         !
         DO ibnd = nb_s, nb_e
            !
            IF (noncolin) THEN
               !
               DO is1 = 1, npol
               DO is2 = 1, npol
                  dns(m1,m2,npol*(is1-1)+is2,na1,na2,1) = &
                  dns(m1,m2,npol*(is1-1)+is2,na1,na2,1) + &
                  wg(ibnd,ik) * &
              dcmplx(proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) &
           * CONJG(dproj_nc(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd)) + &
                   dproj_nc(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) &
           * CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd)))
               ENDDO ! is2
               ENDDO ! is1
               !
            ELSE
               dns(m1,m2,current_spin,na1,na2,1) = &
               dns(m1,m2,current_spin,na1,na2,1) + &
                                     wg(ibnd,ik) * &
                     ((proj(ofstv(na2)+m2,ibnd)) * CONJG(dproj%k(ofstv(na1)+m1,ibnd)) + &
                   (dproj%k(ofstv(na2)+m2,ibnd)) * CONJG(proj(ofstv(na1)+m1,ibnd)))
            ENDIF
         ENDDO !ibnd
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
            DO ibnd = nb_s, nb_e
              IF (noncolin) THEN
                DO is1 = 1, npol
                DO is2 = 1, npol
                   dns(m1,m2,npol*(is1-1)+is2,na1,na2,nnc) = &
                   dns(m1,m2,npol*(is1-1)+is2,na1,na2,nnc) + &
                   wg(ibnd,ik) * dcmplx( kphase *( &
                   proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) &
           * CONJG(dproj_nc(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd)) + &
                   dproj_nc(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) &
           * CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))))
                ENDDO ! is2
                ENDDO ! is1
              ELSE
                dns(m1,m2,current_spin,na1,na2,nnc) = &
                dns(m1,m2,current_spin,na1,na2,nnc) + &
                          kphase*wg(ibnd,ik) * &
                        ((proj(ofstv(na2)+m2,ibnd)) * CONJG(dproj%k(ofstv(na1)+m1,ibnd)) + &
                      (dproj%k(ofstv(na2)+m2,ibnd)) * CONJG(proj(ofstv(na1)+m1,ibnd)))
              ENDIF
            ENDDO ! ibnd
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
10 CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   IF ( nspin == 1 ) dns = 0.5D0 * dns
   !
   IF (noncolin) THEN
     DEALLOCATE(dproj_nc)
   ELSE
     CALL deallocate_bec_type (dproj)
   ENDIF
   !
   RETURN
   !
END SUBROUTINE dndepsilon_k_v
!
SUBROUTINE dprojdepsilon_k_v (spsi, ik, ipol, jpol, nb_s, nb_e, mykey, dproj)
   !-----------------------------------------------------------------------
   ! This routine computes the first derivative of the projection
   ! \( \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle \)
   ! with respect to the strain \( \epsilon(i,j) \).
   !
   ! We remember that: $$ \text{ns}_{I,s,m_1,m_2} = \sum_{k,v}
   ! f_{kv} \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle
   ! \langle\psi_{k,v,s}|S|\phi^{at}_{I,m2}\rangle $$
   !
   USE kinds,                ONLY : DP
   USE cell_base,            ONLY : tpiba
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk, igk_k, ngk
   USE ldaV,                 ONLY : wfcv, nwfcv, ofst, ofstv, ehub_lm
   USE lsda_mod,             ONLY : lsda, nspin, isk
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at, okvan
   USE uspp_param,           ONLY : upf, nhm, nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : becp, calbec
   USE basis,                ONLY : natomwfc, wfcatom, swfcatom
   USE force_mod,            ONLY : eigenval, eigenvect, overlap_inv, at_dy, &
   at_dj
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE noncollin_module,     ONLY : noncolin, npol
   USE control_flags,        ONLY : offload_type

   IMPLICIT NONE
   !
   ! I/O variables 
   !
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx*npol,nbnd)
   !! S|evc>
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: ipol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: jpol
   !! component index 1 to 3
   INTEGER, INTENT(IN) :: nb_s
   !! band number start
   INTEGER, INTENT(IN) :: nb_e
   !! band number end
   INTEGER, INTENT(IN) :: mykey
   !! labels how many times each band appears (mykey=0 first time etc.)  
   COMPLEX(DP), INTENT(OUT) :: dproj(nwfcv,nbnd)
   !
   ! ... local variables
   !
   INTEGER                   :: na, nt
   INTEGER                   :: i, ibnd
   INTEGER                   :: ig, ijkb0, ih, jh, iwf, npw, m1, m2


   REAL(DP)                  :: q, a1, a2

   REAL(DP), PARAMETER       :: eps = 1.0D-8

   COMPLEX (DP), ALLOCATABLE :: &
       dproj0(:,:),      & ! derivative of the projector
       dproj_us(:,:),    & ! USPP contribution to dproj0
       dwfc(:,:),        & ! the derivative of the (ortho-atomic) wavefunction
       doverlap(:,:),    & ! derivative of the overlap matrix
       doverlap_us(:,:), & ! USPP contribution to doverlap
       doverlap_inv(:,:)   ! derivative of (O^{-1/2})_JI (note the transposition)
   !
   REAL (DP), ALLOCATABLE    :: &
       gk(:,:), & ! k+G
       qm1(:),  & ! 1/|k+G|
       a1_temp(:), a2_temp(:)
   COMPLEX (DP) :: temp, temp2
   !
   ! Number of plane waves at the k point with the index ik
   npw = ngk(ik)
   !
   dproj(:,:) = (0.D0, 0.D0)
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx*npol,nwfcv) )
   ALLOCATE (a1_temp(npw), a2_temp(npw))
   dwfc(:,:) = (0.d0, 0.d0)
   !
   DO ig = 1, npw
      !
      gk(1,ig) = (xk(1,ik) + g(1,igk_k(ig,ik))) * tpiba
      gk(2,ig) = (xk(2,ik) + g(2,igk_k(ig,ik))) * tpiba
      gk(3,ig) = (xk(3,ik) + g(3,igk_k(ig,ik))) * tpiba
      !
      q = SQRT(gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2)
      !
      IF ( q > eps ) THEN
         qm1(ig) = 1.D0/q
      ELSE
         qm1(ig) = 0.D0
      ENDIF
      !
      ! (k+G)_jpol
      a1_temp(ig) = -gk(jpol,ig)
      !
      ! (k+G)_ipol * (k+G)_jpol / |k+G|
      a2_temp(ig) = -gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
      !
   ENDDO !ig
   !
   DO na = 1, nat
      nt = ityp(na)
      IF(ehub_lm(nt) >= 1) THEN 
          DO m1 = 1, ehub_lm(nt)*npol
             DO ig = 1, npw
                temp = (0.0_DP, 0.0_DP)
                temp2 = (0.0_DP, 0.0_DP)
                DO m2 = 1, natomwfc
                   temp = temp + overlap_inv(ofst(na)+m1,m2) * &
                       ( at_dy(ig,m2) * a1_temp(ig) + at_dj(ig,m2) * a2_temp(ig) )
                   IF (noncolin) temp2 = temp2 + overlap_inv(ofst(na)+m1,m2) * &
                            ( at_dy(ig+npwx,m2) * a1_temp(ig) + at_dj(ig+npwx,m2) * a2_temp(ig) )
                ENDDO !m2
                dwfc(ig, ofstv(na)+m1) = dwfc(ig, ofstv(na)+m1) + temp
                IF (noncolin) dwfc(ig+npwx, ofstv(na)+m1) = &
                         dwfc(ig+npwx, ofstv(na)+m1) + temp2
             ENDDO !ig
          ENDDO !m1
      ENDIF
   ENDDO !na 
   !
   ! The diagonal term
   IF ( ipol == jpol ) THEN
      dwfc(1:npw,:) = dwfc(1:npw,:) - wfcv(1:npw,:) * 0.5D0
      IF (noncolin) dwfc(1+npwx:npwx+npw,:) = dwfc(1+npwx:npwx+npw,:) - wfcv(1+npwx:npwx+npw,:) * 0.5D0
   ENDIF
   !
   ! 2. Contribution due to the derivative of (O^{-1/2})_JI which
   !    is multiplied by atomic wavefunctions (only for ortho-atomic case)
   !
   !
   ! Compute the derivative dO_IJ/d\epsilon(ipol,jpol)
   !
   ALLOCATE (doverlap(natomwfc,natomwfc))
   ALLOCATE (doverlap_inv(natomwfc,natomwfc))
   doverlap(:,:) = (0.0_DP, 0.0_DP)
   doverlap_inv(:,:) = (0.0_DP, 0.0_DP)
   !
   ! Calculate:
   ! doverlap = < dphi_I/d\epsilon(ipol,jpol) | S | phi_J > 
   !            + < phi_I | S | dphi_J/d\epsilon(ipol,jpol) >
   !
   DO m2 = 1, natomwfc
      DO m1 = 1, natomwfc
         temp = (0.0_DP,0.0_DP)
         temp2 = (0.0_DP,0.0_DP)
         DO ig = 1, npw
            temp = temp + CONJG((at_dy(ig,m1)*a1_temp(ig) + at_dj(ig,m1)*a2_temp(ig))) * swfcatom(ig,m2) &
                        + CONJG(swfcatom(ig,m1)) * (at_dy(ig,m2)*a1_temp(ig) + at_dj(ig,m2)*a2_temp(ig))
            IF (noncolin) THEN
               temp2 = temp2 + CONJG((at_dy(ig+npwx,m1)*a1_temp(ig) + at_dj(ig+npwx,m1)*a2_temp(ig))) * swfcatom(ig+npwx,m2) &
                             + CONJG(swfcatom(ig+npwx,m1)) * (at_dy(ig+npwx,m2)*a1_temp(ig) + at_dj(ig+npwx,m2)*a2_temp(ig))
            ENDIF
         ENDDO
         doverlap(m1,m2) = doverlap(m1,m2) + temp
         IF (noncolin) doverlap(m1,m2) = doverlap(m1,m2) + temp2
      ENDDO
   ENDDO
   !
   IF (ipol.EQ.jpol) THEN
      DO m2 = 1, natomwfc
         DO m1 = 1, natomwfc
            temp = (0.0_DP,0.0_DP)
            temp2 = (0.0_DP,0.0_DP)
            !$acc parallel loop reduction(+:temp)
            DO ig = 1, npw
               temp = temp + CONJG((-wfcatom(ig,m1)*0.5d0)) * swfcatom(ig,m2) &
                           + CONJG(swfcatom(ig,m1)) * (-wfcatom(ig,m2)*0.5d0)
               IF (noncolin) THEN
                  temp2 = temp2 + CONJG((-wfcatom(ig+npwx,m1)*0.5d0)) * swfcatom(ig+npwx,m2) &
                           + CONJG(swfcatom(ig+npwx,m1)) * (-wfcatom(ig+npwx,m2)*0.5d0)
               ENDIF
            ENDDO
            doverlap(m1,m2) = doverlap(m1,m2) + temp
            IF (noncolin) doverlap(m1,m2) = doverlap(m1,m2) + temp2
         ENDDO
      ENDDO
   ENDIF        
   !
   ! Sum over G vectors
   CALL mp_sum( doverlap, intra_bgrp_comm )
   !
   ! USPP term in dO_IJ/d\epsilon(ipol,jpol)
   !
   IF (okvan) THEN
      ! Calculate doverlap_us = < phi_I | dS/d\epsilon(ipol,jpol) | phi_J >
      ALLOCATE (doverlap_us(natomwfc,natomwfc))
      CALL matrix_element_of_dSdepsilon_v (ik, ipol, jpol, &
           natomwfc, wfcatom, natomwfc, wfcatom, doverlap_us, 1, natomwfc, 0, .false.)
      ! Sum up the "normal" and ultrasoft terms
      DO m1 = 1, natomwfc
         DO m2 = 1, natomwfc
            doverlap(m1,m2) = doverlap(m1,m2) + doverlap_us(m1,m2)
         ENDDO
      ENDDO
      DEALLOCATE (doverlap_us)
   ENDIF
   !
   ! Now compute dO^{-1/2}_JI/d\epsilon(ipol,jpol) using dO_IJ/d\epsilon(ipol,jpol)
   ! Note the transposition!
   ! 
   CALL calculate_doverlap_inv (natomwfc, eigenval, eigenvect, &
                                  doverlap, doverlap_inv)
   !
   ! Now compute \sum_J dO^{-1/2}_JI/d\epsilon(ipol,jpol) \phi_J
   ! and add it to another term (see above)
   ! Note, doverlap_inv is d(O^{-1/2}) not transposed. The transposition
   ! of d(O^{-1/2}) is taken into account via a proper usage of the order
   ! of indices in doverlap_inv:
   ! dwfc(ig,ofst(na)+m1) = dwfc(ig,ofst(na)+m1) + wfcatom(ig,m2) * doverlap_inv(m2,ofst(na)+m1)
   ! where m1=1,ehub_lm(nt); m2=1,natomwfc; ig=1,npw
   !
   !DO na = 1, nat
   !   nt = ityp(na)
   !   IF (ehub_lm(nt) >= 1 ) THEN
   !      !DO m1 = 1, ehub_lm(nt)
   !      !   DO m2 = 1, natomwfc
   !      !      DO ig = 1, npw
   !      !         dwfc(ig,ofstv(na)+m1) = dwfc(ig,ofstv(na)+m1) + &
   !      !                       doverlap_inv(ofst(na)+m1,m2) * wfcatom(ig,m2)
   !      !      ENDDO !ig
   !      !   ENDDO !m2
   !      !ENDDO !m1
   !      CALL ZGEMM('N','N', npw, ehub_lm(nt), natomwfc, (1.d0,0.d0), &
   !           wfcatom, npwx, doverlap_inv(:,ofst(na)+1:ofst(na)+ehub_lm(nt)), &
   !           natomwfc, (1.d0,0.d0), dwfc(:,ofstv(na)+1:ofstv(na)+ehub_lm(nt)), npwx)
   !   ENDIF
   !ENDDO !na
   !
   DO na = 1, nat
      nt = ityp(na)
      IF (ehub_lm(nt) >= 1) THEN
         CALL MYZGEMM('N','N', npwx*npol, ehub_lm(nt)*npol, natomwfc, (1.d0,0.d0), &
               wfcatom, npwx*npol, doverlap_inv(:,ofst(na)+1:ofst(na)+ehub_lm(nt)*npol), &
               natomwfc, (1.d0,0.d0), dwfc(:,ofstv(na)+1:ofstv(na)+ehub_lm(nt)*npol), npwx*npol)
      ENDIF
   ENDDO
   !
   DEALLOCATE (doverlap)
   DEALLOCATE (doverlap_inv)
   !
   ! Compute dproj = <dwfc|S|psi> = <dwfc|spsi>
   IF (noncolin) THEN
      CALL ZGEMM('C','N', nwfcv, nbnd, npwx*npol, (1.d0,0.d0), &
            dwfc, npwx*npol, spsi, npwx*npol, (0.d0,0.d0), &
            dproj, nwfcv)   
      CALL mp_sum( dproj, intra_bgrp_comm )
   ELSE   
      CALL calbec( offload_type, npw, dwfc, spsi, dproj )
   ENDIF
   !
   DEALLOCATE ( dwfc, qm1, gk)
   DEALLOCATE(a1_temp, a2_temp)
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\phi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   IF (okvan) THEN
      ALLOCATE(dproj_us(nwfcV,nb_s:nb_e))
      CALL matrix_element_of_dSdepsilon_v (ik, ipol, jpol, &
                         nwfcv, wfcv, nbnd, evc, dproj_us, nb_s, nb_e, mykey, .true.)
      ! dproj + dproj_us
      DO m1 = 1, nwfcv
         dproj(m1,nb_s:nb_e) = dproj(m1,nb_s:nb_e) + dproj_us(m1,:)
      ENDDO
      DEALLOCATE(dproj_us)
   ENDIF
   !
   RETURN
   !
END SUBROUTINE dprojdepsilon_k_v
!
SUBROUTINE matrix_element_of_dSdepsilon_v (ik, ipol, jpol, lA, A, lB, B, A_dS_B, lB_s, lB_e, mykey, flag)
   !
   ! This routine computes the matrix element < A | dS/d\epsilon(ipol,jpol) | B >
   ! Written by I. Timrov (2020)
   !
   ! Compute the term <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE wvfct,                ONLY : npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at, qq_so, okvan
   USE uspp_param,           ONLY : nh, upf
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : calbec
   USE control_flags,        ONLY : offload_type
   USE klist,                ONLY : xk, igk_k, ngk
   USE force_mod,            ONLY : us_dy, us_dj
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE noncollin_module,     ONLY : noncolin, npol
   USE ldav,                 ONLY : ofstv, ehub_lm, ehub_l, counter_lm
   !
   IMPLICIT NONE
   !
   ! Input/Output
   !
   INTEGER, INTENT(IN)      :: ik          ! k point
   INTEGER, INTENT(IN)      :: ipol, jpol  ! Cartesian components
   INTEGER, INTENT(IN)      :: lA, lB, lB_s, lB_e
   ! There is a possibility to parallelize over lB,
   ! where lB_s (start) and lB_e (end)
   LOGICAL, INTENT(IN)      :: flag  ! noncollinear: controlling whether 
                                     ! calculating <phi|dS|PSI> 
                                     ! or          <phi|dS|PHI> (= .true.)   
   COMPLEX(DP), INTENT(IN)  :: A(npwx*npol,lA)
   COMPLEX(DP), INTENT(IN)  :: B(npwx*npol,lB)
   COMPLEX(DP), INTENT(OUT) :: A_dS_B(lA,lB_s:lB_e)
   INTEGER,     INTENT(IN)  :: mykey
   !
   ! Local variables
   !
   INTEGER :: npw, i, nt, na, nb, ih, jh, ig, iA, iB, ijkb0, mU, mD, nt1, nt2, is1
   REAL(DP) :: gvec
   COMPLEX (DP), ALLOCATABLE :: Adbeta(:,:), Abeta(:,:), &
                                dbetaB(:,:), betaB(:,:), &
                                aux(:,:), qq(:,:)
   REAL (DP) :: q, a1, a2
   REAL (DP), PARAMETER :: eps = 1.0d-8
   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   INTEGER :: nh_nt
   REAL (DP), ALLOCATABLE :: a1_temp(:), a2_temp(:)
   !
   IF (.NOT.okvan) RETURN
   !
   !$acc data present_or_copyin(A,B) present_or_copyout(A_ds_B)
   !
   !$acc kernels
   A_dS_B(:,:) = (0.0d0, 0.0d0)
   !$acc end kernels
   npw = ngk(ik)
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE (a1_temp(npw), a2_temp(npw))
   !
   ! Compute k+G and 1/|k+G|
   DO ig = 1, npw
      DO i = 1, 3
         gk(i,ig) = (xk(i,ik) + g(i,igk_k(ig,ik))) * tpiba
      ENDDO
      q = SQRT(gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2)
      IF (q.GT.eps) THEN
         qm1(ig)=1.d0/q
      ELSE
         qm1(ig)=0.d0
      ENDIF
      a1_temp(ig) = -gk(jpol,ig)
      a2_temp(ig) = -gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
   ENDDO
   !
   ijkb0 = 0
   !
   !$acc data copyin(us_dj, qq_at, a1_temp, a2_temp)
   DO nt = 1, ntyp
      !
      ALLOCATE ( Adbeta(lA,npol*nh(nt)) )
      ALLOCATE ( Abeta(lA,npol*nh(nt)) )
      ALLOCATE ( dbetaB(npol*nh(nt),lB) )
      ALLOCATE ( betaB(npol*nh(nt),lB) )
      ALLOCATE ( qq(npol*nh(nt),npol*nh(nt)) )
      !
      nh_nt = nh(nt)
      !$acc data create(Adbeta,Abeta,dbetaB,betaB,qq)
      !
      DO na = 1, nat
         !
         IF ( ityp(na).EQ.nt ) THEN
            !
            IF (noncolin) THEN
               IF ( upf(nt)%has_so ) THEN
                  qq(1:nh(nt),1:nh(nt)) = qq_so(:,:,1,nt)
                  qq(1:nh(nt),1+nh(nt):nh(nt)*npol) = qq_so(:,:,2,nt)
                  qq(1+nh(nt):nh(nt)*npol,1:nh(nt)) = qq_so(:,:,3,nt)
                  qq(1+nh(nt):nh(nt)*npol,1+nh(nt):nh(nt)*npol) = qq_so(:,:,4,nt)
               ELSE
                  qq(1:nh(nt),1:nh(nt)) = qq_at(1:nh(nt),1:nh(nt),na)
                  qq(1:nh(nt),1+nh(nt):nh(nt)*npol) = (0.0,0.0)
                  qq(1+nh(nt):nh(nt)*npol,1:nh(nt)) = (0.0,0.0)
                  qq(1+nh(nt):nh(nt)*npol,1+nh(nt):nh(nt)*npol) = qq_at(1:nh(nt),1:nh(nt),na)
               ENDIF
            ELSE
               !$acc parallel loop collapse(2) 
               DO jh = 1, nh_nt 
                  DO ih = 1, nh_nt
                     qq(ih,jh) = CMPLX(qq_at(ih,jh,na), 0.0d0, kind=DP)
                  ENDDO
               ENDDO
            ENDIF
            !
            ! aux is used as a workspace
            ALLOCATE ( aux(npwx*npol,nh(nt)*npol) )
            !$acc data create(aux)
            aux=(0.0,0.0)
            !
            !$acc parallel loop collapse(2) 
            DO ih = 1, nh_nt !nh(nt)
               ! now we compute the true dbeta function
               DO ig = 1, npw
                  aux(ig,ih) = us_dy(ig,ijkb0+ih) * a1_temp(ig) + us_dj(ig,ijkb0+ih) * a2_temp(ig)
                  IF (noncolin) aux(ig+npwx,ih+nh(nt)) = us_dy(ig,ijkb0+ih) * a1_temp(ig) &
                                                + us_dj(ig,ijkb0+ih) * a2_temp(ig)  
               ENDDO
            ENDDO
            !
            IF (ipol.EQ.jpol) THEN
               !$acc parallel loop collapse(2)
               DO ih = 1, nh_nt !nh(nt)
                  DO ig = 1, npw
                     aux(ig,ih) = aux(ig,ih) - vkb(ig,ijkb0+ih)*0.5d0
                     IF (noncolin) aux(ig+npwx,ih+nh(nt)) = aux(ig+npwx,ih+nh(nt)) &
                                   - vkb(ig,ijkb0+ih)*0.5d0
                  ENDDO
               ENDDO   
            ENDIF 
            IF (noncolin) THEN
               ! Calculate betaB = <dbeta|B>
               ! dbetaB(:,1       : nh(nt))      = spin up
               ! dbetaB(:,1+nh(nt): nh(nt)*npol) = spin down
               dbetaB=(0.0,0.0)
               CALL ZGEMM ('C', 'N', nh(nt)*npol, lB, npwx*npol, (1.0_DP, 0.0_DP), aux, &
                          npwx*npol, B, npwx*npol, (0.0_DP, 0.0_DP), dbetaB, nh(nt)*npol)
               CALL mp_sum( dbetaB(:, 1:lB) , intra_bgrp_comm )               
               !
               ! Calculate Adbeta = <A|dbeta>
               ! Adbeta(:,1       : nh(nt))      = spin up
               ! Adbeta(:,1+nh(nt): nh(nt)*npol) = spin down      
               Adbeta=(0.0,0.0)
               CALL ZGEMM ('C', 'N', lA, nh(nt)*npol, npwx*npol, (1.0_DP, 0.0_DP), A, &
                          npwx*npol, aux, npwx*npol, (0.0_DP, 0.0_DP), Adbeta, lA)
               CALL mp_sum( Adbeta(:, 1:nh(nt)*npol) , intra_bgrp_comm )
            ELSE        
               ! Calculate dbetaB = <dbeta|B> 
               CALL calbec(offload_type, npw, aux, B, dbetaB )
               !
               ! Calculate Adbeta = <A|dbeta>
               CALL calbec(offload_type, npw, A, aux, Adbeta )
            ENDIF   
            !
            ! aux is now used as a work space to store vkb
            !$acc parallel loop collapse(2)
            DO ih = 1, nh_nt  !nh(nt)
               DO ig = 1, npw
                  aux(ig,ih) = vkb(ig,ijkb0+ih)
                  IF (noncolin) aux(ig+npwx,ih+nh(nt)) = vkb(ig,ijkb0+ih)
               ENDDO
            ENDDO
            !
            IF (noncolin) THEN
                ! Calculate Abeta = <A|beta>      
                ! (same as Adbeta)
                !
                Abeta=(0.0,0.0)
                CALL ZGEMM ('C', 'N', lA, nh(nt)*npol, npwx*npol, (1.0_DP, 0.0_DP), A, &
                           npwx*npol, aux, npwx*npol, (0.0_DP, 0.0_DP), Abeta, lA)
                CALL mp_sum( Abeta(:, 1:nh(nt)*npol) , intra_bgrp_comm )
                !
                ! Calculate betaB = <beta|B>
                ! (same as dbetaB)
                !
                betaB=(0.0,0.0)
                CALL ZGEMM ('C', 'N', nh(nt)*npol, lB, npwx*npol, (1.0_DP, 0.0_DP), aux, &
                           npwx*npol, B, npwx*npol, (0.0_DP, 0.0_DP), betaB, nh(nt)*npol)
                CALL mp_sum( betaB(:, 1:lB) , intra_bgrp_comm )
             ELSE
               ! Calculate Abeta = <A|beta>
               CALL calbec( offload_type, npw, A, aux, Abeta )
               !
               ! Calculate betaB = <beta|B>
               CALL calbec( offload_type, npw, aux, B, betaB )
             ENDIF
            !
            !$acc end data
            DEALLOCATE ( aux )
            !
            ALLOCATE ( aux(nh(nt)*npol, lB) )
            !$acc data create(aux)
            aux(:,:)=(0.0,0.0)            
            !
            IF (noncolin) THEN
               ! aux(:, 1:nh(nt))             = \sum_jh qq(1,jh)*dbetaB(1,jh) + qq(2,jh)*dbetaB(2,jh)
               ! aux(:, 1+nh(nt):nh(nt)*npol) = \sum_jh qq(3,jh)*dbetaB(1,jh) + qq(4,jh)*dbetaB(2,jh)
               !
               ! spin up
               CALL ZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1, 1), nh(nt)*npol, &
                          dbetaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1, lB_s), nh(nt)*npol)
               ! spin down   
               CALL ZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1+nh(nt), 1), nh(nt)*npol, &
                          dbetaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1+nh(nt), lB_s), nh(nt)*npol)
            ELSE
               ! Calculate \sum_jh qq_at(ih,jh) * dbetaB(jh)     
               !$acc host_data use_device(qq,dbetaB,aux)
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                            qq, nh(nt), dbetaB(1,lB_s),    nh(nt), (0.0d0,0.0d0), &
                            aux(1,lB_s), nh(nt))
               !$acc end host_data
            ENDIF
            !$acc kernels
            dbetaB(:,:) = aux(:,:)
            !$acc end kernels
            !
            IF (noncolin) THEN
               ! (same as dbetaB)     
               !
               ! spin up 
               CALL ZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1, 1), nh(nt)*npol, &
                          betaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1, lB_s), nh(nt)*npol)
               ! spin down   
               CALL ZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                          qq(1+nh(nt), 1), nh(nt)*npol, &
                          betaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                          aux(1+nh(nt), lB_s), nh(nt)*npol)
            ELSE
               ! Calculate \sum_jh qq_at(ih,jh) * betaB(jh)
               !$acc host_data use_device(qq,betaB,aux)
               CALL MYZGEMM('N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                            qq, nh(nt), betaB(1,lB_s),     nh(nt), (0.0d0,0.0d0), &
                            aux(1,lB_s), nh(nt))
              !$acc end host_data
            ENDIF
            !
            !$acc kernels
            betaB(:,:) = aux(:,:)
            !$acc end kernels
            !
            !$acc end data
            DEALLOCATE ( aux )
            !
            ijkb0 = ijkb0 + nh(nt)
            !
            ! A_dS_B(iA,iB) = \sum_ih [Adbeta(iA,ih) * betapsi(ih,iB) +
            !                          Abeta(iA,ih)  * dbetaB(ih,iB)] 
            ! Only A_dS_B(:,lB_s:lB_e) are calculated
            !
            IF ( mykey == 0 ) THEN
               IF (noncolin) THEN
                  nt1 = nh(nt) + 1
                  IF ( flag ) THEN
                     DO nb = 1, nat
                        nt2 = ityp(nb)
                        IF ( ehub_lm(nt2) > 0 ) THEN
                           mU = ofstv(nb) + 1
                           mD = mU + ehub_lm(nt2)
                           !
                           ! actually, if we consider more than one Hubard channel, the order of wfcatom is
                           ! 121112221111122222.... (1: spin up, 2: spin down) 
                           ! temp
                           DO ih = 1, ehub_lm(nt2)
                             DO jh = lB_s, lB_e
                               DO ig = 1, nh(nt)
                                 DO is1 = 1, npol
                                   A_dS_B(mU-1+ih+counter_lm(ih,nt2)+(2*(ehub_l(ih,nt2)-1)+1)*(is1-1), jh) &
                                 = A_dS_B(mU-1+ih+counter_lm(ih,nt2)+(2*(ehub_l(ih,nt2)-1)+1)*(is1-1), jh) &
                                 + Adbeta(mU-1+ih+counter_lm(ih,nt2)+(2*(ehub_l(ih,nt2)-1)+1)*(is1-1), &
                                   ig+nh(nt)*(is1-1)) * betaB(ig+nh(nt)*(is1-1), jh) &
                                 + Abeta(mU-1+ih+counter_lm(ih,nt2)+(2*(ehub_l(ih,nt2)-1)+1)*(is1-1), &
                                   ig+nh(nt)*(is1-1)) * dbetaB(ig+nh(nt)*(is1-1), jh) 
                                 ENDDO
                               ENDDO
                             ENDDO
                           ENDDO
                           ! temp
                           !! spin up
                           !CALL ZGEMM('N', 'N', ehub_lm(nt2), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                           !            Adbeta(mU,1), lA, &
                           !            betaB(1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                           !            A_dS_B(mU, lB_s), lA)
                           !CALL ZGEMM('N', 'N', ehub_lm(nt2)), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                           !         Abeta(mU,1), lA, &
                           !         dbetaB(1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                           !         A_dS_B(mU, lB_s), lA)
                           !! spin down
                           !CALL ZGEMM('N', 'N', ehub_lm(nt2), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                           !            Adbeta(mD, nt1), lA, &
                           !            betaB(nt1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                           !            A_dS_B(mD,lB_s), lA)
                           !CALL ZGEMM('N', 'N', ehub_lm(nt2), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                           !            Abeta(mD, nt1), lA, &
                           !            dbetaB(nt1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                           !            A_dS_B(mD,lB_s), lA)
                        ENDIF ! ehub_lm
                     ENDDO
                  ELSEIF ( .not. flag ) THEN
                     CALL ZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                                  Adbeta, lA, betaB(1,lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                  A_dS_B(1,lB_s), lA)
                     CALL ZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                                  Abeta, lA, dbetaB(1,lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                                  A_dS_B(1,lB_s), lA)
                  ENDIF ! flag
               ELSE
                  !$acc host_data use_device(Adbeta,betaB,Abeta,dbetaB,A_dS_B)      
                  CALL MYZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                             Adbeta, lA, betaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                             A_dS_B(1,lB_s), lA)
                  CALL MYZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                             Abeta, lA, dbetaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                             A_dS_B(1,lB_s), lA)
                  !$acc end host_data   
                  !
               ENDIF !noncolin
            ENDIF ! mykey
         ENDIF
         !
      ENDDO
      !
      !$acc end data
      DEALLOCATE (dbetaB, betaB, Abeta, Adbeta, qq)
      ! 
   ENDDO
   !$acc end data
   !
   DEALLOCATE (a1_temp, a2_temp)
   DEALLOCATE ( qm1, gk )
   !
   !$acc end data
   !
   RETURN
   !
END SUBROUTINE matrix_element_of_dSdepsilon_v

