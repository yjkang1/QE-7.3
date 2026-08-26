! Following codes are from force_hub.f90 by Iruii Timrov
!
! Copyright (C) 2002-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE force_hub_u ( forceh )
   !
   USE kinds,                ONLY : DP 
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp 
   USE cell_base,            ONLY : at, bg 
   USE basis,                ONLY : natomwfc 
   USE symme,                ONLY : symvector 
   USE io_files,             ONLY : prefix 
   USE wvfct,                ONLY : nbnd, npwx 
   USE control_flags,        ONLY : gamma_only 
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk 
   USE scf,                  ONLY : v, rho 
   USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, me_pool, &
                                    nproc_pool 
   USE mp,                   ONLY : mp_sum 
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, &
                                    deallocate_bec_type 
   USE uspp,                 ONLY : nkb, vkb, ofsbeta 
   USE uspp_param,           ONLY : nh 
   USE wavefunctions,        ONLY : evc 
   USE klist,                ONLY : nks, xk, ngk, igk_k 
   USE io_files,             ONLY : nwordwfc, iunwfc 
   USE buffers,              ONLY : get_buffer 
   USE ldaU,                 ONLY : Hubbard_projectors 
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, acbn0_type,  & 
                                    wfcv, nwfcv, ofstv, copy_v_wfc, ehub_v
   USE uspp_init,            ONLY : init_us_2
   !
   IMPLICIT NONE
   !
   REAL(DP)                 :: forceh(3,nat)
   !
   ! ... local variables
   !
   TYPE(bec_type)           :: proj

   COMPLEX(DP), ALLOCATABLE :: spsi(:,:), wfcatom(:,:) 
   COMPLEX(DP), ALLOCATABLE :: dns(:,:,:,:)

   INTEGER                  :: ics, nnc
   INTEGER                  :: is, ik, ijkb0, ipol, alpha, npw
   INTEGER                  :: m1, m2, na1, na2, nt1, nt2
   INTEGER                  :: n1, n2, n3
   INTEGER                  :: nb_s, nb_e, mykey
   !
   IF ( Hubbard_projectors /= "ortho-atomic" ) CALL errore("force_hub_u", &
                   " eACBN0 for this Hubbard_projectors type not implemented", 1 )
   ! even though Hubbard_projectors is ortho-atomic, I consider atomic calculcation in this SUBROUTINE
   CALL start_clock ('force_hub_u')
   !
   ALLOCATE ( dns(ehub_lm_max,ehub_lm_max,nspin,nat) )
   !
   ALLOCATE ( spsi(npwx,nbnd) )
   ALLOCATE ( wfcatom(npwx,natomwfc) )
   !
   CALL allocate_bec_type(nkb, nbnd, becp)
   CALL allocate_bec_type(nwfcv, nbnd, proj)
   !
   ! poor-man parallelization over bands
   ! - if nproc_pool=1     : nb_s=1, nb_e=nbnd, mykey=0
   ! - if nproc_pool<=nbnd : each processor calculates band nb_s to nb_e; mykey=0
   ! - if nproc_pool>nbnd  : each processor takes care of band na_s=nb_e;
   !   mykey labels how many times each band appears (mykey=0 first time etc.)
   !
   CALL block_distribute(nbnd, me_pool, nproc_pool, nb_s, nb_e, mykey)
   !
   forceh(:,:) = 0.D0
   !
   DO ik = 1, nks
      !
      IF (lsda) current_spin = isk(ik)
      npw = ngk(ik)
      !
      IF (nks > 1) &
         CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
      CALL calbec(npw, vkb, evc, becp)
      CALL s_psi(npwx, npw, nbnd, evc, spsi)
      !
      ! re-calculate atomic wfc - wfcatom is used here as work space
      !
      CALL atomic_wfc(ik, wfcatom)
      CALL copy_v_wfc(wfcatom)
      !
      ! wfcv contains Hubbard-U atomic wavefunctions
      ! proj=<wfcv|S|evc> - no need to read S*wfcv from buffer
      !
      CALL calbec(npw, wfcv, spsi, proj)
      !
      ! now we need the first derivative of proj with respect to tau(alpha,ipol)
      !
      DO alpha = 1, nat ! forces are calculated for atom alpha ...
         !
         ijkb0 = ofsbeta(alpha) ! positions of beta functions for atom alpha
         !
         DO ipol = 1, 3  ! forces are calculated for coordinate ipol ...
            !
            IF ( gamma_only ) THEN
               CALL dndtau_gamma_u(proj%r, spsi, alpha, ijkb0, ipol, ik, &
                                   nb_s, nb_e, mykey, dns)
            ELSE
               CALL dndtau_k_u(proj%k, spsi, alpha, ijkb0, ipol, ik, &
                               nb_s, nb_e, mykey, dns)
            ENDIF
            !
            DO is = 1, nspin
               !
               DO na1 = 1, nat
                  nt1 = ityp(na1)
                  !
                  IF ( ehub_lm(nt1) == 0 ) CYCLE
                  !
                  DO m1 = 1, ehub_lm(nt1)
                  DO m2 = 1, ehub_lm(nt1)
                  !
                  !
                     forceh(ipol,alpha) = forceh(ipol,alpha) &
                                        - DBLE(v%na0(m1,m2,is,na1) &
                                        * dns(m1,m2,is,na1))
                  !                      
                  ENDDO ! m2
                  ENDDO ! m1
                  ! 
               ENDDO ! na1
               !
            ENDDO ! is
            !
         ENDDO ! ipol
         !
      ENDDO ! alpha
      !
   ENDDO ! ik
   !
   CALL mp_sum(forceh, inter_pool_comm)
   !
   CALL deallocate_bec_type(becp)
   CALL deallocate_bec_type(proj)
   !
   DEALLOCATE( dns )
   DEALLOCATE( wfcatom ) 
   DEALLOCATE( spsi ) 
   !
   IF ( nspin == 1 ) forceh(:,:) = 2.D0 * forceh(:,:)
   !
   ! symmetrize
   !
   CALL symvector(nat, forceh)
   !
#if defined(__DEBUG)
   WRITE ( 66,'("Hubbard contribution Begin")' )
   WRITE ( 66,'(3f12.6)' ) forceh(:,:)
   WRITE ( 66,'("Hubbard contribution End")' )
#endif
   !
   CALL stop_clock('force_hub_u')
   !
   RETURN
   !
END SUBROUTINE force_hub_u
!
!-------------------------------------------------------------------------
SUBROUTINE dndtau_k_u &
            (proj, spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dns)
   !----------------------------------------------------------------------
   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : at
   USE klist,                ONLY : xk
   USE lsda_mod,             ONLY : nspin, current_spin
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, nwfcv, ofstv

   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)  :: proj(nwfcv,nbnd)
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx,nbnd)

   INTEGER, INTENT(IN)      :: alpha, jkb0
   INTEGER, INTENT(IN)      :: ipol, ik
   INTEGER, INTENT(IN)      :: nb_s, nb_e, mykey

   COMPLEX(DP), INTENT(OUT) :: dns(ehub_lm_max,ehub_lm_max,nspin,nat)
   !
   ! ... local variables
   !
   COMPLEX(DP), ALLOCATABLE :: dproj(:,:)

   INTEGER                  :: ibnd, is
   INTEGER                  :: m1, m2
   INTEGER                  :: na1, na2, nt1, nt2
   !
   CALL start_clock('dndtau_k_u')
   !
   ALLOCATE ( dproj(nwfcv,nb_s:nb_e) )
   !
   ! Compute the derivative of the generalized occupation matrices 
   ! (the quantities dnsg(m1,m2)) of the atomic orbitals. 
   ! They are complex quantities as well as nsg(m1,m2).
   !
   CALL dprojdtau_k_v(spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !
   dns(:,:,:,:) = (0.D0, 0.D0)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO m1 = 1,  ehub_lm(nt1)
      DO m2 = m1, ehub_lm(nt1)
         !
         DO ibnd = nb_s, nb_e
            dns(m1,m2,current_spin,na1) = &
            dns(m1,m2,current_spin,na1) + &
                                  wg(ibnd,ik) * &
              DBLE((proj(ofstv(na1)+m2,ibnd)) * &
             CONJG(dproj(ofstv(na1)+m1,ibnd)) + &
                  (dproj(ofstv(na1)+m1,ibnd)) * &
              CONJG(proj(ofstv(na1)+m2,ibnd)))
         ENDDO ! ibnd
         !
      ENDDO ! m1
      ENDDO ! m2
      !
   ENDDO ! na1
   !
10 DEALLOCATE ( dproj ) 
   !
   CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   IF ( nspin == 1 ) dns = 0.5D0 * dns
   !
   ! Impose hermiticity of dns_{m1,m2}
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = m1+1, ehub_lm(nt2)
         dns(m2,m1,current_spin,na1) = dns(m1,m2,current_spin,na1)
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na1
   !
   CALL stop_clock('dndtau_k_u')
   !
   RETURN
   !
END SUBROUTINE
!
!-------------------------------------------------------------------------
SUBROUTINE dndtau_gamma_u &
            (rproj, spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dns)
   !----------------------------------------------------------------------
   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : at
   USE klist,                ONLY : xk
   USE lsda_mod,             ONLY : nspin, current_spin
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, nwfcv, ofstv

   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)  :: rproj(nwfcv,nbnd)
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx,nbnd)

   INTEGER, INTENT(IN)      :: alpha, jkb0
   INTEGER, INTENT(IN)      :: ipol, ik
   INTEGER, INTENT(IN)      :: nb_s, nb_e, mykey

   COMPLEX(DP), INTENT(OUT) :: dns(ehub_lm_max,ehub_lm_max,nspin,nat)
   !
   ! ... local variables
   !
   COMPLEX(DP), ALLOCATABLE :: dproj(:,:)

   INTEGER                  :: ibnd, is
   INTEGER                  :: m1, m2
   INTEGER                  :: na1, na2, nt1, nt2
   !
   CALL start_clock('dndtau_gamma_u')
   !
   ALLOCATE ( dproj(nwfcv,nb_s:nb_e) )
   !
   ! Compute the derivative of the generalized occupation matrices 
   ! (the quantities dnsg(m1,m2)) of the atomic orbitals. 
   ! They are complex quantities as well as nsg(m1,m2).
   !
   CALL dprojdtau_gamma_v(spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !
   dns(:,:,:,:) = (0.D0, 0.D0)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO m1 = 1,  ehub_lm(nt1)
      DO m2 = m1, ehub_lm(nt1)
         !
         DO ibnd = nb_s, nb_e
            dns(m1,m2,current_spin,na1) = &
            dns(m1,m2,current_spin,na1) + &
                                  wg(ibnd,ik) * &
                                  rproj(ofstv(na1)+m2,ibnd) * &
                                  dproj(ofstv(na1)+m1,ibnd) + &
                                  dproj(ofstv(na1)+m1,ibnd) * &
                                  rproj(ofstv(na1)+m2,ibnd)
         ENDDO ! ibnd
         !
      ENDDO ! m1
      ENDDO ! m2
      !
   ENDDO ! na1
   !
10 DEALLOCATE ( dproj ) 
   !
   CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   IF ( nspin == 1 ) dns = 0.5D0 * dns
   !
   ! Impose hermiticity of dns_{m1,m2}
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      !
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = m1+1, ehub_lm(nt2)
         dns(m2,m1,current_spin,na1) = dns(m1,m2,current_spin,na1)
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na1
   !
   CALL stop_clock('dndtau_gamma_u')
   !
   RETURN
   !
END SUBROUTINE
!
!-------------------------------------------------------------------------
SUBROUTINE force_hub_v ( forceh )
   !----------------------------------------------------------------------
   ! This routine computes the Hubbard contribution to the force. It gives
   ! as output the product:
   ! $$ \frac{dE_\text{hub}}{dn_{ij}^\alpha}\cdot\frac{dn_{ij}^\alpha} 
   ! {du}(\alpha,\text{ipol}) \ ,$$
   ! which is the force acting on the atom at \(\text{tau_alpha}\)
   ! (in the unit cell) along the direction \(\text{ipol}\).
   ! Note: DFT+U+V force does not support OpenMP.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : at, bg
   USE basis,                ONLY : natomwfc, wfcatom, swfcatom
   USE symme,                ONLY : symvector
   USE io_files,             ONLY : prefix
   USE wvfct,                ONLY : nbnd, npwx
   USE control_flags,        ONLY : gamma_only, offload_type
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE scf,                  ONLY : v, rho
   USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, me_pool, &
                                    nproc_pool
   USE mp,                   ONLY : mp_sum
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, &
                                    deallocate_bec_type
   USE uspp,                 ONLY : nkb, vkb, ofsbeta, okvan
   USE uspp_param,           ONLY : nh
   USE wavefunctions,        ONLY : evc
   USE klist,                ONLY : nks, xk, ngk, igk_k
   USE io_files,             ONLY : nwordwfc, iunwfc
   USE buffers,              ONLY : get_buffer
   !by WIY
   USE mp_bands,             ONLY : use_bgrp_in_hpsi, intra_bgrp_comm
   USE noncollin_module,     ONLY : npol, noncolin
   USE force_mod,            ONLY : eigenval, eigenvect, overlap_inv
   USE ldaU,                 ONLY : Hubbard_projectors 
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, ehub_csize,     &
                                    ehub_nn_vect, ehub_nn_dim, ehub_nnnc, &
                                    wfcv, nwfcv, ofstv, copy_v_wfc, ehub_v
   !by WIY
   USE uspp_init,            ONLY : init_us_2
   !
   IMPLICIT NONE
   !
   REAL(DP)                 :: forceh(3,nat)
   !
   ! ... local variables
   !
   TYPE(bec_type)           :: proj     ! proj(nwfcV,nbnd)
   !COMPLEX(DP), ALLOCATABLE :: proj_nc     ! proj(nwfcV,nbnd)
   COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
   COMPLEX(DP), ALLOCATABLE :: dns(:,:,:,:,:,:)
   ! dns(ehub_lm,ehub_lm,nspin,nat1,nat2,ics)
   INTEGER                  :: ics, nnc
   INTEGER                  :: is, is_prime, ik, ijkb0, ipol, alpha, npw, npw_
   INTEGER                  :: m1, m2, na1, na2, nt1, nt2
   INTEGER                  :: n1, n2, n3
   INTEGER                  :: nb_s, nb_e, mykey
   LOGICAL                  :: save_flag
   !
   CALL start_clock ('force_hub_v')
   !
   save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
   !
   IF ( Hubbard_projectors /= "ortho-atomic" ) CALL errore("force_hub_v", &
                   " forces for this Hubbard_projectors type not implemented", 1 )
   !
   !IF (noncolin) CALL errore ("forceh","Noncollinear case is not supported",1)
   !
   ALLOCATE ( dns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc) )
   !
   ALLOCATE ( spsi(npwx*npol,nbnd) )
   ALLOCATE ( wfcatom(npwx*npol,natomwfc) )
   !! in case Hubbard_projectors == "ortho-atomic"
   ALLOCATE ( swfcatom(npwx*npol,natomwfc))
   ALLOCATE ( eigenval(natomwfc))
   ALLOCATE ( eigenvect(natomwfc, natomwfc))
   ALLOCATE ( overlap_inv(natomwfc,natomwfc))
   !
   IF (noncolin) THEN
    ALLOCATE ( proj%k(nwfcv,nbnd))
    !ALLOCATE ( proj_nc(nwfcv,nbnd))
    !proj_nc(:,:) = (0.0_DP, 0.0_DP)
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
   forceh(:,:) = 0.D0
   !
   ics = ehub_csize + 1
   !
   !
   ! we start a loop on k points
   !
   DO ik = 1, nks
      !
      IF (lsda) current_spin = isk(ik)
      npw = ngk(ik)
      npw_ = npwx * npol
      !
      IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE.)
      ! Compute Spsi = S * psi
      CALL allocate_bec_type(nkb, nbnd, becp)
      CALL calbec(offload_type, npw, vkb, evc, becp)
      CALL s_psi(npwx, npw, nbnd, evc, spsi)
      !CALL s_psi_acc(npwx, npw, nbnd, evc, spsi)
      CALL deallocate_bec_type(becp)
      !
      ! Set up various quantities, in particular wfcv which
      ! contains Hubbard-U ortho-atomic wavefunctions (w.o. ultrasoft S)
      CALL orthowfcv2 (ik)  
      !
      ! wfcv contains Hubbard-U atomic wavefunctions
      ! proj=<wfcv|S|evc> - no need to read S*wfcv from buffer
      !
      IF (noncolin) THEN
         !CALL calbec( offload_type, npw_, wfcv, spsi, proj_nc)
         CALL ZGEMM ('C', 'N', nwfcv, nbnd, npwx*npol, (1.0_DP, 0.0_DP), wfcv, &
                    npwx*npol, spsi, npwx*npol, (0.0_DP, 0.0_DP),  proj%k, nwfcv)
         CALL mp_sum( proj%k( :, 1:nbnd ), intra_bgrp_comm )
      ELSE
         CALL calbec( offload_type, npw, wfcv, spsi, proj)
      ENDIF
      !
      ! now we need the first derivative of proj with respect to tau(alpha,ipol)
      !
      DO alpha = 1, nat ! forces are calculated for atom alpha ...
         !
         ijkb0 = ofsbeta(alpha) ! positions of beta functions for atom alpha
         !
         DO ipol = 1, 3  ! forces are calculated for coordinate ipol ...
            !
            IF ( gamma_only ) THEN
               ! .. we don't support gamma calculation
               !
               !CALL dndtau_gamma_v(proj%r, spsi, alpha, ijkb0, ipol, ik, &
               !                    nb_s, nb_e, mykey, dns)
            ELSE
               ! 
               IF (noncolin) THEN
                  !CALL dndtau_k_v(proj_nc, spsi, alpha, ijkb0, ipol, ik, &
                  CALL dndtau_k_v(proj%k, spsi, alpha, ijkb0, ipol, ik, &
                                  nb_s, nb_e, mykey, dns)
               ELSE
                  CALL dndtau_k_v(proj%k, spsi, alpha, ijkb0, ipol, ik, &
                                  nb_s, nb_e, mykey, dns)
               ENDIF
               ! 
            ENDIF
            !
            DO is = 1, nspin
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
                     IF (noncolin) THEN
                        forceh(ipol,alpha) = forceh(ipol,alpha) &
                                           - DBLE(v%ne0_nc(m1,m2,is,na1,na2,1) &
                                           * dns(m1,m2,is,na1,na2,1))            
                     ELSE
                        forceh(ipol,alpha) = forceh(ipol,alpha) &
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
                     IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
                     !
                     DO m1 = 1, ehub_lm(nt1)
                     DO m2 = 1, ehub_lm(nt2)
                         !
                         IF (noncolin) THEN
                            forceh(ipol,alpha) = forceh(ipol,alpha) &
                                               - DBLE(v%ne0_nc(m1,m2,is,na1,na2,nnc) &
                                               * dns(m1,m2,is,na1,na2,nnc))
                         ELSE
                            forceh(ipol,alpha) = forceh(ipol,alpha) &
                                               - DBLE(v%ne0(m1,m2,is,na1,na2,nnc) &
                                               * dns(m1,m2,is,na1,na2,nnc))
                         ENDIF
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
         ENDDO ! ipol
         !
      ENDDO ! alpha
      !
   ENDDO ! ik
   !
   CALL mp_sum(forceh, inter_pool_comm)
   !
   CALL deallocate_bec_type(proj)
   !
   DEALLOCATE( dns )
   DEALLOCATE( wfcatom ) 
   DEALLOCATE( spsi )
   !
   DEALLOCATE( swfcatom )
   DEALLOCATE( eigenval )
   DEALLOCATE( eigenvect )
   DEALLOCATE( overlap_inv )
   !
   IF ( nspin == 1 ) forceh(:,:) = 2.D0 * forceh(:,:)
   !
   ! symmetrize
   !
   CALL symvector(nat, forceh)
   !
#if defined(__DEBUG)
   WRITE ( 66,'("Hubbard contribution Begin")' )
   WRITE ( 66,'(3f12.6)' ) forceh(:,:)
   WRITE ( 66,'("Hubbard contribution End")' )
#endif
   use_bgrp_in_hpsi = save_flag
   !
   CALL stop_clock('force_hub_v')
   !
   RETURN
   !
END SUBROUTINE force_hub_v
!
!
SUBROUTINE dndtau_k_v ( proj, spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dns )
   !-------------------------------------------------------------------------
   ! This routine computes the derivative of the nsg (generalized occupation
   ! matrix of the DFT+U+V scheme) with respect to the ionic
   ! displacement \(u(\text{alpha,ipol})\) used to obtain the generalized 
   ! Hubbard contribution to the atomic forces.
   !
   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : at
   USE klist,                ONLY : xk
   USE lsda_mod,             ONLY : nspin, current_spin
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE wavefunctions,        ONLY : evc
   USE uspp,                 ONLY : okvan
   USE force_mod,            ONLY : doverlap_inv
   USE basis,                ONLY : natomwfc
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, ehub_csize,     &
                                    ehub_nn_vect, ehub_nn_dim, ehub_nnnc, &
                                    wfcv, nwfcv, ofstv, ehub_l, counter_lm
   USE noncollin_module,     ONLY : npol, noncolin
   USE symm_base,            ONLY : time_reversal


   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)  :: proj(nwfcv,nbnd)
   !! projection 
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx*npol,nbnd)
   !! \(S|\ \text{evc}\rangle\)
   INTEGER, INTENT(IN)      :: alpha
   !! the displaced atom index
   INTEGER, INTENT(IN)      :: jkb0
   !! positions of beta functions for atom alpha
   INTEGER, INTENT(IN)      :: ipol
   !! the component of displacement
   INTEGER, INTENT(IN)      :: ik
   !! k-point index
   INTEGER, INTENT(IN)      :: nb_s, nb_e, mykey
   !! For band parallelization, starting and ending band number
   !! If each band appears more than once,
   !! compute ist contribution only once (i.e. when mykey=0)
   COMPLEX(DP), INTENT(OUT) :: dns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
   !! the derivative of the generalized atomic occupations
   !
   ! ... local variables
   !
   COMPLEX(DP), ALLOCATABLE :: dproj1(:,:,:), dproj_us(:,:)
   !
   COMPLEX(DP)              :: kphase
   !
   REAL(DP)                 :: arg
   !
   INTEGER                  :: ibnd, is1, is2
   INTEGER                  :: ics, nnc
   INTEGER                  :: m1, m2
   INTEGER                  :: na1, na2, nt1, nt2
   INTEGER                  :: n1, n2, n3
   INTEGER                  :: nnnc_count
   !
   CALL start_clock('dndtau_k_v')
   !
   ALLOCATE ( dproj1(nwfcv,nb_s:nb_e,nat) )
   IF (okvan) ALLOCATE ( dproj_us(nwfcv,nb_s:nb_e) )
   !
   ! Compute the derivative of the generalized occupation matrices 
   ! (the quantities dns(m1,m2)) of the atomic orbitals. 
   ! They are complex quantities as well as ns(m1,m2).
   !
   dns(:,:,:,:,:,:) = (0.0_DP, 0.0_DP)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   nnc = 1
   ics = ehub_csize + 1
   !
   !
   ! At this ik, need to compute the phase for each atom because the overlap
   ! matrix has the derivative of ortho-atomic orbital in na2 (dproj2).
   !! For the calculation of phase foctor, I will decide later to make new 
   !! subroutine or in the for-loop.
   ! 
   ! Compute the USPP contribution to dproj1:
   ! <\phi^{at}_{I,m1}|dS/du(alpha,ipol)|\psi_{k,v,s}>
   !
   IF (okvan) THEN
    CALL matrix_element_of_dSdtau_v (alpha, ipol, ik, jkb0, nwfcv, wfcv, nbnd, &
                        evc, dproj_us, nb_s, nb_e, mykey, .false.)
   ENDIF
   !
   ! In the 'ortho-atomic' case calculate d[(O^{-1/2})^T]
   ALLOCATE ( doverlap_inv(natomwfc,natomwfc) )
   CALL calc_doverlap_inv_v (alpha, ipol, ik, jkb0)
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) > 0 ) THEN
         CALL dprojdtau_k_v (spsi, alpha, na1, jkb0, ipol, ik, nb_s, nb_e, mykey, dproj1)
         IF (okvan) THEN
         dproj1(:,:,na1) = dproj1(:,:,na1) + dproj_us(:,:)
         ENDIF
      ENDIF !ehub_lm(nt1)
      !
   ENDDO !na1
   !
   ! Band parallelization. If each band appears more then once
   ! compute its contribution only once (i.e. when mykey=0)
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
      IF (mykey==0) THEN
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
             * CONJG(dproj1(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd,na1)) + &
                     dproj1(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd,na2) &
             * CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd)))
                  ENDDO
                  ENDDO
                  !
               ELSE
                  dns(m1,m2,current_spin,na1,na2,1) = &
                  dns(m1,m2,current_spin,na1,na2,1) + &
                   wg(ibnd,ik) *   DBLE( &
                  (proj(ofstv(na2)+m2,ibnd) * CONJG(dproj1(ofstv(na1)+m1,ibnd,na1)) + &
                   dproj1(ofstv(na2)+m2,ibnd,na2) * CONJG(proj(ofstv(na1)+m1,ibnd))))
               ENDIF
            ENDDO ! ibnd
            !
         ENDDO ! m2
         ENDDO ! m1
         !
      ENDIF ! mykey
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
         IF ( mykey == 0 ) THEN
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
                        dns(m1,m2,npol*(is1-1)+is2,na1,na2,nnc) = &
                        dns(m1,m2,npol*(is1-1)+is2,na1,na2,nnc) + &
                        wg(ibnd,ik) * dcmplx( kphase *( &
                        proj(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd) &
                * CONJG(dproj1(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd,na1)) + &
                        dproj1(ofstv(na2)+m2+counter_lm(m2,nt2)+(2*(ehub_l(m2,nt2)-1)+1)*(is2-1),ibnd,na2) &
                * CONJG(proj(ofstv(na1)+m1+counter_lm(m1,nt1)+(2*(ehub_l(m1,nt1)-1)+1)*(is1-1),ibnd))))
                     ENDDO
                     ENDDO
                  ELSE
                     dns(m1,m2,current_spin,na1,na2,nnc) = &
                     dns(m1,m2,current_spin,na1,na2,nnc) + &
                       wg(ibnd,ik) * DBLE( kphase * &
                    (proj(ofstv(na2)+m2,ibnd) * CONJG(dproj1(ofstv(na1)+m1,ibnd,na1)) + &
                     dproj1(ofstv(na2)+m2,ibnd,na2) * CONJG(proj(ofstv(na1)+m1,ibnd))))
                  ENDIF
                  !
               ENDDO ! ibnd
               !
            ENDDO ! m2
            ENDDO ! m1
            !
         ENDIF ! mykey
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
   DEALLOCATE ( dproj1 ) 
   IF (ALLOCATED(doverlap_inv)) DEALLOCATE ( doverlap_inv )
   IF (okvan)     DEALLOCATE ( dproj_us )
   !
   CALL mp_sum(dns, intra_pool_comm)   
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   !DO is2 = 1, ehub_nnnc
   !  DO na1 = 1, nat
   !  DO na2 = 1, nat
   !    nt1 = ityp(na1)
   !    nt2 = ityp(na2)
   !    DO m1 = 1, ehub_lm(nt1)
   !    DO m2 = 1, ehub_lm(nt2)
   !      IF(noncolin) THEN
   !        DO is1 = 1, nspin
   !         WRITE(11,'(6I4, f17.12)') m1, m2, is1, na1, na2, is2, ABS(dns(m1,m2,is1,na1,na2,is2))
   !        ENDDO
   !      ELSE
   !         WRITE(11,'(6I4, f17.12)') m1, m2, current_spin, na1, na2, is2, ABS(dns(m1,m2,current_spin,na1,na2,is2))
   !      ENDIF
   !    ENDDO
   !    ENDDO
   !  ENDDO
   !  ENDDO
   !ENDDO
   ! temp
   IF ( nspin == 1 ) dns = 0.5D0 * dns
   !
   ! Impose hermiticity of dns_{m1,m2}
   !
!   DO nnnc_count = 1, ehub_nnnc
!    DO na1 = 1, nat
!       nt1 = ityp(na1)
!       !
!       IF ( ehub_lm(nt1) == 0 ) CYCLE
!       !
!    DO na2 = na1, nat
!       nt2 = ityp(na2)
!       !
!       IF ( ehub_lm(nt2) == 0 ) CYCLE
!       !
!       IF ( na1 /= na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
!       !
!       DO m1 = 1, ehub_lm(nt1)
!       DO m2 = 1, ehub_lm(nt2)
!          dns(m1,m2,current_spin,na1,na2,nnnc_count) = (dns(m1,m2,current_spin,na1,na2,nnnc_count)   + &
!                                         CONJG(dns(m2,m1,current_spin,na2,na1,nnnc_count))) * 0.5D0
!          dns(m2,m1,current_spin,na2,na1,nnnc_count) = &
!                                         CONJG(dns(m1,m2,current_spin,na1,na2,nnnc_count))
!       ENDDO ! m2
!       ENDDO ! m1
!       !
!    ENDDO ! na2
!    ENDDO ! na1
!   ENDDO ! nnnc_count
   !
   CALL stop_clock('dndtau_k_v')
   !
   RETURN
   !
END SUBROUTINE dndtau_k_v
!
!
SUBROUTINE dndtau_gamma_v (rproj, spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dns)
   !--------------------------------------------------------------------------
   ! This routine computes the derivative of the nsg (generalized occupation
   ! matrix of the DFT+U+V scheme) with respect to the ionic
   ! displacement \(u(\text{alpha,ipol})\) used to obtain the generalized 
   ! Hubbard contribution to the atomic forces.
   !
   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : at
   USE klist,                ONLY : xk
   USE lsda_mod,             ONLY : nspin, current_spin
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE ldaV,                 ONLY : ehub_lm, ehub_lm_max, ehub_csize,     &
                                    ehub_nn_vect, ehub_nn_dim, ehub_nnnc, &
                                    nwfcv, ofstv

   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx,nbnd)

   REAL(DP), INTENT(IN)     :: rproj(nwfcv,nbnd)

   INTEGER, INTENT(IN)      :: alpha, jkb0
   INTEGER, INTENT(IN)      :: ipol, ik
   INTEGER, INTENT(IN)      :: nb_s, nb_e, mykey

   COMPLEX(DP), INTENT(OUT) :: dns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
   !
   ! ... local variables
   !
   REAL(DP), ALLOCATABLE    :: dproj(:,:)

   COMPLEX(DP)              :: kphase

   REAL(DP)                 :: arg

   INTEGER                  :: ibnd, is
   INTEGER                  :: ics, nnc
   INTEGER                  :: m1, m2
   INTEGER                  :: na1, na2, nt1, nt2
   INTEGER                  :: n1, n2, n3
   !
   CALL start_clock( 'dndtau_v' )
   !
   ALLOCATE ( dproj(nwfcv,nb_s:nb_e) )
   !
   ! Compute the derivative of the generalized occupation matrices
   ! (the quantities dnsg(m1,m2)) of the atomic orbitals.
   ! They are complex quantities as well as nsg(m1,m2).
   !
   CALL dprojdtau_gamma_v (spsi, alpha, jkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !
   dns(:,:,:,:,:,:) = (0.D0, 0.D0)
   !
   ! Band parallelization. If each band appears more than once
   ! compute its contribution only once (i.e. when mykey=0)
   !
   IF ( mykey /= 0 ) GO TO 10
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
      DO m1 = 1, ehub_lm(nt1)
      DO m2 = 1, ehub_lm(nt2)
         !
         DO ibnd = nb_s, nb_e
            dns(m1,m2,current_spin,na1,na2,1) = &
            dns(m1,m2,current_spin,na1,na2,1) + &
                                  wg(ibnd,ik) * &
                   (rproj(ofstv(na2)+m2,ibnd) * &
                    dproj(ofstv(na1)+m1,ibnd) + &
                    dproj(ofstv(na2)+m2,ibnd) * &
                    rproj(ofstv(na1)+m1,ibnd))
         ENDDO ! ibnd
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
               dns(m1,m2,current_spin,na1,na2,nnc) = &
               dns(m1,m2,current_spin,na1,na2,nnc) + &
                          DBLE(kphase)*wg(ibnd,ik) * &
                        (rproj(ofstv(na2)+m2,ibnd) * &
                         dproj(ofstv(na1)+m1,ibnd) + &
                         dproj(ofstv(na2)+m2,ibnd) * &
                         rproj(ofstv(na1)+m1,ibnd))
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
10 DEALLOCATE ( dproj )
   !
   CALL mp_sum(dns, intra_pool_comm)
   !
   ! In nspin=1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dnsg of one spin component
   !
   IF ( nspin == 1 ) dns = 0.5D0 * dns
   !
   ! Impose hermiticity of dns_{m1,m2}
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
         dns(m1,m2,current_spin,na1,na2,1) = (dns(m1,m2,current_spin,na1,na2,1)   + &
                                        CONJG(dns(m2,m1,current_spin,na2,na1,1))) * 0.5D0
         dns(m1,m2,current_spin,na1,na2,1) =                                        &
                                        CONJG(dns(m2,m1,current_spin,na2,na1,1))
      ENDDO ! m2
      ENDDO ! m1
      !
   ENDDO ! na2
   ENDDO ! na1
   !
   CALL stop_clock('dndtau_v')
   !
   RETURN
   !
END SUBROUTINE dndtau_gamma_v
!
SUBROUTINE dprojdtau_k_v(spsi, alpha, na, ijkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !-----------------------------------------------------------------------------
   ! This routine computes the first derivative of the projection
   ! \(\langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle\) with respect to 
   ! the atomic displacement \(u(\text{alpha,ipol})\). We remind that:
   ! $$ \text{ns}_{I,s,m1,m2} = \sum_{k,v}
   !    f_{kv} \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle
   !           \langle\psi_{k,v,s}|S|\phi^{at}_{I,m2}\rangle $$
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE klist,                ONLY : xk, ngk, igk_k
   USE wvfct,                ONLY : nbnd, npwx, wg
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, &
                                    deallocate_bec_type 
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE basis,                ONLY : natomwfc, wfcatom, swfcatom
   USE force_mod,            ONLY : eigenval, eigenvect, overlap_inv, doverlap_inv
   USE ldaV,                 ONLY : counter_lm, ehub_l, ehub_lm, wfcv, nwfcv, ofst, ofstv
   ! ofst: index of atomic wfcs used as projectors
   ! ofstv: For each atom, index of its projectors (among projectors only)
   USE noncollin_module,     ONLY : npol, noncolin
   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)   :: spsi(npwx*npol,nbnd)
   !! \(S\ |\text{evc}\rangle\)
   INTEGER, INTENT(IN)       :: alpha
   !! the displaced atom
   INTEGER, INTENT(IN)       :: na
   !! the atom for which the force is computed
   INTEGER, INTENT(IN)       :: ijkb0
   !! the postion of beta functions for atom alpha
   INTEGER, INTENT(IN)       :: ipol
   !! the component of displacement
   INTEGER, INTENT(IN)       :: ik
   !! k-point index
   INTEGER, INTENT(IN)       :: nb_s
   !! starting band number (for band parallelization)
   INTEGER, INTENT(IN)       :: nb_e
   !! ending band number (for band parallelization)
   INTEGER, INTENT(IN)       :: mykey
   !! If each band appears more than once, compute
   !! its contribution only once (i.e. when \(\text{mykey}=0)\)
   COMPLEX(DP), INTENT(OUT)  :: dproj(nwfcv,nb_s:nb_e,na)
   !! derivative of projection
   !
   ! ... local variables
   !
   REAL(DP)                  :: gvec, xki
   INTEGER                   :: nh_nt

   INTEGER                   :: npw, npw_, nt, ig, m1, m2, ibnd, m_start, m_end
   ! need to set local variable properly to calculate dproj

   COMPLEX (DP), ALLOCATABLE :: dproj0(:,:), dwfc(:,:)
   !! dproj0(:,:)   ! derivative of the projector    
   !! dwfc(:,:)     ! the derivative of the ortho-atomic wavefunction 
   !
   CALL start_clock('dprojdtau')
   !
   nt = ityp(na)
   npw = ngk(ik)
   npw_ = npwx*npol
   xki = xk(ipol,ik)
   nh_nt = nh(nt)
   !
   !$acc kernels
   dproj(:,:,na) = (0.d0, 0.d0)
   !$acc end kernels
   !
   !!! ORTHO-ATOMIC CASE !!!
   !
   ! Compute the derivative of the ortho-atomic wfc 'na' when displacing atom
   ! 'alpha'.
   ! Note, this derivative is different from zero not only when na = alpha but
   ! also when na /= alpha, i.e. when we displace a non-Hubbard atom this will
   ! give a non-zero contribution to the derivative of the ortho-atomic wfc na.
   ! This is so due to the definition of the ortho-atomic wfc:
   ! \phi_ortho_I = \sum_J O^{-1/2}_JI \phi_J
   ! Note: parallelization here is over plane waves, not over bands!
   !
   ALLOCATE( dwfc(npwx*npol,ehub_lm(nt)*npol) )
   !$acc data create(dwfc) present_or_copyin(wfcatom,overlap_inv)
   !
   !$acc kernels
   dwfc(:,:) = (0.d0,0.d0)
   !$acc end kernels
   !
   ! Determine how many atomic wavefunctions there are for atom 'alpha' and
   ! dertermine their position in the list of all atomic wavefunctions of all
   ! atoms
   CALL natomwfc_per_atom_v(alpha, m_start, m_end)
   !
   ! 1. Derivative of the atomic wavefunctions (the only one which is different
   ! from zero)) times O^-0.5 transposed.
   ! NOTE: overlap_inv is already transposed (it is O^{-1/2}_JI),
   !       hence we obtain \sum_J O^{-1/2}_JI \dphi_J/d\tau(alpha,ipol)
   DO ig = 1, npw
      !
      ! DFT+U: In the expression of dwfc we don't need (k+G) but just G; k always
      ! multiplies the underived quantity and gives an opposite contribution
      ! in c.c. term because the sign of the imaginary unit.
      ! DFT+U+V: the k-point coordinate is needed, i.e. (k+G) instead of just G
      !
      gvec = (g(ipol,igk_k(ig,ik)) + xk(ipol,ik)) * tpiba 
      DO m1 = 1, ehub_lm(nt)*npol
         DO m2 = m_start, m_end
         dwfc(ig,m1) = dwfc(ig,m1) + (0.D0,-1.D0) * gvec * &
                    overlap_inv(ofst(na)+m1,m2) *  wfcatom(ig,m2)
         IF (noncolin) THEN
         dwfc(ig+npwx,m1) = dwfc(ig+npwx,m1) + (0.D0,-1.D0) * gvec * &
                    overlap_inv(ofst(na)+m1,m2) *  wfcatom(ig+npwx,m2)
         ENDIF
         ENDDO !m2
      ENDDO !m1
      !
   ENDDO !ig
   !
   ! 2. Contribution due to the derivative of (O^{-1/2})_JI which is multiplied
   ! by atomic wavefunctions
   !
   ! Now compute \sum_J dO^{-1/2}_JI/d\tau(alpha,ipol) \phi_J
   ! and add it to another term (see above)
   ! Note, doverlap_inv is d(O^{-1/2}) not transposed. The transposition
   ! of d(O^{-1/2}) is taken into account via a proper usage of the order
   ! of indices in doverlap_inv:
   ! dwfc(ig,m1) = dwfc(ig,m1) + wfcatom(ig,m2) * doverlap_inv(m2,ofst(na)+m1)
   ! where m1=1, ehub_lm(nt); m2=1,natomwfc, ig=1,npw
   !
   !DO ig = 1, npw
   !   DO m1 = 1, ehub_lm(nt)
   !      DO m2 = 1, natomwfc
   !         dwfc(ig,m1) = dwfc(ig,m1) + &
   !             doverlap_inv(ofst(na)+m1,m2) * wfcatom(ig,m2)
   !      ENDDO !m2
   !   ENDDO !m1
   !ENDDO !ig
   !$acc host_data use_device(wfcatom,doverlap_inv,dwfc)
   CALL MYZGEMM( 'N','N', npwx*npol, ehub_lm(nt)*npol, natomwfc, (1.d0,0.d0), &
                 wfcatom, npwx*npol, doverlap_inv(:,ofst(na)+1:ofst(na)+ehub_lm(nt)*npol), &
                 natomwfc, (1.d0,0.d0), dwfc, npwx*npol )
   !$acc end host_data
   !   CALL ZGEMM('N','N', npw, ehub_lm(nt), natomwfc, (1.d0,0.d0), &
   !               wfcatom, npwx, doverlap_inv(:,ofst(na)+1:ofst(na)+ehub_lm(nt)), &
   !               natomwfc, (1.d0,0.d0), dwfc, npwx)
   !
   ! 3. Final step: compute dproj0 = <dwfc|spsi>
   !
   ALLOCATE ( dproj0(ehub_lm(nt)*npol,nbnd) )
   dproj0(:,:) = (0.0d0, 0.0d0)
   !$acc data create(dproj0)
   !
   !$acc host_data use_device(dwfc,spsi,dproj0)
   CALL MYZGEMM( 'C','N',ehub_lm(nt)*npol, nbnd, npwx*npol, (1.d0,0.d0), &
                 dwfc, npwx*npol, spsi, npwx*npol,  (0.d0,0.d0), &
                 dproj0, ehub_lm(nt)*npol )         
   !IF (noncolin) THEN
   !   CALL ZGEMM('C','N',ehub_lm(nt), nbnd, npw_, (1.D0,0.D0), &
   !           dwfc, npw_, spsi, npw_, (0.D0,0.D0), dproj0, &
   !           ehub_lm(nt)) 
   !ELSE
   !   CALL ZGEMM('C','N',ehub_lm(nt), nbnd, npw, (1.D0,0.D0), &
   !           dwfc, npwx, spsi, npwx, (0.D0,0.D0), dproj0, &
   !           ehub_lm(nt)) 
   !ENDIF
   !
   CALL mp_sum( dproj0, intra_bgrp_comm )
   !$acc end host_data
   !
   ! copy to dproj results for the bands treated by this processor
   !
   IF (mykey==0) THEN
      DO ibnd = nb_s, nb_e
         DO m1 = 1, ehub_lm(nt)*npol
            dproj(ofstv(na)+m1,ibnd,na) = dproj0(m1,ibnd)
            !IF (noncolin) dproj(ofstv(na)+m1+ehub_lm(nt),ibnd,na) = dproj0(m1+ehub_lm(nt),ibnd)
         ENDDO
      ENDDO 
   ENDIF
   !
   !$acc end data
   !$acc end data
   DEALLOCATE ( dproj0 ) 
   DEALLOCATE ( dwfc   )
   !
   CALL stop_clock('dprojdtau')
   !
   RETURN
   !
END SUBROUTINE dprojdtau_k_v
!
!------------------------------------------------------------------------
SUBROUTINE natomwfc_per_atom_v(alpha, m_start, m_end)
   !-----------------------------------------------------------------------
   !! This routine determines the starting (m_start) and the last (m_end)
   !! index for all atomic wavefunctions of a given atom 'alpha'
   !! when referring to the total list of all atomic wavefunctions of all atoms.
   !
   USE ions_base,    ONLY : nat, ityp
   USE uspp_param,   ONLY : upf
   USE io_global,    ONLY : stdout
   USE noncollin_module, ONLY : noncolin
   USE kinds,            ONLY : DP
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)  :: alpha
   INTEGER, INTENT(OUT) :: m_start
   INTEGER, INTENT(OUT) :: m_end
   !
   ! Local variables
   !
   INTEGER :: counter, l, na, nt, nb
   !
   counter = 0
   m_start = 0
   m_end   = 0
   DO na = 1, nat
      IF (na == alpha) m_start = counter + 1
      nt = ityp(na)
      DO nb = 1, upf(nt)%nwfc
         IF (upf(nt)%oc(nb) >= 0.d0) THEN
            l = upf(nt)%lchi(nb)
            !
            IF (noncolin) THEN
               IF (upf(nt)%has_so) THEN
                  ! j = l-1/2, degeneracy 2*l
                  counter = counter + 2*l
                  ! j = l+1/2, degeneracy 2*l+2
                  IF (ABS( upf(nt)%jchi(nb)-l-0.5_DP ) < 1.D-6) &
                  & counter = counter + 2
                  !
               ELSE
                  counter = counter + 2*( 2*l + 1 )
               ENDIF
            ELSE
                counter = counter + 2*l + 1
            ENDIF
         ENDIF
      ENDDO
      IF (na == alpha) THEN
         m_end = counter
         GO TO 11
      ENDIF
   ENDDO
   !
11 CONTINUE
   !
   IF (m_start==0 .OR. m_end==0) CALL errore("natomwfc_per_atom_v", &
                                   "m_start=0 or m_end=0",1)
   IF (m_start > m_end) CALL errore("natomwfc_per_atom_v", &
                                   "m_start > m_end",1)
   !
   RETURN
   !
END SUBROUTINE natomwfc_per_atom_v
! by WIY
! by WIY
!--------------------------------------------------------------------
SUBROUTINE calc_doverlap_inv_v (alpha, ipol, ik, ijkb0)
   !-----------------------------------------------------------------
   ! This routine computes the derivative of O^{-1/2} transposed 
   !
   USE kinds,          ONLY : DP
   USE cell_base,      ONLY : tpiba
   USE gvect,          ONLY : g
   USE uspp,           ONLY : okvan
   USE klist,          ONLY : xk, ngk, igk_k
   USE basis,          ONLY : natomwfc, wfcatom, swfcatom
   USE force_mod,      ONLY : eigenval, eigenvect, overlap_inv, doverlap_inv
   USE mp_bands,       ONLY : intra_bgrp_comm
   USE mp,             ONLY : mp_sum
   USE ldaU,           ONLY : Hubbard_projectors
   USE noncollin_module, ONLY : npol, noncolin
   USE wvfct,            ONLY : nbnd, npwx, wg
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: alpha
   !! the displaced atom
   INTEGER, INTENT(IN) :: ipol
   !! the component of displacement
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   INTEGER, INTENT(IN) :: ijkb0
   !! position of beta functions for atom alpha
   !
   INTEGER :: ig, m1, m2, npw, m_start, m_end
   REAL(DP) :: gvec, xki
   COMPLEX(DP) :: temp
   COMPLEX(DP), ALLOCATABLE :: doverlap(:,:), doverlap_us(:,:)
   !! derivative of the overlap matrix  
   !
   CALL start_clock( 'calc_doverlap_inv_v' )
   !
   IF (Hubbard_projectors.NE."ortho-atomic") RETURN
   !
   xki = xk(ipol,ik)
   !
   ALLOCATE (doverlap(natomwfc,natomwfc))
   !
   !$acc data present_or_copyin(wfcatom,swfcatom) 
   !
   !$acc kernels
   doverlap_inv(:,:) = (0.0d0,0.0d0)
   !$acc end kernels
   doverlap(:,:) = (0.0d0, 0.0d0)
   !
   npw = ngk(ik)
   !
   ! Determine how many atomic wafefunctions there are for atom 'alpha'
   ! and determine their position in the list of all atomic 
   ! wavefunctions of all atoms
   CALL natomwfc_per_atom_v(alpha, m_start, m_end)
   !
   ! Compute the derivative dO_IJ/d\tau(alpha,ipol)
   ! Calculate < dphi_I/d\tau(alpha,ipol) | S | phi_J >
   DO m1 = m_start, m_end
      DO m2 = 1, natomwfc
         temp = (0.d0,0.d0)
         DO ig = 1, npw
            ! (k+G) * 2pi/a
            gvec = (g(ipol,igk_k(ig,ik)) + xki) * tpiba
            temp = temp + (0.d0,1.d0) * gvec * CONJG(wfcatom(ig,m1))  *&
                          swfcatom(ig,m2)
            IF (noncolin) THEN
            temp = temp + (0.d0,1.d0) * gvec * CONJG(wfcatom(ig+npwx,m1)) *&
                          swfcatom(ig+npwx,m2)
            ENDIF
         ENDDO
         doverlap(m1,m2) = temp
      ENDDO
   ENDDO
   ! Calculate < phi_I | S | dphi_J/d\tau(alpha,ipol) >
   DO m1 = 1, natomwfc
      DO m2 = m_start, m_end
         temp = (0.d0,0.d0)
         !$acc parallel loop reduction(+:temp)
         DO ig = 1, npw
            ! (k+G) * 2pi/a
            gvec = (g(ipol,igk_k(ig,ik)) + xki) * tpiba
            temp = temp + (0.d0,-1.d0) * gvec * CONJG(swfcatom(ig,m1))  *&
                          wfcatom(ig,m2)
            IF (noncolin) THEN
            temp = temp + (0.d0,-1.d0) * gvec * CONJG(swfcatom(ig+npwx,m1)) *&
                          wfcatom(ig+npwx,m2)
            ENDIF
         ENDDO
         doverlap(m1,m2) = doverlap(m1,m2) + temp
      ENDDO
   ENDDO
   !
   ! Sum over G vectors
   CALL mp_sum( doverlap, intra_bgrp_comm )
   !
   ! Add the USPP term in dO_IJ/d\tau(alpha,ipol):
   ! < phi_I | dS/d\tau(alpha,ipol) | phi_J >
   !
   IF (okvan) THEN
      ! ... Calculate doverlap_us = < phi_I | dS/d\tau(alpha,ipol) | phi_J >
      ALLOCATE(doverlap_us(natomwfc,natomwfc))
      !$acc data create(doverlap_us)
      CALL matrix_element_of_dSdtau_v (alpha, ipol, ik, ijkb0, &
               natomwfc, wfcatom, natomwfc, wfcatom, doverlap_us, 1, natomwfc, 0, .TRUE.)
      !$acc kernels
      doverlap(:,:) = doverlap(:,:) + doverlap_us(:,:)
      !$acc end kernels
      !$acc end data 
      DEALLOCATE(doverlap_us)
   ENDIF
   !
   ! Now compute dO^{-1/2}_JI/d\tau(alpha,ipol) using dO_IJ/d\tau(alpha,ipol)
   ! Note the transposition!
   !
   CALL calculate_doverlap_inv (natomwfc, eigenval, eigenvect, &
                                doverlap, doverlap_inv)
   !
   !$acc end data
   !$acc end data
   DEALLOCATE (doverlap)
   !
   CALL stop_clock('calc_doverlap_inv_v')
   !
   RETURN
   !
END SUBROUTINE calc_doverlap_inv_v
! by WIY
! by WIY
SUBROUTINE matrix_element_of_dSdtau_v( alpha, ipol, ik, ijkb0, lA, A, &
                                       lB, B, A_dS_B, lB_s, lB_e, mykey, flag )
   !--------------------------------------------------------------------
   !! This routine computes the matrix element \(\langle A | 
   !! dS/d\tau(\alpha,\text{ipol}) | B \rangle\).  
   !! Written by I. Timrov (2020).
   !
   USE ldav,                 ONLY : ehub_lm, ofstv, ehub_l, counter_lm
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE wvfct,                ONLY : npwx, wg
   USE uspp,                 ONLY : nkb, vkb, qq_at, qq_so, okvan
   USE uspp_param,           ONLY : nh, upf
   USE klist,                ONLY : igk_k, ngk
   USE becmod,               ONLY : calbec
   USE gvect,                ONLY : g
   USE control_flags,        ONLY : offload_type
   USE noncollin_module,     ONLY : noncolin, npol, lspinorb
   USE mp,                   ONLY : mp_sum   
   USE mp_bands,             ONLY : intra_bgrp_comm
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: alpha
   !! the displaced atom
   INTEGER, INTENT(IN) :: ipol
   !! the component of displacement
   INTEGER, INTENT(IN) :: ik
   !! the k point
   INTEGER, INTENT(IN) :: ijkb0
   !! position of beta functions for atom alpha 
   INTEGER, INTENT(IN) :: lA
   INTEGER, INTENT(IN) :: lB
   !! There is a possibility to parallelize over lB
   INTEGER, INTENT(IN) :: lB_s
   !! lB start
   INTEGER, INTENT(IN) :: lB_e
   !! lB end
   LOGICAL, INTENT(IN) :: flag   ! noncollinear: controlling whether 
                                 ! calculating <phi|dS|PSI> 
                                 ! or          <phi|dS|PHI> (= .true.)
   COMPLEX(DP), INTENT(IN)  :: A(npwx*npol,lA)
   COMPLEX(DP), INTENT(IN)  :: B(npwx*npol,lB)
   COMPLEX(DP), INTENT(OUT) :: A_dS_B(lA,lB_s:lB_e)
   INTEGER, INTENT(IN) :: mykey
   !
   ! ... local variables
   !
   INTEGER :: npw, nt, ih, jh, ig, iA, iB, mU, mD, nt1, na, nh_nt, nt2, is1
   REAL(DP) :: gvec
   COMPLEX(DP), ALLOCATABLE :: Adbeta(:,:), Abeta(:,:), dbetaB(:,:), &
                               betaB(:,:), aux(:,:), qq(:,:)
   !
   A_dS_B(:,:) = (0.0d0, 0.0d0)
   !
   IF (.NOT.okvan) RETURN
   !
   !$acc data present_or_copyin(A,B) present_or_copyout(A_dS_B)
   !
   nt = ityp(alpha)
   npw = ngk(ik)
   nh_nt = nh(nt)
   !
   ALLOCATE( Adbeta(lA,npol*nh(nt)) )
   ALLOCATE( Abeta(lA,npol*nh(nt))  )
   ALLOCATE( dbetaB(npol*nh(nt),lB) )
   ALLOCATE( betaB(npol*nh(nt),lB)  )
   ALLOCATE( qq(npol*nh(nt),npol*nh(nt)) )
   !$acc data create(Adbeta,Abeta,dbetaB,betaB,qq)
   !
   !$acc parallel loop collapse(2) present(qq_at)
   DO jh = 1, nh_nt
      DO ih = 1, nh_nt
         IF (noncolin) THEN
            IF ( upf(nt)%has_so ) THEN     
               qq(ih,jh) =             qq_so(ih,jh,1,nt)
               qq(ih,jh+nh_nt) =       qq_so(ih,jh,2,nt)
               qq(ih+nh_nt,jh) =       qq_so(ih,jh,3,nt)
               qq(ih+nh_nt,jh+nh_nt) = qq_so(ih,jh,4,nt)
            ELSE
               qq(ih,jh) =             CMPLX(qq_at(ih,jh,alpha), 0.0d0, kind=DP)
               qq(ih,jh+nh_nt) =       (0.0,0.0)
               qq(ih+nh_nt,jh) =       (0.0,0.0)
               qq(ih+nh_nt,jh+nh_nt) = CMPLX(qq_at(ih,jh,alpha), 0.0d0, kind=DP)     
            ENDIF   
         ELSE
            qq(ih,jh) = CMPLX(qq_at(ih,jh,alpha), 0.0d0, kind=DP)
         ENDIF
      ENDDO
   ENDDO
   !
   ! ... aux is used as a workspace
   ALLOCATE( aux(npwx*npol,nh(nt)*npol) )
   !$acc data create(aux)
   !
   !$acc parallel loop collapse(2) 
   DO ih = 1, nh_nt*npol
      DO ig = 1, npwx*npol
         aux(ig,ih) = (0.0d0, 0.0d0)
      ENDDO
   ENDDO
   !
   ! ... Beta function
! !omp parallel do default(shared) private(ig,ih)
   !$acc parallel loop collapse(2) present(vkb(:,:))
   DO ih = 1, nh_nt
      DO ig = 1, npw
         aux(ig,ih) = vkb(ig,ijkb0+ih)
         IF (noncolin) aux(ig+npwx,ih+nh_nt) = vkb(ig,ijkb0+ih)
      ENDDO
   ENDDO
! !omp end parallel do
   !
   IF (noncolin) THEN
      ! Calculate Abeta = <A|beta>
      ! Abeta(:,1       : nh(nt))      = spin up
      ! Abeta(:,1+nh(nt): nh(nt)*npol) = spin down      
      !
      Abeta=(0.0,0.0)
      CALL ZGEMM ('C', 'N', lA, nh_nt*npol, npwx*npol, (1.0_DP, 0.0_DP), A, &
                 npwx*npol, aux, npwx*npol, (0.0_DP, 0.0_DP), Abeta, lA)
      CALL mp_sum( Abeta(:, 1:nh_nt*npol) , intra_bgrp_comm )
      !
      ! Calculate betaB = <beta|B>
      ! betaB(:,1       : nh(nt))      = spin up
      ! betaB(:,1+nh(nt): nh(nt)*npol) = spin down
      !      
      betaB=(0.0,0.0)
      CALL ZGEMM ('C', 'N', nh_nt*npol, lB, npwx*npol, (1.0_DP, 0.0_DP), aux, &
                 npwx*npol, B, npwx*npol, (0.0_DP, 0.0_DP), betaB, nh_nt*npol)
      CALL mp_sum( betaB(:, 1:lB) , intra_bgrp_comm )
   ELSE 
      ! ... Calculate Abeta = <A|beta>
      CALL calbec( offload_type, npw, A, aux, Abeta )
      ! ... Calculate betaB = <beta|B>
      CALL calbec( offload_type, npw, aux, B, betaB )
   ENDIF
   !
   ! ... Calculate the derivative of the beta function
! !omp parallel do default(shared) private(ig,ih)
   !
   !$acc parallel loop
   DO ih = 1, nh_nt
      DO ig = 1, npw
         gvec = g(ipol,igk_k(ig,ik)) * tpiba
         aux(ig,ih) = (0.d0,-1.d0) * aux(ig,ih) * gvec
         IF (noncolin) aux(ig+npwx,ih+nh(nt)) = &
                    (0.d0,-1.d0) * aux(ig+npwx,ih+nh(nt)) * gvec
      ENDDO
   ENDDO
! !omp end parallel do
   !
   IF (noncolin) THEN
      ! Calculate Adbeta = <A|dbeta>      
      ! (same as Abeta)
      !
      Adbeta=(0.0,0.0)
      CALL ZGEMM ('C', 'N', lA, nh(nt)*npol, npwx*npol, (1.0_DP, 0.0_DP), A, &
               npwx*npol, aux, npwx*npol, (0.0_DP, 0.0_DP), Adbeta, lA)
      CALL mp_sum( Adbeta(:, 1:nh(nt)*npol) , intra_bgrp_comm )
      !
      ! Calculate dbetaB = <dbeta|B>
      ! (same as betaB)
      !
      dbetaB=(0.0,0.0)
      CALL ZGEMM ('C', 'N', nh(nt)*npol, lB, npwx*npol, (1.0_DP, 0.0_DP), aux, &
               npwx*npol, B, npwx*npol, (0.0_DP, 0.0_DP), dbetaB, nh(nt)*npol)
      CALL mp_sum( dbetaB(:, 1:lB) , intra_bgrp_comm )
   ELSE
      ! ... Calculate Abeta = <A|beta>
      CALL calbec( offload_type, npw, A, aux, Adbeta )
      ! ... Calculate betaB = <beta|B>
      CALL calbec( offload_type, npw, aux, B, dbetaB )
      !
   ENDIF
   !$acc end data
   DEALLOCATE( aux )
   ALLOCATE( aux(nh(nt)*npol,lB) )
   !$acc data create(aux)
   !
   IF (noncolin) THEN    
      aux(:,:) = (0.0,0.0)
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
      ! ... Calculate \sum_jh qq_at(ih,jh) * dbetaB(jh)
      !$acc host_data use_device(qq,dbetaB,aux)
      CALL MYZGEMM( 'N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                  qq, nh(nt), dbetaB(1,lB_s),    nh(nt), (0.0d0,0.0d0), &
                  aux(1,lB_s), nh(nt) )
      !$acc end host_data
   ENDIF
   !
   !$acc kernels
   dbetaB(:,:) = aux(:,:)
   !$acc end kernels
   !
   IF (noncolin) THEN
      aux(:,:) = (0.0,0.0)
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
      ! ... Calculate \sum_jh qq_at(ih,jh) * betaB(jh)
      !$acc host_data use_device(qq,betaB,aux)
      CALL MYZGEMM( 'N', 'N', nh(nt), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                  qq, nh(nt), betaB(1,lB_s),     nh(nt), (0.0d0,0.0d0), &
                  aux(1,lB_s), nh(nt) )
      !$acc end host_data
   ENDIF
   !
   !$acc kernels
   betaB(:,:) = aux(:,:)
   !$acc end kernels
   !
   !$acc end data
   DEALLOCATE( aux )
   !
   ! ... A_dS_B(iA,iB) = \sum_ih [Adbeta(iA,ih) * betaB(ih,iB) +
   ! ...                          Abeta(iA,ih)  * dbetaB(ih,iB)] 
   ! ... Only A_dS_B(:,lB_s:lB_e) are calculated
   !
   IF ( mykey == 0 ) THEN
      IF (noncolin) THEN
         nt1 = nh(nt) + 1
         IF ( .NOT.flag ) THEN
            DO na = 1, nat
               nt2 = ityp(na)
               IF ( ehub_lm(nt2) > 0 ) THEN
                  mU = ofstv(na) + 1
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
                  !         Adbeta(mU,1), lA, &
                  !         betaB(1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                  !         A_dS_B(mU, lB_s), lA)
                  !CALL ZGEMM('N', 'N', ehub_lm(nt2), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                  !         Abeta(mU,1), lA, &
                  !         dbetaB(1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                  !         A_dS_B(mU, lB_s), lA)        
                  !! spin down
                  !CALL ZGEMM('N', 'N', ehub_lm(nt2), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                  !            Adbeta(mD, nt1), lA, &
                  !            betaB(nt1, lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                  !            A_dS_B(mD,lB_s), lA)
                  !CALL ZGEMM('N', 'N', ehub_lm(nt2), lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                  !            Abeta(mD, nt1), lA, &
                  !            dbetaB(nt1, lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                  !            A_dS_B(mD,lB_s), lA)   
               ENDIF    
            ENDDO 
         ELSEIF ( flag ) THEN
            CALL ZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                        Adbeta, lA, betaB(1,lB_s), nh(nt)*npol, (0.0d0,0.0d0), &
                        A_dS_B(1,lB_s), lA)
            CALL ZGEMM('N', 'N', lA, lB_e-lB_s+1, nh(nt)*npol, (1.0d0,0.0d0), &
                        Abeta, lA, dbetaB(1,lB_s), nh(nt)*npol, (1.0d0,0.0d0), &
                        A_dS_B(1,lB_s), lA)
         ENDIF        
      ELSE
         !$acc host_data use_device(Adbeta,betaB,Abeta,dbetaB,A_dS_B)
         CALL MYZGEMM( 'N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                     Adbeta, lA, betaB(1,lB_s), nh(nt), (0.0d0,0.0d0), &
                     A_dS_B(1,lB_s), lA )
         CALL MYZGEMM( 'N', 'N', lA, lB_e-lB_s+1, nh(nt), (1.0d0,0.0d0), &
                     Abeta, lA, dbetaB(1,lB_s), nh(nt), (1.0d0,0.0d0), &
                     A_dS_B(1,lB_s), lA )
         !$acc end host_data
      ENDIF
   ENDIF
   !
   !$acc end data
   DEALLOCATE( Abeta  )
   DEALLOCATE( Adbeta )
   DEALLOCATE( dbetaB )
   DEALLOCATE( betaB  )
   DEALLOCATE( qq     )
   !
   !$acc end data
   !
   RETURN
   !
END SUBROUTINE matrix_element_of_dSdtau_v
! by WIY
!
!
SUBROUTINE dprojdtau_gamma_v(spsi, alpha, ijkb0, ipol, ik, nb_s, nb_e, mykey, dproj)
   !-----------------------------------------------------------------------
   ! This routine is the gamma version of \(\texttt{dprojdtau_k}\).
   ! It computes the first derivative of the projection
   ! \(\langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle\) with respect to 
   ! the atomic displacement \(u(\text{alpha,ipol})\). We remind that:
   ! $$ \text{ns}_{I,s,m1,m2} = \sum_{k,v}
   !    f_{kv} \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle
   !           \langle\psi_{k,v,s}|S|\phi^{at}_{I,m2}\rangle $$
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk, ngk, igk_k
   USE wvfct,                ONLY : nbnd, npwx,  wg
   USE uspp,                 ONLY : nkb, vkb, qq_at
   USE uspp_param,           ONLY : nh
   USE wavefunctions,        ONLY : evc
   USE becmod,               ONLY : bec_type, becp, calbec
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE ldaV,                 ONLY : ehub_lm, wfcv, nwfcv, ofstv

   IMPLICIT NONE
   !
   ! I/O variables
   !
   COMPLEX(DP), INTENT(IN)  :: spsi(npwx,nbnd)

   INTEGER, INTENT(IN)      :: alpha, ijkb0
   INTEGER, INTENT(IN)      :: ipol, ik
   INTEGER, INTENT(IN)      :: nb_s, nb_e, mykey

   REAL(DP), INTENT(OUT)    :: dproj(nwfcv,nb_s:nb_e)
   !
   ! ... local variables
   !
   REAL(DP)                 :: gvec

   INTEGER                  :: nt
   INTEGER                  :: ih, jh, ig, npw
   INTEGER                  :: m, ibnd

   COMPLEX(DP), ALLOCATABLE :: dwfc(:,:), dbeta(:,:)

   REAL(DP), ALLOCATABLE    :: dproj0(:,:), betapsi(:,:), dbetapsi(:,:), &
                               wfatbeta(:,:), wfatdbeta(:,:), bproj(:,:)
   !      dwfc(npwx,ldim),     ! the derivative of the atomic wavefunction
   !      dbeta(npwx,nhm),     ! the derivative of the beta function
   !      betapsi(nhm,nbnd),   ! <beta|evc>
   !      dbetapsi(nhm,nbnd),  ! <dbeta|evc>
   !      wfatbeta(nwfcv,nhm), ! <wfcv|beta>
   !      wfatdbeta(nwfcv,nhm) ! <wfcv|dbeta>
   !
   CALL start_clock('dprojdtau')
   !
   nt = ityp(alpha)
   npw = ngk(ik)
   !
   dproj(:,:) = 0.0_dp
   !
   ! First the derivatives of the atomic wfc and the beta are computed
   ! Note: parallelization here is over plane waves, not over bands!
   !
   ALLOCATE ( dproj0(ehub_lm(nt),nbnd) )
   ALLOCATE ( dwfc(npwx,ehub_lm(nt)) )
   !
   dproj0(:,:) =  0.D0
   dwfc(:,:)   = (0.D0,0.D0)
   !
   DO ig = 1, npw
      !
      ! In the expression of dwfc we don't need (k+G) but just G; k always
      ! multiplies the underived quantity and gives an opposite contribution
      ! in c.c. term because the sign of the imaginary unit. But in any case,
      ! here we consider the situation when k = 0.
      !
      gvec = g(ipol,igk_k(ig,ik)) * tpiba
      !
      DO m = 1, ehub_lm(nt)
         dwfc(ig,m) = (0.D0,-1.D0) * gvec * wfcv(ig,ofstv(alpha)+m)
      ENDDO
      !
   ENDDO
   !
   ! there is no G=0 term
   CALL DGEMM('T', 'N', ehub_lm(nt), nbnd, 2*npw, 2.0_dp, &
              dwfc, 2*npwx, spsi, 2*npwx, 0.0_dp, dproj0, &
              ehub_lm(nt))
   !
   DEALLOCATE ( dwfc )
   !
   CALL mp_sum(dproj0, intra_bgrp_comm)
   !
   ! copy to dproj results for the bands treated by this processor
   !
   DO m = 1, ehub_lm(nt)
      dproj(ofstv(alpha)+m,:) = dproj0(m,nb_s:nb_e)
   ENDDO
   !
   DEALLOCATE ( dproj0 ) 
   !
   !
   ALLOCATE ( dbetapsi(nh(nt),nbnd) ) 
   ALLOCATE ( wfatdbeta(nwfcv,nh(nt)) )
   ALLOCATE ( wfatbeta(nwfcv,nh(nt)) )
   ALLOCATE(dbeta(npwx,nh(nt)))
   !
   DO ih = 1, nh(nt)
      DO ig = 1, npw
         dbeta(ig,ih) = vkb(ig,ijkb0+ih)
      ENDDO
   ENDDO
   !
   CALL calbec(npw, wfcv, dbeta, wfatbeta)
   !
   DO ih = 1, nh(nt)
      DO ig = 1, npw
         gvec = g(ipol,igk_k(ig,ik)) * tpiba
         dbeta(ig,ih) = (0.D0,-1.D0) * dbeta(ig,ih) * gvec
      ENDDO
   ENDDO
   !
   CALL calbec(npw, dbeta, evc, dbetapsi)
   CALL calbec(npw, wfcv, dbeta, wfatdbeta)
   !
   DEALLOCATE ( dbeta )
   !
   ! calculate \sum_j qq(i,j)*dbetapsi(j)
   ! betapsi is used here as work space
   !
   ALLOCATE ( betapsi(nh(nt), nbnd) )
   !
   betapsi(:,:) = (0.0_dp, 0.0_dp)
   !
   ! here starts band parallelization
   DO ih = 1, nh(nt)
      DO ibnd = nb_s, nb_e
         DO jh = 1, nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                             qq_at(ih,jh,alpha) * &
                              dbetapsi(jh,ibnd)
         ENDDO
      ENDDO
   ENDDO
   !
   dbetapsi(:,:) = betapsi(:,:)
   !
   ! calculate \sum_j qq(i,j)*betapsi(j)
   !
   betapsi(:,:) = (0.0_dp, 0.0_dp)
   !
   DO ih = 1, nh(nt)
      DO ibnd = nb_s, nb_e
         DO jh = 1, nh(nt)
            betapsi(ih,ibnd) = betapsi(ih,ibnd) + &
                             qq_at(ih,jh,alpha) * &
                          becp%r(ijkb0+jh,ibnd)
         ENDDO
      ENDDO
   ENDDO
   !
   ! dproj(iwf,ibnd) = \sum_ih wfatdbeta(iwf,ih)*betapsi(ih,ibnd) +
   !                           wfatbeta(iwf,ih)*dbetapsi(ih,ibnd) 
   !
   IF ( mykey == 0 .AND. nh(nt) > 0 ) THEN
      CALL DGEMM('N', 'N', nwfcv, nb_e-nb_s+1, nh(nt), 1.0_dp, &
           wfatdbeta, nwfcv, betapsi(1,nb_s), nh(nt), 1.0_dp,  &
           dproj(1,nb_s), nwfcv)
      CALL DGEMM('N', 'N', nwfcv, nb_e-nb_s+1, nh(nt), 1.0_dp, &
           wfatbeta, nwfcv, dbetapsi(1,nb_s), nh(nt), 1.0_dp,  &
           dproj(1,nb_s), nwfcv)
   ENDIF
   !
   ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
   !
   DEALLOCATE ( betapsi )
   DEALLOCATE ( wfatbeta )
   DEALLOCATE ( wfatdbeta )
   DEALLOCATE ( dbetapsi )
   !
   CALL stop_clock('dprojdtau')
   !
   RETURN
   !
END SUBROUTINE dprojdtau_gamma_v
