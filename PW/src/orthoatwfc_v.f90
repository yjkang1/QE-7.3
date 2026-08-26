! Refer to SUBROUTINE orthoUwfc in orhoatwfc.f90
! by WIY 
SUBROUTINE orthowfcv

  USE kinds,            ONLY : DP
  USE buffers,          ONLY : get_buffer, save_buffer
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunhubv, nwordwfcv
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc, swfcatom
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE ldaU,             ONLY : Hubbard_projectors
  USE ldaV,             ONLY : wfcv, nwfcv, copy_v_wfc
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type,   &
                               deallocate_bec_type, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
  USE uspp_init,        ONLY : init_us_2

  IMPLICIT NONE

  INTEGER :: ik, npw

  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag

  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)

  ! by WIY
  ! Since eACBN0 code uses ortho-atomic wfcs only, we impose Hubbard_projectors to 'ortho-atomic'
  Hubbard_projectors = 'ortho-atomic'
  ! by WIY

  IF ( Hubbard_projectors=="ortho-atomic" ) THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.
     WRITE ( stdout,'(/5x,a,/)') 'Atomic wfcs used for DFT+U+V projectors are orthogonalized'
     IF ( gamma_only ) CALL errore('orthoatwfcv', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE ( stdout,'(/5x,"Hubbard_projectors = ",a)') Hubbard_projectors
     CALL errore ("orthoatwfcv"," This type of Hubbard projectors is not valid",1)
  ENDIF
  !
  ALLOCATE ( wfcatom(npwx*npol, natomwfc), swfcatom(npwx*npol, natomwfc) )
  !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
  !
  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type(nkb, natomwfc, becp)
  !
  DO ik = 1, nks
     !
     IF ( noncolin ) THEN
       CALL atomic_wfc_v_nc(ik, wfcatom)
       !CALL atomic_wfc_nc_updown_v(ik, wfcatom)
     ELSE
       CALL atomic_wfc(ik, wfcatom)
     ENDIF
     !
     npw = ngk(ik)
     !
     CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
     CALL calbec(npw, vkb, wfcatom, becp)
     CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     IF ( orthogonalize_wfc ) &
        CALL ortho_swfc(npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE.)
     !
     ! copy S * atomic wavefunctions with Hubbard term only in wfcv
     ! (this is used during the self-consistent solution of Kohn-Sham equations)
     ! save to unit iunhubv
     ! 
     CALL copy_v_wfc(swfcatom, noncolin)
     !
     IF ( nks > 1 ) &
        CALL save_buffer(wfcv, nwordwfcv, iunhubv, ik)
     !
  ENDDO
  !
  DEALLOCATE ( wfcatom, swfcatom )
  !
  CALL deallocate_bec_type(becp)
  !
  use_bgrp_in_hpsi = save_flag
  !
  RETURN

END SUBROUTINE orthowfcv
!
!-----------------------------------------------------------------------
SUBROUTINE orthowfcv2 (ik)
  !
  ! For a given k point "ik", this routine computes (ortho-)atomic wavefunctions 
  ! having an associated Hubbard U term * S, for DFT+U(+V) calculations. 
  ! Also without S (this is then used to computed Hubbard forces 
  ! and stresses). 
  ! wfcatom and swfcatom must be allocated on input.
  ! Beta functions vkb must be already computed before.
  ! 
  ! please refer to orthoUwfc_k in orthoatwfc.f90
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunhubv, nwordwfcv
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc, wfcatom, swfcatom
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE ldaU,             ONLY : Hubbard_projectors
  USE ldaV,             ONLY : wfcv, nwfcv, copy_v_wfc
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik ! the k point under consideration
  !
  INTEGER :: npw
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)

  IF ( Hubbard_projectors=="ortho-atomic" ) THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.
     IF ( gamma_only ) CALL errore('orthoatwfcv2', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE ( stdout,'(/5x,"Hubbard_projectors = ",a)') Hubbard_projectors
     CALL errore ("orthoatwfcv2"," This type of Hubbard projectors is not valid",1)
  ENDIF
  !
  ! Compute atomic wfc at this k (phi)
  IF (noncolin) THEN
     CALL atomic_wfc_v_nc (ik, wfcatom)
  ELSE
     CALL atomic_wfc (ik, wfcatom)
  ENDIF
  !
  IF (Hubbard_projectors=="ortho-atomic") THEN
     ALLOCATE(aux(npwx*npol,natomwfc))
     ! Copy atomic wfcs (phi)
     aux(:,:) = wfcatom(:,:)
  ENDIF
  !
  IF (orthogonalize_wfc) THEN
     !
     ! Number of plane waves at this k point
     npw = ngk(ik)
     !
     ! Allocate the array becp = <beta|wfcatom>
     CALL allocate_bec_type (nkb,natomwfc, becp)
     CALL calbec (npw, vkb, wfcatom, becp)
     ! Calculate swfcatom = S * phi
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     CALL deallocate_bec_type (becp)
     !  
     ! Compute the overlap matrix
     ! On the output: wfcatom = O^{-1/2} \phi (no ultrasoft S)
     CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .TRUE. )
     !
  ENDIF
  !
  ! Copy (ortho-)atomic wavefunctions with Hubbard U term only 
  ! in wfcU (no ultrasoft S).
  CALL copy_v_wfc (wfcatom, noncolin)
  !
  IF (Hubbard_projectors=="ortho-atomic") THEN
     ! Copy atomic wfcs
     wfcatom(:,:) = aux(:,:)
     DEALLOCATE(aux)
  ENDIF
  !
  RETURN
  !   

END SUBROUTINE orthowfcv2
!
! by WIY
