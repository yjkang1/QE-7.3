SUBROUTINE offset_atom_wfc_v( hubbard_only, offset, counter )
  !
  USE kinds,            ONLY : DP
  USE uspp_param,       ONLY : upf
  USE ions_base,        ONLY : nat, ityp
  USE ldaV,             ONLY : ehub_l_choice
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN)  :: hubbard_only 
  INTEGER, INTENT(OUT) :: offset(nat), counter
  !
  INTEGER  :: na, nt, nw, l
  LOGICAL  :: hubbard_wfc
  CHARACTER(LEN=2) :: s
  !
  offset(:) = -99
  counter   =   0
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     ! WIY
     WRITE(s,'(i2)') nt
     IF ( ( ANY(ehub_l_choice(nt,:) >= 1) ) .AND. upf(nt)%nwfc < 1 ) THEN
        CALL errore('offset_atom_wfc_v', 'no atomic wavefunctions in &
               &pseudopotential file for species #' // s // new_line('a') // &
               &'use a pseudopotential file with atomic wavefunctions!', 1)
     ENDIF
     ! WIY
     !
     DO nw = 1, upf(nt)%nwfc
        !
        IF ( upf(nt)%oc(nw) >= 0.D0 ) THEN
           !
           l = upf(nt)%lchi(nw)
           !
           ! WIY
           IF ( noncolin ) THEN
              !
              hubbard_wfc = .FALSE.
              !
              IF ( ehub_l_choice(nt,l+1) == nw .OR. &
                 ( ehub_l_choice(nt,l+1) == nw+1 .AND. l /= 0)) hubbard_wfc = .TRUE.
              !
              IF (upf(nt)%has_so ) THEN
                 !
                 ! offset to be set at the first occurrence of required l
                 IF ( hubbard_wfc .AND. offset(na) == -99 ) &
                  & offset(na) = counter
                 !
                 IF ( .NOT. hubbard_only .OR. hubbard_wfc ) THEN
                    ! j = l-1/2, degeneracy 2l
                    counter = counter + 2*l
                    ! j = l+1/2, degeneracy 2*l+2
                    IF (ABS( upf(nt)%jchi(nw)-l-0.5_DP ) < 1.D-6) &
                    &  counter = counter + 2
                    !
                 ENDIF
                 !
              ELSE
                 !
                 IF ( ehub_l_choice(nt,l+1) == nw .AND. offset(na) == -99 ) &
                  & offset(na) = counter
                 !
                 IF ( .NOT. hubbard_only .OR. ehub_l_choice(nt,l+1) == nw ) &
                  & counter = counter + 2*( 2*l + 1 )
                 !
              ENDIF
              !
           ELSE
              !
              IF ( ehub_l_choice(nt,nw) == 1 .AND. offset(na) == -99 ) &
               & offset(na) = counter
              !
              IF ( .NOT. hubbard_only .OR. ehub_l_choice(nt,nw) == 1 ) &
               & counter = counter + 2*l + 1
              !
           ENDIF
           ! WIY
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN

END SUBROUTINE offset_atom_wfc_v
