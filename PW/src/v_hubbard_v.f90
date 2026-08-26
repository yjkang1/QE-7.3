SUBROUTINE v_hubbard_v( ns, v_hub, ehubene ) 

  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp
  USE ldaV,          ONLY : ehub_lm_max, ehub_lm,   &
                            ehub_v, ehub_nn_vect,   &
                            ehub_nn_dim, ehub_nnnc, &
                            ehub_csize
  USE lsda_mod,      ONLY : nspin
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout

  IMPLICIT NONE

  REAL(DP), INTENT(IN)  :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
  REAL(DP), INTENT(OUT) :: v_hub(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
  REAL(DP), INTENT(OUT) :: ehubene

  INTEGER :: is
  INTEGER :: na1, na2, nt1, nt2, m1, m2
  INTEGER :: n1, n2, n3, nnc, ics

  v_hub = 0.0_DP
  ehubene = 0.0_DP
  !
  ics = ehub_csize + 1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( ehub_lm(nt1) == 0 ) CYCLE
     !
     DO is = 1, nspin
        !
        DO m1 = 1, ehub_lm(nt1)
           !
           ehubene = ehubene + 0.5_DP*ehub_v(na1,na1,m1,m1,1)*ns(m1,m1,is,na1,na1,1)
           v_hub(m1,m1,is,na1,na1,1) = v_hub(m1,m1,is,na1,na1,1) + 0.5_DP*ehub_v(na1,na1,m1,m1,1)
           !
           DO m2 = 1, ehub_lm(nt1)
              !
              ehubene = ehubene - &
                        0.5_DP*ehub_v(na1,na1,m1,m2,1)*ns(m1,m2,is,na1,na1,1)*ns(m2,m1,is,na1,na1,1)
              v_hub(m1,m2,is,na1,na1,1) = v_hub(m1,m2,is,na1,na1,1) - ehub_v(na1,na1,m1,m2,1)*ns(m2,m1,is,na1,na1,1)
              !
           ENDDO ! m2
           !
        ENDDO ! m1
        !
     ENDDO ! is
     !
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
     DO is = 1, nspin
        !
        DO m1 = 1, ehub_lm(nt1)
        DO m2 = 1, ehub_lm(nt2)
           !
           ehubene = ehubene - &
                     0.5_DP*ehub_v(na1,na2,m1,m2,1)*ns(m1,m2,is,na1,na2,1)*ns(m2,m1,is,na2,na1,1)
           v_hub(m1,m2,is,na1,na2,1) = v_hub(m1,m2,is,na1,na2,1) - ehub_v(na1,na2,m1,m2,1)*ns(m2,m1,is,na2,na1,1)
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! is
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
     IF( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
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
        DO is = 1, nspin
           !
           DO m1 = 1, ehub_lm(nt1)
           DO m2 = 1, ehub_lm(nt2)
              ! NOTE: ns* = ns
              ehubene = ehubene - &
                        0.5_DP*ehub_v(na1,na2,m1,m2,nnc)*ns(m1,m2,is,na1,na2,nnc)*ns(m1,m2,is,na1,na2,nnc)
              v_hub(m1,m2,is,na1,na2,nnc) = v_hub(m1,m2,is,na1,na2,nnc) - ehub_v(na1,na2,m1,m2,nnc)*ns(m1,m2,is,na1,na2,nnc)
              !
           ENDDO ! m2
           ENDDO ! m1
           !
        ENDDO ! is
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
  IF ( nspin == 1 ) ehubene = 2.0_DP * ehubene
  !
  IF ( iverbosity > 0 ) THEN
     WRITE ( stdout,             * ) ''
     WRITE ( stdout,             * ) '    ========================='
     WRITE ( stdout, '(a,f12.9,a)' ) '     E(+U+V) = ', ehubene, ' Ry'
     WRITE ( stdout,             * ) '    ========================='
  ENDIF
  !
  RETURN

END SUBROUTINE v_hubbard_v
!
SUBROUTINE v_hubbard_a0( ns, v_hub, ehubene )

  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp
  USE ldaV,          ONLY : ehub_lm_max, ehub_lm, ehub_u
  USE lsda_mod,      ONLY : nspin
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout

  IMPLICIT NONE

  REAL(DP), INTENT(IN)  :: ns(ehub_lm_max,ehub_lm_max,nspin,nat)
  REAL(DP), INTENT(OUT) :: v_hub(ehub_lm_max,ehub_lm_max,nspin,nat)
  REAL(DP), INTENT(OUT) :: ehubene

  INTEGER :: is
  INTEGER :: na, nt, m1, m2

  v_hub = 0.0_DP
  ehubene = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( ehub_lm(nt) == 0 ) CYCLE
     !
     DO is = 1, nspin
        !
        DO m1 = 1, ehub_lm(nt)
           !
           ehubene = ehubene + 0.5_DP*ehub_u(na,m1,m1)*ns(m1,m1,is,na)
           v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + 0.5_DP*ehub_u(na,m1,m1)
           !
           DO m2 = 1, ehub_lm(nt)
              !
              ehubene = ehubene - &
                        0.5_DP*ehub_u(na,m1,m2)*ns(m1,m2,is,na)*ns(m2,m1,is,na)
              v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - ehub_u(na,m1,m2)*ns(m2,m1,is,na)
              !
           ENDDO ! m2
           !
        ENDDO ! m1
        !
     ENDDO ! is
     !
  ENDDO ! na
  !
  IF ( nspin == 1 ) ehubene = 2.0_DP * ehubene
  !
  IF ( iverbosity > 0 ) THEN
     WRITE ( stdout,             * ) ''
     WRITE ( stdout,             * ) '    ======================='
     WRITE ( stdout, '(a,f12.9,a)' ) '     E(+U) = ', ehubene, ' Ry'
     WRITE ( stdout,             * ) '    ======================='
  ENDIF
  !
  RETURN

END SUBROUTINE v_hubbard_a0
!
! WIY
SUBROUTINE v_hubbard_v_nc( ns, v_hub, ehubene ) 
   !
   !! Computes the extended Hubbard potential and energy (noncollinear case)
   !
   USE kinds,              ONLY : DP
   USE ions_base,          ONLY : nat, ityp
   USE noncollin_module,   ONLY : noncolin
   USE ldaV,               ONLY : ehub_lm_max, ehub_lm,   &
                                  ehub_v, ehub_nn_vect,   &
                                  ehub_nn_dim, ehub_nnnc, &
                                  ehub_csize
   USE lsda_mod,           ONLY : nspin
   USE control_flags,      ONLY : iverbosity
   USE io_global,          ONLY : stdout
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN)  :: ns(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
   !! Occupation matrix
   COMPLEX(DP), INTENT(OUT) :: v_hub(ehub_lm_max,ehub_lm_max,nspin,nat,nat,ehub_nnnc)
   !! Hubbard potential
   REAL(DP), INTENT(OUT) :: ehubene
   !! Hubbard energy
   !
   ! ... local variables
   !
   INTEGER :: is, is_prime
   INTEGER :: na1, na2, nt1, nt2, m1, m2
   INTEGER :: n1, n2, n3, nnc, ics
   !REAL(DP) :: mx, my, mz, mag2
   !COMPLEX(DP) :: n_tot
   !
   v_hub = 0.0_DP
   ehubene = 0.0_DP
   !
   ics = ehub_csize + 1
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      !
      IF ( ehub_lm(nt1) == 0 ) CYCLE
      ! u_matrix in plus_u_full.f90 means <m,m''|V_ee|m',m'''> = \sum(ak*Fk)
      ! In simplied formular, u_matrix can be replaced with the averaged Hubbard parameters.
      !
      DO is = 1, nspin ! is = 1 -> upup, 2 -> updown, 3 -> downup, 4 -> downdown
         ! Following from Physical Review B 102, 155117 (2020)
         !
         IF (is == 2) THEN 
            is_prime = 3
         ELSEIF (is == 3) THEN 
            is_prime = 2
         ELSE
            is_prime = is
         ENDIF
         !
         DO m1 = 1, ehub_lm(nt1)
            !
            IF ( is_prime == is ) THEN
               ehubene = ehubene + 0.5_DP*ehub_v(na1,na1,m1,m1,1)*DBLE(ns(m1,m1,is,na1,na1,1))
               v_hub(m1,m1,is,na1,na1,1) = v_hub(m1,m1,is,na1,na1,1) + 0.5_DP*ehub_v(na1,na1,m1,m1,1)
            ENDIF
            !
            DO m2 = 1, ehub_lm(nt1)
               !
               ehubene = ehubene - 0.5_DP*ehub_v(na1,na1,m1,m2,1) &
                                    *DBLE(ns(m1,m2,is,na1,na1,1)*ns(m2,m1,is_prime,na1,na1,1))
               v_hub(m1,m2,is,na1,na1,1) = v_hub(m1,m2,is,na1,na1,1) &
                                          & - ehub_v(na1,na1,m1,m2,1)*ns(m2,m1,is_prime,na1,na1,1)
               !
            ENDDO ! m2
            !
         ENDDO ! m1
         !
      ENDDO ! is
      !
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
      DO is = 1, nspin
         !
         IF (is == 2) THEN 
            is_prime = 3
         ELSEIF (is == 3) THEN 
            is_prime = 2
         ELSE
            is_prime = is
         ENDIF
         !
         DO m1 = 1, ehub_lm(nt1)
         DO m2 = 1, ehub_lm(nt2)
            !
            ehubene = ehubene - 0.5_DP*ehub_v(na1,na2,m1,m2,1) &
                              & *DBLE(ns(m1,m2,is,na1,na2,1)*ns(m2,m1,is_prime,na2,na1,1))
            v_hub(m1,m2,is,na1,na2,1) = v_hub(m1,m2,is,na1,na2,1) &
                                       & - ehub_v(na1,na2,m1,m2,1)*ns(m2,m1,is_prime,na2,na1,1)
            !
         ENDDO ! m2
         ENDDO ! m1
         !
      ENDDO ! is
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
      IF( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
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
         DO is = 1, nspin
            !
            DO m1 = 1, ehub_lm(nt1)
            DO m2 = 1, ehub_lm(nt2)
               ehubene = ehubene - 0.5_DP*ehub_v(na1,na2,m1,m2,nnc) & 
                                 & *DBLE(ns(m1,m2,is,na1,na2,nnc)*CONJG(ns(m1,m2,is,na1,na2,nnc)))
               v_hub(m1,m2,is,na1,na2,nnc) = v_hub(m1,m2,is,na1,na2,nnc) &
                                            & - ehub_v(na1,na2,m1,m2,nnc)*CONJG(ns(m1,m2,is,na1,na2,nnc))
               !
            ENDDO ! m2
            ENDDO ! m1
            !
         ENDDO ! is
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
   IF ( iverbosity > 0 ) THEN
      WRITE ( stdout,             * ) ''
      WRITE ( stdout,             * ) '    ========================='
      WRITE ( stdout, '(a,f12.9,a)' ) '     E(+U+V) = ', ehubene, ' Ry'
      WRITE ( stdout,             * ) '    ========================='
   ENDIF
   !
   RETURN
 
 END SUBROUTINE v_hubbard_v_nc
 !WIY
