      MODULE hp_acbn0_eri_mod

      USE ISO_C_BINDING, ONLY : C_INT, C_DOUBLE

      INTERFACE
         REAL(C_DOUBLE) FUNCTION hp_acbn0_coulomb_repulsion               &
                                (xa,ya,za,norma,la,ma,na,alphaa, &
                                 xb,yb,zb,normb,lb,mb,nb,alphab, &
                                 xc,yc,zc,normc,lc,mc,nc,alphac, &
                                 xd,yd,zd,normd,ld,md,nd,alphad) &
                                 BIND(C, name="hp_acbn0_coulomb_repulsion")

            IMPORT

            IMPLICIT NONE

            REAL(C_DOUBLE), VALUE :: xa, ya, za, norma, alphaa
            REAL(C_DOUBLE), VALUE :: xb, yb, zb, normb, alphab
            REAL(C_DOUBLE), VALUE :: xc, yc, zc, normc, alphac
            REAL(C_DOUBLE), VALUE :: xd, yd, zd, normd, alphad
            INTEGER(C_INT), VALUE :: la, ma, na
            INTEGER(C_INT), VALUE :: lb, mb, nb
            INTEGER(C_INT), VALUE :: lc, mc, nc
            INTEGER(C_INT), VALUE :: ld, md, nd

         END FUNCTION hp_acbn0_coulomb_repulsion
      END INTERFACE

      END MODULE hp_acbn0_eri_mod

      FUNCTION hp_acbn0_compute_eri(lm1, lm2, lm3, lm4, x1, x2, x3, x4, c1, c2, c3, c4)

      USE hp_acbn0_eri_mod

      IMPLICIT NONE

      INTEGER,  PARAMETER  :: DP = selected_real_kind(14,200)

      REAL(DP), INTENT(IN) :: x1(3), x2(3), x3(3), x4(3)
      REAL(DP), INTENT(IN) :: c1(6), c2(6), c3(6), c4(6)

      INTEGER,  INTENT(IN) :: lm1, lm2, lm3, lm4

      REAL(DP)             :: hp_acbn0_compute_eri

      REAL(DP)             :: norm1, norm2, norm3, norm4
      REAL(DP)             :: norm1_tmp, norm2_tmp, &
                              norm3_tmp, norm4_tmp

      INTEGER              :: l1(3), l2(3), l3(3), l4(3)
      INTEGER              :: m1(3), m2(3), m3(3), m4(3)
      INTEGER              :: n1(3), n2(3), n3(3), n4(3)

      INTEGER              :: ia1, ia2, ia3, ia4
      INTEGER              :: ilmn1, ilmn2, ilmn3, ilmn4

      hp_acbn0_compute_eri = 0.0_DP
      CALL hp_acbn0_expnorm(lm1, l1, m1, n1, norm1_tmp)
      CALL hp_acbn0_expnorm(lm2, l2, m2, n2, norm2_tmp)
      CALL hp_acbn0_expnorm(lm3, l3, m3, n3, norm3_tmp)
      CALL hp_acbn0_expnorm(lm4, l4, m4, n4, norm4_tmp)

      DO ia1 = 1, 3
      DO ilmn1 = 1, 3

         ! normalization constant of phi_1
         IF (lm1 /= 5 .AND. lm1 /= 8) THEN
            IF (ilmn1 > 1) CYCLE
            norm1 = norm1_tmp * c1(ia1)
         ELSEIF (lm1 == 5) THEN
            IF (ilmn1 == 1) THEN
               norm1 = norm1_tmp * c1(ia1) *  2.0_DP
            ELSE
               norm1 = norm1_tmp * c1(ia1) * -1.0_DP
            ENDIF
         ELSEIF (lm1 == 8) THEN
            IF (ilmn1 == 1) THEN
               norm1 = norm1_tmp * c1(ia1) *  1.0_DP
            ELSEIF (ilmn1 == 2) THEN
               norm1 = norm1_tmp * c1(ia1) * -1.0_DP
            ELSE
               CYCLE
            ENDIF
         ENDIF

         DO ia2 = 1, 3
         DO ilmn2 = 1, 3

            ! normalization constant of phi_2
            IF (lm2 /= 5 .AND. lm2 /= 8) THEN
               IF (ilmn2 > 1) CYCLE
               norm2 = norm2_tmp * c2(ia2)
            ELSEIF (lm2 == 5) THEN
               IF (ilmn2 == 1) THEN
                  norm2 = norm2_tmp * c2(ia2) *  2.0_DP
               ELSE
                  norm2 = norm2_tmp * c2(ia2) * -1.0_DP
               ENDIF
            ELSEIF (lm2 == 8) THEN
               IF (ilmn2 == 1) THEN
                  norm2 = norm2_tmp * c2(ia2) *  1.0_DP
               ELSEIF (ilmn2 == 2) THEN
                  norm2 = norm2_tmp * c2(ia2) * -1.0_DP
               ELSE
                  CYCLE
               ENDIF
            ENDIF

            ! normalization constant of phi_3
            DO ia3 = 1, 3
            DO ilmn3 = 1, 3

               IF (lm3 /= 5 .AND. lm3 /= 8) THEN
                  IF (ilmn3 > 1) CYCLE
                  norm3 = norm3_tmp * c3(ia3)
               ELSEIF (lm3 == 5) THEN
                  IF (ilmn3 == 1) THEN
                     norm3 = norm3_tmp * c3(ia3) *  2.0_DP
                  ELSE
                     norm3 = norm3_tmp * c3(ia3) * -1.0_DP
                  ENDIF
               ELSEIF (lm3 == 8) THEN
                  IF (ilmn3 == 1) THEN
                     norm3 = norm3_tmp * c3(ia3) *  1.0_DP
                  ELSEIF (ilmn3 == 2) THEN
                     norm3 = norm3_tmp * c3(ia3) * -1.0_DP
                  ELSE
                     CYCLE
                  ENDIF
               ENDIF

               ! normalization constant of phi_4
               DO ia4 = 1, 3
               DO ilmn4 = 1, 3

                  IF (lm4 /= 5 .AND. lm4 /= 8) THEN
                     IF (ilmn4 > 1) CYCLE
                     norm4 = norm4_tmp * c4(ia4)
                  ELSEIF (lm4 == 5) THEN
                     IF (ilmn4 == 1) THEN
                        norm4 = norm4_tmp * c4(ia4) *  2.0_DP
                     ELSE
                        norm4 = norm4_tmp * c4(ia4) * -1.0_DP
                     ENDIF
                  ELSEIF (lm4 == 8) THEN
                     IF (ilmn4 == 1) THEN
                        norm4 = norm4_tmp * c4(ia4) *  1.0_DP
                     ELSEIF (ilmn4 == 2) THEN
                        norm4 = norm4_tmp * c4(ia4) * -1.0_DP
                     ELSE
                        CYCLE
                     ENDIF
                  ENDIF

                  hp_acbn0_compute_eri = hp_acbn0_compute_eri + hp_acbn0_coulomb_repulsion(            &
                                x1(1), x1(2), x1(3), norm1,                 &
                                l1(ilmn1), m1(ilmn1), n1(ilmn1), c1(ia1+3), &
                                x2(1), x2(2), x2(3), norm2,                 &
                                l2(ilmn2), m2(ilmn2), n2(ilmn2), c2(ia2+3), &
                                x3(1), x3(2), x3(3), norm3,                 &
                                l3(ilmn3), m3(ilmn3), n3(ilmn3), c3(ia3+3), &
                                x4(1), x4(2), x4(3), norm4,                 &
                                l4(ilmn4), m4(ilmn4), n4(ilmn4), c4(ia4+3))

               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      ! Hartree to Rydberg
      hp_acbn0_compute_eri = 2.0_DP * hp_acbn0_compute_eri

      RETURN

      END FUNCTION hp_acbn0_compute_eri

      SUBROUTINE hp_acbn0_expnorm(lm, l, m, n, norm)

      IMPLICIT NONE

      INTEGER,  PARAMETER  :: DP = selected_real_kind(14,200)

      REAL(DP), PARAMETER  :: PI = 3.14159265358979323846_DP

      INTEGER, INTENT(IN)  :: lm
      INTEGER, INTENT(OUT) :: l(3), m(3), n(3)
      REAL(DP),INTENT(OUT) :: norm


      l(:) = 0 ; m(:) = 0 ; n(:) = 0

      ! NOTE: (2*0+1)!! = 1, (2*1+1)!! = 3, (2*2+1)!! = 15
      SELECT CASE(lm)
         CASE(1)
            l(1) = 0 ; m(1) = 0 ; n(1) = 0 ! 1
            norm  = 0.25_DP * sqrt( 1.0_DP/(PI*0.25_DP))
         CASE(2)
            l(1) = 0 ; m(1) = 0 ; n(1) = 1 ! z
            norm  = 0.25_DP * sqrt( 3.0_DP/(PI*0.25_DP))
         CASE(3)
            l(1) = 1 ; m(1) = 0 ; n(1) = 0 ! x
            norm  = 0.25_DP * sqrt( 3.0_DP/(PI*0.25_DP))
         CASE(4)
            l(1) = 0 ; m(1) = 1 ; n(1) = 0 ! y
            norm  = 0.25_DP * sqrt( 3.0_DP/(PI*0.25_DP))
         CASE(5)
            l(1) = 0 ; m(1) = 0 ; n(1) = 2 ! 2z^2 - x^2 - y^2
            l(2) = 2 ; m(2) = 0 ; n(2) = 0 ! 2z^2 - x^2 - y^2
            l(3) = 0 ; m(3) = 2 ; n(3) = 0 ! 2z^2 - x^2 - y^2
            norm  = 0.25_DP * sqrt(15.0_DP/(PI*3.00_DP))
         CASE(6)
            l(1) = 1 ; m(1) = 0 ; n(1) = 1 ! xz
            norm  = 0.25_DP * sqrt(15.0_DP/(PI*0.25_DP))
         CASE(7)
            l(1) = 0 ; m(1) = 1 ; n(1) = 1 ! yz
            norm  = 0.25_DP * sqrt(15.0_DP/(PI*0.25_DP))
         CASE(8)
            l(1) = 2 ; m(1) = 0 ; n(1) = 0 ! x^2 - y^2
            l(2) = 0 ; m(2) = 2 ; n(2) = 0 ! x^2 - y^2
            norm  = 0.25_DP * sqrt(15.0_DP/(PI*1.00_DP))
         CASE(9)
            l(1) = 1 ; m(1) = 1 ; n(1) = 0 ! xy
            norm  = 0.25_DP * sqrt(15.0_DP/(PI*0.25_DP))
      END SELECT

      END SUBROUTINE hp_acbn0_expnorm
