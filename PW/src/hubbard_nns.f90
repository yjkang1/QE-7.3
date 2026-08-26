! similar to the routine in Modules/dist.f90
SUBROUTINE hubbard_nns()

      USE kinds,     ONLY : dp
      USE constants, ONLY : bohr_radius_angs
      USE cell_base, ONLY : at, alat
      USE ions_base, ONLY : nat, tau
      USE ldaV,      ONLY : ehub_nn_vect, ehub_nn_distance, &
                            ehub_nn_leng, ehub_nn_dim,      &
                            ehub_nn_leng_all, ehub_csize

      REAL(DP) :: rmax, dtau(3)
      REAL(DP) :: r(3), r2

      INTEGER  :: ics
      INTEGER  :: na, nb
      INTEGER  :: n1, n2, n3, ix

      rmax = ehub_nn_distance / bohr_radius_angs / alat
      !
      ehub_nn_vect(:,:,:,:,:) = 0
      ehub_nn_leng(:,:,:,:,:) = 1000.0_DP
      !
      IF ( rmax == 0.0_DP ) RETURN
      !
      ics = ehub_csize + 1
      !
      DO na = 1, nat
      DO nb = 1, nat
         !
         dtau(:) = tau(:, na) - tau(:, nb)
         !
         DO n1 = -1*ehub_csize, ehub_csize
         DO n2 = -1*ehub_csize, ehub_csize
         DO n3 = -1*ehub_csize, ehub_csize
            !
            r2 = 0.0_DP
            !
            DO ix = 1, 3
               r(ix) = n1 * at(ix,1) + n2 * at(ix,2) + n3 * at(ix,3) - dtau(ix)
               r2 = r2 + r(ix) * r(ix)
            ENDDO
            !
            ehub_nn_leng_all(n1+ics,n2+ics,n3+ics,na,nb) = SQRT(r2)
            !
            IF ( r2 <= rmax**2 .AND. ABS(r2) > 1.0D-8 ) THEN
               ehub_nn_vect(n1+ics,n2+ics,n3+ics,na,nb) = 1
               ehub_nn_leng(n1+ics,n2+ics,n3+ics,na,nb) = SQRT(r2)
            ENDIF
            !
         ENDDO ! n3
         ENDDO ! n2
         ENDDO ! n1
         !
      ENDDO ! nb
      ENDDO ! na
      !
      ehub_nn_dim(:) = -1000
      !
      DO n1 = -1*ehub_csize, ehub_csize
      DO n2 = -1*ehub_csize, ehub_csize
      DO n3 = -1*ehub_csize, ehub_csize
         !
         DO na = 1, nat
         DO nb = 1, nat
            !
            IF ( ehub_nn_vect(n1+ics,n2+ics,n3+ics,na,nb) == 1 ) THEN
               IF (ABS(n1) > ehub_nn_dim(1)) ehub_nn_dim(1) = ABS(n1)
               IF (ABS(n2) > ehub_nn_dim(2)) ehub_nn_dim(2) = ABS(n2)
               IF (ABS(n3) > ehub_nn_dim(3)) ehub_nn_dim(3) = ABS(n3)
            ENDIF
            !
         ENDDO
         ENDDO
         !
      ENDDO
      ENDDO
      ENDDO
      !
      RETURN
      !
END SUBROUTINE hubbard_nns


SUBROUTINE hubbard_symmetry(atom1, atom2, opatom1, opatom2, n, nnc, sym)
      !-----------------------------------------------------------------------
      !
      !  This subroutine acts the symmetry operatorion number p_sym over the
      !  atomic position of atom at1 in the reference unit cell and over the
      !  atomic position of atom at2 in the supercell.
      !  Then it identifies, in the reference unit cell,
      !  the equivalent atom to the first 'rotation' above, storing the result
      !  in rat1 and the equivalent atom to the second 'rotation', storing the
      !  result in rat2 (temporarily).
      !  The difference dx = O(p_sym)at1 - rat1 is also subtracted from
      !  O(p_sym)at2 to find the final position of the second atom.
      !  To determine in which unit cell the second atom ended up, we
      !  compute the difference between its final position and the position
      !  of its equivalent in the reference unit cell.
      !  Finally, we store in rat2 the position, in the full list of atoms in
      !  the supercell, the rotated and translated second atom corresponds to.
      !
      !  slightly modified version of symonpair.f90
      !
      USE kinds
      USE symm_base, ONLY : s, ft, nsym 
      USE ions_base, ONLY : nat, ityp, tau
      USE cell_base, ONLY : bg
      USE ldaV,      ONLY : ehub_nn_vect, ehub_nn_dim
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: atom1, atom2
      INTEGER, INTENT(IN)  :: n(3)
      INTEGER, INTENT(IN)  :: sym
      INTEGER, INTENT(OUT) :: opatom1, opatom2
      INTEGER, INTENT(OUT) :: nnc
      !
      ! Local variables
      !
      INTEGER :: i, j
      INTEGER :: n1, n2, n3
      INTEGER :: na
      INTEGER :: opn(3)
      !
      REAL(DP) :: ss(3,3)
      REAL(DP) :: xd(3,nat), x1(3), x2(3)
      REAL(DP) :: diff, dx(3)
      !
      REAL(DP), PARAMETER :: eps = 5.0D-5
      !
      ! Convert the symmetry matrix from integer to real type
      ! for a given symmetry operation p_sym
      !
      DO i = 1, 3
         DO j = 1, 3
            ss(i,j) = DBLE(s(i,j,sym))
         ENDDO
      ENDDO
      !
      ! Compute the coordinates of each atom
      ! in the basis of the direct lattice vectors
      !
      DO na = 1, nat
         DO i = 1, 3
            xd(i,na) = bg(1,i)*tau(1,na) &
                     + bg(2,i)*tau(2,na) &
                     + bg(3,i)*tau(3,na)
         ENDDO
      ENDDO
      !
      ! Apply the symmetry operation: S*atom - f
      !
      DO i = 1, 3
         !
         x1(i) = 0.D0
         x2(i) = 0.D0
         !
         ! Rotate atom1 and atom2
         !
         DO j = 1, 3
            x1(i) = x1(i) + (xd(j,atom1)             ) * ss(j,i)
            x2(i) = x2(i) + (xd(j,atom2) + DBLE(n(j))) * ss(j,i)
         ENDDO
         !
         ! Subtract vectors of fractional translations
         !
         x1(i) = x1(i) - ft(i,sym)
         x2(i) = x2(i) - ft(i,sym)
         !
      ENDDO
      !
      ! Atom equivalent to x2 : irt(sym,atom2)
      !
      diff = 1.D0
      na = 1
      !
      DO WHILE ( na <= nat .AND. diff > eps )
         !
         IF ( ityp(na) == ityp(atom2) ) THEN
            diff = 0.D0
            DO i = 1, 3
               dx(i) = x2(i) - xd(i,na)
               diff = diff + ABS(dx(i)-NINT(dx(i)))
            ENDDO
         ELSE
            diff = 1.D0
         ENDIF
         !
         na = na + 1
         !
      ENDDO
      !
      IF ( diff > eps ) THEN
         WRITE(*,*) "diff > 0, diff= ", diff, "atom1= ", atom1, "atom2= ", atom2
         CALL errore('hubbard_symmetry', 'No atom equivalent to x2', 1)
      ENDIF
      !
      opatom2 = na - 1
      !
      ! Atom equivalent to x1 : irt(sym,atom1)
      !
      diff = 1.D0
      na = 1
      !
      DO WHILE ( na <= nat .AND. diff > eps )
         !
         IF ( ityp(na) == ityp(atom1) ) THEN
            diff = 0.D0
            DO i = 1, 3
               dx(i) = x1(i) - xd(i,na)
               diff = diff + ABS(dx(i)-NINT(dx(i)))
            ENDDO
         ELSE
            diff = 1.D0
         ENDIF
         !
         na = na + 1
         !
      ENDDO
      !
      IF ( diff > eps ) THEN
        WRITE(*,*) "diff > 0, diff= ", diff, "atom1= ", atom1, "atom2= ", atom2
        CALL errore('hubbard_symmetry', 'No atom equivalent to x1', 1)
      ENDIF
      !
      opatom1 = na - 1
      !
      DO i = 1, 3
         x2(i) = x2(i) - dx(i)
         opn(i) = NINT(x2(i) - xd(i,opatom2))
      ENDDO
      !
      CALL hubbard_get_nnc(opn,nnc)
      !
      RETURN
      !
END SUBROUTINE hubbard_symmetry


SUBROUTINE hubbard_get_nnc(n, nnc)

      USE ions_base, ONLY : nat
      USE ldaV,      ONLY : ehub_nn_vect, ehub_nn_dim, ehub_csize

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n(3)
      INTEGER, INTENT(OUT) :: nnc

      INTEGER              :: ics
      INTEGER              :: n1, n2, n3, opnnc

      ics = ehub_csize + 1
      !
      IF ( n(1) == 0 .AND. n(2) == 0 .AND. n(3) == 0 ) THEN
         nnc = 1
      ELSE
         opnnc = 2
         !
         DO n1 = -1*ehub_nn_dim(1), ehub_nn_dim(1)
         DO n2 = -1*ehub_nn_dim(2), ehub_nn_dim(2)
         DO n3 = -1*ehub_nn_dim(3), ehub_nn_dim(3)
            !
            IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
            !
            IF ( MAXVAL(ehub_nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) THEN
               IF (n(1) == n1 .AND. n(2) == n2 .AND. n(3) == n3) THEN
                  nnc = opnnc
               ELSE
                  opnnc = opnnc + 1
               ENDIF
            ENDIF
            !
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      RETURN

END SUBROUTINE hubbard_get_nnc
