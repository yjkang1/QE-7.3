!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
!  eACBN0.f90  --  ACBN0 / eACBN0 Hubbard parameters inside hp.x
!-----------------------------------------------------------------------
!
!  Driven by the input tag  acbn0  of &INPUTHP:
!
!     acbn0 = 0 (default)    DFPT linear response  -> official hp.x, untouched
!     acbn0 = 1              on-site U             -> ACBN0
!     acbn0 = 2              on-site U + V         -> eACBN0        (stage S4)
!
!  When acbn0 > 0 hp.x computes the ACBN0 parameters from the ground state in
!  outdir and stops; it never enters the DFPT q-loop.  The point is that both
!  estimators - ACBN0 and linear response - can then be produced from ONE
!  ground state by one binary, which is what makes them comparable.
!
!-----------------------------------------------------------------------
!  STAGE S3 - ldaU bookkeeping
!-----------------------------------------------------------------------
!
!  S1/S2 read the Hubbard manifold out of MODULE ldaV, which only exists if
!  pw.x was run with the fork's lacbn0 = .true. and wrote an <ehub> block to
!  the XML.  That made the S0/S1/S2 comparison exact, but it also meant hp.x
!  could only do ACBN0 on a ground state produced by the patched pw.x.
!
!  S3 moves the bookkeeping to MODULE ldaU - official QE's own - so ACBN0 runs on
!  ANY ortho-atomic DFT+U ground state, including the very one the DFPT branch
!  of hp.x uses.  Nothing about the physics changes; only where the indices
!  come from.
!
!  The translation, as of stage V4 (background channels):
!
!      ldaV                        ldaU                        this file
!      ------------------------    ------------------------    -----------
!      wfcv                        wfcU                        wfcU
!      nwfcv                       nwfcU                       nwfcU
!      ofstv(na) + m               offsetU / offsetU_back(1)   ofs(m,na)
!      iunhubv / nwordwfcv         iunhub / nwordwfcU          idem
!      ehub_lm(nt)                 ldim_u(nt)                  lm(nt)
!      ehub_lm_max                 ldmx_tot                    lm_max
!      ehub_l(m,nt)                ll(m,nt) + 1                lch(m,nt)
!      orthowfcv()                 orthoUwfc(.FALSE.)          idem
!      copy_v_wfc                  copy_U_wfc                  (inside those)
!      ehub_u(na,m1,m2)            -                           hub_u argument
!      ehub_l_choice(nt,l+1)       -                           acbn0_l_choice (ldaU_hp)
!      remove_ehub_u(nt,l+1)       -                           acbn0_rm_ehub_u (ldaU_hp, 1 = remove)
!
!  Stages S3 and V3 collapsed ehub_l(m,nt) to the scalar Hubbard_l(nt)+1,
!  which is correct only while a type carries ONE channel, and multi-channel
!  types were refused in hp_acbn0_check.  Stage V4 removes that restriction.
!  ldaU does hold the same information - ll(m,nt) labels the concatenated
!  manifold exactly as ehub_l does, and ldim_u(nt) is its length - it merely
!  stores the PROJECTORS in up to three separate blocks (offsetU,
!  offsetU_back, offsetU_back1) instead of contiguously.  MODULE
!  hp_acbn0_man rebuilds the contiguous view as ofs(m,na); see it for the
!  layout and for why the XML round-trip is safe.
!
!  Three things S2 needed and S3 does not, all of them ldaV artefacts:
!
!    - remove_ehub_u / remove_ehub_v.  Uninitialised in hp.x (set only in
!      pw.x's dftU_iosys_v1, and not round-tripped through the XML), and
!      new_ns_a0 multiplies the final U by them twice.  That was Bug 2 of the
!      S1 run: a silently scaled or zeroed U.  ldaU has no such mask, so the
!      bug class is gone rather than worked around.
!    - ehub_mixing.  A self-consistency knob for pw.x's ACBN0 loop; hp.x is
!      one-shot, so it was forced to 1.0.  Now simply absent.
!    - the ehub_u shape guard.  ehub_u was a module array whose XML-stored
!      dimension could disagree with the recomputed one.  hub_u is now a local
!      allocated here, so it cannot disagree with anything.
!
!  What stays, because both are properties of hp.x's read path and not of the
!  bookkeeping:
!
!    - CALL d_matrix.  read_file only builds d1/d2/d3 under okpaw
!      (read_file_new.f90:348-351) and hp_acbn0_symm_n4 needs them.
!    - reopening iunwfc.  read_file closes it at read_file_new.f90:62 with
!      close_buffer(...,'KEEP'); 'KEEP' preserves the FILES but the call also
!      deregisters the unit, so get_buffer falls through to davcio on a closed
!      unit (Bug 1 of the S1 run).  Official hp.x never trips this because it
!      reads wavefunctions through iuwfc (units_lr), opened in hp_openfil_q.
!
!  Unit note: hp.x also uses iunhub, in hp_openfil_q.f90:76 and
!  hp_close_q.f90:39.  There is no collision because acbn0 > 0 returns from
!  hp_main before the q-loop opens it, but this file must therefore close it
!  itself.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
MODULE hp_acbn0_nn
  !---------------------------------------------------------------------
  !
  ! STAGE V3 - the neighbour bookkeeping for the eACBN0 (+V) branch, owned by
  ! hp.x instead of by MODULE ldaV.
  !
  ! Everything the +V kernels need from ldaV that is NOT the Hubbard manifold
  ! reduces to one table:
  !
  !     ldaV                V3 (here)     meaning
  !     ehub_csize          nn_csize      half-width of the cell grid
  !     ehub_nn_vect        nn_vect       1 if (na2 in cell n) neighbours na1
  !     ehub_nn_leng        nn_leng       the corresponding distance
  !     ehub_nn_dim         nn_dim        largest |n| actually used, per axis
  !     ehub_nnnc           nnnc          number of distinct occupied cells
  !
  ! The manifold and the projectors come from official QE's MODULE ldaU
  ! (Hubbard_l, Hubbard_lmax, offsetU, wfcU, nwfcU), exactly as the +U branch
  ! has done since S3.  So the +V kernels no longer touch ldaV at all, and
  ! acbn0 = 2 runs on an ordinary DFT+U ground state.
  !
  ! WHY A CUTOFF AND NOT ldaU's neighood/at_sc.  Lee & Son define the
  ! inter-site set by a cutoff distance d_IJ (PRR 2, 043410 (2020)); the fork
  ! implements exactly that, in hubbard_nns(), from tau/at/alat alone - no
  ! bookkeeping of either flavour is involved.  ldaU's neighood/at_sc is not
  ! an alternative definition, it is how official QE stores a pair list the user
  ! typed into the HUBBARD card.  Taking the cutoff directly is therefore both
  ! closer to the method and more useful: hp.x can now sweep d_IJ without
  ! re-running pw.x, which the fork cannot do.
  !
  ! hp_acbn0_build_nn reproduces hubbard_nns() (PW/src/hubbard_nns.f90:2) and
  ! the nnnc count in input.f90:2144-2180 verbatim, so that for the same
  ! cutoff the table is identical to the one init_lda_plus_v would have built.
  ! That equivalence is what the V3a checkpoint tests.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  !
  INTEGER, PARAMETER :: nn_csize = 10       ! = ehub_csize
  !
  INTEGER,  ALLOCATABLE :: nn_vect(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: nn_leng(:,:,:,:,:)
  INTEGER               :: nn_dim(3)
  INTEGER               :: nnnc
  REAL(DP)              :: nn_distance             ! Angstrom
  !
CONTAINS
  !
  !-------------------------------------------------------------------
  SUBROUTINE hp_acbn0_build_nn( dist )
    !-----------------------------------------------------------------
    !
    ! Build the neighbour table for a cutoff of "dist" Angstrom.
    ! Line-for-line equivalent to hubbard_nns() plus the nnnc count that
    ! PW/src/input.f90:2144 performs, with the ldaV arrays replaced by ours.
    !
    USE io_global, ONLY : stdout
    USE constants, ONLY : bohr_radius_angs
    USE cell_base, ONLY : at, alat
    USE ions_base, ONLY : nat, tau
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: dist
    !
    REAL(DP) :: rmax, dtau(3), r(3), r2
    INTEGER  :: ics, na, nb, n1, n2, n3, ix, counter
    !
    nn_distance = dist
    ics  = nn_csize + 1
    rmax = dist / bohr_radius_angs / alat
    !
    IF ( .NOT.ALLOCATED(nn_vect) ) &
       ALLOCATE( nn_vect(2*nn_csize+1,2*nn_csize+1,2*nn_csize+1,nat,nat) )
    IF ( .NOT.ALLOCATED(nn_leng) ) &
       ALLOCATE( nn_leng(2*nn_csize+1,2*nn_csize+1,2*nn_csize+1,nat,nat) )
    !
    nn_vect(:,:,:,:,:) = 0
    nn_leng(:,:,:,:,:) = 1000.0_DP
    nn_dim(:)          = -1000
    nnnc               = 1          ! index 1 is always the home cell
    !
    IF ( rmax == 0.0_DP ) THEN
       ! +U only: no inter-site partners, one cell.  This is the route
       ! Howto2.pdf recommends for a pure DFT+U calculation.
       nn_dim(:) = 0
       RETURN
    ENDIF
    !
    DO na = 1, nat
    DO nb = 1, nat
       dtau(:) = tau(:,na) - tau(:,nb)
       DO n1 = -nn_csize, nn_csize
       DO n2 = -nn_csize, nn_csize
       DO n3 = -nn_csize, nn_csize
          r2 = 0.0_DP
          DO ix = 1, 3
             r(ix) = n1*at(ix,1) + n2*at(ix,2) + n3*at(ix,3) - dtau(ix)
             r2 = r2 + r(ix)*r(ix)
          ENDDO
          IF ( r2 <= rmax**2 .AND. ABS(r2) > 1.0D-8 ) THEN
             nn_vect(n1+ics,n2+ics,n3+ics,na,nb) = 1
             nn_leng(n1+ics,n2+ics,n3+ics,na,nb) = SQRT(r2)
          ENDIF
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    !
    ! nn_dim: the largest |n| that carries any pair, per axis.
    !
    DO n1 = -nn_csize, nn_csize
    DO n2 = -nn_csize, nn_csize
    DO n3 = -nn_csize, nn_csize
       IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) THEN
          IF (ABS(n1) > nn_dim(1)) nn_dim(1) = ABS(n1)
          IF (ABS(n2) > nn_dim(2)) nn_dim(2) = ABS(n2)
          IF (ABS(n3) > nn_dim(3)) nn_dim(3) = ABS(n3)
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !
    IF ( MAXVAL(nn_dim) > 10 ) CALL errore( 'hp_acbn0_build_nn', &
       'reduce nn_distance (nn_dim > 10)', 1 )
    IF ( MAXVAL(nn_dim) < 0 ) THEN
       nn_dim(:) = 0
       RETURN
    ENDIF
    !
    ! nnnc: number of distinct cells carrying at least one pair, with index 1
    ! reserved for the home cell whether or not it carries one.  Same count
    ! input.f90:2144-2180 performs.
    !
    counter = 0
    nnnc    = 0
    DO n1 = -nn_csize, nn_csize
    DO n2 = -nn_csize, nn_csize
    DO n3 = -nn_csize, nn_csize
       IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) THEN
          IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) counter = 1
          nnnc = nnnc + 1
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    IF ( counter == 0 ) nnnc = nnnc + 1
    !
    ! Report, in the same form pw.x uses, so the two are comparable by eye.
    !
    WRITE( stdout,'(5x,a)' ) '-----------------------------------------------------------'
    WRITE( stdout,'(5x,a)' ) 'List of the nearest neighbors for inter-site interactions V'
    WRITE( stdout,'(5x,a)' ) '-----------------------------------------------------------'
    DO n1 = -nn_dim(1), nn_dim(1)
    DO n2 = -nn_dim(2), nn_dim(2)
    DO n3 = -nn_dim(3), nn_dim(3)
       DO na = 1, nat
       DO nb = 1, nat
          IF ( nn_vect(n1+ics,n2+ics,n3+ics,na,nb) == 0 ) CYCLE
          WRITE( stdout,'(5x,a,i2,a,i2,a,3i6)' ) &
                'atom 1: ', na, ',  atom 2: ', nb, ',', n1, n2, n3
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    WRITE( stdout,'(/5x,"nn_distance = ",f6.3," A     nnnc = ",i4,/)') dist, nnnc
    !
    RETURN
    !
  END SUBROUTINE hp_acbn0_build_nn
  !
  !-------------------------------------------------------------------
  SUBROUTINE hp_acbn0_free_nn()
    IMPLICIT NONE
    IF ( ALLOCATED(nn_vect) ) DEALLOCATE( nn_vect )
    IF ( ALLOCATED(nn_leng) ) DEALLOCATE( nn_leng )
    RETURN
  END SUBROUTINE hp_acbn0_free_nn
  !
  SUBROUTINE hp_acbn0_symmetry(atom1, atom2, opatom1, opatom2, n, nnc, sym)
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
         CALL errore('hp_acbn0_symmetry', 'No atom equivalent to x2', 1)
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
        CALL errore('hp_acbn0_symmetry', 'No atom equivalent to x1', 1)
      ENDIF
      !
      opatom1 = na - 1
      !
      DO i = 1, 3
         x2(i) = x2(i) - dx(i)
         opn(i) = NINT(x2(i) - xd(i,opatom2))
      ENDDO
      !
      CALL hp_acbn0_get_nnc(opn,nnc)
      !
      RETURN
      !
  END SUBROUTINE hp_acbn0_symmetry


  SUBROUTINE hp_acbn0_get_nnc(n, nnc)

      USE ions_base, ONLY : nat

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n(3)
      INTEGER, INTENT(OUT) :: nnc

      INTEGER              :: ics
      INTEGER              :: n1, n2, n3, opnnc

      ics = nn_csize + 1
      !
      ! BUGFIX vs stock PW/src/hubbard_nns.f90:241.  Stock leaves nnc UNDEFINED
      ! whenever the requested cell carries no neighbour pair: the enumeration
      ! below only assigns inside "IF (MAXVAL(nn_vect(...)) == 1)", so an empty
      ! cell falls through untouched.  hubbard_symmetry asks for the index of a
      ! ROTATED cell, which is not guaranteed to be one of the occupied ones,
      ! and its caller then indexes n(nnc) with whatever was on the stack.
      ! With bounds checking:
      !     Subscript #1 of the array N has value 1059660741 ... upper bound 19
      ! Returning 0 lets the callers skip such a pair, which is the correct
      ! contribution: a rotated image outside the neighbour set contributes
      ! nothing to the symmetrised occupation matrix.
      !
      nnc = 0
      !
      IF ( n(1) == 0 .AND. n(2) == 0 .AND. n(3) == 0 ) THEN
         nnc = 1
      ELSE
         opnnc = 2
         !
         DO n1 = -1*nn_dim(1), nn_dim(1)
         DO n2 = -1*nn_dim(2), nn_dim(2)
         DO n3 = -1*nn_dim(3), nn_dim(3)
            !
            IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
            !
            IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) THEN
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

  END SUBROUTINE hp_acbn0_get_nnc

END MODULE hp_acbn0_nn

!-----------------------------------------------------------------------
MODULE hp_acbn0_man
  !---------------------------------------------------------------------
  !
  ! STAGE V4 - the Hubbard MANIFOLD in the shape the transplanted kernels
  ! expect, background channels included.
  !
  ! The fork's MODULE ldaV stores a type's manifold as ONE concatenated
  ! block: ehub_lm(nt) states, ehub_l(m,nt) giving the l+1 of each state, and
  ! the projectors laid out contiguously in wfcv from ofstv(na).  Every
  ! kernel in PW/src/new_ns_v.f90 is written against that shape.
  !
  ! Official ldaU holds the same information split across up to three pieces -
  ! the first manifold at offsetU(na), background ones at offsetU_back(na)
  ! and offsetU_back1(na).  Crucially it ALREADY labels them the way the
  ! kernels want:
  !
  !     ldaV                ldaU                    here
  !     ehub_lm(nt)         ldim_u(nt)              lm(nt)
  !     ehub_lm_max         ldmx_tot                lm_max
  !     ehub_l(m,nt)        ll(m,nt) + 1            lch(m,nt)
  !     ofstv(na) + m       (three-way, below)      ofs(m,na)
  !
  ! ll(m,nt) is built by init_hubbard (PW/src/ldaU.f90:290-320) and runs over
  ! the concatenated manifold in the same order ldaV uses - Hubbard_l first,
  ! then Hubbard_l2, then Hubbard_l3.  So the only thing genuinely missing is
  ! the map from a manifold index to a wfcU column, since ldaU does not store
  ! the manifold contiguously.  That map is ofs(m,na), and it follows
  ! force_hub.f90:1060-1071 and copy_U_wfc (ldaU.f90:604-620) exactly:
  !
  !     m <=   2*l+1          ->  ofs(m,na)
  !     m <= 2*(l+l2+1)       ->  offsetU_back(na)  + m -   (2*l+1)
  !     otherwise, backall    ->  offsetU_back1(na) + m - 2*(l+l2+1)
  !
  ! This lifts stage V3's one restriction - a single Hubbard channel per
  ! atomic type - so GaAs (Ga 4s+4p) runs on the ldaU branch.
  !
  ! ON THE XML.  Hubbard_l2/l3/U2/backall DO round-trip through the schema
  ! (qexsd_copy.f90:434-441, called at pw_restart_new.f90:1322) and
  ! init_hubbard is re-run on hp.x's read path, so ll, ldim_u, ldmx_tot,
  ! offsetU_back and offsetU_back1 are all valid here.  This was checked
  ! against the same criterion that found bugs 1-3: is it assigned in a file
  ! hp.x actually executes?
  !
  ! A background channel is ACTIVE only when Hubbard_U2 /= 0 (ldaU.f90:249),
  ! so a HUBBARD card that wants one must give it a nonzero U2.  0.0001 eV is
  ! enough and is numerically nothing - ACBN0 recomputes U from scratch and
  ! never reads the value.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  !
  INTEGER              :: lm_max = 0   ! = ehub_lm_max
  INTEGER, ALLOCATABLE :: lm(:)        ! (ntyp)         = ehub_lm
  INTEGER, ALLOCATABLE :: lch(:,:)     ! (lm_max,ntyp)  = ehub_l
  INTEGER, ALLOCATABLE :: ofs(:,:)     ! (lm_max,nat)   = ofstv(na) + m
  !
CONTAINS
  !
  !-------------------------------------------------------------------
  SUBROUTINE hp_acbn0_build_man()
    !-----------------------------------------------------------------
    !
    USE ions_base, ONLY : nat, ityp, ntyp => nsp
    USE io_global, ONLY : stdout
    USE ldaU,      ONLY : Hubbard_l, Hubbard_l2, Hubbard_l3, backall,     &
                          is_hubbard, is_hubbard_back, ll, ldim_u,        &
                          ldmx_tot, offsetU, offsetU_back, offsetU_back1, &
                          nwfcU, lda_plus_u_kind
    !
    IMPLICIT NONE
    !
    INTEGER            :: na, nt, m, n0, n0b
    CHARACTER(LEN=256) :: msg
    !
    CALL hp_acbn0_free_man()
    !
    ! Liechtenstein's formulation never builds ll / ldim_u / ldmx_tot
    ! (ldaU.f90:323-364), so there is no manifold table to read.
    !
    IF ( lda_plus_u_kind == 1 ) CALL errore( 'hp_acbn0_build_man', &
       'acbn0 needs lda_plus_u_kind = 0 or 2; the full (Liechtenstein) '// &
       'formulation does not build the ll manifold table', 1 )
    !
    IF ( .NOT.ALLOCATED(ll) .OR. .NOT.ALLOCATED(ldim_u) ) &
       CALL errore( 'hp_acbn0_build_man', &
       'MODULE ldaU carries no manifold table: the pw.x run had no '// &
       'HUBBARD card', 1 )
    !
    lm_max = ldmx_tot
    !
    IF ( lm_max < 1 ) CALL errore( 'hp_acbn0_build_man', &
       'ldmx_tot < 1: the ground state declares no Hubbard states', 1 )
    !
    ALLOCATE ( lm (ntyp) )        ; lm  = 0
    ALLOCATE ( lch(lm_max,ntyp) ) ; lch = 0
    ALLOCATE ( ofs(lm_max,nat) )  ; ofs = 0
    !
    ! ... the manifold, per type
    !
    DO nt = 1, ntyp
       !
       IF ( .NOT.( is_hubbard(nt) .OR. is_hubbard_back(nt) ) ) CYCLE
       !
       lm(nt) = ldim_u(nt)
       !
       DO m = 1, lm(nt)
          ! ll is -1 on slots that carry no state, which maps onto the same
          ! 0 sentinel every "IF (lm(nt) == 0) CYCLE" in the kernels tests.
          IF ( ll(m,nt) >= 0 ) lch(m,nt) = ll(m,nt) + 1
       ENDDO
       !
    ENDDO
    !
    ! ... the manifold index -> wfcU column map, per atom
    !
    DO na = 1, nat
       !
       nt = ityp(na)
       IF ( lm(nt) == 0 ) CYCLE
       !
       n0 = 0
       !
       IF ( is_hubbard(nt) ) THEN
          DO m = 1, 2*Hubbard_l(nt) + 1
             ofs(m,na) = offsetU(na) + m
          ENDDO
          n0 = 2*Hubbard_l(nt) + 1
       ENDIF
       !
       IF ( is_hubbard_back(nt) ) THEN
          !
          IF ( .NOT.ALLOCATED(offsetU_back) ) CALL errore( 'hp_acbn0_build_man', &
             'is_hubbard_back is set but offsetU_back was never allocated', na )
          !
          DO m = 1, 2*Hubbard_l2(nt) + 1
             ofs(n0+m,na) = offsetU_back(na) + m
          ENDDO
          n0b = n0 + 2*Hubbard_l2(nt) + 1
          !
          IF ( backall(nt) ) THEN
             IF ( .NOT.ALLOCATED(offsetU_back1) ) CALL errore( 'hp_acbn0_build_man', &
                'backall is set but offsetU_back1 was never allocated', na )
             DO m = 1, 2*Hubbard_l3(nt) + 1
                ofs(n0b+m,na) = offsetU_back1(na) + m
             ENDDO
          ENDIF
          !
       ENDIF
       !
       ! Every state in the manifold must have landed on a real wfcU column.
       ! If ldim_u and the offsets ever disagree this fires here, loudly,
       ! instead of silently indexing column 0 or past nwfcU - which is the
       ! failure mode bugs 3 and 4 had.
       !
       DO m = 1, lm(nt)
          IF ( ofs(m,na) < 1 .OR. ofs(m,na) > nwfcU ) THEN
             WRITE( msg, '("atom ",i4," type ",i3," state ",i3," of ",i3, &
                   &" maps to wfcU column ",i6," of ",i6)' ) &
                   na, nt, m, lm(nt), ofs(m,na), nwfcU
             CALL errore( 'hp_acbn0_build_man', TRIM(msg), na )
          ENDIF
       ENDDO
       !
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE hp_acbn0_build_man
  !
  !-------------------------------------------------------------------
  SUBROUTINE hp_acbn0_free_man()
    !-----------------------------------------------------------------
    IMPLICIT NONE
    IF ( ALLOCATED(lm)  ) DEALLOCATE( lm  )
    IF ( ALLOCATED(lch) ) DEALLOCATE( lch )
    IF ( ALLOCATED(ofs) ) DEALLOCATE( ofs )
    lm_max = 0
    RETURN
  END SUBROUTINE hp_acbn0_free_man
  !
  !-------------------------------------------------------------------
  FUNCTION hp_acbn0_chanlab( nt, c ) RESULT( lab )
    !-----------------------------------------------------------------
    !
    ! Printable label of ONE channel of type nt: c is the lch value
    ! (= l+1), the result is e.g. '4s' or '3d'.  The principal quantum
    ! number comes from whichever of Hubbard_l / l2 / l3 carries that l.
    !
    USE ldaU, ONLY : Hubbard_n, Hubbard_l, Hubbard_n2, Hubbard_l2, &
                     Hubbard_n3, Hubbard_l3, backall
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nt, c
    CHARACTER(LEN=2)    :: lab
    CHARACTER(LEN=1)    :: spdf(4) = ['s','p','d','f']
    INTEGER             :: n
    !
    lab = '??'
    IF ( c < 1 .OR. c > 4 ) RETURN
    n = -1
    IF ( Hubbard_l(nt)  == c-1 ) THEN
       n = Hubbard_n(nt)
    ELSEIF ( Hubbard_l2(nt) == c-1 ) THEN
       n = Hubbard_n2(nt)
    ELSEIF ( backall(nt) .AND. Hubbard_l3(nt) == c-1 ) THEN
       n = Hubbard_n3(nt)
    ENDIF
    IF ( n >= 1 .AND. n <= 9 ) THEN
       WRITE( lab, '(i1,a1)' ) n, spdf(c)
    ELSE
       lab = ' ' // spdf(c)
    ENDIF
    RETURN
  END FUNCTION hp_acbn0_chanlab
  !
  !-------------------------------------------------------------------
  FUNCTION hp_acbn0_manlab( nt ) RESULT( lab )
    !-----------------------------------------------------------------
    !
    ! Printable label for a type's manifold: "d" for a single channel,
    ! "sp" when a background channel is present.  Under stage V3 this was
    ! always one character, because a type had exactly one channel.
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nt
    CHARACTER(LEN=4)    :: lab
    CHARACTER(LEN=1)    :: spdf(4) = ['s','p','d','f']
    INTEGER             :: m
    !
    lab = ' '
    IF ( .NOT.ALLOCATED(lm) ) RETURN
    DO m = 1, lm(nt)
       IF ( lch(m,nt) < 1 .OR. lch(m,nt) > 4 ) CYCLE
       IF ( INDEX(TRIM(lab), spdf(lch(m,nt))) == 0 ) &
          lab = TRIM(lab)//spdf(lch(m,nt))
    ENDDO
    !
    RETURN
    !
  END FUNCTION hp_acbn0_manlab
  !
  !
END MODULE hp_acbn0_man
!-----------------------------------------------------------------------


SUBROUTINE hp_acbn0_u()
  !---------------------------------------------------------------------
  !
  ! Driver.  Everything hp.x needs that read_file did not already do.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunhub, nwordwfcU, iunwfc, nwordwfc
  USE buffers,          ONLY : open_buffer, close_buffer
  USE ions_base,        ONLY : nat
  USE lsda_mod,         ONLY : nspin
  USE wvfct,            ONLY : npwx, nbnd
  USE control_flags,    ONLY : io_level
  USE noncollin_module, ONLY : npol
  USE symm_base,        ONLY : d1, d2, d3
  USE ldaU,             ONLY : Hubbard_lmax, nwfcU
  USE ldaU_hp,          ONLY : acbn0
  USE hp_acbn0_man,         ONLY : lm_max, hp_acbn0_build_man, &
                                   hp_acbn0_free_man
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: ns(:,:,:,:)      ! occupation matrix
  REAL(DP), ALLOCATABLE :: hub_u(:,:,:)     ! hub_u(na,m1,m2), in Ry
  !
  LOGICAL               :: exst
  !
  CALL start_clock('hp_acbn0')
  !
  WRITE( stdout, '(/5x,70("="))')
  WRITE( stdout, '(5x,"ACBN0 Hubbard parameters (acbn0 = ",i1,")")') acbn0
  WRITE( stdout, '(5x,"pseudo-hybrid estimator - NOT linear response")')
  WRITE( stdout, '(5x,70("="),/)')
  !
  CALL hp_acbn0_check()
  !
  ! Build the manifold table (lm / lch / ofs).  This is what makes the
  ! kernels multi-channel: see MODULE hp_acbn0_man.
  !
  CALL hp_acbn0_build_man()
  !
  !
  ! d1/d2/d3 are the Wigner rotation matrices used by hp_acbn0_symm_n4.
  ! read_file builds them only under okpaw, and our pseudos are norm-conserving.
  !
  CALL d_matrix( d1, d2, d3 )
  !
  ! Reopen iunwfc - see the Bug 1 note in the file header.  These are
  ! read_file's own two lines (read_file_new.f90:59-62).
  !
  nwordwfc = nbnd * npwx * npol
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )
  IF ( io_level /= 0 .AND. .NOT.exst ) CALL errore( 'hp_acbn0_u', &
     'the KS wavefunction files (prefix.wfc*) are missing from outdir - &
     &rerun pw.x, or set disk_io to at least "low"', 1 )
  !
  ! Projector buffer.  In pw.x this is opened by openfil.f90; hp.x's own
  ! open of iunhub lives in hp_openfil_q, which acbn0 > 0 never reaches.
  !
  nwordwfcU = npwx * nwfcU * npol
  CALL open_buffer( iunhub, 'hub', nwordwfcU, io_level, exst )
  !
  ! Fill wfcU = O^{-1/2} S phi on the Hubbard manifold, for every k point.
  !
  CALL orthoUwfc( .FALSE. )
  !
  ALLOCATE( ns(lm_max,lm_max,nspin,nat) )   ; ns(:,:,:,:)    = 0.0_DP
  ALLOCATE( hub_u(nat,lm_max,lm_max) )      ; hub_u(:,:,:)   = 0.0_DP
  !
  CALL hp_acbn0_ns_a0( ns, hub_u )
  !
  CALL hp_acbn0_write( ns, hub_u )
  CALL hp_acbn0_dump_u( ns, hub_u )
  !
  DEALLOCATE( ns, hub_u )
  CALL hp_acbn0_free_man()
  !
  CALL close_buffer( iunhub, 'delete' )   ! ours, scratch
  CALL close_buffer( iunwfc, 'KEEP'   )   ! the user's wavefunctions - NEVER delete
  !
  CALL stop_clock('hp_acbn0')
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_u

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_check()
  !---------------------------------------------------------------------
  !
  ! Everything the ACBN0 kernels silently assume.  Each of these fails as a
  ! wrong number rather than as a crash if it is not caught here.
  !
  USE io_global,        ONLY : stdout
  USE ions_base,        ONLY : ntyp => nsp, atm
  USE lsda_mod,         ONLY : nspin
  USE control_flags,    ONLY : gamma_only
  USE ldaU,             ONLY : lda_plus_u, Hubbard_projectors, Hubbard_lmax, &
                               Hubbard_l, is_hubbard, is_hubbard_back, nwfcU, &
                               Hubbard_l2, Hubbard_l3, backall
  !
  IMPLICIT NONE
  !
  INTEGER :: nt, nath
  !
  IF ( .NOT.lda_plus_u ) CALL errore( 'hp_acbn0_check', &
     'acbn0 > 0 needs a Hubbard manifold, but none was selected: set &
     &acbn0_l_choice(nt,b) = 1 in &INPUTHP (the pw.x ground state itself must &
     &be a plain non-Hubbard SCF).', 1 )
  !
  ! ortho-atomic is a requirement, not a preference: ACBN0's occupation
  ! ratio is derived for an orthonormal projector basis.  Note the change
  ! from S1/S2, where orthowfcv silently FORCED Hubbard_projectors to
  ! ortho-atomic whatever the input said; orthoUwfc honours the input, so
  ! the contract has to be enforced here instead.
  !
  IF ( TRIM(Hubbard_projectors) /= 'ortho-atomic' ) CALL errore( 'hp_acbn0_check', &
     'ACBN0 requires ortho-atomic Hubbard projectors; rerun pw.x with &
     &HUBBARD (ortho-atomic).', 1 )
  !
  IF ( Hubbard_lmax < 0 .OR. nwfcU <= 0 ) CALL errore( 'hp_acbn0_check', &
     'no Hubbard manifold was found (Hubbard_lmax < 0 or nwfcU = 0).', 1 )
  !
  ! Multi-channel types are supported from stage V4 on: MODULE hp_acbn0_man
  ! rebuilds the concatenated manifold ldaV used out of ldaU's ll(m,nt) plus
  ! the offsetU / offsetU_back / offsetU_back1 blocks.  What must still hold
  ! is that ldaU actually laid the background out, which it only does when
  ! the channel is ACTIVE - i.e. Hubbard_U2 /= 0 (ldaU.f90:249).
  !
  DO nt = 1, ntyp
     IF ( is_hubbard_back(nt) .AND. Hubbard_l2(nt) < 0 ) &
        CALL errore( 'hp_acbn0_check', &
        'a background Hubbard channel is active (Hubbard_U2 /= 0) but &
        &Hubbard_l2 was never set for this type.', nt )
     IF ( backall(nt) .AND. Hubbard_l3(nt) < 0 ) &
        CALL errore( 'hp_acbn0_check', &
        'backall is set for this type but Hubbard_l3 was never set.', nt )
  ENDDO
  !
  ! f manifolds: hp_acbn0_pao_3g carries 3-Gaussian PAO fits for s, p and d only, and
  ! returns c = 0 without complaining for anything else - which would make
  ! every eri vanish and U come out as 0/0.
  !
  DO nt = 1, ntyp
     IF ( is_hubbard(nt) .AND. Hubbard_l(nt) > 2 ) CALL errore( 'hp_acbn0_check', &
        'ACBN0 is implemented for s, p and d manifolds only: hp_acbn0_pao_3g has no &
        &f-shell Gaussian fit and returns zero coefficients silently.', nt )
     IF ( is_hubbard_back(nt) ) THEN
        IF ( Hubbard_l2(nt) > 2 ) CALL errore( 'hp_acbn0_check', &
           'background channel Hubbard_l2 is an f manifold; hp_acbn0_pao_3g has no &
           &f-shell fit.', nt )
        IF ( backall(nt) .AND. Hubbard_l3(nt) > 2 ) CALL errore( 'hp_acbn0_check', &
           'background channel Hubbard_l3 is an f manifold; hp_acbn0_pao_3g has no &
           &f-shell fit.', nt )
     ENDIF
  ENDDO
  !
  ! nspin: the kernels hard-code spin indices 1 and 2 in the A0/HF accumulators.
  !
  IF ( nspin /= 2 ) CALL errore( 'hp_acbn0_check', &
     'ACBN0 in hp.x is implemented for collinear spin-polarised only &
     &(nspin = 2).', 1 )
  !
  IF ( gamma_only ) CALL errore( 'hp_acbn0_check', &
     'gamma_only is not supported by the ACBN0 kernels', 1 )
  !
  ! Report the manifold that was actually picked up, so the run is auditable
  ! from stdout alone.
  !
  nath = 0
  WRITE( stdout, '(5x,"Hubbard manifold taken from MODULE ldaU:")')
  DO nt = 1, ntyp
     IF ( .NOT.is_hubbard(nt) ) CYCLE
     nath = nath + 1
     WRITE( stdout, '(7x,"type ",i2,"  ",a6,"  l = ",i1)', ADVANCE='NO') &
           nt, TRIM(atm(nt)), Hubbard_l(nt)
     IF ( is_hubbard_back(nt) ) THEN
        WRITE( stdout, '(" + ",i1)', ADVANCE='NO') Hubbard_l2(nt)
        IF ( backall(nt) ) WRITE( stdout, '(" + ",i1)', ADVANCE='NO') Hubbard_l3(nt)
     ENDIF
     WRITE( stdout, * )
  ENDDO
  WRITE( stdout, '(7x,"nwfcU = ",i5,/)') nwfcU
  !
  IF ( nath == 0 ) CALL errore( 'hp_acbn0_check', &
     'lda_plus_u is set but no atomic type is Hubbard-active.', 1 )
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_check

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_write( ns, hub_u )
  !---------------------------------------------------------------------
  !
  ! Write the ACBN0 Hubbard parameters to eACBN0_hub_dev.dat (and to stdout).
  !
  ! Note on the "scalar U" column: hub_u(na,m1,m2) is an (m1,m2) matrix, but
  ! the kernel builds it as
  !     hub_u(na,m1,m2) = A0%U(l,l,na)/HF%U(l,l,na) - J-term
  ! i.e. it depends on m1,m2 ONLY through their l channel.  With one channel
  ! per type - which hp_acbn0_check now enforces - every element is therefore
  ! the same number, and that number is the scalar U to compare against the
  ! DFPT value in prefix.Hubbard_parameters.dat.  The full matrix is printed
  ! anyway, so that the claim is checkable rather than asserted.
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : rytoev
  USE io_global,  ONLY : stdout, ionode
  USE io_files,   ONLY : prefix, tmp_dir
  USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp
  USE lsda_mod,   ONLY : nspin
  USE ldaU,       ONLY : Hubbard_l, Hubbard_lmax, is_hubbard
  USE ldaU_hp,    ONLY : acbn0
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch, hp_acbn0_manlab, hp_acbn0_chanlab
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: ns   (lm_max,lm_max,nspin,nat)
  REAL(DP), INTENT(IN) :: hub_u(nat,lm_max,lm_max)
  !
  INTEGER              :: na, nt, m1, m2, is, iunit
  REAL(DP)             :: utr, occ, e_on
  CHARACTER(LEN=1)     :: spdf(4) = ['s','p','d','f']
  CHARACTER(LEN=256)   :: filename
  !
  INTEGER, EXTERNAL    :: find_free_unit
  !
  IF ( .NOT.ionode ) RETURN
  !
  ! Written into outdir (tmp_dir), next to prefix.save, not the launch dir.
  !
  filename = TRIM(tmp_dir) // 'eACBN0_hub_dev.dat'
  iunit = find_free_unit()
  OPEN( iunit, file = TRIM(filename), form = 'formatted', status = 'unknown' )
  !
  WRITE( iunit, '(/2x,79("="))')
  WRITE( iunit, '(/25x,"Hubbard parameters from ACBN0")')
  WRITE( iunit, '(/22x,"computed by hp.x  (acbn0 = ",i1,")")') acbn0
  WRITE( iunit, '(/2x,79("="),/)')
  WRITE( iunit, '(2x,"prefix     : ",a)') TRIM(prefix)
  WRITE( iunit, '(2x,"projectors : ortho-atomic")')
  WRITE( iunit, '(2x,"bookkeeping: MODULE ldaU (official QE)")')
  WRITE( iunit, '(2x,"mixing     : 1.000  (one-shot; no SCF history)")')
  WRITE( iunit, '(2x,"estimator  : pseudo-hybrid ACBN0 ratio, NOT linear response")')
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"Summary")')
  WRITE( iunit, '(2x,79("-"))')
  WRITE( iunit, '(/4x,"site n.  type  label  channel   occupation   U_ACBN0 (eV)")')
  WRITE( iunit, '( 4x,"(one row per channel; occupation and U of that channel only)")')
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     !
     ! trace of the occupation matrix, summed over spin
     !
     occ = 0.0_DP
     DO is = 1, nspin
        DO m1 = 1, lm(nt)
           occ = occ + ns(m1,m1,is,na)
        ENDDO
     ENDDO
     !
     ! One row per channel.  U is constant inside an (l1,l2) block (it is
     ! built per l, see the header note), so the diagonal block's value is
     ! that channel's U.  With several channels the (1,1) element alone
     ! would show only the first channel - or 0 when its on-site U is masked
     ! by acbn0_rm_ehub_u.  Cross blocks are in the full matrix below.
     !
     DO m1 = 1, lm(nt)
        IF ( ANY( lch(1:m1-1,nt) == lch(m1,nt) ) ) CYCLE     ! first m of this channel
        occ = 0.0_DP
        DO is = 1, nspin
           DO m2 = 1, lm(nt)
              IF ( lch(m2,nt) == lch(m1,nt) ) occ = occ + ns(m2,m2,is,na)
           ENDDO
        ENDDO
        utr = hub_u(na,m1,m1) * rytoev
        WRITE( iunit, '(6x,i3,4x,i3,3x,a6,4x,a2,6x,f9.4,4x,f10.4)') &
              na, nt, TRIM(atm(nt)), hp_acbn0_chanlab(nt,lch(m1,nt)), occ, utr
     ENDDO
  ENDDO
  !
  CALL hp_acbn0_energy_u( ns, hub_u, e_on )
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"One-shot Hubbard energy  (NOT self-consistent)")')
  WRITE( iunit, '(2x,79("-"))')
  WRITE( iunit, '(/4x,"E(+U) : ",f16.9," Ry  = ",f12.6," eV")') e_on, e_on*rytoev
  WRITE( iunit, '( 4x,"Evaluated with the one-shot U above and the occupations it was computed")')
  WRITE( iunit, '( 4x,"from (v_hubbard_a0.f90 loop); not a total-energy contribution.")')
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"Full U matrix  U[na, m1, m2]  (eV)")')
  WRITE( iunit, '(2x,79("-"))')
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"   (",i2," x",i2,")")') &
           na, TRIM(atm(nt)), lm(nt), lm(nt)
     DO m1 = 1, lm(nt)
        WRITE( iunit, '(6x,16f10.4)') (hub_u(na,m1,m2)*rytoev, m2 = 1, lm(nt))
     ENDDO
  ENDDO
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"Occupation matrix  ns[na, m1, m2, spin]")')
  WRITE( iunit, '(2x,79("-"))')
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     DO is = 1, nspin
        WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"  spin ",i1)') na, TRIM(atm(nt)), is
        DO m1 = 1, lm(nt)
           WRITE( iunit, '(6x,16f10.4)') (ns(m1,m2,is,na), m2 = 1, lm(nt))
        ENDDO
     ENDDO
  ENDDO
  !
  WRITE( iunit, '(/2x,79("="),/)')
  CLOSE( iunit )
  !
  WRITE( stdout, '(/5x,"ACBN0 Hubbard parameters written to ",a,/)') TRIM(filename)
  !
  ! Echo the summary to stdout as well, so it is greppable from the job output.
  !
  WRITE( stdout, '(5x,"site n.  type  label  channel   U_ACBN0 (eV)")')
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     DO m1 = 1, lm(nt)
        IF ( ANY( lch(1:m1-1,nt) == lch(m1,nt) ) ) CYCLE
        WRITE( stdout, '(7x,i3,4x,i3,3x,a6,4x,a2,6x,f10.4)') &
              na, nt, TRIM(atm(nt)), hp_acbn0_chanlab(nt,lch(m1,nt)), hub_u(na,m1,m1)*rytoev
     ENDDO
  ENDDO
  WRITE( stdout, * )
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_write
SUBROUTINE hp_acbn0_ns_a0(ns, hub_u)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE cell_base,            ONLY : at, alat
  USE klist,                ONLY : nks, ngk, xk
  USE ldaU,                 ONLY : Hubbard_l, Hubbard_lmax, is_hubbard, &
                                   wfcU, nwfcU
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, wg
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch, ofs
  USE ldaU_hp,              ONLY : acbn0_rm_ehub_u

  IMPLICIT NONE

  REAL(DP),    INTENT(OUT) :: ns   (lm_max,lm_max,nspin,nat)
  REAL(DP),    INTENT(OUT) :: hub_u(nat,lm_max,lm_max)

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

  REAL(DP),    ALLOCATABLE :: hub_u_(:,:,:)

  REAL(DP), EXTERNAL       :: hp_acbn0_compute_eri
  REAL(DP)                 :: c1(6), c2(6), c3(6), c4(6)

  INTEGER                  :: i, ik, ibnd, is, npw
  INTEGER                  :: na, nt, l, m
  INTEGER                  :: na1, na2
  INTEGER                  :: nt1, nt2
  INTEGER                  :: m1, m2, m3, m4
  INTEGER                  :: l1, l2, l3, l4
  !
  !

  CALL start_clock('acbn0_ns_a0')
  !
  !
  !
  CALL allocate_bec_type(nwfcU, nbnd, proj)
  !
  ALLOCATE ( n%O(lm_max,lm_max,nspin,nat) ) ; n%O = 0.0_DP
  ALLOCATE ( p%O(lm_max,lm_max,nspin,nat) ) ; p%O = 0.0_DP
  ALLOCATE ( A0%U(4,4,nat), A0%J(4,4,nat) )
  ALLOCATE ( HF%U(4,4,nat), HF%J(4,4,nat) )
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcU | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( nks > 1 ) CALL get_buffer(wfcU, nwordwfcU, iunhub, ik)
     CALL calbec(npw, wfcU, evc, proj)
     !
     ! calculate the occupation matrices n
     !
     DO na = 1, nat
        nt = ityp(na)
        !
        IF ( lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt)
        DO m2 = 1, lm(nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n%O(m1,m2,current_spin,na) = &
                 n%O(m1,m2,current_spin,na) + &
                                wg(ibnd,ik) * &
                  proj%r(ofs(m2,na),ibnd) * &
                  proj%r(ofs(m1,na),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 n%O(m1,m2,current_spin,na) = &
                 n%O(m1,m2,current_spin,na) + &
                                wg(ibnd,ik) * &
             DBLE(proj%k(ofs(m2,na),ibnd) * &
            CONJG(proj%k(ofs(m1,na),ibnd)))
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
        IF ( lm(nt) == 0 ) CYCLE
        !
        DO m = 1, lm(nt)
           !
           l = lch(m,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt) + &
                         proj%r(ofs(m,na),ibnd) * &
                         proj%r(ofs(m,na),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt) + &
                    DBLE(proj%k(ofs(m,na),ibnd) * &
                   CONJG(proj%k(ofs(m,na),ibnd)))
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
        IF ( lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt)
           !
           l1 = lch(m1,nt)
           !
        DO m2 = 1, lm(nt)
           !
           l2 = lch(m2,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 p%O(m1,m2,current_spin,na) = &
                 p%O(m1,m2,current_spin,na) + &
                         num(ibnd)%O(l1,l2,   &
                           current_spin,nt) * &
                                wg(ibnd,ik) * &
                  proj%r(ofs(m2,na),ibnd) * &
                  proj%r(ofs(m1,na),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 p%O(m1,m2,current_spin,na) = &
                 p%O(m1,m2,current_spin,na) + &
                         num(ibnd)%O(l1,l2,   &
                           current_spin,nt) * &
                                wg(ibnd,ik) * &
             DBLE(proj%k(ofs(m2,na),ibnd) * &
            CONJG(proj%k(ofs(m1,na),ibnd)))
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
  CALL hp_acbn0_symm_n4(n)
  !
  ns = 0.0_DP
  !
  ns = n%O
  !
  DEALLOCATE (n%O)
  !
  IF ( nspin == 1 ) p%O = 0.5_DP * p%O
  !
  CALL hp_acbn0_symm_n4(p)
  !
  ! calculate ACBN0 U & J
  !
  A0%U = 0.0_DP ; A0%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, lm(nt)
        CALL hp_acbn0_orberi(lch(m1,nt), l1)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m1,nt), c1)
        !
     l2 = 0
     DO m2 = 1, lm(nt)
        CALL hp_acbn0_orberi(lch(m2,nt), l2)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m2,nt), c2)
        !
     l3 = 0
     DO m3 = 1, lm(nt)
        IF ( lch(m3,nt) /= lch(m1,nt) ) CYCLE
        CALL hp_acbn0_orberi(lch(m3,nt), l3)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m3,nt), c3)
        !
     l4 = 0
     DO m4 = 1, lm(nt)
        IF ( lch(m4,nt) /= lch(m2,nt) ) CYCLE
        CALL hp_acbn0_orberi(lch(m4,nt), l4)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m4,nt), c4)
        !
        IF ( lch(m1,nt) == 1 .AND. lch(m2,nt) == 1 ) THEN
              A0%U(lch(m1,nt),lch(m2,nt),na)  = &
              A0%U(lch(m1,nt),lch(m2,nt),na)  + &
                  (p%O(m3,m1,1,na)*p%O(m4,m2,2,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,1,na)) * &
                          hp_acbn0_compute_eri(l3,l1,l4,l2,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                                      c3,c1,c4,c2)
        ELSE
              A0%U(lch(m1,nt),lch(m2,nt),na)  = &
              A0%U(lch(m1,nt),lch(m2,nt),na)  + &
                  (p%O(m3,m1,1,na)*p%O(m4,m2,1,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,2,na)  + &
                   p%O(m3,m1,1,na)*p%O(m4,m2,2,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,1,na)) * &
                          hp_acbn0_compute_eri(l3,l1,l4,l2,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                    tau(:,na)*alat,tau(:,na)*alat,    &
                                      c3,c1,c4,c2)
              A0%J(lch(m1,nt),lch(m2,nt),na)  = &
              A0%J(lch(m1,nt),lch(m2,nt),na)  + &
                  (p%O(m3,m1,1,na)*p%O(m4,m2,1,na)  + &
                   p%O(m3,m1,2,na)*p%O(m4,m2,2,na)) * &
                          hp_acbn0_compute_eri(l3,l2,l4,l1,    &
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
     IF ( lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt)
     DO m2 = 1, lm(nt)
        !
        IF ( m1 == m2 ) THEN
             HF%U(lch(m1,nt),lch(m2,nt),na) = &
             HF%U(lch(m1,nt),lch(m2,nt),na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,2,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,1,na)
        ELSE
             HF%U(lch(m1,nt),lch(m2,nt),na) = &
             HF%U(lch(m1,nt),lch(m2,nt),na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,1,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,2,na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,2,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,1,na)
 
             HF%J(lch(m1,nt),lch(m2,nt),na) = &
             HF%J(lch(m1,nt),lch(m2,nt),na) + &
                    ns(m1,m1,1,na)*ns(m2,m2,1,na) + &
                    ns(m1,m1,2,na)*ns(m2,m2,2,na)
        ENDIF
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  ALLOCATE (hub_u_(nat,lm_max,lm_max)) ; hub_u_ = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt)
     DO m2 = 1, lm(nt)
        !
        IF ( lch(m1,nt) == 1 .AND. lch(m2,nt) == 1 ) THEN
           !
           ! both indices in an s channel: no exchange term, ACBN0 is the
           ! direct ratio alone.  With background channels m1 and m2 can sit
           ! in DIFFERENT channels, so both must be tested - collapsing this
           ! to a single test was only valid while a type had one channel.
           !
           hub_u_(na,m1,m2) = A0%U(lch(m1,nt),lch(m2,nt),na) /  &
                              HF%U(lch(m1,nt),lch(m2,nt),na)
        ELSE
           hub_u_(na,m1,m2) = A0%U(lch(m1,nt),lch(m2,nt),na) /  &
                              HF%U(lch(m1,nt),lch(m2,nt),na) -  &
                              A0%J(lch(m1,nt),lch(m2,nt),na) /  &
                              HF%J(lch(m1,nt),lch(m2,nt),na)
        ENDIF
        !
        ! One-shot (no mixing).  acbn0_rm_ehub_u is Lee & Son's
        ! remove_ehub_u with the flag inverted (1 = remove), applied exactly
        ! as new_ns_v.f90:1362 does - a product over both channel indices,
        ! on the on-site block only, AFTER the ratio is formed.  The channel
        ! stays in the manifold, so the N^IJ weights and the ERI sums above
        ! are untouched.
        !
        hub_u(na,m1,m2) = hub_u_(na,m1,m2)                              &
                        * DBLE( 1 - acbn0_rm_ehub_u(nt,lch(m1,nt)) )    &
                        * DBLE( 1 - acbn0_rm_ehub_u(nt,lch(m2,nt)) )
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( lm(nt1) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt1)
     DO m2 = m1, lm(nt1)
        !
        hub_u(na1,m1,m2) = 0.5D0 * (hub_u(na1,m1,m2) + &
                                     hub_u(na1,m2,m1))
        !
        hub_u(na1,m2,m1) = hub_u(na1,m1,m2)
        !
     ENDDO
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE ( hub_u_ )
  DEALLOCATE ( A0%U, A0%J )
  DEALLOCATE ( HF%U, HF%J )

  CALL stop_clock('acbn0_ns_a0')

  RETURN

END SUBROUTINE hp_acbn0_ns_a0
!
SUBROUTINE hp_acbn0_symm_n4( n )

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE ldaU,                 ONLY : Hubbard_l, Hubbard_lmax, is_hubbard
  USE lsda_mod,             ONLY : nspin
  USE symm_base,            ONLY : d1, d2, d3, nsym, irt
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch

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

  !

  ALLOCATE (n%T(lm_max,lm_max,nspin,nat)) ; n%T = 0.0_DP

  DO nt = 1, ntyp
     !
     ! BUGFIX vs PW/src/new_ns_v.f90:2993 (Finding #1 of logs/new_ns_v_review.md).
     ! Without this guard a type carrying no Hubbard manifold has lch(nt)=0,
     ! and the next line writes l_separation(0,nt) - one element BEFORE the array.
     ! Dormant for NiO (all three types are Hubbard-active); fires the moment a
     ! spectator species is added.
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     temp = lch(1,nt)
     l_separation(lch(1,nt),nt) = 1
     !
     DO m = 2, lm(nt)
        !
        IF ( lch(m,nt) /= temp ) THEN
           temp = lch(m,nt)
           l_separation(lch(m,nt),nt) = m
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
        IF ( lm(nt1) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
           m1_ = m1 - l_separation(lch(m1,nt1),nt1) + 1
           !
        DO m2 = 1, lm(nt1)
           m2_ = m2 - l_separation(lch(m2,nt1),nt1) + 1
           !
           DO isym = 1, nsym
              !
              na2 = irt(isym, na1)
              !
              DO m3 = 1, lm(nt1)
                 IF ( lch(m1,nt1) /= lch(m3,nt1) ) CYCLE
                 m3_ = m3 - l_separation(lch(m3,nt1),nt1) + 1
                 !
              DO m4 = 1, lm(nt1)
                 IF ( lch(m2,nt1) /= lch(m4,nt1) ) CYCLE
                 m4_ = m4 - l_separation(lch(m4,nt1),nt1) + 1
                 !
                 IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)
  
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d1(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)

                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 4 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d2(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d3(m4_,m2_,isym) / &
                                                   (nsym)

                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 1 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 2 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d1(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 3 ) THEN
                    n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) + &
                                         d3(m3_,m1_,isym) * &
                                        n%O(m3,m4,is,na2) * &
                                         d2(m4_,m2_,isym) / &
                                                   (nsym)
                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 4 ) THEN
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
     IF ( lm(nt1) == 0 ) CYCLE
     !
     DO is = 1, nspin
        !
        DO m1 = 1, lm(nt1)
        DO m2 = m1, lm(nt1)
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

END SUBROUTINE hp_acbn0_symm_n4
!
SUBROUTINE hp_acbn0_orberi(lp1, eri_l)

  ! lp1 = 1: s orbital [eri_l = 1(s)]
  ! lp1 = 2: p orbital [eri_l = 2(pz), 3(px), 4(py)]
  ! lp1 = 3: d orbital [eri_l = 5(dz^2-r^2), 6(dxz), 7(dyz), 8(dx^2-y^2),
  ! 9(dxy)]

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: lp1
  INTEGER, INTENT(INOUT) :: eri_l

  IF (lp1 == 1) THEN
     eri_l = 1
  ELSEIF (lp1 == 2) THEN
     IF (eri_l == 2) THEN
        eri_l = 3
     ELSEIF (eri_l == 3) THEN
        eri_l = 4
     ELSE
        eri_l = 2
     ENDIF
  ELSEIF (lp1 == 3) THEN
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

END SUBROUTINE hp_acbn0_orberi
!-----------------------------------------------------------------------
!  Dispatcher
!-----------------------------------------------------------------------
!
!    acbn0 = 1   on-site U        -> hp_acbn0_u   (S3)
!    acbn0 = 2   on-site U + V    -> hp_acbn0_v   (V3)
!
!  Both branches now read the manifold and the projectors from official QE's
!  MODULE ldaU, so both run on an ordinary ortho-atomic DFT+U ground state.
!  They remain separate because the occupation matrix has a different rank:
!  ns(m1,m2,spin,na) on-site versus ns(m1,m2,spin,na1,na2,nnc) inter-site.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0()
  !---------------------------------------------------------------------
  !
  USE ldaU_hp, ONLY : acbn0
  !
  IMPLICIT NONE
  !
  IF ( acbn0 == 1 ) THEN
     CALL hp_acbn0_u()
  ELSEIF ( acbn0 == 2 ) THEN
     CALL hp_acbn0_v()
  ELSE
     CALL errore( 'hp_acbn0', 'acbn0 must be 1 or 2 here', ABS(acbn0) )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE hp_acbn0

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_v()
  !---------------------------------------------------------------------
  !
  ! STAGE V3 - on-site U AND inter-site V, on ldaU bookkeeping.
  !
  ! The manifold and the projectors come from official QE's MODULE ldaU exactly
  ! as the +U branch has since S3; the neighbour table is built here, by
  ! MODULE hp_acbn0_nn, from geometry and a cutoff.  MODULE ldaV is no longer
  ! read by any kernel.
  !
  ! The consequence is the one that matters: acbn0 = 2 now runs on an ordinary
  ! ortho-atomic DFT+U ground state - the same one hp.x's DFPT branch consumes
  ! - so U_ACBN0, V_ACBN0 and U_DFPT all come from a single pw.x run.
  !
  ! The cutoff d_IJ (Lee & Son, PRR 2, 043410 (2020)) is taken from the
  ! &INPUTHP tag nn_distance, and from nowhere else.  Earlier stages fell back
  ! to ehub_nn_distance in the XML, which exists only for a ground state made
  ! by the Lee-Son pw.x; that fallback is gone, because this file must build
  ! and run against an unpatched QE 7.3 where MODULE ldaV does not exist.
  ! hp_readin refuses acbn0 = 2 without nn_distance.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunhub, nwordwfcU, iunwfc, nwordwfc
  USE buffers,          ONLY : open_buffer, close_buffer
  USE ions_base,        ONLY : nat
  USE lsda_mod,         ONLY : nspin
  USE wvfct,            ONLY : npwx, nbnd
  USE control_flags,    ONLY : io_level
  USE noncollin_module, ONLY : npol
  USE symm_base,        ONLY : d1, d2, d3
  USE ldaU,             ONLY : Hubbard_lmax, nwfcU
  USE ldaU_hp,          ONLY : nn_distance
  USE hp_acbn0_nn,      ONLY : hp_acbn0_build_nn, hp_acbn0_free_nn, nnnc
  USE hp_acbn0_man,         ONLY : lm_max, hp_acbn0_build_man, &
                                   hp_acbn0_free_man
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: ns(:,:,:,:,:,:)   ! ns(m1,m2,spin,na1,na2,nnc)
  REAL(DP), ALLOCATABLE :: hub_v(:,:,:,:,:)  ! hub_v(na1,na2,m1,m2,nnc), Ry
  !
  REAL(DP)              :: dist
  LOGICAL               :: exst
  !
  CALL start_clock('hp_acbn0')
  !
  WRITE( stdout, '(/5x,70("="))')
  WRITE( stdout, '(5x,"eACBN0 Hubbard parameters (acbn0 = 2): on-site U + intersite V")')
  WRITE( stdout, '(5x,"pseudo-hybrid estimator - NOT linear response")')
  WRITE( stdout, '(5x,"stage V3: ldaU bookkeeping, own kernels, own neighbour table")')
  WRITE( stdout, '(5x,70("="),/)')
  !
  CALL hp_acbn0_check_v()
  !
  ! Build the manifold table (lm / lch / ofs).  This is what makes the
  ! kernels multi-channel: see MODULE hp_acbn0_man.
  !
  CALL hp_acbn0_build_man()
  !
  !
  ! Inter-site cutoff.  hp_readin has already refused the case where neither
  ! source is available.
  !
  IF ( nn_distance < 0.0_DP ) CALL errore('hp_acbn0_v', &
     'acbn0 = 2 needs the inter-site cutoff: set nn_distance (Angstrom) in &INPUTHP', 1)
  !
  dist = nn_distance
  WRITE( stdout, '(5x,"inter-site cutoff from &INPUTHP: nn_distance = ",f6.3," A")') dist
  !
  CALL hp_acbn0_build_nn( dist )
  !
  ! Wigner d-matrices: read_file builds them only under okpaw
  ! (read_file_new.f90:348-351) and symm_n4/_n4m/_n5 all need them.
  !
  CALL d_matrix( d1, d2, d3 )
  !
  ! Bug 1: read_file closes iunwfc at read_file_new.f90:62 with
  ! close_buffer(...,'KEEP').  'KEEP' preserves the FILES but the call also
  ! deregisters the unit, so get_buffer falls through to davcio on a closed
  ! unit.  Same two lines read_file itself uses.
  !
  nwordwfc = nbnd * npwx * npol
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )
  IF ( io_level /= 0 .AND. .NOT.exst ) CALL errore( 'hp_acbn0_v', &
     'the KS wavefunction files (prefix.wfc*) are missing from outdir - &
     &rerun pw.x, or set disk_io to at least "low"', 1 )
  !
  nwordwfcU = npwx * nwfcU * npol
  CALL open_buffer( iunhub, 'hub', nwordwfcU, io_level, exst )
  !
  ! wfcU = O^{-1/2} S phi on the Hubbard manifold, for every k point.
  !
  CALL orthoUwfc( .FALSE. )
  !
  ALLOCATE( ns(lm_max,lm_max,nspin,nat,nat,nnnc) )   ; ns(:,:,:,:,:,:)  = 0.0_DP
  ALLOCATE( hub_v(nat,nat,lm_max,lm_max,nnnc) )      ; hub_v(:,:,:,:,:) = 0.0_DP
  !
  CALL hp_acbn0_ns_e0_1( ns, hub_v )
  !
  ! pw.x-format tables first, so V0/V1/V2/V3 stay comparable by diff, then
  ! the readable file.
  !
  CALL hp_acbn0_print_uv( hub_v )
  CALL hp_acbn0_write_v( ns, hub_v )
  CALL hp_acbn0_dump_v( ns, hub_v )
  !
  DEALLOCATE( ns, hub_v )
  CALL hp_acbn0_free_nn()
  CALL hp_acbn0_free_man()
  !
  CALL close_buffer( iunhub, 'delete' )   ! ours, scratch
  CALL close_buffer( iunwfc, 'KEEP'   )   ! the user's wavefunctions - NEVER delete
  !
  CALL stop_clock('hp_acbn0')
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_v

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_check_v()
  !---------------------------------------------------------------------
  !
  ! Preconditions of the +V branch.  Identical in spirit to hp_acbn0_check:
  ! each of these fails as a wrong number rather than a crash if not caught.
  !
  USE io_global,        ONLY : stdout
  USE ions_base,        ONLY : ntyp => nsp, atm
  USE lsda_mod,         ONLY : nspin
  USE control_flags,    ONLY : gamma_only
  USE ldaU,             ONLY : lda_plus_u, Hubbard_projectors, Hubbard_lmax, &
                               Hubbard_l, is_hubbard, is_hubbard_back, nwfcU, &
                               Hubbard_l2, Hubbard_l3, backall
  !
  IMPLICIT NONE
  !
  INTEGER :: nt, nath
  !
  IF ( .NOT.lda_plus_u ) CALL errore( 'hp_acbn0_check_v', &
     'acbn0 > 0 needs a Hubbard manifold, but none was selected: set &
     &acbn0_l_choice(nt,b) = 1 in &INPUTHP (the pw.x ground state itself must &
     &be a plain non-Hubbard SCF).', 1 )
  !
  ! ortho-atomic is a requirement, not a preference: the ACBN0 ratio is
  ! derived for an orthonormal projector basis.  orthoUwfc honours the input
  ! value rather than forcing it, so the contract is enforced here.
  !
  IF ( TRIM(Hubbard_projectors) /= 'ortho-atomic' ) CALL errore( 'hp_acbn0_check_v', &
     'eACBN0 requires ortho-atomic Hubbard projectors; rerun pw.x with &
     &HUBBARD (ortho-atomic).', 1 )
  !
  IF ( Hubbard_lmax < 0 .OR. nwfcU <= 0 ) CALL errore( 'hp_acbn0_check_v', &
     'no Hubbard manifold was found (Hubbard_lmax < 0 or nwfcU = 0).', 1 )
  !
  ! Multi-channel types are supported from stage V4 on: MODULE hp_acbn0_man
  ! rebuilds the concatenated manifold ldaV used out of ldaU's ll(m,nt) plus
  ! the offsetU / offsetU_back / offsetU_back1 blocks.  What must still hold
  ! is that ldaU actually laid the background out, which it only does when
  ! the channel is ACTIVE - i.e. Hubbard_U2 /= 0 (ldaU.f90:249).
  !
  DO nt = 1, ntyp
     IF ( is_hubbard_back(nt) .AND. Hubbard_l2(nt) < 0 ) &
        CALL errore( 'hp_acbn0_check_v', &
        'a background Hubbard channel is active (Hubbard_U2 /= 0) but &
        &Hubbard_l2 was never set for this type.', nt )
     IF ( backall(nt) .AND. Hubbard_l3(nt) < 0 ) &
        CALL errore( 'hp_acbn0_check_v', &
        'backall is set for this type but Hubbard_l3 was never set.', nt )
  ENDDO
  !
  ! f manifolds: hp_acbn0_pao_3g carries 3-Gaussian PAO fits for s, p and d only and
  ! returns c = 0 without complaining for anything else, which would make
  ! every eri vanish and U come out as 0/0.
  !
  DO nt = 1, ntyp
     IF ( is_hubbard(nt) .AND. Hubbard_l(nt) > 2 ) CALL errore( 'hp_acbn0_check_v', &
        'eACBN0 is implemented for s, p and d manifolds only: hp_acbn0_pao_3g has no &
        &f-shell Gaussian fit and returns zero coefficients silently.', nt )
     IF ( is_hubbard_back(nt) ) THEN
        IF ( Hubbard_l2(nt) > 2 ) CALL errore( 'hp_acbn0_check_v', &
           'background channel Hubbard_l2 is an f manifold; hp_acbn0_pao_3g has no &
           &f-shell fit.', nt )
        IF ( backall(nt) .AND. Hubbard_l3(nt) > 2 ) CALL errore( 'hp_acbn0_check_v', &
           'background channel Hubbard_l3 is an f manifold; hp_acbn0_pao_3g has no &
           &f-shell fit.', nt )
     ENDIF
  ENDDO
  !
  ! nspin: the E0/HF assembly hard-codes spin slots 1 and 2 while the
  ! occupation arrays are dimensioned (...,nspin,...), so nspin = 1 reads past
  ! the end of the array, silently.
  !
  IF ( nspin /= 2 ) CALL errore( 'hp_acbn0_check_v', &
     'eACBN0 in hp.x is implemented for collinear spin-polarised only &
     &(nspin = 2).  For a non-magnetic system use nspin = 2 with &
     &starting_magnetization = 0.0.', 1 )
  !
  IF ( gamma_only ) CALL errore( 'hp_acbn0_check_v', &
     'gamma_only is not supported by the eACBN0 kernels', 1 )
  !
  nath = 0
  WRITE( stdout, '(5x,"Hubbard manifold taken from MODULE ldaU:")')
  DO nt = 1, ntyp
     IF ( .NOT.is_hubbard(nt) ) CYCLE
     nath = nath + 1
     WRITE( stdout, '(7x,"type ",i2,"  ",a6,"  l = ",i1)', ADVANCE='NO') &
           nt, TRIM(atm(nt)), Hubbard_l(nt)
     IF ( is_hubbard_back(nt) ) THEN
        WRITE( stdout, '(" + ",i1)', ADVANCE='NO') Hubbard_l2(nt)
        IF ( backall(nt) ) WRITE( stdout, '(" + ",i1)', ADVANCE='NO') Hubbard_l3(nt)
     ENDIF
     WRITE( stdout, * )
  ENDDO
  WRITE( stdout, '(7x,"nwfcU = ",i5,/)') nwfcU
  !
  IF ( nath == 0 ) CALL errore( 'hp_acbn0_check_v', &
     'lda_plus_u is set but no atomic type is Hubbard-active.', 1 )
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_check_v

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_energy_v( ns, hub_v, e_on, e_int )
  !---------------------------------------------------------------------
  !
  ! One-shot extended Hubbard energy, evaluated with the U and V just
  ! computed and the occupations they were computed from.  Loop for loop
  ! the same as PW/src/v_hubbard_v.f90 (Lee & Son), minus the potential:
  !
  !   E_on  = sum_{I,s,m}  U_mm/2 n_mm - sum_{I,s,m1,m2} U_m1m2/2 n_m1m2 n_m2m1
  !   E_int = - sum_{I/=J,s,m1,m2} V_m1m2/2 n^IJ_m1m2 n^JI_m2m1
  !
  ! NOT a physical energy: U, V and the density are not self-consistent
  ! here (no feedback loop).  It is the number v_hubbard_v would print as
  ! E(+U+V) on the first iteration of a run started from this density with
  ! ehub_mixing = 1, and exists to validate the functional port.
  !
  USE kinds,       ONLY : DP
  USE ions_base,   ONLY : nat, ityp
  USE lsda_mod,    ONLY : nspin
  USE hp_acbn0_nn, ONLY : nnnc, nn_vect, nn_dim, nn_csize
  USE hp_acbn0_man, ONLY : lm_max, lm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns   (lm_max,lm_max,nspin,nat,nat,nnnc)
  REAL(DP), INTENT(IN)  :: hub_v(nat,nat,lm_max,lm_max,nnnc)
  REAL(DP), INTENT(OUT) :: e_on, e_int
  !
  INTEGER :: is, na1, na2, nt1, nt2, m1, m2, n1, n2, n3, nnc, ics
  !
  e_on = 0.0_DP ; e_int = 0.0_DP
  ics  = nn_csize + 1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( lm(nt1) == 0 ) CYCLE
     DO is = 1, nspin
        DO m1 = 1, lm(nt1)
           e_on = e_on + 0.5_DP*hub_v(na1,na1,m1,m1,1)*ns(m1,m1,is,na1,na1,1)
           DO m2 = 1, lm(nt1)
              e_on = e_on - 0.5_DP*hub_v(na1,na1,m1,m2,1) &
                          * ns(m1,m2,is,na1,na1,1)*ns(m2,m1,is,na1,na1,1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( lm(nt1) == 0 ) CYCLE
     DO na2 = 1, nat
        nt2 = ityp(na2)
        IF ( lm(nt2) == 0 ) CYCLE
        IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        DO is = 1, nspin
           DO m1 = 1, lm(nt1)
              DO m2 = 1, lm(nt2)
                 e_int = e_int - 0.5_DP*hub_v(na1,na2,m1,m2,1) &
                               * ns(m1,m2,is,na1,na2,1)*ns(m2,m1,is,na2,na1,1)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  nnc = 2
  DO n1 = -nn_dim(1), nn_dim(1)
  DO n2 = -nn_dim(2), nn_dim(2)
  DO n3 = -nn_dim(3), nn_dim(3)
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     DO na1 = 1, nat
        nt1 = ityp(na1)
        IF ( lm(nt1) == 0 ) CYCLE
        DO na2 = 1, nat
           nt2 = ityp(na2)
           IF ( lm(nt2) == 0 ) CYCLE
           IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           DO is = 1, nspin
              DO m1 = 1, lm(nt1)
                 DO m2 = 1, lm(nt2)
                    ! NOTE: ns* = ns  (v_hubbard_v.f90)
                    e_int = e_int - 0.5_DP*hub_v(na1,na2,m1,m2,nnc) &
                                  * ns(m1,m2,is,na1,na2,nnc)*ns(m1,m2,is,na1,na2,nnc)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
  ENDDO
  ENDDO
  ENDDO
  !
  IF ( nspin == 1 ) THEN
     e_on = 2.0_DP*e_on ; e_int = 2.0_DP*e_int
  ENDIF
  !
  RETURN
END SUBROUTINE hp_acbn0_energy_v

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_energy_u( ns, hub_u, e_on )
  !---------------------------------------------------------------------
  !
  ! On-site-only counterpart of hp_acbn0_energy_v, mirroring v_hubbard_a0.
  !
  USE kinds,        ONLY : DP
  USE ions_base,    ONLY : nat, ityp
  USE lsda_mod,     ONLY : nspin
  USE hp_acbn0_man, ONLY : lm_max, lm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns   (lm_max,lm_max,nspin,nat)
  REAL(DP), INTENT(IN)  :: hub_u(nat,lm_max,lm_max)
  REAL(DP), INTENT(OUT) :: e_on
  !
  INTEGER :: is, na, nt, m1, m2
  !
  e_on = 0.0_DP
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     DO is = 1, nspin
        DO m1 = 1, lm(nt)
           e_on = e_on + 0.5_DP*hub_u(na,m1,m1)*ns(m1,m1,is,na)
           DO m2 = 1, lm(nt)
              e_on = e_on - 0.5_DP*hub_u(na,m1,m2)*ns(m1,m2,is,na)*ns(m2,m1,is,na)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  IF ( nspin == 1 ) e_on = 2.0_DP*e_on
  !
  RETURN
END SUBROUTINE hp_acbn0_energy_u

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_dump_v( ns, hub_v )
  !---------------------------------------------------------------------
  !
  ! Machine-comparable dumps in EXACTLY the formats PW/src/io_rho_xml.f90
  ! of the Lee & Son patch writes into prefix.save/, so the two codes can
  ! be diffed file against file:
  !
  !   ehub_uv.txt         list-directed dump of ehub_v(nat,nat,9,9,nnnc), Ry
  !                       (the forked code allocates the orbital dims as 9,
  !                        input.f90:2188 - padded here to match)
  !   ehub_uv0.txt        indexed:  na1 na2 m1 m2 n1 n2 n3 nnc  V(Ry)
  !   occupation_uv.txt   list-directed dump of ne0(lm,lm,nspin,nat,nat,nnnc)
  !   occupation_uv0.txt  indexed, nonzero entries only
  !
  ! Written into outdir/prefix.save/, where the forked code puts them.
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : ionode, stdout
  USE io_files,    ONLY : tmp_dir, prefix
  USE ions_base,   ONLY : nat, ityp
  USE lsda_mod,    ONLY : nspin
  USE hp_acbn0_nn, ONLY : nnnc, nn_vect, nn_dim, nn_csize
  USE hp_acbn0_man, ONLY : lm_max, lm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: ns   (lm_max,lm_max,nspin,nat,nat,nnnc)
  REAL(DP), INTENT(IN) :: hub_v(nat,nat,lm_max,lm_max,nnnc)
  !
  REAL(DP), ALLOCATABLE :: v9(:,:,:,:,:)
  INTEGER  :: iu, is, na1, na2, nt1, nt2, m1, m2, n1, n2, n3, nnc, ics
  CHARACTER(LEN=256) :: dirname
  INTEGER, EXTERNAL  :: find_free_unit
  !
  IF ( .NOT.ionode ) RETURN
  !
  dirname = TRIM(tmp_dir) // TRIM(prefix) // '.save/'
  ics = nn_csize + 1
  !
  ! ---- ehub_uv.txt --------------------------------------------------------
  ALLOCATE( v9(nat,nat,9,9,nnnc) ) ; v9 = 0.0_DP
  v9(:,:,1:lm_max,1:lm_max,:) = hub_v(:,:,1:lm_max,1:lm_max,:)
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'ehub_uv.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, * ) v9
  CLOSE( iu )
  DEALLOCATE( v9 )
  !
  ! ---- ehub_uv0.txt -------------------------------------------------------
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'ehub_uv0.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, '(a)' ) '     # V[atom i, atom j, #orb(i), #orb(j), l, m, n, nnc]'
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( lm(nt1) == 0 ) CYCLE
     DO na2 = 1, nat
        nt2 = ityp(na2)
        IF ( lm(nt2) == 0 ) CYCLE
        IF ( na1 /= na2 .AND. nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        DO m1 = 1, lm(nt1)
           DO m2 = 1, lm(nt2)
              WRITE( iu, '(8i4, 1f25.15)' ) na1, na2, m1, m2, 0, 0, 0, 1, hub_v(na1,na2,m1,m2,1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  nnc = 2
  DO n1 = -nn_dim(1), nn_dim(1)
  DO n2 = -nn_dim(2), nn_dim(2)
  DO n3 = -nn_dim(3), nn_dim(3)
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     DO na1 = 1, nat
        nt1 = ityp(na1)
        IF ( lm(nt1) == 0 ) CYCLE
        DO na2 = 1, nat
           nt2 = ityp(na2)
           IF ( lm(nt2) == 0 ) CYCLE
           IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           DO m1 = 1, lm(nt1)
              DO m2 = 1, lm(nt2)
                 WRITE( iu, '(8i4, 1f25.15)' ) na1, na2, m1, m2, n1, n2, n3, nnc, hub_v(na1,na2,m1,m2,nnc)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
  ENDDO
  ENDDO
  ENDDO
  CLOSE( iu )
  !
  ! ---- occupation_uv.txt --------------------------------------------------
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'occupation_uv.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, * ) ns
  CLOSE( iu )
  !
  ! ---- occupation_uv0.txt -------------------------------------------------
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'occupation_uv0.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, '(a)' ) '# ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc'
  WRITE( iu, '(A, I3, 5I4)' ) '#', lm_max, lm_max, nspin, nat, nat, nnnc
  DO m1 = 1, lm_max
     DO m2 = 1, lm_max
        DO is = 1, nspin
           DO na1 = 1, nat
              DO na2 = 1, nat
                 DO nnc = 1, nnnc
                    IF ( ns(m1,m2,is,na1,na2,nnc) /= 0.0_DP ) &
                       WRITE( iu, '(6i4, 1f25.15)' ) m1, m2, is, na1, na2, nnc, ns(m1,m2,is,na1,na2,nnc)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  CLOSE( iu )
  !
  WRITE( stdout, '(5x,"eACBN0 dumps (forked-code formats) written to ",a)') TRIM(dirname)
  WRITE( stdout, '(7x,"ehub_uv.txt  ehub_uv0.txt  occupation_uv.txt  occupation_uv0.txt")')
  !
  RETURN
END SUBROUTINE hp_acbn0_dump_v

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_dump_u( ns, hub_u )
  !---------------------------------------------------------------------
  !
  ! acbn0 = 1 counterpart of hp_acbn0_dump_v.  The forked code (acbn0_type
  ! = 1) dumps ehub_u(nat,9,9) into ehub_uv.txt, an eV matrix listing into
  ! ehub_u0.txt, and na0(lm,lm,nspin,nat) into occupation_uv.txt.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : rytoev
  USE io_global,   ONLY : ionode, stdout
  USE io_files,    ONLY : tmp_dir, prefix
  USE ions_base,   ONLY : nat, ityp
  USE lsda_mod,    ONLY : nspin
  USE hp_acbn0_man, ONLY : lm_max, lm, lch
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: ns   (lm_max,lm_max,nspin,nat)
  REAL(DP), INTENT(IN) :: hub_u(nat,lm_max,lm_max)
  !
  REAL(DP), ALLOCATABLE :: u9(:,:,:)
  INTEGER  :: iu, is, na, nt, m1, m2, lmx
  CHARACTER(LEN=256) :: dirname
  INTEGER, EXTERNAL  :: find_free_unit
  !
  IF ( .NOT.ionode ) RETURN
  !
  dirname = TRIM(tmp_dir) // TRIM(prefix) // '.save/'
  !
  ALLOCATE( u9(nat,9,9) ) ; u9 = 0.0_DP
  u9(:,1:lm_max,1:lm_max) = hub_u(:,1:lm_max,1:lm_max)
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'ehub_uv.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, * ) u9
  CLOSE( iu )
  DEALLOCATE( u9 )
  !
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'ehub_u0.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, '(a)' ) '                          '
  WRITE( iu, '(a)' ) '     ---------------------'
  WRITE( iu, '(a)' ) '     U[atom i, #orb, #orb]'
  WRITE( iu, '(a)' ) '     ---------------------'
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     lmx = 2*MAXVAL(lch(1:lm(nt),nt)) - 1
     WRITE( iu, '(3(a,i2),a)' ) '     [ ', na, ', ', lmx, ', ', lmx, '  ]'
     DO m1 = 1, lm(nt)
        WRITE( iu, '(4x,16f9.4)' ) (hub_u(na,m1,m2)*rytoev, m2 = 1, lm(nt))
     ENDDO
     WRITE( iu, * ) ''
  ENDDO
  CLOSE( iu )
  !
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'occupation_uv.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, * ) ns
  CLOSE( iu )
  !
  iu = find_free_unit()
  OPEN( iu, FILE = TRIM(dirname)//'occupation_uv0.txt', FORM = 'FORMATTED', STATUS = 'unknown' )
  WRITE( iu, '(a)' ) '# ehub_lm_max, ehub_lm_max, nspin, nat'
  WRITE( iu, '(A, I3, 3I4)' ) '#', lm_max, lm_max, nspin, nat
  DO m1 = 1, lm_max
     DO m2 = 1, lm_max
        DO is = 1, nspin
           DO na = 1, nat
              IF ( ns(m1,m2,is,na) /= 0.0_DP ) &
                 WRITE( iu, '(4i4, 1f25.15)' ) m1, m2, is, na, ns(m1,m2,is,na)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  CLOSE( iu )
  !
  WRITE( stdout, '(5x,"ACBN0 dumps (forked-code formats) written to ",a)') TRIM(dirname)
  WRITE( stdout, '(7x,"ehub_uv.txt  ehub_u0.txt  occupation_uv.txt  occupation_uv0.txt")')
  !
  RETURN
END SUBROUTINE hp_acbn0_dump_u

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_write_v( ns, hub_v )
  !---------------------------------------------------------------------
  !
  ! Write the eACBN0 parameters to eACBN0_hub_dev.dat.
  !
  ! The on-site U is the na1 = na2, cell [0 0 0] block of hub_v, not a
  ! separate array.  Howto2.pdf: "A = B and l = m = n = 0 >>> on-site U
  ! matrix; Others >>> inter-site V matrix".
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : rytoev, bohr_radius_angs
  USE cell_base,   ONLY : alat
  USE io_global,   ONLY : stdout, ionode
  USE io_files,    ONLY : prefix, tmp_dir
  USE ions_base,   ONLY : nat, ityp, atm, ntyp => nsp
  USE lsda_mod,    ONLY : nspin
  USE ldaU,        ONLY : Hubbard_l, Hubbard_lmax, is_hubbard
  USE hp_acbn0_nn, ONLY : nnnc, nn_distance, nn_vect, nn_leng, nn_dim, nn_csize
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch, hp_acbn0_manlab, hp_acbn0_chanlab
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: ns   (lm_max,lm_max,nspin,nat,nat,nnnc)
  REAL(DP), INTENT(IN) :: hub_v(nat,nat,lm_max,lm_max,nnnc)
  !
  INTEGER              :: na, na1, na2, nt, nt1, nt2, m1, m2, is, nnc, iunit
  INTEGER              :: npair, ics, n1, n2, n3, thisnnc
  REAL(DP)             :: utr, vtr, occ, e_on, e_int
  CHARACTER(LEN=1)     :: spdf(4) = ['s','p','d','f']
  CHARACTER(LEN=256)   :: filename
  !
  REAL(DP), PARAMETER  :: veps = 1.0D-8   ! Ry; below this a pair carries no V
  !
  INTEGER, EXTERNAL    :: find_free_unit
  !
  IF ( .NOT.ionode ) RETURN
  !
  !
  ics = nn_csize + 1
  !
  ! Written into outdir (tmp_dir), next to prefix.save, not the launch dir.
  !
  filename = TRIM(tmp_dir) // 'eACBN0_hub_dev.dat'
  iunit = find_free_unit()
  OPEN( iunit, file = TRIM(filename), form = 'formatted', status = 'unknown' )
  !
  WRITE( iunit, '(/2x,79("="))')
  WRITE( iunit, '(/24x,"Hubbard parameters from eACBN0")')
  WRITE( iunit, '(/24x,"computed by hp.x  (acbn0 = 2)")')
  WRITE( iunit, '(/2x,79("="),/)')
  WRITE( iunit, '(2x,"prefix     : ",a)') TRIM(prefix)
  WRITE( iunit, '(2x,"projectors : ortho-atomic")')
  WRITE( iunit, '(2x,"bookkeeping: MODULE ldaU (official QE) + own neighbour table")')
  WRITE( iunit, '(2x,"mixing     : 1.000  (one-shot; no SCF history)")')
  WRITE( iunit, '(2x,"neighbours : cutoff d_IJ = ",f6.3," A,  nnnc = ",i4)') &
        nn_distance, nnnc
  WRITE( iunit, '(2x,"estimator  : pseudo-hybrid eACBN0 ratio, NOT linear response")')
  !
  ! ---- on-site U ------------------------------------------------------
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"On-site U   = hub_v[na, na, m1, m2, cell 0]")')
  WRITE( iunit, '(2x,79("-"))')
  WRITE( iunit, '(/4x,"site n.  type  label  channel   occupation    U (eV)")')
  WRITE( iunit, '( 4x,"(one row per channel; occupation and U of that channel only)")')
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     ! One row per channel (see the same loop in hp_acbn0_write).
     DO m1 = 1, lm(nt)
        IF ( ANY( lch(1:m1-1,nt) == lch(m1,nt) ) ) CYCLE
        occ = 0.0_DP
        DO is = 1, nspin
           DO m2 = 1, lm(nt)
              IF ( lch(m2,nt) == lch(m1,nt) ) occ = occ + ns(m2,m2,is,na,na,1)
           ENDDO
        ENDDO
        utr = hub_v(na,na,m1,m1,1) * rytoev
        WRITE( iunit, '(6x,i3,4x,i3,3x,a6,4x,a2,6x,f9.4,4x,f10.4)') &
              na, nt, TRIM(atm(nt)), hp_acbn0_chanlab(nt,lch(m1,nt)), occ, utr
     ENDDO
  ENDDO
  !
  ! ---- inter-site V, with the cell vector and distance -----------------
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"Inter-site V")')
  WRITE( iunit, '(2x,79("-"))')
  WRITE( iunit, '(/4x,"na1  type1     na2  type2      cell [n1 n2 n3]", &
        &"    d (A)    block      V (eV)")')
  WRITE( iunit, '( 4x,"(one row per (channel of na1, channel of na2) block; ", &
        &"V is constant inside a block)")')
  !
  npair = 0
  DO n1 = -nn_dim(1), nn_dim(1)
  DO n2 = -nn_dim(2), nn_dim(2)
  DO n3 = -nn_dim(3), nn_dim(3)
     CALL hp_acbn0_get_nnc_wrap( (/n1,n2,n3/), thisnnc )
     DO na1 = 1, nat
        nt1 = ityp(na1)
        IF ( lm(nt1) == 0 ) CYCLE
        DO na2 = 1, nat
           nt2 = ityp(na2)
           IF ( lm(nt2) == 0 ) CYCLE
           IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 .AND. na1 == na2 ) CYCLE
           IF ( MAXVAL(ABS(hub_v(na1,na2,1:lm(nt1),1:lm(nt2),thisnnc))) < veps ) CYCLE
           npair = npair + 1
           DO m1 = 1, lm(nt1)
              IF ( ANY( lch(1:m1-1,nt1) == lch(m1,nt1) ) ) CYCLE
              DO m2 = 1, lm(nt2)
                 IF ( ANY( lch(1:m2-1,nt2) == lch(m2,nt2) ) ) CYCLE
                 vtr = hub_v(na1,na2,m1,m2,thisnnc)
                 WRITE( iunit, '(4x,i3,2x,a6,4x,i3,2x,a6,4x,"[",3i3," ]",4x,f7.4,3x,a2,"-",a2,4x,f10.4)') &
                       na1, TRIM(atm(nt1)), na2, TRIM(atm(nt2)), n1, n2, n3, &
                       nn_leng(n1+ics,n2+ics,n3+ics,na1,na2)*alat*bohr_radius_angs, &
                       hp_acbn0_chanlab(nt1,lch(m1,nt1)), hp_acbn0_chanlab(nt2,lch(m2,nt2)), &
                       vtr*rytoev
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  ENDDO
  ENDDO
  !
  IF ( npair == 0 ) WRITE( iunit, '(/4x,"(none - check nn_distance)")')
  WRITE( iunit, '(/4x,i5," pair(s) with |V| > ",es8.1," Ry")') npair, veps
  !
  ! ---- one-shot Hubbard energy ------------------------------------------
  !
  CALL hp_acbn0_energy_v( ns, hub_v, e_on, e_int )
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"One-shot extended Hubbard energy  (NOT self-consistent, see below)")')
  WRITE( iunit, '(2x,79("-"))')
  WRITE( iunit, '(/4x,"on-site   (U)  : ",f16.9," Ry  = ",f12.6," eV")') e_on,        e_on*rytoev
  WRITE( iunit, '( 4x,"inter-site (V) : ",f16.9," Ry  = ",f12.6," eV")') e_int,       e_int*rytoev
  WRITE( iunit, '( 4x,"E(+U+V)        : ",f16.9," Ry  = ",f12.6," eV")') e_on+e_int, (e_on+e_int)*rytoev
  WRITE( iunit, '(/4x,"Evaluated with the one-shot U, V above and the occupations they were")')
  WRITE( iunit, '( 4x,"computed from; U, V and the density are not mutually self-consistent,")')
  WRITE( iunit, '( 4x,"so this is not a total-energy contribution.  It equals what the forked")')
  WRITE( iunit, '( 4x,"pw.x prints as E(+U+V) on the first iteration of a run started from this")')
  WRITE( iunit, '( 4x,"density with ehub_mixing = 1 (v_hubbard_v.f90).")')
  !
  ! ---- full on-site matrices ------------------------------------------
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"Full on-site U matrix  (eV)")')
  WRITE( iunit, '(2x,79("-"))')
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"   (",i2," x",i2,")")') &
           na, TRIM(atm(nt)), lm(nt), lm(nt)
     DO m1 = 1, lm(nt)
        WRITE( iunit, '(6x,16f10.4)') (hub_v(na,na,m1,m2,1)*rytoev, m2 = 1, lm(nt))
     ENDDO
  ENDDO
  !
  ! ---- on-site occupation matrices ------------------------------------
  !
  WRITE( iunit, '(/2x,79("-"))')
  WRITE( iunit, '(2x,"On-site occupation matrix  ns[na, m1, m2, spin]")')
  WRITE( iunit, '(2x,79("-"))')
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     DO is = 1, nspin
        WRITE( iunit, '(/4x,"atom ",i3,"  type ",a,"  spin ",i1)') na, TRIM(atm(nt)), is
        DO m1 = 1, lm(nt)
           WRITE( iunit, '(6x,16f10.4)') (ns(m1,m2,is,na,na,1), m2 = 1, lm(nt))
        ENDDO
     ENDDO
  ENDDO
  !
  WRITE( iunit, '(/2x,79("="),/)')
  CLOSE( iunit )
  !
  WRITE( stdout, '(/5x,"eACBN0 Hubbard parameters written to ",a,/)') TRIM(filename)
  !
  WRITE( stdout, '(5x,"site n.  type  label  channel      U (eV)")')
  DO na = 1, nat
     nt = ityp(na)
     IF ( lm(nt) == 0 ) CYCLE
     DO m1 = 1, lm(nt)
        IF ( ANY( lch(1:m1-1,nt) == lch(m1,nt) ) ) CYCLE
        WRITE( stdout, '(7x,i3,4x,i3,3x,a6,4x,a2,6x,f10.4)') &
              na, nt, TRIM(atm(nt)), hp_acbn0_chanlab(nt,lch(m1,nt)), hub_v(na,na,m1,m1,1)*rytoev
     ENDDO
  ENDDO
  WRITE( stdout, '(/5x,i5," inter-site pair(s) with non-zero V; see ",a)') &
        npair, TRIM(filename)
  WRITE( stdout, '(/5x,"one-shot E(+U+V) = ",f14.9," Ry   (on-site ",f14.9,",  inter-site ",f14.9,")")') &
        e_on+e_int, e_on, e_int
  WRITE( stdout, '(5x,"(not self-consistent: U, V and the density do not feed back into each other)",/)')
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_write_v

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_print_uv( hub_v )
  !---------------------------------------------------------------------
  !
  ! Print the U and V tables in the EXACT format pw.x uses (write_ehub_uv,
  ! PW/src/ldaV.f90:617).  Not decoration: it is what makes every V0/V1/V2/V3
  ! comparison a plain diff instead of an eyeball check.  The loop order and
  ! the nnc bookkeeping are copied line for line, with ldaV's arrays replaced
  ! by ldaU + MODULE hp_acbn0_nn.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : rytoev
  USE io_global,   ONLY : stdout, ionode
  USE ions_base,   ONLY : nat, ityp, ntyp => nsp
  USE ldaU,        ONLY : Hubbard_l, Hubbard_lmax, is_hubbard
  USE hp_acbn0_nn, ONLY : nn_vect, nn_dim, nn_csize, nnnc
  USE hp_acbn0_man,         ONLY : lm_max, lm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: hub_v(nat,nat,lm_max,lm_max,nnnc)
  !
  INTEGER :: na1, na2, nt, nt1, nt2, m1, m2, n1, n2, n3, nnc, ics
  !
  IF ( .NOT.ionode ) RETURN
  !
  !
  WRITE ( stdout,'(a)' ) '                                        '
  WRITE ( stdout,'(a)' ) '     -----------------------------------'
  WRITE ( stdout,'(a)' ) '     V[atom i, atom j, #orb(i), #orb(j)]'
  WRITE ( stdout,'(a)' ) '     -----------------------------------'
  !
  ics = nn_csize + 1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( lm(nt1) == 0 ) CYCLE
  DO na2 = 1, nat
     nt2 = ityp(na2)
     IF ( lm(nt2) == 0 ) CYCLE
     IF ( na1 /= na2 .AND. nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     WRITE ( stdout,'(4(a,i2),a)' ) &
           '     [ ', na1, ', ', na2, ', ', lm(nt1), ', ', lm(nt2), '  ]'
     WRITE ( stdout,'(a,i2,2i6,a)' ) '     [ ', 0, 0, 0, '  ]'
     DO m1 = 1, lm(nt1)
        WRITE ( stdout,'(4x,16f9.4)' ) &
              (hub_v(na1,na2,m1,m2,1) * rytoev, m2 = 1, lm(nt2))
     ENDDO
     WRITE ( stdout,* ) ''
  ENDDO ! na2
  ENDDO ! na1
  !
  nnc = 2
  DO n1 = -1*nn_dim(1), nn_dim(1)
  DO n2 = -1*nn_dim(2), nn_dim(2)
  DO n3 = -1*nn_dim(3), nn_dim(3)
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     DO na1 = 1, nat
        nt1 = ityp(na1)
        IF ( lm(nt1) == 0 ) CYCLE
     DO na2 = 1, nat
        nt2 = ityp(na2)
        IF ( lm(nt2) == 0 ) CYCLE
        IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
        WRITE ( stdout,'(4(a,i2),a)' ) &
              '     [ ', na1, ', ', na2, ', ', lm(nt1), ', ', lm(nt2), '  ]'
        WRITE ( stdout,'(a,i2,2i6,a)' ) '     [ ', n1, n2, n3, '  ]'
        DO m1 = 1, lm(nt1)
           WRITE ( stdout,'(4x,16f9.4)' ) &
                 (hub_v(na1,na2,m1,m2,nnc) * rytoev, m2 = 1, lm(nt2))
        ENDDO
        WRITE ( stdout,* ) ''
     ENDDO ! na2
     ENDDO ! na1
     IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
  ENDDO ! n3
  ENDDO ! n2
  ENDDO ! n1
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_print_uv

!-----------------------------------------------------------------------
SUBROUTINE hp_acbn0_get_nnc_wrap( n, nnc )
  !---------------------------------------------------------------------
  ! Thin wrapper so the writer, which is an external subroutine, can reach
  ! the module procedure without USE-ing it into a scope that already has a
  ! local named nnc.
  USE hp_acbn0_nn, ONLY : hp_acbn0_get_nnc
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: n(3)
  INTEGER, INTENT(OUT) :: nnc
  CALL hp_acbn0_get_nnc( n, nnc )
  RETURN
END SUBROUTINE hp_acbn0_get_nnc_wrap

!-----------------------------------------------------------------------
!  V3 - transplanted eACBN0 (U+V) kernels, on ldaU bookkeeping
!-----------------------------------------------------------------------
!
!  Same four kernels V2 validated against pw.x, with the ldaV bookkeeping
!  replaced.  symmetrize_n4 is no longer duplicated: V2 needed a ldaV copy
!  because the S3 one had already been converted, and after V3 they are the
!  same routine, so hp_acbn0_symm_n4 above serves both branches.
!
!      hp_acbn0_ns_e0_1   <- new_ns_e0_1     636-1808
!      hp_acbn0_ns_v      <- new_ns_v         10- 232
!      hp_acbn0_symm_n4m  <- symmetrize_n4m 3168-3399  (+bugfix)
!      hp_acbn0_symm_n5   <- symmetrize_n5  3400-3779  (+bugfix)
!
SUBROUTINE hp_acbn0_ns_e0_1(ns, hub_v)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi, rytoev
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE cell_base,            ONLY : at, alat
  USE klist,                ONLY : nks, ngk, xk
  USE ldaU,                 ONLY : Hubbard_l, Hubbard_lmax, is_hubbard, &
                                   wfcU, nwfcU
  USE hp_acbn0_nn,          ONLY : nn_vect, nn_leng, nn_dim, nn_csize, nnnc, &
                                   hp_acbn0_get_nnc
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch, ofs
  USE ldaU_hp,              ONLY : acbn0_rm_ehub_u

  IMPLICIT NONE

  REAL(DP), INTENT(OUT)    :: ns(lm_max,lm_max,nspin,nat,nat,nnnc)
  REAL(DP), INTENT(OUT)    :: hub_v(nat,nat,lm_max,lm_max,nnnc)

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

  TYPE(hubbard_ene_4)      :: E0(nnnc)
  TYPE(hubbard_ene_4)      :: HF(nnnc)

  TYPE(occmat_5)           :: n(nnnc)
  TYPE(occmat_5)           :: p(nnnc)
  TYPE(occmat_5)           :: num(nbnd)

  TYPE(occmat_4)           :: dm0
  TYPE(occmat_4)           :: dm1(ntyp,4)

  COMPLEX(DP)              :: kphase

  REAL(DP), ALLOCATABLE    :: hub_v_(:,:,:,:)

  REAL(DP)                 :: arg
  REAL(DP)                 :: ntau(3,nat)
  REAL(DP), EXTERNAL       :: hp_acbn0_compute_eri
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
  !
  !

  CALL start_clock('acbn0_ns_e0_1')
  !
  !
  !
  ! The stop_acbn0 short-circuit that stood here is SCF-convergence
  ! bookkeeping for pw.x's ACBN0 loop: when the V's have stopped moving it
  ! skips the expensive kernel and just refreshes the occupations.  hp.x is
  ! one-shot, so there is no previous iteration to compare against and the
  ! branch is unreachable.  Removing it also removes the only readers of
  ! ehub_v_diff and ehub_conv_thr, neither of which hp.x can initialise.

  CALL allocate_bec_type(nwfcU, nbnd, proj)
  !
  ALLOCATE ( dm0%O(lm_max,lm_max,nspin,nat) ) ; dm0%O = 0.0_DP
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        ALLOCATE ( dm1(nt,l)%O(lm_max,lm_max,nspin,nat) ) ; dm1(nt,l)%O = 0.0_DP
     ENDDO
  ENDDO
  !
  DO i = 1, nnnc
     ALLOCATE ( n(i)%O(lm_max,lm_max,nspin,nat,nat) ) ; n(i)%O = 0.0_DP
     ALLOCATE ( p(i)%O(lm_max,lm_max,nspin,nat,nat) ) ; p(i)%O = 0.0_DP
     ALLOCATE ( E0(i)%V(4,4,nat,nat), E0(i)%J(4,4,nat,nat) )
     ALLOCATE ( HF(i)%V(4,4,nat,nat), HF(i)%J(4,4,nat,nat) )
  ENDDO

  ics = nn_csize + 1
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcU | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( nks > 1 ) CALL get_buffer(wfcU, nwordwfcU, iunhub, ik)
     CALL calbec(npw, wfcU, evc, proj)
     !
     ! calculate the occupation matrices n
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. nn_vect(ics,ics,ics,na1,na2 ) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
        DO m2 = 1, lm(nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                         proj%r(ofs(m2,na2),ibnd) * &
                         proj%r(ofs(m1,na1),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                    DBLE(proj%k(ofs(m2,na2),ibnd) * &
                   CONJG(proj%k(ofs(m1,na1),ibnd)))
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
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
        DO m2 = 1, lm(nt2)
           n(1)%O(m2,m1,current_spin,na2,na1) = n(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     nnc = 2
     !
     DO n1 = -1*nn_dim(1), nn_dim(1)
     DO n2 = -1*nn_dim(2), nn_dim(2)
     DO n3 = -1*nn_dim(3), nn_dim(3)
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
           IF ( lm(nt1) == 0 ) CYCLE
           !
        DO na2 = 1, nat
           nt2 = ityp(na2)
           !
           IF ( lm(nt2) == 0 ) CYCLE
           !
           IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           !
           DO m1 = 1, lm(nt1)
           DO m2 = 1, lm(nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                             wg(ibnd,ik) * &
                              proj%r(ofs(m2,na2),ibnd) * &
                              proj%r(ofs(m1,na1),ibnd)
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                             wg(ibnd,ik) * &
                                           DBLE((kphase) * &
                              proj%k(ofs(m2,na2),ibnd) * &
                        CONJG(proj%k(ofs(m1,na1),ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
           ENDDO ! m2
           ENDDO ! m1
           !
        ENDDO ! na2
        ENDDO ! na1
        !
        IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
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
        IF ( lm(nt) == 0 ) CYCLE
        !
        DO m = 1, lm(nt)
           !
           l = lch(m,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt,nt) + &
                            proj%r(ofs(m,na),ibnd) * &
                            proj%r(ofs(m,na),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 num(ibnd)%O(l,l,current_spin,nt,nt) = &
                 num(ibnd)%O(l,l,current_spin,nt,nt) + &
                       DBLE(proj%k(ofs(m,na),ibnd) * &
                      CONJG(proj%k(ofs(m,na),ibnd)))
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
        IF ( lm(nt) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt)
           !
           l1 = lch(m1,nt)
           !
        DO m2 = 1, lm(nt)
           !
           l2 = lch(m2,nt)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 dm0%O(m1,m2,current_spin,na) = &
                 dm0%O(m1,m2,current_spin,na) + &
                           num(ibnd)%O(l1,l2,   &
                          current_spin,nt,nt) * &
                                  wg(ibnd,ik) * &
                    proj%r(ofs(m2,na),ibnd) * &
                    proj%r(ofs(m1,na),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                 dm0%O(m1,m2,current_spin,na) = &
                 dm0%O(m1,m2,current_spin,na) + &
                           num(ibnd)%O(l1,l2,   &
                          current_spin,nt,nt) * &
                                  wg(ibnd,ik) * &
               DBLE(proj%k(ofs(m2,na),ibnd) * &
              CONJG(proj%k(ofs(m1,na),ibnd)))
              ENDDO ! ibnd
           ENDIF
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na

     DO nt = 1, ntyp
        !
        IF ( lm(nt) == 0 ) CYCLE
        !
        DO l = 1, 4
           !
           DO na1 = 1, nat
              nt1 = ityp(na1)
              !
              IF ( lm(nt1) == 0 ) CYCLE
              !
              DO m1 = 1, lm(nt1)
              DO m2 = 1, lm(nt1)
                 !
                 IF ( gamma_only ) THEN
                    DO ibnd = 1, nbnd
                       dm1(nt,l)%O(m1,m2,current_spin,na1) = &
                       dm1(nt,l)%O(m1,m2,current_spin,na1) + &
                             num(ibnd)%O(lch(m1,nt1),l,   &
                                      current_spin,nt1,nt) * &
                                               wg(ibnd,ik) * &
                                proj%r(ofs(m2,na1),ibnd) * &
                                proj%r(ofs(m1,na1),ibnd)
                    ENDDO ! ibnd
                 ELSE
                    DO ibnd = 1, nbnd
                       dm1(nt,l)%O(m1,m2,current_spin,na1) = &
                       dm1(nt,l)%O(m1,m2,current_spin,na1) + &
                             num(ibnd)%O(lch(m1,nt1),l,   &
                                      current_spin,nt1,nt) * &
                                               wg(ibnd,ik) * &
                           DBLE(proj%k(ofs(m2,na1),ibnd) * &
                          CONJG(proj%k(ofs(m1,na1),ibnd)))
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
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
           !
           l1 = lch(m1,nt1)
           !
        DO m2 = 1, lm(nt2)
           !
           l2 = lch(m2,nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                   p(1)%O(m1,m2,current_spin,na1,na2) = &
                   p(1)%O(m1,m2,current_spin,na1,na2) + &
                                   num(ibnd)%O(l1,l2,   &
                                current_spin,nt1,nt2) * &
                                          wg(ibnd,ik) * &
                           proj%r(ofs(m2,na2),ibnd) * &
                           proj%r(ofs(m1,na1),ibnd)
              ENDDO ! ibnd
           ELSE
              DO ibnd = 1, nbnd
                   p(1)%O(m1,m2,current_spin,na1,na2) = &
                   p(1)%O(m1,m2,current_spin,na1,na2) + &
                                   num(ibnd)%O(l1,l2,   &
                                current_spin,nt1,nt2) * &
                                          wg(ibnd,ik) * &
                      DBLE(proj%k(ofs(m2,na2),ibnd) * &
                     CONJG(proj%k(ofs(m1,na1),ibnd)))
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
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
        DO m2 = 1, lm(nt2)
           p(1)%O(m2,m1,current_spin,na2,na1) = p(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     nnc = 2
     !
     DO n1 = -1*nn_dim(1), nn_dim(1)
     DO n2 = -1*nn_dim(2), nn_dim(2)
     DO n3 = -1*nn_dim(3), nn_dim(3)
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
           IF ( lm(nt1) == 0 ) CYCLE
           !
        DO na2 = 1, nat
           nt2 = ityp(na2)
           !
           IF ( lm(nt2) == 0 ) CYCLE
           !
           IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           !
           DO m1 = 1, lm(nt1)
              !
              l1 = lch(m1,nt1)
              !
           DO m2 = 1, lm(nt2)
              !
              l2 = lch(m2,nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    p(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    p(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                      num(ibnd)%O(l1,l2,   &
                                   current_spin,nt1,nt2) * &
                                             wg(ibnd,ik) * &
                              proj%r(ofs(m2,na2),ibnd) * &
                              proj%r(ofs(m1,na1),ibnd)
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    p(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    p(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                      num(ibnd)%O(l1,l2,   &
                                   current_spin,nt1,nt2) * &
                                             wg(ibnd,ik) * &
                                           DBLE((kphase) * &
                              proj%k(ofs(m2,na2),ibnd) * &
                        CONJG(proj%k(ofs(m1,na1),ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
           ENDDO ! m2
           ENDDO ! m1
           !
        ENDDO ! na2
        ENDDO ! na1
        !
        IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
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
  DO i = 1, nnnc
     CALL mp_sum(n(i)%O, inter_pool_comm)
     CALL mp_sum(p(i)%O, inter_pool_comm)
  ENDDO

  DO i = 1, nnnc
     IF ( nspin == 1 ) n(i)%O = 0.5_DP * n(i)%O
  ENDDO
  !
  CALL hp_acbn0_symm_n5(n)
  !
  ns = 0.0_DP
  !
  DO i = 1, nnnc
     ns(:,:,:,:,:,i) = n(i)%O(:,:,:,:,:)
     DEALLOCATE ( n(i)%O )
  ENDDO
  !
  IF ( nspin == 1 ) dm0%O = 0.5_DP * dm0%O
  !
  CALL hp_acbn0_symm_n4(dm0)
  !
  DO nt = 1, ntyp
     DO l = 1, 4
        IF ( nspin == 1 ) dm1(nt,l)%O = 0.5_DP * dm1(nt,l)%O
     ENDDO
  ENDDO
  !
  CALL hp_acbn0_symm_n4m(dm1)
  !
  DO i = 1, nnnc
     IF ( nspin == 1 ) p(i)%O = 0.5_DP * p(i)%O
  ENDDO
  !
  CALL hp_acbn0_symm_n5(p)
  !
  ! on-site interaction: U_ij^II
  !
  E0(1)%V = 0.0_DP ; E0(1)%J = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, lm(nt)
        CALL hp_acbn0_orberi(lch(m1,nt), l1)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m1,nt), c1)
        !
     l2 = 0
     DO m2 = 1, lm(nt)
        CALL hp_acbn0_orberi(lch(m2,nt), l2)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m2,nt), c2)
        !
     l3 = 0
     DO m3 = 1, lm(nt)
        IF ( lch(m3,nt) /= lch(m1,nt) ) CYCLE
        CALL hp_acbn0_orberi(lch(m3,nt), l3)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m3,nt), c3)
        !
     l4 = 0
     DO m4 = 1, lm(nt)
        IF ( lch(m4,nt) /= lch(m2,nt) ) CYCLE
        CALL hp_acbn0_orberi(lch(m4,nt), l4)
        CALL hp_acbn0_pao_3g(upf(nt)%psd, lch(m4,nt), c4)
        !
        IF ( lch(m1,nt) == 1 .AND. lch(m2,nt) == 1 ) THEN
           E0(1)%V(lch(m1,nt),lch(m2,nt),na,na)  = &
           E0(1)%V(lch(m1,nt),lch(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,1,na)) * &
                             hp_acbn0_compute_eri(l3,l1,l4,l2,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c1,c4,c2)
        ELSE
           E0(1)%V(lch(m1,nt),lch(m2,nt),na,na)  = &
           E0(1)%V(lch(m1,nt),lch(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,1,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,1,na)*dm0%O(m4,m2,2,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,1,na)) * &
                             hp_acbn0_compute_eri(l3,l1,l4,l2,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                       tau(:,na)*alat,tau(:,na)*alat,    &
                                         c3,c1,c4,c2)

           E0(1)%J(lch(m1,nt),lch(m2,nt),na,na)  = &
           E0(1)%J(lch(m1,nt),lch(m2,nt),na,na)  + &
                 (dm0%O(m3,m1,1,na)*dm0%O(m4,m2,1,na)  + &
                  dm0%O(m3,m1,2,na)*dm0%O(m4,m2,2,na)) * &
                             hp_acbn0_compute_eri(l3,l2,l4,l1,    &
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
     IF ( lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt)
     DO m2 = 1, lm(nt)
        !
        IF ( m1 == m2 ) THEN
            HF(1)%V(lch(m1,nt),lch(m2,nt),na,na) = &
            HF(1)%V(lch(m1,nt),lch(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,1,na,na,1)
        ELSE
            HF(1)%V(lch(m1,nt),lch(m2,nt),na,na) = &
            HF(1)%V(lch(m1,nt),lch(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,2,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,1,na,na,1)

            HF(1)%J(lch(m1,nt),lch(m2,nt),na,na) = &
            HF(1)%J(lch(m1,nt),lch(m2,nt),na,na) + &
             ns(m1,m1,1,na,na,1) * ns(m2,m2,1,na,na,1) + &
             ns(m1,m1,2,na,na,1) * ns(m2,m2,2,na,na,1)
        ENDIF
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na
  !
  ALLOCATE ( hub_v_(nat,nat,lm_max,lm_max) ) ; hub_v_ = 0.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt)
     DO m2 = 1, lm(nt)
        !
        IF ( lch(m1,nt) == 1 .AND. lch(m2,nt) == 1 ) THEN
           !
           ! both indices in an s channel: no exchange term.  See the same
           ! test in hp_acbn0_ns_a0 for why both must be checked.
           !
           hub_v_(na,na,m1,m2) = E0(1)%V(lch(m1,nt),lch(m2,nt),na,na) / &
                                 HF(1)%V(lch(m1,nt),lch(m2,nt),na,na)
        ELSE
           hub_v_(na,na,m1,m2) = E0(1)%V(lch(m1,nt),lch(m2,nt),na,na)  / &
                                 HF(1)%V(lch(m1,nt),lch(m2,nt),na,na)  - &
                                 E0(1)%J(lch(m1,nt),lch(m2,nt),na,na)  / &
                                 HF(1)%J(lch(m1,nt),lch(m2,nt),na,na)
        ENDIF
        !
        ! One-shot (no mixing).  acbn0_rm_ehub_u is Lee & Son's
        ! remove_ehub_u with the flag inverted (1 = remove), applied exactly
        ! as new_ns_v.f90:1362-1371 does - a product over both channel
        ! indices, on the on-site block only, AFTER the ratio is formed.
        ! V (every na1 /= na2 or nnc > 1 block) is never masked; the channel
        ! stays in the manifold, so the N^IJ weights and the ERI sums above
        ! are untouched.
        !
        hub_v(na,na,m1,m2,1) = hub_v_(na,na,m1,m2)                          &
                             * DBLE( 1 - acbn0_rm_ehub_u(nt,lch(m1,nt)) )    &
                             * DBLE( 1 - acbn0_rm_ehub_u(nt,lch(m2,nt)) )
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
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     l1 = 0
     DO m1 = 1, lm(nt1)
        CALL hp_acbn0_orberi(lch(m1,nt1), l1)
        CALL hp_acbn0_pao_3g(upf(nt1)%psd, lch(m1,nt1), c1)
        !
     l2 = 0
     DO m2 = 1, lm(nt2)
        CALL hp_acbn0_orberi(lch(m2,nt2), l2)
        CALL hp_acbn0_pao_3g(upf(nt2)%psd, lch(m2,nt2), c2)
        !
     l3 = 0
     DO m3 = 1, lm(nt1)
        IF ( lch(m3,nt1) /= lch(m1,nt1) ) CYCLE
        CALL hp_acbn0_orberi(lch(m3,nt1), l3)
        CALL hp_acbn0_pao_3g(upf(nt1)%psd, lch(m3,nt1), c3)
        !
     l4 = 0
     DO m4 = 1, lm(nt2)
        IF ( lch(m4,nt2) /= lch(m2,nt2) ) CYCLE
        CALL hp_acbn0_orberi(lch(m4,nt2), l4)
        CALL hp_acbn0_pao_3g(upf(nt2)%psd, lch(m4,nt2), c4)
        !
        E0(1)%V(lch(m1,nt1),lch(m2,nt2),na1,na2)  = &
        E0(1)%V(lch(m1,nt1),lch(m2,nt2),na1,na2)  + &
               (dm1(nt2,lch(m2,nt2))%O(m3,m1,1,na1)  * &
                dm1(nt1,lch(m1,nt1))%O(m4,m2,1,na2)  + &
                dm1(nt2,lch(m2,nt2))%O(m3,m1,2,na1)  * &
                dm1(nt1,lch(m1,nt1))%O(m4,m2,2,na2)  + &
                dm1(nt2,lch(m2,nt2))%O(m3,m1,1,na1)  * &
                dm1(nt1,lch(m1,nt1))%O(m4,m2,2,na2)  + &
                dm1(nt2,lch(m2,nt2))%O(m3,m1,2,na1)  * &
                dm1(nt1,lch(m1,nt1))%O(m4,m2,1,na2)  - &
                               p(1)%O(m3,m2,1,na1,na2)  * &
                               p(1)%O(m4,m1,1,na2,na1)  - &
                               p(1)%O(m3,m2,2,na1,na2)  * &
                               p(1)%O(m4,m1,2,na2,na1)) * &
                              hp_acbn0_compute_eri(l3,l1,l4,l2,    &
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
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt1)
     DO m2 = 1, lm(nt2)
        !
        E0(1)%V(lch(m2,nt2),lch(m1,nt1),na2,na1) = E0(1)%V(lch(m1,nt1),lch(m2,nt2),na1,na2)
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
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     IF ( nn_vect(ics,ics,ics,na1,na2) == 0) CYCLE
     !
     DO m1 = 1, lm(nt1)
     DO m2 = 1, lm(nt2)
        !
        HF(1)%V(lch(m1,nt1),lch(m2,nt2),na1,na2) = &
        HF(1)%V(lch(m1,nt1),lch(m2,nt2),na1,na2) + &
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
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt1)
     DO m2 = 1, lm(nt2)
        HF(1)%V(lch(m2,nt2),lch(m1,nt1),na2,na1) = HF(1)%V(lch(m1,nt1),lch(m2,nt2),na1,na2)
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = 1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt1)
        !
        l1 = lch(m1,nt1)
        !
     DO m2 = 1, lm(nt2)
        !
        l2 = lch(m2,nt2)
        !
        hub_v_(na1,na2,m1,m2) = &
               (E0(1)%V(l1,l2,na1,na2) / HF(1)%V(l1,l2,na1,na2)) / 2.0_DP
        !
        hub_v(na1,na2,m1,m2,1) = hub_v_(na1,na2,m1,m2)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  ! ... out of the unit cell
  !
  DO i = 2, nnnc
     E0(i)%V = 0.0_DP
  ENDDO
  !
  nnc = 2
  !
  DO n1 = -1*nn_dim(1), nn_dim(1)
  DO n2 = -1*nn_dim(2), nn_dim(2)
  DO n3 = -1*nn_dim(3), nn_dim(3)
     !
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = 1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
        !
        DO ix = 1, 3
           ntau(ix,na2) = n1 * at(ix,1) + n2 * at(ix,2) + n3 * at(ix,3) + tau(ix,na2)
        ENDDO
        !
        l1 = 0
        DO m1 = 1, lm(nt1)
           CALL hp_acbn0_orberi(lch(m1,nt1), l1)
           CALL hp_acbn0_pao_3g(upf(nt1)%psd, lch(m1,nt1), c1)
           !
        l2 = 0
        DO m2 = 1, lm(nt2)
           CALL hp_acbn0_orberi(lch(m2,nt2), l2)
           CALL hp_acbn0_pao_3g(upf(nt2)%psd, lch(m2,nt2), c2)
           !
        l3 = 0
        DO m3 = 1, lm(nt1)
           IF ( lch(m3,nt1) /= lch(m1,nt1) ) CYCLE
           CALL hp_acbn0_orberi(lch(m3,nt1), l3)
           CALL hp_acbn0_pao_3g(upf(nt1)%psd, lch(m3,nt1), c3)
           !
        l4 = 0
        DO m4 = 1, lm(nt2)
           IF ( lch(m4,nt2) /= lch(m2,nt2) ) CYCLE
           CALL hp_acbn0_orberi(lch(m4,nt2), l4)
           CALL hp_acbn0_pao_3g(upf(nt2)%psd, lch(m4,nt2), c4)
           !
           E0(nnc)%V(lch(m1,nt1),lch(m2,nt2),na1,na2)  = &
           E0(nnc)%V(lch(m1,nt1),lch(m2,nt2),na1,na2)  + &
                    (dm1(nt2,lch(m2,nt2))%O(m3,m1,1,na1)  * &
                     dm1(nt1,lch(m1,nt1))%O(m4,m2,1,na2)  + &
                     dm1(nt2,lch(m2,nt2))%O(m3,m1,2,na1)  * &
                     dm1(nt1,lch(m1,nt1))%O(m4,m2,2,na2)  + &
                     dm1(nt2,lch(m2,nt2))%O(m3,m1,1,na1)  * &
                     dm1(nt1,lch(m1,nt1))%O(m4,m2,2,na2)  + &
                     dm1(nt2,lch(m2,nt2))%O(m3,m1,2,na1)  * &
                     dm1(nt1,lch(m1,nt1))%O(m4,m2,1,na2)  - &
                                  p(nnc)%O(m3,m2,1,na1,na2)  * &
                                  p(nnc)%O(m1,m4,1,na1,na2)  - &
                                  p(nnc)%O(m3,m2,2,na1,na2)  * &
                                  p(nnc)%O(m1,m4,2,na1,na2)) * &
                                   hp_acbn0_compute_eri(l3,l1,l4,l2,    &
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
     IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
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
  DO i = 2, nnnc
     HF(i)%V = 0.0_DP
  ENDDO
  !
  nnc = 2
  !
  DO n1 = -1*nn_dim(1), nn_dim(1)
  DO n2 = -1*nn_dim(2), nn_dim(2)
  DO n3 = -1*nn_dim(3), nn_dim(3)
     !
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = 1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
           !
           l1 = lch(m1,nt1)
           !
        DO m2 = 1, lm(nt2)
           !
           l2 = lch(m2,nt2)
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
     IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
     !
  ENDDO ! n3
  ENDDO ! n2
  ENDDO ! n1
  !
  nnc = 2
  !
  DO n1 = -1*nn_dim(1), nn_dim(1)
  DO n2 = -1*nn_dim(2), nn_dim(2)
  DO n3 = -1*nn_dim(3), nn_dim(3)
     !
     IF ( n1 == 0 .AND. n2 == 0 .AND. n3 == 0 ) CYCLE
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = 1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
           !
           l1 = lch(m1,nt1)
           !
        DO m2 = 1, lm(nt2)
           !
           l2 = lch(m2,nt2)
           !
           hub_v_(na1,na2,m1,m2) = &
                  (E0(nnc)%V(l1,l2,na1,na2) / HF(nnc)%V(l1,l2,na1,na2)) / 2.0_DP
           !
           hub_v(na1,na2,m1,m2,nnc) = hub_v_(na1,na2,m1,m2)
           !
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
     !
  ENDDO ! n3
  ENDDO ! n2
  ENDDO ! n1
  !
  DEALLOCATE ( hub_v_ )
  !
  DO i = 1, nnnc
     DEALLOCATE ( E0(i)%V, E0(i)%J )
     DEALLOCATE ( HF(i)%V, HF(i)%J )
  ENDDO
  !
  ! Hermiticity of Hubbard parameters U & V
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     !
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt1)
     DO m2 = 1, lm(nt2)
        hub_v(na1,na2,m1,m2,1) = 0.5D0 * (hub_v(na1,na2,m1,m2,1) + &
                                           hub_v(na2,na1,m2,m1,1))
        !
        hub_v(na2,na1,m2,m1,1) = hub_v(na1,na2,m1,m2,1)
        !
     ENDDO ! m2
     ENDDO ! m1
     !
  ENDDO ! na2
  ENDDO ! na1
  !
  DO na1 = 1, nat
  DO na2 = 1, nat
     IF ( na1 /= na2 .AND. nn_vect(ics,ics,ics,na1,na2) == 0 ) &
        hub_v(na1,na2,:,:,1) = 0.0_DP
  ENDDO ! na2
  ENDDO ! na1

  CALL stop_clock('acbn0_ns_e0_1')

  RETURN

END SUBROUTINE hp_acbn0_ns_e0_1

SUBROUTINE hp_acbn0_ns_v(ns)

  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at
  USE klist,                ONLY : nks, ngk, xk
  USE ldaU,                 ONLY : Hubbard_projectors, Hubbard_l, Hubbard_lmax, &
                                   is_hubbard, wfcU, nwfcU
  USE hp_acbn0_nn,          ONLY : nn_vect, nn_dim, nn_csize, nnnc
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec,  &
                                   allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : upf
  USE hp_acbn0_man,         ONLY : lm_max, lm, ofs

  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: ns(lm_max,lm_max,nspin,nat,nat,nnnc)

  TYPE occmat_5
     REAL(DP), ALLOCATABLE :: O(:,:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:,:)
  END TYPE occmat_5

  TYPE(bec_type)        :: proj

  COMPlEX(DP)           :: kphase

  REAL(DP)              :: arg

  TYPE(occmat_5)        :: n(nnnc)

  INTEGER               :: ik, ibnd, is, i, npw
  INTEGER               :: na1, na2, nt1, nt2
  INTEGER               :: m1, m2, m3, m4
  INTEGER               :: n1, n2, n3, nnc, ics
  INTEGER               :: nt

  CALL start_clock('acbn0_ns_v')
  !
  !
  !
  CALL allocate_bec_type(nwfcU, nbnd, proj)
  !
  DO i = 1, nnnc
     ALLOCATE ( n(i)%O(lm_max,lm_max,nspin,nat,nat) ) ; n(i)%O = 0.0_DP
  ENDDO
  !
  ics = nn_csize + 1
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk(ik)
     !
     ! projection < wfcU | evc >
     !
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)     
     IF ( Hubbard_projectors == 'pseudo' ) THEN
        CALL errore('hp_acbn0_ns_v', 'Hubbard_projectors = pseudo is not supported', 1)
     ELSE
        IF ( nks > 1 ) CALL get_buffer(wfcU, nwordwfcU, iunhub, ik)
        CALL calbec(npw, wfcU, evc, proj)
     ENDIF
     !
     ! calculate the occupation matrices n
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( na1 /= na2 .AND. nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
        DO m2 = 1, lm(nt2)
           !
           IF ( gamma_only ) THEN
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                         proj%r(ofs(m2,na2),ibnd) * &
                         proj%r(ofs(m1,na1),ibnd)
              ENDDO
           ELSE
              DO ibnd = 1, nbnd
                 n(1)%O(m1,m2,current_spin,na1,na2) = &
                 n(1)%O(m1,m2,current_spin,na1,na2) + &
                                        wg(ibnd,ik) * &
                    DBLE(proj%k(ofs(m2,na2),ibnd) * &
                   CONJG(proj%k(ofs(m1,na1),ibnd)))
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
        IF ( lm(nt1) == 0 ) CYCLE
        !
     DO na2 = na1, nat
        nt2 = ityp(na2)
        !
        IF ( lm(nt2) == 0 ) CYCLE
        !
        IF ( nn_vect(ics,ics,ics,na1,na2) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
        DO m2 = 1, lm(nt2)
           n(1)%O(m2,m1,current_spin,na2,na1) = n(1)%O(m1,m2,current_spin,na1,na2)
        ENDDO ! m2
        ENDDO ! m1
        !
     ENDDO ! na2
     ENDDO ! na1
     !
     nnc = 2
     !
     DO n1 = -1*nn_dim(1), nn_dim(1)
     DO n2 = -1*nn_dim(2), nn_dim(2)
     DO n3 = -1*nn_dim(3), nn_dim(3)
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
           IF ( lm(nt1) == 0 ) CYCLE
           !
        DO na2 = 1, nat
           nt2 = ityp(na2)
           !
           IF ( lm(nt2) == 0 ) CYCLE
           !
           IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
           !
           DO m1 = 1, lm(nt1)
           DO m2 = 1, lm(nt2)
              !
              IF ( gamma_only ) THEN
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                 DBLE(kphase*wg(ibnd,ik) * &
                              proj%r(ofs(m2,na2),ibnd) * &
                              proj%r(ofs(m1,na1),ibnd))
                 ENDDO ! ibnd
              ELSE
                 DO ibnd = 1, nbnd
                    n(nnc)%O(m1,m2,current_spin,na1,na2) = &
                    n(nnc)%O(m1,m2,current_spin,na1,na2) + &
                                 DBLE(kphase*wg(ibnd,ik) * &
                              proj%k(ofs(m2,na2),ibnd) * &
                        CONJG(proj%k(ofs(m1,na1),ibnd)))
                 ENDDO ! ibnd
              ENDIF
              !
           ENDDO ! m2
           ENDDO ! m1
           !
        ENDDO ! na2
        ENDDO ! na1
        !
        IF ( MAXVAL(nn_vect(n1+ics,n2+ics,n3+ics,:,:)) == 1 ) nnc = nnc + 1
        !
     ENDDO ! n3
     ENDDO ! n2
     ENDDO ! n1
     !
  ENDDO ! ik
  !
  CALL deallocate_bec_type(proj)
  !
  DO i = 1, nnnc
     CALL mp_sum(n(i)%O, inter_pool_comm)
  ENDDO
  !
  DO i = 1, nnnc
     IF ( nspin == 1 ) n(i)%O = 0.5_DP * n(i)%O
  ENDDO
  !
  CALL hp_acbn0_symm_n5(n)
  !
  ns = 0.0_DP
  !
  DO i = 1, nnnc
     ns(:,:,:,:,:,i) = n(i)%O(:,:,:,:,:)
     DEALLOCATE ( n(i)%O )
  ENDDO
  !
  CALL stop_clock('acbn0_ns_v')
  !
  RETURN
  !
END SUBROUTINE hp_acbn0_ns_v

SUBROUTINE hp_acbn0_symm_n4m( n )

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE ldaU,                 ONLY : Hubbard_l, Hubbard_lmax, is_hubbard
  USE hp_acbn0_nn,          ONLY : nn_vect, nn_dim, nn_csize, nnnc, &
                                   hp_acbn0_symmetry, hp_acbn0_get_nnc
  USE lsda_mod,             ONLY : nspin
  USE symm_base,            ONLY : d1, d2, d3, nsym, irt
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch

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

  !

  DO i = 1, ntyp
     DO l = 1, 4
        ALLOCATE (n(i,l)%T(lm_max,lm_max,nspin,nat)) ; n(i,l)%T = 0.0_DP
     ENDDO
  ENDDO
  !
  DO nt = 1, ntyp
     !
     ! BUGFIX vs PW/src/new_ns_v.f90 (same defect in symmetrize_n4, _n4m and
     ! _n5).  Without this guard a type carrying no Hubbard manifold has
     ! lch(nt) = 0, and the next line writes l_separation(0,nt) - one
     ! element BEFORE the array.  Dormant whenever every type is
     ! Hubbard-active, which is why it has never been hit; fires the moment a
     ! spectator species is added.
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     temp = lch(1,nt)
     l_separation(lch(1,nt),nt) = 1
     !
     DO m = 2, lm(nt)
        !
        IF ( lch(m,nt) /= temp ) THEN
           temp = lch(m,nt)
           l_separation(lch(m,nt),nt) = m
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  DO nt = 1, ntyp
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     DO l = 1, 4
        !
        DO is = 1, nspin
           !
           DO na1 = 1, nat
              nt1 = ityp(na1)
              !
              IF ( lm(nt1) == 0 ) CYCLE
              !
              DO m1 = 1, lm(nt1)
                 m1_ = m1 - l_separation(lch(m1,nt1),nt1) + 1
                 !
              DO m2 = 1, lm(nt1)
                 m2_ = m2 - l_separation(lch(m2,nt1),nt1) + 1
                 !
                 DO isym = 1, nsym
                    !
                    na2 = irt(isym, na1)
                    !
                    DO m3 = 1, lm(nt1)
                       IF ( lch(m1,nt1) /= lch(m3,nt1) ) CYCLE
                       m3_ = m3 - l_separation(lch(m3,nt1),nt1) + 1
                       !
                    DO m4 = 1, lm(nt1)
                       IF ( lch(m2,nt1) /= lch(m4,nt1) ) CYCLE
                       m4_ = m4 - l_separation(lch(m4,nt1),nt1) + 1
                       !
                       IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)
      
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d1(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)

                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d2(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)

                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 1 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 2 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d1(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 3 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d2(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 4 ) THEN
                          n(nt,l)%T(m1,m2,is,na1) = n(nt,l)%T(m1,m2,is,na1) + &
                                                           d3(m3_,m1_,isym) * &
                                                    n(nt,l)%O(m3,m4,is,na2) * &
                                                           d3(m4_,m2_,isym) / &
                                                                     (nsym)
                       ELSE
                          CALL errore('hp_acbn0_ns_v','angular momentum not implemented',1)
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
     IF ( lm(nt) == 0 ) CYCLE
     !
     DO l = 1, 4
        !
        DO na1 = 1, nat
           nt1 = ityp(na1)
           !
           IF ( lm(nt1) == 0 ) CYCLE
           !
           DO is = 1, nspin
              !
              DO m1 = 1, lm(nt1)
              DO m2 = m1, lm(nt1)
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

END SUBROUTINE hp_acbn0_symm_n4m

SUBROUTINE hp_acbn0_symm_n5( n )

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE ldaU,                 ONLY : Hubbard_l, Hubbard_lmax, is_hubbard
  USE hp_acbn0_nn,          ONLY : nn_vect, nn_dim, nn_csize, nnnc, &
                                   hp_acbn0_symmetry, hp_acbn0_get_nnc
  USE lsda_mod,             ONLY : nspin
  USE symm_base,            ONLY : d1, d2, d3, nsym, irt
  USE hp_acbn0_man,         ONLY : lm_max, lm, lch

  IMPLICIT NONE

  TYPE occmat_5
     REAL(DP), ALLOCATABLE :: O(:,:,:,:,:)
     REAL(DP), ALLOCATABLE :: T(:,:,:,:,:)
  END TYPE occmat_5

  TYPE(occmat_5), INTENT(INOUT) :: n(nnnc)

  INTEGER                       :: is, isym, i, m, nt
  INTEGER                       :: na1, na2, na3, na4
  INTEGER                       :: nt1, nt2, nt3, nt4
  INTEGER                       :: m1, m2, m3, m4
  INTEGER                       :: m1_, m2_, m3_, m4_
  INTEGER                       :: n0(3), n1, n2, n3
  INTEGER                       :: nnc0, nnc, ics
  INTEGER                       :: temp, l_separation(5,ntyp)

  !

  DO i = 1, nnnc
     ALLOCATE (n(i)%T(lm_max,lm_max,nspin,nat,nat)) ; n(i)%T = 0.0_DP
  ENDDO
  !
  DO nt = 1, ntyp
     !
     ! BUGFIX vs PW/src/new_ns_v.f90 (same defect in symmetrize_n4, _n4m and
     ! _n5).  Without this guard a type carrying no Hubbard manifold has
     ! lch(nt) = 0, and the next line writes l_separation(0,nt) - one
     ! element BEFORE the array.  Dormant whenever every type is
     ! Hubbard-active, which is why it has never been hit; fires the moment a
     ! spectator species is added.
     !
     IF ( lm(nt) == 0 ) CYCLE
     !
     temp = lch(1,nt)
     l_separation(lch(1,nt),nt) = 1
     !
     DO m = 2, lm(nt)
        !
        IF ( lch(m,nt) /= temp ) THEN
           temp = lch(m,nt)
           l_separation(lch(m,nt),nt) = m
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ics = nn_csize + 1
  !
  DO is = 1, nspin
     !
     DO na1 = 1, nat
        nt1 = ityp(na1)
        !
        IF ( lm(nt1) == 0 ) CYCLE
        !
        DO m1 = 1, lm(nt1)
           m1_ = m1 - l_separation(lch(m1,nt1),nt1) + 1
           !
        DO m2 = 1, lm(nt1)
           m2_ = m2 - l_separation(lch(m2,nt1),nt1) + 1
           !
           DO isym = 1, nsym
              !
              na2 = irt(isym, na1)
              !
              DO m3 = 1, lm(nt1)
                 IF ( lch(m1,nt1) /= lch(m3,nt1) ) CYCLE
                 m3_ = m3 - l_separation(lch(m3,nt1),nt1) + 1
                 !
              DO m4 = 1, lm(nt1)
                 IF ( lch(m2,nt1) /= lch(m4,nt1) ) CYCLE
                 m4_ = m4 - l_separation(lch(m4,nt1),nt1) + 1
                 !
                 IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)

                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d1(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)

                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d2(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)

                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 1 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 2 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d1(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 3 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d2(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt1) == 4 ) THEN
                    n(1)%T(m1,m2,is,na1,na1) = n(1)%T(m1,m2,is,na1,na1) + &
                                                       d3(m3_,m1_,isym) * &
                                               n(1)%O(m3,m4,is,na2,na2) * &
                                                       d3(m4_,m2_,isym) / &
                                                                 (nsym)
                 ELSE
                    CALL errore('hp_acbn0_ns_v','angular momentum not implemented',1)
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
        IF ( lm(nt1) == 0 ) CYCLE
        !
        DO n1 = -1*nn_dim(1), nn_dim(1)
        DO n2 = -1*nn_dim(2), nn_dim(2)
        DO n3 = -1*nn_dim(3), nn_dim(3)
           !
           n0(1) = n1 ; n0(2) = n2 ; n0(3) = n3
           !
           CALL hp_acbn0_get_nnc(n0,nnc0)
           !
           DO na2 = 1, nat
              nt2 = ityp(na2)
              !
              IF ( lm(nt2) == 0 ) CYCLE
              !
              IF ( nn_vect(n1+ics,n2+ics,n3+ics,na1,na2) == 0 ) CYCLE
              !
              DO m1 = 1, lm(nt1)
                 m1_ = m1 - l_separation(lch(m1,nt1),nt1) + 1
                 !
              DO m2 = 1, lm(nt2)
                 m2_ = m2 - l_separation(lch(m2,nt2),nt2) + 1
                 !
                 DO isym = 1, nsym
                    !
                    CALL hp_acbn0_symmetry(na1,na2,na3,na4,n0,nnc,isym)
                    !
                    ! See the bugfix note in hp_acbn0_get_nnc: nnc = 0 means
                    ! this symmetry maps the pair outside the neighbour set,
                    ! so it contributes nothing.
                    !
                    IF ( nnc == 0 ) CYCLE
                    !
                    DO m3 = 1, lm(nt1)
                       IF ( lch(m1,nt1) /= lch(m3,nt1) ) CYCLE
                       m3_ = m3 - l_separation(lch(m3,nt1),nt1) + 1
                       !
                    DO m4 = 1, lm(nt2)
                       IF ( lch(m2,nt2) /= lch(m4,nt2) ) CYCLE
                       m4_ = m4 - l_separation(lch(m4,nt2),nt2) + 1
                       !
                       IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 1 .AND. lch(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)

                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 2 .AND. lch(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d1(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)

                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                            (nsym)
                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 3 .AND. lch(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d2(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)

                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt2) == 1 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt2) == 2 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d1(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt2) == 3 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d2(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE IF ( lch(m1,nt1) == 4 .AND. lch(m2,nt2) == 4 ) THEN
                          n(nnc0)%T(m1,m2,is,na1,na2) = n(nnc0)%T(m1,m2,is,na1,na2) + &
                                                                   d3(m3_,m1_,isym) * &
                                                        n(nnc )%O(m3,m4,is,na3,na4) * &
                                                                   d3(m4_,m2_,isym) / &
                                                                             (nsym)
                       ELSE
                          CALL errore('hp_acbn0_ns_v','angular momentum not implemented',1)
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
     IF ( lm(nt1) == 0 ) CYCLE
     !
  DO na2 = na1, nat
     nt2 = ityp(na2)
     !
     IF ( lm(nt2) == 0 ) CYCLE
     !
     DO m1 = 1, lm(nt1)
     DO m2 = 1, lm(nt2)
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
        IF ( na1 /= na2 .AND. nn_vect(ics,ics,ics,na1,na2) == 0 ) &
           n(1)%T(:,:,is,na1,na2) = 0.0_DP
     ENDDO
     !
  ENDDO
  ENDDO
  !
  DO i = 1, nnnc
     n(i)%O = n(i)%T
     DEALLOCATE (n(i)%T)
  ENDDO

END SUBROUTINE hp_acbn0_symm_n5
