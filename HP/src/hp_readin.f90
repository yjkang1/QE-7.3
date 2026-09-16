!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE hp_readin()
  !----------------------------------------------------------------------------
  !
  !  This routine reads the input parameters for HP and reads
  !  the data produced by PWscf.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : meta_ionode, meta_ionode_id
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE check_stop,       ONLY : max_seconds
  USE io_files,         ONLY : tmp_dir, prefix, create_directory
  USE control_flags,    ONLY : iverbosity
  USE control_lr,       ONLY : ethr_nscf, lrpa
  USE ldaU_hp,          ONLY : conv_thr_chi, thresh_init, find_atpert, skip_atom, skip_type, &
                               equiv_type, background, compute_hp, tmp_dir_hp,         &
                               perturb_only_atom, sum_pertq, determine_num_pert_only,  &
                               skip_equivalence_q, tmp_dir_save, niter_max, dist_thr,  &
                               disable_type_analysis, docc_thr, num_neigh, lmin, rmax, &
                               nmix, nq1, nq2, nq3, alpha_mix, start_q, last_q, maxter, &
                               determine_q_mesh_only, acbn0, nn_distance,       &
                         acbn0_l_choice, acbn0_rm_ehub_u
  USE paw_variables,    ONLY : okpaw
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER             :: ios, &    ! integer variable for I/O control
                         iter      ! counter on iterations
  CHARACTER (LEN=256) :: outdir
  CHARACTER(LEN=6)    :: int_to_char
  !
  NAMELIST / INPUTHP / prefix, outdir, nq1, nq2, nq3, skip_equivalence_q, dist_thr,   &
                         conv_thr_chi, skip_atom, skip_type, equiv_type, iverbosity,  &
                         background, thresh_init, find_atpert, max_seconds, rmax,     &
                         niter_max, alpha_mix, nmix, compute_hp, perturb_only_atom,   &
                         start_q, last_q, sum_pertq, ethr_nscf, num_neigh, lmin,      &
                         determine_num_pert_only, disable_type_analysis, docc_thr,    &
                         determine_q_mesh_only, acbn0, nn_distance,       &
                         acbn0_l_choice, acbn0_rm_ehub_u
  !
  ! Note: meta_ionode is a single processor that reads the input
  !       Data read from input is subsequently broadcast to all processors
  !       from meta_ionode_id (using the default communicator world_comm)
  !
  IF (meta_ionode) CALL input_from_file ()
  !
  ! Set default values for variables in namelist
  !
  prefix             = 'pwscf'
  conv_thr_chi       = 1.D-5
  thresh_init        = 1.D-14
  ethr_nscf          = 1.D-11
  docc_thr           = 5.D-5
  dist_thr           = 6.D-4  ! same as eps_dist in PW/src/ldaU.f90
  rmax               = 100.D0
  skip_atom(:)       = .FALSE.
  skip_type(:)       = .FALSE.
  perturb_only_atom(:)    = .FALSE.
  skip_equivalence_q      = .FALSE.
  determine_num_pert_only = .FALSE.
  determine_q_mesh_only   = .FALSE.
  disable_type_analysis   = .FALSE.
  equiv_type(:)      = 0
  find_atpert        = 1
  acbn0              = 0
  nn_distance        = -1.0D0
  acbn0_l_choice(:,:)   = 0
  acbn0_rm_ehub_u(:,:)  = 0
  background         = 'no'
  compute_hp         = .FALSE.
  sum_pertq          = .FALSE.
  num_neigh          = 6
  lmin               = 2
  nq1                = 1
  nq2                = 1
  nq3                = 1
  start_q            = 1
  last_q             = -1
  iverbosity         = 1
  niter_max          = 100     
  alpha_mix(:)       = 0.D0   
  alpha_mix(1)       = 0.3D0  
  nmix               = 4     
  max_seconds        = 1.E+7_DP
  lrpa               = .FALSE.   ! Needed in dv_of_drho
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  ! Reading the namelist inputhp
  !
  IF (meta_ionode) READ( 5, INPUTHP, IOSTAT = ios )
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore( 'hp_readin', 'reading inputhp namelist', ABS( ios ) )
  !
  ! Setup tmp_dir
  !
  tmp_dir = trimcheck (outdir)
  !
  ! Broadcast input parameters over all processors
  !
  CALL hp_bcast_input()
  !
  ! Here we finished the reading of the input file.
  ! Now allocate space for pwscf variables, read and check them.
  !
  tmp_dir_save = tmp_dir
  tmp_dir_hp = TRIM (tmp_dir) // 'HP' //'/'
  CALL create_directory(tmp_dir_hp)
  !
  ! Read various data produced by PWscf.
  ! In particular, read the unperturbed occupation matrices
  ! via calling the routine read_rho.
  ! read_file calls init_tab_atwfc which initializes tab_atwfc.
  !
  CALL read_file()
  !
  ! ACBN0 only: if the pw.x run carried no HUBBARD card, MODULE ldaU is empty
  ! and the transplanted kernels have no manifold to project onto.  With
  ! acbn0 > 0 the manifold is declared here from acbn0_l_choice and stock
  ! QE's own init_hubbard builds everything from it.  See the routine below.
  !
  CALL acbn0_manifold_setup()
  !
  ! Compute the trace of the occupation matrix and then compute
  ! the occupations and the magnetization
  !
  CALL hp_ns_trace()
  !
  ! Make sure all the features used in the PWscf calculation 
  ! are actually supported by HP.
  !
  CALL input_sanity()
  !
  RETURN
  !
CONTAINS
  !
!-----------------------------------------------------------------------
SUBROUTINE acbn0_manifold_setup()
  !---------------------------------------------------------------------
  !
  ! Declare the Hubbard manifold from acbn0_l_choice in &INPUTHP and let
  ! stock QE build everything else from it.
  !
  ! WHY THIS EXISTS.  The transplanted ACBN0 kernels read the manifold out of
  ! MODULE ldaU (as lm / lch / ofs, see MODULE hp_acbn0_man).  ldaU is filled
  ! by init_hubbard, which read_file calls ONLY when the XML says lda_plus_u
  ! (PW/src/read_file_new.f90:318-320) - i.e. only when the pw.x input carried
  ! a HUBBARD card.  A plain LDA/PBE ground state therefore arrives with
  ! is_hubbard all .false. and nwfcU = 0.  Here the declaration is made on the
  ! hp.x side instead - which is what Lee & Son do with ehub_l_choice in
  ! &system, one program earlier - and the pw.x ground state stays a plain
  ! SCF, exactly as in their workflow.
  !
  ! CHANNEL SELECTION mirrors ldaV.f90 (init_lda_plus_v): acbn0_l_choice(nt,b)
  ! = 1 selects angular momentum l = b-1 of type nt; of the PP_CHI channels
  ! with that l the one with the highest principal quantum number (the valence
  ! one) is taken; selected channels are then ordered as they appear in the
  ! pseudopotential.  That ordering is what fixes the row/column layout of the
  ! printed U and V matrices, so it must match ldaV for the two codes to be
  ! compared block by block.
  !
  ! The first selected channel becomes stock QE's main manifold (Hubbard_l),
  ! the second its background (Hubbard_l2), the third Hubbard_l3 with backall.
  ! That is storage only: MODULE hp_acbn0_man rebuilds the single contiguous
  ! manifold ldaV uses.
  !
  ! WHAT IT DOES NOT DO.  It does not reimplement init_hubbard.  It sets the
  ! handful of module variables init_hubbard reads, then CALLS it, so ldim_u,
  ! ldmx_tot, ll, Hubbard_lmax, hubbard_occ, offsetU, offsetU_back(1) and
  ! nwfcU are built by stock QE exactly as they would have been.
  !
  ! The projectors are unaffected by any of this: orthoUwfc builds
  ! O^-1/2 S phi from the pseudo-atomic orbitals in the .upf, which lives in
  ! the .save.  Whether pw.x applied a Hubbard potential is irrelevant to them.
  !
  USE ions_base,        ONLY : ntyp => nsp
  USE ldaU,             ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_projectors, &
                               Hubbard_n, Hubbard_l, Hubbard_n2, Hubbard_l2,    &
                               Hubbard_n3, Hubbard_l3, Hubbard_U, Hubbard_U2,   &
                               nwfcU, wfcU, backall, init_hubbard
  USE uspp_param,       ONLY : upf
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : noncolin, npol
  USE wvfct,            ONLY : npwx
  USE io_global,        ONLY : stdout
  !
  IMPLICIT NONE
  !
  ! Placeholder occupancy, in Ry.  Never enters the physics: ACBN0 recomputes U
  ! from the occupations and the two-electron integrals and never reads
  ! Hubbard_U.  It exists only because is_hubbard is a "/= 0" test.
  !
  REAL(DP), PARAMETER :: U_FLAG = 1.0D-10
  !
  INTEGER, PARAMETER  :: lmax_acbn0 = 2        ! s, p, d: hp_acbn0_pao_3g has no f fits
  !
  INTEGER :: nt, nw, l, b, n, nsel, ios_n, best(0:5), nchar
  LOGICAL :: given
  CHARACTER(LEN=2) :: lab
  CHARACTER(LEN=1) :: spdf(0:5) = ['s','p','d','f','g','h']
  !
  given = ANY( acbn0_l_choice(1:ntyp,:) /= 0 )
  !
  IF ( acbn0 == 0 ) THEN
     IF ( given ) CALL errore( 'acbn0_manifold_setup', &
        'acbn0_l_choice is meaningful only for acbn0 = 1 or 2', 1 )
     RETURN
  ENDIF
  !
  ! acbn0 > 0: the manifold is declared HERE, and only here.  The pw.x ground
  ! state must be a plain (non-Hubbard) SCF.  The former alternative - a
  ! HUBBARD card with a placeholder U in the pw.x input, manifold read back
  ! from the XML - is no longer accepted: it made pw.x apply a Hubbard
  ! potential during the SCF, split a second channel into stock QE's
  ! background bookkeeping, and tripped the ortho-atomic force guard in
  ! dprojdtau_k whenever tprnfor was set.
  !
  IF ( lda_plus_u ) CALL errore( 'acbn0_manifold_setup', &
     'acbn0 > 0 needs a plain (non-Hubbard) pw.x ground state, but this one &
     &carries a HUBBARD card.  Remove the HUBBARD card from the pw.x input, &
     &re-run the SCF, and select the manifold with acbn0_l_choice in &INPUTHP.', 1 )
  !
  IF ( .NOT.given ) CALL errore( 'acbn0_manifold_setup', &
     'acbn0 > 0 needs the Hubbard manifold selected in &INPUTHP: set &
     &acbn0_l_choice(nt,b) = 1 (b = 1 s, 2 p, 3 d) for at least one type.', 1 )
  !
  DO nt = 1, ntyp
     DO b = 1, 6
        IF ( acbn0_l_choice(nt,b) /= 0 .AND. acbn0_l_choice(nt,b) /= 1 ) &
           CALL errore( 'acbn0_manifold_setup', &
           'acbn0_l_choice must be 0 or 1', nt )
        IF ( acbn0_l_choice(nt,b) == 1 .AND. b-1 > lmax_acbn0 ) &
           CALL errore( 'acbn0_manifold_setup', &
           'ACBN0 is implemented for s, p and d only (acbn0_l_choice b = 1, 2, 3)', nt )
        IF ( acbn0_rm_ehub_u(nt,b) /= 0 .AND. acbn0_rm_ehub_u(nt,b) /= 1 ) &
           CALL errore( 'acbn0_manifold_setup', &
           'acbn0_rm_ehub_u must be 0 or 1', nt )
        IF ( acbn0_rm_ehub_u(nt,b) == 1 .AND. acbn0_l_choice(nt,b) /= 1 ) &
           CALL errore( 'acbn0_manifold_setup', &
           'acbn0_rm_ehub_u set for a channel that is not in the manifold &
           &(acbn0_l_choice is 0 there)', nt )
     ENDDO
  ENDDO
  !
  DO nt = 1, ntyp
     !
     ! For every l, the PP_CHI channel with the highest principal quantum
     ! number - the valence one.  ldaV does the same with ICHAR(els(1:1)).
     !
     best(:) = 0
     nchar   = 0
     DO l = 0, lmax_acbn0
        nchar = 0
        DO nw = upf(nt)%nwfc, 1, -1
           IF ( upf(nt)%lchi(nw) /= l )      CYCLE
           IF ( upf(nt)%oc(nw)   <  0.0_DP ) CYCLE
           lab = upf(nt)%els(nw)
           IF ( nchar < ICHAR(lab(1:1)) ) THEN
              nchar   = ICHAR(lab(1:1))
              best(l) = nw
           ENDIF
        ENDDO
        IF ( acbn0_l_choice(nt,l+1) == 1 .AND. best(l) == 0 ) &
           CALL errore( 'acbn0_manifold_setup', &
           'acbn0_l_choice selects an l that the pseudopotential of this type &
           &has no (occupied, labelled) atomic wavefunction for', nt )
     ENDDO
     !
     ! Walk the pseudopotential in file order and hand the selected channels to
     ! stock QE in that order: main, then background, then background-1.
     !
     nsel = 0
     DO nw = 1, upf(nt)%nwfc
        l = upf(nt)%lchi(nw)
        IF ( l > lmax_acbn0 )                CYCLE
        IF ( acbn0_l_choice(nt,l+1) /= 1 )   CYCLE
        IF ( nw /= best(l) )                 CYCLE
        !
        lab = upf(nt)%els(nw)
        READ( lab(1:1), '(i1)', IOSTAT = ios_n ) n
        IF ( ios_n /= 0 .OR. n < 1 ) CALL errore( 'acbn0_manifold_setup', &
           'the PP_CHI label of the selected channel does not start with the &
           &principal quantum number (expected e.g. "3D")', nt )
        !
        nsel = nsel + 1
        SELECT CASE ( nsel )
        CASE ( 1 )
           Hubbard_n(nt)  = n ;  Hubbard_l(nt)  = l ;  Hubbard_U(nt)  = U_FLAG
        CASE ( 2 )
           Hubbard_n2(nt) = n ;  Hubbard_l2(nt) = l ;  Hubbard_U2(nt) = U_FLAG
        CASE ( 3 )
           Hubbard_n3(nt) = n ;  Hubbard_l3(nt) = l ;  backall(nt)    = .TRUE.
        CASE DEFAULT
           CALL errore( 'acbn0_manifold_setup', &
              'more than three channels selected for one type', nt )
        END SELECT
     ENDDO
     !
  ENDDO
  !
  Hubbard_projectors = 'ortho-atomic'
  lda_plus_u_kind    = 0
  lda_plus_u         = .TRUE.
  !
  ! Stock QE builds the rest: ldim_u, ldmx_tot, ll, hubbard_occ, offsetU,
  ! offsetU_back(1) and nwfcU (PW/src/ldaU.f90:202-527).
  !
  CALL init_hubbard( upf(1:ntyp)%psd, nspin, noncolin )
  !
  IF ( nwfcU <= 0 ) CALL errore( 'acbn0_manifold_setup', &
     'init_hubbard found no Hubbard states.', 1 )
  !
  ! allocate_wfc skipped wfcU on the read path because lda_plus_u was .false.
  ! at that point (PW/src/allocate_wfc.f90:34).  It is ours to allocate.
  !
  IF ( ALLOCATED(wfcU) ) DEALLOCATE( wfcU )
  ALLOCATE( wfcU(npwx*npol, nwfcU) )
  !
  WRITE( stdout, '(/5x,"Hubbard manifold selected in &INPUTHP (acbn0_l_choice)")')
  WRITE( stdout, '(5x,"(the pw.x run carried no HUBBARD card)")')
  DO nt = 1, ntyp
     IF ( ALL( acbn0_l_choice(nt,:) == 0 ) ) CYCLE
     WRITE( stdout, '(7x,"type ",i2,":")', ADVANCE='NO') nt
     WRITE( stdout, '(3x,i1,a1)', ADVANCE='NO') Hubbard_n(nt), spdf(Hubbard_l(nt))
     IF ( Hubbard_l2(nt) >= 0 .AND. Hubbard_U2(nt) /= 0.0_DP ) &
        WRITE( stdout, '("  + ",i1,a1)', ADVANCE='NO') Hubbard_n2(nt), spdf(Hubbard_l2(nt))
     IF ( backall(nt) ) &
        WRITE( stdout, '("  + ",i1,a1)', ADVANCE='NO') Hubbard_n3(nt), spdf(Hubbard_l3(nt))
     IF ( ANY( acbn0_rm_ehub_u(nt,:) == 1 ) ) THEN
        WRITE( stdout, '("     on-site U removed for:")', ADVANCE='NO')
        DO b = 1, lmax_acbn0+1
           IF ( acbn0_rm_ehub_u(nt,b) == 1 ) WRITE( stdout, '(1x,a1)', ADVANCE='NO') spdf(b-1)
        ENDDO
     ENDIF
     WRITE( stdout, * )
  ENDDO
  WRITE( stdout, '(7x,"projectors: ortho-atomic     nwfcU = ",i5,/)') nwfcU
  !
  RETURN
  !
END SUBROUTINE acbn0_manifold_setup

SUBROUTINE input_sanity()
  !-------------------------------------------------------------------------- 
  ! 
  ! This subroutine aims to gather all of the input sanity checks
  ! (features enabled in PWscf which are unsupported in HP).
  !
  USE klist,            ONLY : lgauss, ltetra, two_fermi_energies
  USE control_flags,    ONLY : gamma_only, tqr
  USE fixed_occ,        ONLY : tfixed_occ
  USE cellmd,           ONLY : lmovecell
  USE noncollin_module, ONLY : i_cons, noncolin
  USE mp_bands,         ONLY : nbgrp
  USE xc_lib,           ONLY : xclib_dft_is
  USE ldaU,             ONLY : lda_plus_u, Hubbard_projectors, lda_plus_u_kind, Hubbard_J0, &
                               Hubbard_V, is_hubbard_back
  !
  IMPLICIT NONE
  !
  IF (acbn0 < 0 .OR. acbn0 > 2) CALL errore ('hp_readin', &
     ' acbn0 must be 0 (DFPT, default), 1 (ACBN0 U) or 2 (eACBN0 U+V) ',1)
  !
  IF (conv_thr_chi <= 0.D0) CALL errore ('hp_readin', ' Wrong conv_thr_chi ',1)
  !
  IF (nq1.LE.0 .OR. nq2.LE.0 .OR. nq3.LE.0) &
       CALL errore('hp_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  IF (start_q <= 0 ) CALL errore('hp_readin', ' Wrong start_q ',1)
  !
  IF (compute_hp .AND. ANY(perturb_only_atom(:))) &
     CALL errore ('hp_readin', 'compute_hp and perturb_only_atom are not allowed to be true together', 1)
  !
  IF ( ANY(Hubbard_V(:,:,2).NE.0.d0) .OR. &
       ANY(Hubbard_V(:,:,3).NE.0.d0) .OR. &
       ANY(Hubbard_V(:,:,4).NE.0.d0) ) &
     CALL errore ('hp_readin', 'The HP code does not support DFT+U+V with the background', 1)
  !
  ! Background channels are a DFPT limitation, not an ACBN0 one: the eACBN0
  ! kernels index the manifold through ll(m,nt) and handle several channels
  ! per type from stage V4 on (see HP/src/eACBN0.f90, MODULE hp_acbn0_man).
  IF (ANY(is_hubbard_back(:)) .AND. acbn0 == 0) CALL errore ("hp_readin", &
          &" Two (or more) Hubbard channels per atomic type is not implemented", 1)
  !
  IF (ANY(Hubbard_J0(:).NE.0.d0)) &
     CALL errore ('hp_readin', 'Hubbard_J0 /= 0 is not allowed.', 1)
  !
  IF (determine_q_mesh_only .AND. .NOT.ANY(perturb_only_atom(:))) &
     CALL errore ('hp_readin', 'determine_q_mesh_only can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
  !
  IF (sum_pertq .AND. .NOT.ANY(perturb_only_atom(:))) &
     CALL errore ('hp_readin', 'sum_pertq can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
  !
  IF (niter_max.LT.1 .OR. niter_max.GT.maxter) &
     CALL errore ('hp_readin', ' Wrong niter_max ', 1)
  !
  DO iter = 1, niter_max
     IF ( alpha_mix(iter).LT.0.D0 .OR. alpha_mix(iter).GT.1.D0 ) &
        CALL errore ('hp_readin', ' Wrong alpha_mix ', iter)
  ENDDO
  !
  IF (num_neigh.LT.1) CALL errore('hp_readin','Not allowed value of num_neigh',1)
  !
  IF (lmin.LT.0 .OR. lmin.GT.3) CALL errore('hp_readin','Not allowed value of lmin',1)
  !
  IF (nmix.LT.1) CALL errore ('hp_readin', ' Wrong nmix ', 1) 
  !
  ! Three of the guards below are preconditions of the DFPT linear-response
  ! branch, not of the code as a whole.  The ACBN0 branch (acbn0 > 0) reads a
  ! converged ground state and evaluates a functional of it; it solves no
  ! Sternheimer equation, so it does not inherit them.  Each is relaxed for
  ! exactly the cases where it does not apply, and the ACBN0 preconditions are
  ! then enforced separately in hp_acbn0_check / hp_acbn0_check_v.
  !
  ! ltetra: DFPT needs smearing or fixed occupations to build the response.
  ! ACBN0 only reads wg (weights x occupations), which the tetrahedron method
  ! supplies perfectly well.  The authors' own reference input
  ! (example/MnO.acbn0.in) uses occupations = 'tetrahedra_opt'.
  !
  IF (ltetra .AND. acbn0 == 0) CALL errore ('hp_readin', &
     'HP with tetrahedra is not supported', 1)
  !
  IF (gamma_only) CALL errore('hp_readin',&
     & 'Cannot start from pw.x data file using Gamma-point tricks',1)
  !
  ! lda_plus_u: required by DFPT and by acbn0 = 1, which projects with ldaU's
  ! wfcU / offsetU / Hubbard_l.  Since stage V3 acbn0 = 2 reads the manifold
  ! from ldaU as well, so it needs a HUBBARD card too; the exemption below is
  ! left over from stage V2 and is backstopped by hp_acbn0_check_v.
  !
  IF (.NOT.lda_plus_u .AND. acbn0 /= 2) CALL errore('hp_readin',&
     & 'The HP code can be used only on top of DFT+Hubbard (i.e. when the HUBBARD card is used in pw.x)',1)
  !
  ! acbn0 = 2 needs an inter-site cutoff d_IJ.  Since V3 the +V branch owns
  ! its own neighbour table, so the cutoff comes from &INPUTHP and from
  ! nowhere else.  The old fallback to ehub_nn_distance in the XML is gone:
  ! that field exists only for a ground state made by the Lee-Son pw.x, and
  ! this code must build against an unpatched QE 7.3 where ldaV is absent.
  !
  IF (acbn0 == 2 .AND. nn_distance < 0.0D0) &
     CALL errore('hp_readin', 'acbn0 = 2 needs an inter-site cutoff: set &
     &nn_distance (in Angstrom) in &INPUTHP.',1)
  !
  IF (nn_distance > 20.0D0) CALL errore('hp_readin', &
     ' nn_distance is unreasonably large (> 20 Angstrom) ',1)
  !
  ! lda_plus_u_kind = 1 is the Liechtenstein DFT+U+J formulation.  It changes
  ! the DFPT response and the ldaU projectors.  The acbn0 /= 2 exemption
  ! below is stage-V2 residue, backstopped by hp_acbn0_check_v.
  !
  IF (lda_plus_u_kind.EQ.1 .AND. acbn0 /= 2) CALL errore("hp_readin", &
     & ' The HP code does not support the Liechtenstein formulation of DFT+U',1)
  !
  IF (Hubbard_projectors.NE."atomic" .AND. Hubbard_projectors.NE."ortho-atomic" &
      .AND. acbn0 /= 2) &
     CALL errore("hp_readin", &
     " The HP code for this Hubbard_projectors type is not implemented",1)
  !   
  IF (lmovecell) CALL errore('hp_readin','The HP code is not working after vc-relax',1)
  !
  IF (nbgrp > 1) CALL errore('hp_readin', &
     & 'band parallelization is not implemented in HP',1)
  !
  IF (i_cons /= 0) CALL errore('hp_readin',&
     & 'The HP code with constrained magnetization is not yet available',1)
  !
  IF (two_fermi_energies .AND. (ltetra .OR. lgauss)) CALL errore('hp_readin', &
     & 'The HP code with two Fermi energies is not available for metals',1)
  !
  IF (tqr) CALL errore('hp_readin',&
     & 'The HP code with Q in real space is not supported',1)
  !
  IF (tfixed_occ) CALL errore('hp_readin', &
     & 'The HP code with arbitrary occupations not tested',1)
  !
  IF ( xclib_dft_is('meta') ) CALL errore('hp_readin',&
     'The HP code with meta-GGA functionals is not yet available',1)
  !
  IF ( xclib_dft_is('hybrid') ) CALL errore('hp_readin',&
     'The HP code with hybrid functionals is not yet available',1)
  !
  RETURN
  !
END SUBROUTINE input_sanity
  !
END SUBROUTINE hp_readin
