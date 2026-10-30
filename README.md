# qe-7.3_thiswork -- ACBN0 / eACBN0 Hubbard parameters inside hp.x

Quantum ESPRESSO 7.3 with ONE addition: `hp.x` can evaluate the on-site Hubbard
U (ACBN0) and the inter-site V (extended ACBN0, eACBN0) of Lee and Son as a
one-shot post-processing step on a plain DFT ground state, next to the official
DFPT linear-response route.  No source file outside `HP/src/` differs from
official QE 7.3; the only other differences are a two-line `readme.txt` and
build scaffolding (section 3).  The official QE README is kept as
`README_QE.md`.

This document covers

1. the three code trees this project talks about and the words used for them
2. what dev hp.x does, and what it deliberately does not do
3. differences from official QE 7.3
4. differences from Lee and Son's patched pw.x
5. the `&INPUTHP` tags, their meaning, and the mapping to Lee and Son's tags
6. preconditions and limitations
7. output files, field by field
8. build and run
9. the `example/` directory
10. known rough edges in the current code
11. references and licences

Papers: ACBN0, Agapito, Curtarolo, Buongiorno Nardelli, Phys. Rev. X 5, 011006
(2015).  eACBN0 (on-site U plus inter-site V), Lee and Son, Phys. Rev. Research
2, 043410 (2020).  DFPT Hubbard parameters, Timrov, Marzari, Cococcioni,
Phys. Rev. B 98, 085127 (2018), 103, 045141 (2021); HP code, Comput. Phys.
Commun. 279, 108455 (2022).

---------------------------------------------------------------------------

## 1. Terminology: the three trees

| name         | tree                                | what it is                                                                                   |
|--------------|-------------------------------------|----------------------------------------------------------------------------------------------|
| **official** | `~/projects/QE/qe-7.3_official`     | pristine QE 7.3 (branch `official`). Hubbard U from DFPT via `hp.x`. No ACBN0 anywhere.        |
| **stock**    | `~/projects/QE/qe-7.3_eACBN0`       | official + Lee and Son's `qe-7.3_ehub_uv.diff` = branch `pw/eACBN0`. U and V recomputed inside the `pw.x` SCF loop. Never touches `HP/`. The working tree there is checked out on `pw/eACBN0-rawprint`, which adds Yujin's `eACBN0_hub.dat` writer and `eacbn0_debug` on top; pristine Lee and Son is `git show origin/pw/eACBN0:<file>`. |
| **dev**      | `~/projects/QE/qe-7.3_thiswork`     | this tree (branch `dev`, HEAD `a1ad49f`). official + the stock ACBN0/eACBN0 kernels transplanted into `HP/src` as a one-shot evaluator. |

All three are branches of the same remote (`yjkang1/qe-7.3`).  "stock" always
means Lee and Son's pw.x; it never means "unmodified upstream".

---------------------------------------------------------------------------

## 2. What dev hp.x does

The `&INPUTHP` tag `acbn0` selects the estimator:

```
acbn0 = 0   DFPT linear response (default). Official hp.x, code path untouched.
acbn0 = 1   on-site U      from the ACBN0 pseudo-hybrid functional   (one-shot)
acbn0 = 2   on-site U + V  from the eACBN0 functional                (one-shot)
```

Workflow for `acbn0 > 0`:

```
pw.x  scf.in      plain LDA/PBE SCF, nspin = 2, NO HUBBARD card, disk_io >= 'low'
hp.x  hp.in       acbn0 = 1 or 2, prefix/outdir of that .save,
                  manifold chosen with acbn0_ehub_l_choice(type, l+1) = 1,
                  (acbn0 = 2) cutoff nn_distance in Angstrom
  -> outdir/eACBN0_hub.dat                       human-readable U, V, energies, occupations
  -> outdir/<prefix>.save/ehub_uv.txt, ...       dumps in Lee and Son's file formats
  -> stdout                                      summary tables
```

`hp.x` reads the converged wavefunctions, builds ortho-atomic projectors
`wfcU = O^-1/2 S phi` with official QE's own `orthoUwfc`, forms the occupation
matrices `n^I_{mm'}` and the ACBN0-renormalised density matrices `P`,
evaluates the two-electron integrals with Rys quadrature, and prints

```
U^I(l1,l2) = A0_U / HF_U  -  A0_J / HF_J             (J term dropped for an s-s block)
  m1 in l1, m2 in l2;  m3 runs over the channel of m1, m4 over the channel of m2
  A0_U = sum [sum_{s,s'} P_{m3m1,s} P_{m4m2,s'}] (m3 m1 | m4 m2)    s-s block: only s /= s' terms
  A0_J = sum [sum_{s}    P_{m3m1,s} P_{m4m2,s} ] (m3 m2 | m4 m1)
  HF_U = sum_{m1,m2} [sum_{s,s'} n_{m1m1,s} n_{m2m2,s'}]            m1 = m2: only s /= s' terms
  HF_J = sum_{m1/=m2} sum_{s} n_{m1m1,s} n_{m2m2,s}
V^IJ = the inter-site analogue for every pair (I in cell 0, J in cell n) with |R_J + n - R_I| <= d_IJ
```

The Coulomb integrals `(..|..)` are NOT taken over the pseudopotential's own
`PP_CHI` orbitals.  They use a fixed per-element table of 3-Gaussian fits
(`hp_acbn0_pao_3g`, keyed on the element symbol and l), exactly as stock does.
They are therefore the same for every pseudopotential of a given element;
only the projectors and occupations depend on the pseudopotential.

(Ry internally, eV in the printed tables.)  Within one (l1, l2) channel block U is a
single number, which is why the summary tables print one U per channel.

exactly as the stock kernels do on every SCF step.  Then it stops.  The DFPT
machinery (`hp_init`, `hp_summary`, the q loop, `hp_postproc`) is never
entered when `acbn0 > 0`.

What it deliberately does NOT do:

- **No self-consistency.**  The stock code recomputes U and V every SCF
  iteration, mixes them, feeds them back into the Hamiltonian, and converges
  density and parameters together.  Dev evaluates the functional ONCE on a
  ground state that was produced WITHOUT any Hubbard potential.  The numbers
  are U_0 = F[rho_DFT].  They reproduce the stock *raw, unmixed* values of the
  first ACBN0 call from that density, not the converged ones.
- **No Hubbard potential, forces, stress, or band structure.**  `hp.x` writes
  parameters; it does not apply them.  Nothing in this tree can run DFT+U+V
  with the eACBN0 parameters, and nothing writes a `HUBBARD` card for official
  pw.x.
- **No mixing knob, no convergence threshold, no restart files to read.**
  Those are SCF-loop concepts.

Why it exists anyway: it puts the ACBN0 estimator in the same binary as the
DFPT estimator, on top of unmodified PW/, so the two can be compared without
building two different pw.x codes.  Note one caveat the code comments gloss
over: the two routes need DIFFERENT pw.x ground states.  `acbn0 = 0` requires
a `HUBBARD` card in the pw.x input (official rule), while `acbn0 > 0` refuses
one (section 5.3).  Sharing one DFT+U ground state between the two routes is
an open project TODO, not implemented in `a1ad49f`.  The earlier route
(placeholder U in a pw.x `HUBBARD` card, manifold read back from the XML) was
dropped because pw.x then applied a Hubbard potential during the SCF, pushed a
second channel into official QE's background bookkeeping, and tripped the
ortho-atomic force guard in `dprojdtau_k` whenever `tprnfor` was set.

To APPLY one-shot parameters, the only consumer is stock pw.x: run it on a
copy of the same `.save` with `lda_plus_v = .true.`, `lacbn0 = .false.`,
`read_ehub_uv_file = .true.` and the manifold and cutoff matching the hp.x run
(section 7.3).  This path has not been tested.

---------------------------------------------------------------------------

## 3. Differences from official QE 7.3

`git diff --stat origin/official..dev`: 15 paths, 8514 insertions, 10
deletions.  Every source change is in `HP/src/`:

| file                       | change                                                                                 |
|----------------------------|----------------------------------------------------------------------------------------|
| `HP/src/eACBN0.f90`        | NEW, 4689 lines. Modules `hp_acbn0_nn` (neighbour table), `hp_acbn0_man` (manifold view); drivers `hp_acbn0`, `hp_acbn0_u`, `hp_acbn0_v`; checks; the two kernels `hp_acbn0_ns_a0` (U) and `hp_acbn0_ns_e0_1` (U+V); symmetrisers; writers. |
| `HP/src/hp_acbn0_pao.f90`  | NEW. Stock `pao_3g.f90` renamed: 3-Gaussian fits of the s, p, d pseudo-atomic orbitals per element. |
| `HP/src/hp_acbn0_eri.f90`  | NEW. Stock `compute_eri.f90` renamed: Fortran side of the two-electron integrals, Hartree to Ry factor. |
| `HP/src/hp_acbn0_crys.c/.h`| NEW. Stock `crys.c/.h` renamed: Rys-quadrature ERIs from PyQuante (modified BSD, see section 11). |
| `HP/src/hp_readin.f90`     | +294/-8. Four new namelist tags and defaults; new `acbn0_manifold_setup`; `hp_ns_trace` skipped for `acbn0 /= 0`; five DFPT guards conditioned on `acbn0`; three new checks (`acbn0` in 0..2; `nn_distance` required for `acbn0 = 2`; `nn_distance <= 20`). |
| `HP/src/hpcom.f90`         | +35. Declarations of `acbn0`, `nn_distance`, `acbn0_ehub_l_choice`, `acbn0_remove_ehub_u` in `MODULE ldaU_hp`. |
| `HP/src/hp_main.f90`       | +22/-1. `IF (acbn0 > 0) THEN CALL hp_acbn0(); GO TO 104` before `hp_init`. |
| `HP/src/hp_bcast_input.f90`| +6/-1. Broadcast of the four new tags. |
| `HP/src/hp_print_clock.f90`| +2. Clocks `hp_acbn0`, `acbn0_ns_a0`. |
| `HP/src/Makefile`          | +4. The four new objects in `HPOBJS`. |
| `readme.txt`               | 2-line note, superseded by this file. |
| `external/devxlib`, `external/mbd`, `install/libcuda_devxlib` | submodule gitlinks at the recorded hashes and an empty build stamp; not code. See section 8. |

`PW/`, `Modules/`, `LR_Modules/`, `PP/`, `UtilXlib/`, `upflib/`, `HP/Doc`,
`HP/examples`, `HP/CMakeLists.txt` are untouched.

**The DFPT route (`acbn0 = 0`) is behaviourally identical to official hp.x.**
With an input valid for official hp.x the dev binary differs only by: four
broadcasts of default values; one `ANY()` over an all-zero array in
`acbn0_manifold_setup` followed by `RETURN`; three checks on new tags at
their defaults; five guards whose extra `.AND. acbn0 == 0` / `.AND. acbn0 /= 2`
clause is true; two `IF` tests on `acbn0` (the one around `hp_ns_trace`, true,
and the one in `hp_main`, false); two `print_clock` calls on clocks that were
never started, which print nothing.  No default of any pre-existing tag
changed.  Numbers, stdout and files of a DFPT run are the official ones.

Guards relaxed for `acbn0 > 0` only:

- tetrahedra occupations are accepted (DFPT needs smearing; ACBN0 only reads
  weights);
- more than one Hubbard channel per type is accepted (the ACBN0 kernels index
  the manifold through `ll(m,nt)`; DFPT does not).

The other three conditioned guards (`lda_plus_u`, `lda_plus_u_kind`,
`Hubbard_projectors`, exempted for `acbn0 = 2`) are dead code: the manifold
setup forces those variables before the guards run (section 10).

Two side effects of the footprint:

- The CMake build does not know the new sources (`HP/CMakeLists.txt` is
  official); `hp_main` would fail to link.  Build with `make`.
- `HP/Doc/INPUT_HP.def` documents none of the new tags.  Section 5 of this file
  is the only documentation.

---------------------------------------------------------------------------

## 4. Differences from Lee and Son's stock pw.x

Stock patches 60 files, 59 of them in `PW/src`, `Modules`, `PP`, `UtilXlib`
plus `external/.gitignore` (15 new files:
`ldaV.f90`, `new_ns_v.f90`, `hubbard_nns.f90`, `pao_3g.f90`, `compute_eri.f90`,
`crys.c/.h`, `v_hubbard_v.f90`, `force_hub_v.f90`, `stres_hub_v.f90`, ...) and
never touches `HP/`.  Dev patches 11 source files, all in `HP/src`, and never
touches `PW/`.

### 4.1 Where and when the parameters are computed

| aspect                       | stock (pw.x)                                                                 | dev (hp.x)                                                                      |
|------------------------------|------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| when                         | every SCF iteration (`new_ns_a0` / `new_ns_e0_1` from `sum_band`)             | once, after reading a converged plain ground state                               |
| feedback                     | U, V enter the Hamiltonian via `v_hubbard_v`; density and parameters converge together | none. U, V are printed and the program stops                                      |
| mixing                       | `U_n = (1 - m) U_{n-1} + m U_raw`, `m = ehub_mixing` (default 0.7) weights the RAW term | none (the .dat prints `mixing     : 1.000  (one-shot; no SCF history)` only so the file lines up with stock's) |
| convergence                  | `dr2 < conv_thr` AND `|dE_hub| < ehub_conv_thr`                                | not applicable                                                                   |
| ground state                 | SCF with `lda_plus_v`, `lacbn0`, manifold in `&system`                        | plain SCF from official pw.x, NO HUBBARD card; manifold in `&INPUTHP`            |
| forces / stress with U, V    | yes (`force_hub_v`, `stres_hub_v`)                                            | no                                                                               |
| restart with fixed U, V      | `read_ehub_uv_file`, `read_ehub_ns_file`                                      | no (dev only writes those files)                                                 |
| noncollinear / SOC           | only for `acbn0_type = 2` (`new_ns_e0_nc_1`); the type 1 and 3 noncollinear calls are commented out | no, `nspin = 2` required                                                         |
| `acbn0_type = 3`             | anisotropic / orbital-dependent U+V (`new_ns_e0_2`), undocumented in the stock README | not implemented, `acbn0 = 3` is refused                                          |
| manifold bookkeeping         | `MODULE ldaV` (own copy of projectors, offsets, `ehub_l`)                     | official `MODULE ldaU` via `init_hubbard`; contiguous view rebuilt in `hp_acbn0_man` |
| channels per type            | several                                                                      | several, up to three (official main + two background blocks)                     |
| neighbour table              | `hubbard_nns` in pw.x from `ehub_nn_distance`                                 | `hp_acbn0_build_nn` in hp.x from `nn_distance`; same algorithm, so d_IJ can be swept without re-running pw.x |
| per-iteration history file   | `outdir/eACBN0_hub.dat`, one block per kernel call (Yujin's addition on branch `pw/eACBN0-rawprint`, not pristine Lee and Son) | `outdir/eACBN0_hub.dat`, one block                                                |

### 4.2 Kernel identity

The arithmetic is a verbatim port.  After renaming the `ldaV` identifiers
(`ehub_lm -> lm`, `ehub_l -> lch`, `ofstv(na)+m -> ofs(m,na)`, `wfcv -> wfcU`,
`ehub_u -> hub_u`, ...) and `diff -w`, `hp_acbn0_ns_a0` equals stock
`new_ns_a0`, `hp_acbn0_ns_e0_1` equals stock `new_ns_e0_1`, the three
symmetrisers equal `symmetrize_n4/n4m/n5`, and `hp_acbn0_pao.f90`,
`hp_acbn0_eri.f90`, `hp_acbn0_crys.c/.h` equal their stock originals line for
line once the `hp_acbn0_` prefix on their own identifiers is stripped
(`ERI`, `compute_ERI`, `coulomb_repulsion`, `contr_coulomb`,
`exponents_and_normalization_constants`, the `crys.h` include).  The only edits:

- the mixing statements are replaced by the raw assignment;
- the `stop_acbn0` short-circuit is removed;
- the `remove_ehub_u` mask is applied as a separate statement after the ratio
  (same product over both channel indices, same result);
- the stdout/dump writers are called from the driver instead of the kernel;
- two guards added that stock lacks: `IF (lm(nt) == 0) CYCLE` before the
  symmetrisers write `l_separation` (stock writes index 0 for a non-Hubbard
  type), and `nnc = 0` for a symmetry image outside the neighbour set (stock
  leaves it undefined);
- the neighbour-cell count uses `r^2 > 1e-8` everywhere; stock's count in
  `input.f90` uses `1e-10` while its `hubbard_nns` uses `1e-8`.  Irrelevant
  unless two atoms sit within 1e-4 alat of each other.

Numerical validation (`/gpfs/yjkang/project-01/{1,2,3,5}-*`, Ni, GaAs, ZnO,
CdO, CrI3, PBE and LDA, U and U+V, 20 cases): dev hp.x run on a copy of the
stock converged `.save` reproduces the stock RAW (unmixed) U and V of the last
SCF iteration to the last printed digit of the `ehub_uv.txt` dumps,
at most 1.4e-14 eV = 1e-15 Ry (`2-print_wo_mixing`), and the per-iteration
`.dat` tables at printed precision (`3-hub_dat_files`).

The gap to the stock MIXED values is the mixing lag
`(1-m)/m * (U_N - U_{N-1})` of the last stock iteration, so it depends on how
settled U was when the stock SCF stopped.  It is about 1e-9 Ry (1e-8 eV) for
Ni and GaAs (`m` = 0.7 / 0.8, `ehub_conv_thr = 1e-8`) and up to 2.7e-5 Ry
(0.37 meV) for CdO (`m = 0.3`, `ehub_conv_thr = 1e-6`).  Those are the
third-to-fourth-decimal deviations shown on 2026-09-17; against the RAW values
they vanish.  For the same reason `ehub_u0.txt` and the one-shot `E(+U)` match
stock exactly only where the lag is below printed precision (Ni, GaAs); for
ZnO, CdO and CrI3 they differ at the 1e-7 to 1e-5 Ry level.

### 4.3 A convention trap in the stock README

Lee and Son's README states `U_new = (1 - ehub_mixing) dU + U_old`.  The code
(`new_ns_v.f90`, every mixing site) does `U_n = (1 - m) U_{n-1} + m U_raw`:
`ehub_mixing = 1` is the raw one-shot value and `0` freezes U.  Dev follows the
code.

---------------------------------------------------------------------------

## 5. Input: the `&INPUTHP` tags

### 5.1 New tags

| tag                          | type, default               | meaning                                                                                                           |
|------------------------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| `acbn0`                      | integer, 0                  | 0 DFPT (official), 1 ACBN0 on-site U, 2 eACBN0 on-site U + inter-site V. Anything else: error.                      |
| `nn_distance`                | real, Angstrom, -1 (unset)  | inter-site cutoff d_IJ. REQUIRED for `acbn0 = 2`. Unused for 0 and 1, but the `> 20` check still runs and aborts for any `acbn0`. `0.0` is legal and gives no pairs (U only through the e0 kernel). |
| `acbn0_ehub_l_choice(a,b)`   | integer array (10,6), 0/1, all 0 | put channel `b` of atomic type `a` into the Hubbard manifold. `a` = position in `ATOMIC_SPECIES`, `b = l + 1` (1 s, 2 p, 3 d). REQUIRED for `acbn0 > 0` (at least one 1), forbidden for `acbn0 = 0`. |
| `acbn0_remove_ehub_u(a,b)`   | integer array (10,6), 0/1, all 1 | 1 keep (default), 0 REMOVE the on-site U of channel `b` of type `a`. Lee and Son's `remove_ehub_u` convention. May be 0 only where `acbn0_ehub_l_choice(a,b) = 1`. |

The array shape is `(ntypx = 10, 6)`; an index `a > 10` or `b > 6` fails at
the namelist read.

`acbn0_ehub_l_choice` in detail (routine `acbn0_manifold_setup`, `hp_readin.f90`):

- For each selected `l`, the PP_CHI channel of that `l` with `oc >= 0` and the
  highest principal quantum number (first character of the `PP_CHI` label) is
  taken.  That is the valence channel, the same choice stock makes with
  `ehub_l_choice`.  A stock pseudo whose highest-n channel has `oc < 0` would
  be silently dropped by stock; dev takes the next occupied one instead.
- Selected channels are ordered as they appear in the pseudopotential file.
  This fixes the row and column layout of every printed U and V matrix and
  matches stock, so the two codes' matrices compare block by block.  Within a
  channel the m index is QE's real-spherical-harmonic order: p = pz, px, py;
  d = dz2, dxz, dyz, dx2-y2, dxy.  The GaAs 4x4 blocks are 4s, 4pz, 4px, 4py.
- The first selected channel becomes official QE's `Hubbard_l`, the second
  `Hubbard_l2` (background), the third `Hubbard_l3` with `backall`.  A
  placeholder `Hubbard_U = 1e-10 Ry` is set only to make `is_hubbard` true; it
  never enters any number.  Then official `init_hubbard` builds `ldim_u`,
  `ll`, `offsetU`, `nwfcU`, and `wfcU` is allocated.
- `Hubbard_projectors` is forced to `'ortho-atomic'`, `lda_plus_u_kind` to 0.
- Consequences inherited from official QE: the `PP_CHI` label must have the
  form `<n><letter>` (e.g. `3D`), and a channel with ZERO occupation in the
  pseudopotential is refused by `offset_atom_wfc` ("Hubbard manifold with zero
  occupations is not allowed").  Stock has no such refusal, so e.g. a GBRV Ga
  with `4P` at `oc = 0.00` runs in stock but not in dev.
- Only s, p, d (`b = 1..3`).  `b = 4` (f) is refused because the Gaussian fits
  have no f shell; stock's `b = 5` semicore-d slot has no dev equivalent.

`acbn0_remove_ehub_u` in detail: applied in one place per route, after the
ACBN0 ratio is formed and before the (m1,m2)/(m2,m1) averaging:

```
hub_u(na,m1,m2) = hub_u_(na,m1,m2) * mask(nt, lch(m1,nt)) * mask(nt, lch(m2,nt))
```

so with the mask 0 on channel c the whole row block and column block touching
c vanish, including the cross blocks with a kept channel.  The channel stays in
the manifold: projectors, occupations, renormalisation weights N^IJ, ERI sums
and every inter-site V are unaffected.  Identical to stock.

### 5.2 Official tags under `acbn0 > 0`

Still meaningful: `prefix`, `outdir` (where the `.save` is read, where
`eACBN0_hub.dat` is written).  `hp.x` also creates `outdir/HP/`, empty on this
route.

Read but ignored: `iverbosity` (no effect in `eACBN0.f90`; the examples set
it to 1 out of habit), `max_seconds` (never polled on this route), and every
DFPT tag: `nq1..3`, `skip_equivalence_q`, `conv_thr_chi`,
`find_atpert`, `skip_type`, `skip_atom`, `equiv_type`, `perturb_only_atom`,
`start_q`, `last_q`, `sum_pertq`, `compute_hp`, `determine_num_pert_only`,
`determine_q_mesh_only`, `docc_thr`, `thresh_init`, `ethr_nscf`, `niter_max`,
`alpha_mix`, `nmix`, `num_neigh`, `lmin`, `rmax`, `dist_thr`, `background`,
`disable_type_analysis`.  The range checks of `conv_thr_chi`, `nq1..3`,
`start_q`, `niter_max`, `alpha_mix`, `num_neigh`, `lmin`, `nmix` and the
`compute_hp` / `perturb_only_atom` / `determine_q_mesh_only` / `sum_pertq`
consistency checks still run, so a nonsense value there still aborts.  The
other listed tags are not validated at all.

### 5.3 Minimal inputs

```
 &INPUTHP                              ! on-site U, Ni 3d
   prefix = 'Ni',  outdir = './out'
   acbn0  = 1
   acbn0_ehub_l_choice(1,3) = 1
 /

 &INPUTHP                              ! U + V, GaAs s+p on both atoms, s on-site U removed
   prefix = 'GaAs',  outdir = './out'
   acbn0  = 2
   nn_distance = 3.0                   ! Angstrom, mandatory for acbn0 = 2
   acbn0_ehub_l_choice(1,1) = 1,  acbn0_ehub_l_choice(1,2) = 1
   acbn0_ehub_l_choice(2,1) = 1,  acbn0_ehub_l_choice(2,2) = 1
   acbn0_remove_ehub_u(1,1) = 0,  acbn0_remove_ehub_u(2,1) = 0
 /
```

The pw.x input is an ordinary official SCF with `nspin = 2` (use
`starting_magnetization = 0.0` for a non-magnetic system), `disk_io = 'low'`
or higher, and NO `HUBBARD` card.  `occupations = 'fixed'`, smearing and
tetrahedra are all accepted; the smearing in `example/GaAs` copies the stock
reference input so the numbers are comparable, and `'fixed'` has not been
validated against stock.  If the `.save` carries a `HUBBARD` card, hp.x stops
with this single line:

```
acbn0 > 0 needs a plain (non-Hubbard) pw.x ground state, but this one carries a HUBBARD card.  Remove the HUBBARD card from the pw.x input, re-run the SCF, and select the manifold with acbn0_ehub_l_choice in &INPUTHP.
```

A stock ground state (`lda_plus_v` / `lacbn0` in `&system`, no `HUBBARD`
card) counts as a plain SCF and is accepted.  Translate `ehub_l_choice` to
`acbn0_ehub_l_choice`, `ehub_nn_distance` to `nn_distance`, `remove_ehub_u` to
`acbn0_remove_ehub_u`.  WARNING: hp.x then overwrites stock's own dumps in that
`prefix.save` (section 7.3).  Copy the `.save` first, or rename stock's files
(the validation runs renamed them `*.ref.txt`).

### 5.4 Lee and Son's `&system` tags and their dev equivalents

| stock tag (type, default)             | what stock does with it                                                                 | dev                                                                           |
|---------------------------------------|-----------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
| `lda_plus_v` (L, .false.)             | master switch of the ldaV machinery: manifold, neighbours, V potential                  | implicit in `acbn0 > 0`. No fixed-U/V SCF mode exists in dev.                 |
| `lacbn0` (L, .false.)                 | recompute U, V from ACBN0 every SCF step; adds the Hubbard-energy convergence test       | implicit in `acbn0 > 0` (one-shot)                                            |
| `acbn0_type` (I, 1)                   | 1 U only, 2 U + V, 3 anisotropic U + V                                                  | `acbn0 = 1`, `acbn0 = 2`; type 3 not implemented                              |
| `ehub_nn_distance` (R, 3.0 A)         | cutoff d_IJ for the neighbour list                                                      | `nn_distance` (no default; mandatory for `acbn0 = 2`)                          |
| `ehub_l_choice(a,b)` (I(10,6), 0)     | manifold, `b = l + 1`; `b = 5` semicore d                                               | `acbn0_ehub_l_choice(a,b)`, `b = 1..3` only                                   |
| `remove_ehub_u(a,b)` (I(10,6), 1)     | 1 keep / 0 remove on-site U, product mask after the ratio                               | `acbn0_remove_ehub_u(a,b)`, same convention, same formula                     |
| `remove_ehub_v(a,b)` (I(10,6), 1)     | read, copied, never consumed by any routine                                             | none (dead in stock too)                                                      |
| `ehub_mixing` (R, 0.7)                | linear mixing, weight on the RAW term                                                   | none (one-shot)                                                               |
| `ehub_conv_thr` (R, 1e-8 Ry)          | SCF stops only when `dr2 < conv_thr` AND the Hubbard-energy change is below this         | none                                                                          |
| `ehub_pot_fix` (L, .false.)           | hold the occupation matrices of the potential at the input ones during the SCF           | none                                                                          |
| `read_ehub_uv_file`, `read_ehub_ns_file` (L) | read `ehub_uv.txt` / `occupation_uv.txt` from `prefix.save`                       | none; dev only writes these files                                             |
| `stop_acbn0` (L, .false.)             | freeze U, V once their change is below `sqrt(ehub_conv_thr)`; acts for `acbn0_type` 2 and 3 only (no effect on a U-only run) | none                                                                          |
| `eacbn0_debug` (I, 0)                 | 1 also writes `eACBN0_hub_raw.dat` and `ehub_uv_raw.txt` (Yujin's rawprint branch only) | none; dev values are raw by construction                                     |
| `ehub_u`, `ehub_v` (R arrays, eV)     | starting / fixed U and V                                                                | none                                                                          |

---------------------------------------------------------------------------

## 6. Preconditions and limitations (as enforced by the code)

Enforced with an error message:

- `acbn0` in 0..2.
- `acbn0 > 0`: no `HUBBARD` card in the pw.x run; at least one
  `acbn0_ehub_l_choice = 1`; values 0/1 only; `b <= 3`; a selected `l` must
  exist as a labelled `PP_CHI` channel with `oc >= 0`, and official
  `offset_atom_wfc` then refuses `oc = 0`; `acbn0_remove_ehub_u = 0` only
  inside the manifold.
- `nspin = 2` (collinear, spin-polarised).  `nspin = 1` and noncollinear are
  refused; the kernels hard-code two spin slots.
- not `gamma_only`.
- `acbn0 = 2`: `nn_distance >= 0`.  Any `acbn0`: `nn_distance <= 20`.
- plus the official ground-state guards that apply to every hp.x run: not
  after vc-relax, `nbgrp = 1`, no constrained magnetisation, no
  `two_fermi_energies` with smearing or tetrahedra, no `tqr`, no
  `occupations = 'from_input'`, no meta-GGA, no hybrid functional.

Required but NOT caught cleanly: the collected wavefunctions must be in
`<prefix>.save` (`disk_io` at least `'low'` in pw.x).  hp.x has its own check
for missing `prefix.wfc*` files, but it cannot fire: `read_file` has already
created those files.  A `disk_io = 'nowf'` ground state fails later, with
read_file's "Wavefunctions in collected format not available" notice followed
by a `davcio` read error.

Accepted for `acbn0 > 0` although DFPT refuses them: tetrahedra occupations;
several channels per type.

Not checked, and not validated: PAW and ultrasoft pseudopotentials (every
validation run used norm-conserving PseudoDojo or GBRV sets with the stock code
as reference).  Do not read "not refused" as "supported".

Not implemented: f manifolds; `acbn0_type = 3` (anisotropic); noncollinear /
SOC; self-consistency; forces and stress; a `HUBBARD` card writer for official
pw.x.

---------------------------------------------------------------------------

## 7. Output

### 7.1 stdout

In order, for `acbn0 = 2` (the `acbn0 = 1` run omits the neighbour and V parts):

1. The official HP banner ("Calculation of Hubbard parameters using the HP
   code based on DFPT ... Timrov ...").  This is official boilerplate printed
   before the route is known; no DFPT is done.
2. `read_file` output ("Reading xml data from directory", "Reading collected,
   re-writing distributed wavefunctions").
3. The manifold declared from `&INPUTHP`:
   ```
   Hubbard manifold selected in &INPUTHP (acbn0_ehub_l_choice)
   (the pw.x run carried no HUBBARD card)
     type  1:   4s  + 4p     on-site U removed for: s
     type  2:   4s  + 4p     on-site U removed for: s
     projectors: ortho-atomic     nwfcU =     8
   ```
4. Route banner (`eACBN0 Hubbard parameters (acbn0 = 2) ... NOT linear
   response`; the third line `stage V3: ...` is a leftover development label).
5. `Hubbard manifold taken from MODULE ldaU:` per type.
6. `inter-site cutoff from &INPUTHP: nn_distance = ...`, then the neighbour
   list in stock's format, `atom 1: na, atom 2: nb, n1 n2 n3` meaning atom nb
   displaced by lattice vector n lies within d_IJ of atom na in the home cell,
   then `nn_distance = ... nnnc = N`.  `nnnc` is the number of DISTINCT cells
   carrying at least one pair, home cell included (GaAs d = 3.0 A: 7; fcc Ni:
   13); it is not the number of pairs.
7. `V[atom i, atom j, #orb(i), #orb(j)]` tables in eV, stock's exact format:
   `[ na1, na2, lm1, lm2 ]` then `[ n1 n2 n3 ]` then lm1 x lm2 numbers.  Rows
   and columns follow the manifold order, i.e. the order of the selected
   channels in the pseudopotential file (4s then 4p for the example
   pseudopotentials; Ga.upf lists 3D first, so d+s+p on Ga would put d first).
   `na1 = na2`, cell 0 is the on-site U block.
8. `eACBN0 Hubbard parameters written to <outdir>eACBN0_hub.dat`, a per-channel
   U table, `N inter-site pair(s) with non-zero V`, `one-shot E(+U+V) = ...
   Ry (on-site ..., inter-site ...)` and the line `(not self-consistent: U, V
   and the density do not feed back into each other)`.
9. `eACBN0 dumps (forked-code formats) written to <outdir><prefix>.save/`.
10. Clocks.  `hp_acbn0` is printed; `acbn0_ns_a0` only for `acbn0 = 1` (the
    U+V kernel clock `acbn0_ns_e0_1` is started but not in the print list).

### 7.2 `outdir/eACBN0_hub.dat`

Written by the ionode into the `outdir` of `&INPUTHP`, next to `<prefix>.save`,
`STATUS = 'unknown'`: a previous file of that name is overwritten silently, and
the `acbn0 = 1` and `acbn0 = 2` runs use the SAME name.  Rename between runs,
as `example/*/run.sh` does.

Sections, `acbn0 = 2`:

| section                                    | content                                                                                                          |
|--------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| header                                     | prefix; `projectors : ortho-atomic`; `bookkeeping: MODULE ldaU (official QE) + own neighbour table`; `mixing : 1.000 (one-shot; no SCF history)`; `neighbours : cutoff d_IJ = x A, nnnc = N`; `estimator : pseudo-hybrid eACBN0 ratio, NOT linear response` |
| `On-site U = hub_v[na, na, m1, m2, cell 0]` | one row per (site, channel): site, type, label, channel (`4s`), occupation = spin-summed trace of the on-site occupation matrix over that channel, U in eV. U is constant inside an (l1,l2) block, so one number per channel is exact. A masked channel prints 0.0000. |
| `Inter-site V`                             | one row per (na1, na2, cell, channel block): atoms and labels, `cell [n1 n2 n3]` of na2, pair distance d in A, block label `4s-4p`, V in eV. Both directions are listed; the (na2,na1) block is the transpose. Pairs whose whole block is below 1e-8 Ry are skipped. Footer `N pair(s) with |V| > 1.0E-08 Ry` counts directed (na1,na2,cell) triples. |
| `One-shot extended Hubbard energy`         | on-site (U), inter-site (V), and E(+U+V) in Ry and eV, with the formula of stock `v_hubbard_v` evaluated with the one-shot U, V and the occupations they came from. It equals what stock pw.x prints as `E(+U+V)` on the first iteration from this density with `ehub_mixing = 1`. It is NOT a total-energy contribution: nothing was iterated. |
| `Full on-site U matrix (eV)`               | the lm x lm matrix per atom                                                                                       |
| `On-site occupation matrix ns[na, m1, m2, spin]` | per atom and spin, dimensionless. `-0.0000` is a numerical zero of order 1e-16 or smaller printed with f10.4. |

`acbn0 = 1` has the same layout without the neighbour and V parts: header,
`Summary` (occupation and `U_ACBN0 (eV)` per channel), `One-shot Hubbard
energy E(+U)`, full U matrix, occupation matrix.

The inter-site occupation matrices and the full V matrices beyond the
per-block representative are not in this file; they are in the dumps below and
in the stdout tables.

### 7.3 Dumps in `outdir/<prefix>.save/` (stock file formats)

| file                  | route  | content                                                                                                         | unit  |
|-----------------------|--------|-----------------------------------------------------------------------------------------------------------------|-------|
| `ehub_uv.txt`         | 1, 2   | list-directed dump of U padded to `(nat,9,9)` (route 1) or V padded to `(nat,nat,9,9,nnnc)` (route 2); exactly the array stock's `write_scf` dumps and `read_ehub_uv` reads | Ry    |
| `ehub_u0.txt`         | 1      | human-readable `U[atom i, #orb, #orb]` blocks, stock format (byte-identical to stock's on the same density). The `[ na, L, L ]` header uses `L = 2 l_max + 1`, so it says 3 for a 4x4 s+p block; inherited quirk. | eV    |
| `ehub_uv0.txt`        | 2      | `na1 na2 m1 m2 n1 n2 n3 nnc V` rows, home cell (nnc = 1, U included) first. The header's `l, m, n` are the cell vector, not angular momenta (stock wording). | Ry    |
| `occupation_uv.txt`   | 1, 2   | list-directed dump of `ns` (`(lm,lm,nspin,nat)` or `(lm,lm,nspin,nat,nat,nnnc)`), machine precision              | -     |
| `occupation_uv0.txt`  | 1, 2   | indexed list of the non-zero `ns` entries. For route 1 stock writes an EMPTY file here; dev writes the listing.   | -     |

Both routes write `ehub_uv.txt`, `occupation_uv.txt` and `occupation_uv0.txt`
under the same names and overwrite silently; only `ehub_u0.txt` (route 1) and
`ehub_uv0.txt` (route 2) are route-specific.  After `example/*/run.sh` the
`.save` holds the route-2 dumps plus a stale route-1 `ehub_u0.txt`.  The same
names are the ones stock's `write_scf` uses, so running hp.x on a stock
`.save` destroys stock's dumps unless they are copied first.

Stock pw.x can read `ehub_uv.txt` and `occupation_uv.txt` back when its input
flags `read_ehub_uv_file` / `read_ehub_ns_file` are set (routines
`read_ehub_uv`, `read_ehub_ns` in `io_rho_xml.f90`), if and only if the shapes agree:
same geometry, same cutoff (`ehub_nn_distance = nn_distance`), same channel
selection in the same order.  Dev never reads any of these files; it reads
only the `.save` named in its own `&INPUTHP`.

Side effects of a run: distributed `<prefix>.wfc*` files, one per MPI process,
are rewritten into outdir by `read_file` (32 files, about 470 MB for the GaAs
example on 32 processes), and an empty `outdir/HP/` is created.

---------------------------------------------------------------------------

## 8. Build and run

Toolchain actually used: Intel 19.1.3.304 (`ifort`, Intel MPI `mpiifort`,
`mpicc` resolving to gcc with Intel MPI paths), MKL for BLAS/LAPACK and FFT
(`-D__DFTI -D__MPI`), no OpenMP, no ScaLAPACK, no GPU.  `make.inc` and
`HP/src/make.depend` are gitignored, so a fresh clone needs `configure`.

```
source ~/.config/qe/env.sh                       # module purge; module load intel/19.1.3.304
cd ~/projects/QE/qe-7.3_thiswork
./configure MPIF90=mpiifort F90=ifort CC=mpicc
make update_submodule=: pw hp                    # see note on submodules
cp HP/src/hp.x bin/hp.x.$(git rev-parse --short HEAD)   # pin the binary to its source commit
```

Submodules: the `libcuda_devxlib` and `libmbd` targets of
`install/extlibs_makefile` call the `update_submodule` macro of
`install/install_utils`, which runs `git submodule update` whenever the tree
has a `.git`.  In THIS working copy both directories are populated at the
hashes recorded in `external/submodule_commit_hash_records` (devxlib
`a6b89ef`, mbd `82005cb`) but the submodules are not initialised in
`.git/config`, so the update step is bypassed with `update_submodule=:` (a
shell no-op).  In a FRESH clone the two directories are empty gitlinks; fill
them first:

```
git submodule update --init -- external/devxlib external/mbd   # or clone the .gitmodules URLs at a6b89ef / 82005cb
```

Two more things to know:

- `install/libcuda_devxlib` is an empty stamp file that got committed.  In a
  fresh clone `make` will consider devxlib already built and skip it; delete
  the stamp or build `external/devxlib` by hand before `make`.
- `external/devxlib/make.inc` is a symlink to the tree's pre-2026-09-28 path
  and dangles.  The make-driven rebuild re-creates it (`ln -fs`); only a
  by-hand `make` inside `external/devxlib` trips over it, so fix the link
  first in that case.

The C file is compiled by `make.inc`'s generic `.c.o` rule, i.e. by whatever
`mpicc` is.  Fortran calls it through `BIND(C)` with `VALUE` arguments, so the
gcc/ifort mix is ABI-safe.

Pinned binaries in `bin/` (named `<code>.x.<short hash of the source commit>`):

| binary            | md5                                | note                                                  |
|-------------------|------------------------------------|-------------------------------------------------------|
| `bin/hp.x.a1ad49f`| `a3d8daac790938795e5d63632f7c0024` | dev hp.x, this README's code state                     |
| `bin/pw.x.e2801ae`| `0a59d7478de1e56ea82ff69f37906841` | pw.x from this tree; PW/ and Modules/ are official     |

`bin/pw.x` and `bin/hp.x` are make's symlinks into `PW/src` and `HP/src` and
change on every rebuild; job scripts should name the hash copies.

Running (this cluster: every QE run is an `sbatch` job from `/gpfs/yjkang`,
never on the login node):

```
module purge; module load intel/19.1.3.304
export MKL_DEBUG_CPU_TYPE=5 MKL_CBWR=AUTO
srun --mpi=pmi2 --kill-on-bad-exit=1 $BIN/pw.x.e2801ae -npool 32 -inp scf.in   > scf.out
srun --mpi=pmi2 --kill-on-bad-exit=1 $BIN/hp.x.a1ad49f -npool 32 -inp hp_U.in  > hp_U.out
mv out/eACBN0_hub.dat eACBN0_hub_U.dat
srun --mpi=pmi2 --kill-on-bad-exit=1 $BIN/hp.x.a1ad49f -npool 32 -inp hp_UV.in > hp_UV.out
mv out/eACBN0_hub.dat eACBN0_hub_UV.dat
```

Parallelisation of the ACBN0 route: the k-point loop runs over the local pool
and the occupation sums are reduced over `inter_pool_comm`; the ERI and
Gaussian-fit parts are serial per rank.  Use `-npool` equal to the number of
ranks.  The route is sub-second for GaAs (s+p) and about 9 s for Ni (d, more
ERIs) on 32 ranks.

---------------------------------------------------------------------------

## 9. `example/`

```
example/pseudo/   Ga.upf As.upf Ni.upf (PseudoDojo NC SR PBE standard v0.4, md5 in README.md)
example/GaAs/     scf.in scf.out  hp_U.in hp_U.out eACBN0_hub_U.dat  hp_UV.in hp_UV.out eACBN0_hub_UV.dat  run.sh
example/Ni/       same nine files
example/README.md reference numbers and what they mean
```

GaAs: zinc blende, s+p manifold on both atoms, 4s on-site U masked
(the manifold choice of the stock GaAs example).  Ni: fcc, 3d, one atom per
cell.  The non-magnetic start copies Lee and Son's `example1/Ni.scf.in`
(`starting_magnetization = 0.0`) so the numbers compare with stock; the two
spin blocks of `ns` are therefore identical.  A ferromagnetic start gives
spin-split occupations and a different U.  Both produced with
`bin/pw.x.e2801ae` and `bin/hp.x.a1ad49f` on 2026-10-30.  Headline numbers in
eV, next to stock self-consistent runs with the same cell, cutoff, k-mesh,
smearing, pseudopotentials and manifold (`/gpfs/yjkang/project-01/3-hub_dat_files`):

| quantity             | dev one-shot (`acbn0 = 1, 2`) | stock SC, U-only run | stock SC, U+V run |
|----------------------|-------------------------------|----------------------|-------------------|
| U(Ga 4p)             | 0.3418                        | 0.3259               | 0.3460            |
| U(As 4p)             | 2.0105                        | 2.0672               | 1.9794            |
| V Ga-As 4s-4s        | 0.8682                        | -                    | 0.8230            |
| V Ga-As 4s-4p        | 1.3418                        | -                    | 1.2971            |
| V Ga-As 4p-4s        | 0.8281                        | -                    | 0.8042            |
| V Ga-As 4p-4p        | 1.6677                        | -                    | 1.6770            |
| U(Ni 3d)             | 16.4018                       | 16.8702              | 16.7860           |
| V Ni-Ni 3d-3d (12 images, 2.4247 A) | 2.5890         | -                    | 2.6496            |

The on-site U of the `acbn0 = 2` run equals the `acbn0 = 1` run exactly for
the same ground state; it does not depend on `nn_distance`.  In stock, the
U-only and U+V self-consistent runs give different U because V feeds back
into the density.  One-shot and self-consistent values differ by 1 to 5 %
here, and self-consistency also moves the occupations (Ga 4p 1.8937 to 1.9443,
Ni 3d 8.8346 to 9.3270 in the U+V runs).

`run.sh` refuses to overwrite `scf.out`, `hp_U.out`, `hp_UV.out` or `out/`.
To rerun, copy `GaAs/` or `Ni/` next to `pseudo/` in a fresh location under
`/gpfs`, delete those files, and `sbatch run.sh` from inside that directory.

`example/.gitignore` contains `!*.out` because the top-level `.gitignore`
ignores every `*.out`.

---------------------------------------------------------------------------

## 10. Known rough edges (code state a1ad49f)

Behaviour is as documented above; the following are cosmetic or historical
and worth knowing before reading the source:

- The DFPT banner is printed for every run, including `acbn0 > 0`.  The route
  banner and the `.dat` header say "NOT linear response".
- `stage V3: ...` is printed in the `acbn0 = 2` banner; the file header calls
  the same code `S4` and the multi-channel support `V4`.  The stage narrative
  in the `eACBN0.f90` header (S1 ... V4) is development history, not current
  behaviour.
- Several comments and error texts predate `acbn0_manifold_setup` and say the
  opposite of what the code does: `hp_readin.f90` input_sanity says `acbn0 = 2`
  "needs a HUBBARD card too"; `hp_acbn0_check` / `hp_acbn0_check_v` advise
  "rerun pw.x with HUBBARD (ortho-atomic)"; `hp_acbn0_build_man` blames "the
  pw.x run had no HUBBARD card"; `hp_acbn0_man` tells the user to give a
  background channel a nonzero U2 in the HUBBARD card; `hp_main.f90` and
  `hp_readin.f90` cite `example/MnO.acbn0.in` and
  `example/GaAs.acbn0.in_noncollienar`, which are files of the stock patch
  repository, not of this tree.  All of these guards are unreachable: the
  manifold setup forces `lda_plus_u`, `lda_plus_u_kind = 0` and
  `'ortho-atomic'` before they run.  The `acbn0 /= 2` exemptions in
  input_sanity are dead code for the same reason.
- The `eACBN0.f90` header says of `remove_ehub_u` that "ldaU has no such
  mask, so the bug class is gone rather than worked around".  The mask exists:
  since `d204215` (as `remove_u`), as `acbn0_rm_ehub_u` with the opposite
  convention since `38a7c8e`, and as `acbn0_remove_ehub_u` since `a1ad49f`.
  The header's own translation table lists it.
- The header says ACBN0 runs on "ANY ortho-atomic DFT+U ground state,
  including the very one the DFPT branch of hp.x uses", and the `hp_acbn0_v`
  comment says "an ordinary ortho-atomic DFT+U ground state - the same one
  hp.x's DFPT branch consumes".  The code refuses a DFT+U ground state for
  `acbn0 > 0` and requires one for `acbn0 = 0`.
- `hp_acbn0_write`'s header note says one channel per type is enforced; the
  body already handles several.
- The `.dat` prose cites `v_hubbard_a0.f90`; the routine lives in stock
  `v_hubbard_v.f90`.
- `hp_acbn0_ns_v` (stock `new_ns_v`, the occupation routine of a fixed-U/V
  stock run, also reached through the `stop_acbn0` short-circuit) is ported
  but never called.
- `nspin == 1` and `gamma_only` branches inside the kernels are dead code
  behind the `nspin = 2` / no-gamma checks.
- The old `&INPUTHP` names `acbn0_l_choice` and `acbn0_rm_ehub_u` (commits
  before `a1ad49f`, opposite mask convention) no longer parse; the namelist
  read fails.  Old job scripts under `/gpfs/yjkang/project-01/{1,2,3}-*` use
  them and are a record, not rerunnable.
- `HP/Doc/INPUT_HP.def` and `HP/CMakeLists.txt` were not updated.

---------------------------------------------------------------------------

## 11. References and licences

- L. A. Agapito, S. Curtarolo, M. Buongiorno Nardelli, Phys. Rev. X 5, 011006
  (2015).  ACBN0.
- S.-H. Lee, Y.-W. Son, Phys. Rev. Research 2, 043410 (2020).  eACBN0, the
  d_IJ cutoff, the stock patch.  Patch repository: W. Yang, S.-H. Lee,
  Y.-W. Son, GPL v2 or later; local copy
  `~/projects/project-01-transplanting-eACBN0-from-PW-to-HP/kias_dft.u.v`.
- I. Timrov, N. Marzari, M. Cococcioni, Phys. Rev. B 98, 085127 (2018);
  103, 045141 (2021); Comput. Phys. Commun. 279, 108455 (2022).  DFPT route.
- Quantum ESPRESSO: P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502
  (2009); 29, 465901 (2017); J. Chem. Phys. 152, 154105 (2020).  GPL v2.

`HP/src/hp_acbn0_crys.c` and `.h` are PyQuante's `crys.c/.h` (Rys quadrature
after Augspurger, Bernholdt, Dykstra, J. Comp. Chem. 11, 972 (1990)),
Copyright (c) 2004 Richard P. Muller, modified BSD licence, carried over from
the stock patch with only symbol renames.  They are the one non-GPL component
of `HP/src`.
