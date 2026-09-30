# eACBN0 in hp.x — validation report

**3 September 2026** · web version with diagrams: <https://claude.ai/code/artifact/c7998bb7-d39c-426b-8752-387f33cb43f5>

Transplanting the DFT+U+V (eACBN0) Hubbard-parameter kernels out of the QE 7.3
fork's `pw.x` and into `hp.x`, so that one binary produces both the
pseudo-hybrid and the linear-response estimate of *U* from a single ground
state.

Method reference: **S. Lee and Y-W. Son, Phys. Rev. Res. 2, 043410 (2020)**
(named in `example/Howto2.pdf`). ACBN0 itself: Agapito, Curtarolo &
Buongiorno Nardelli, PRX 5, 011006 (2015).

---

## 1. Headline

Both branches now run on an ordinary DFT+U ground state, and both reproduce the
original authors' own MnO reference output **exactly** — every U, every V,
every cell vector, 168 lines / 51 blocks, empty diff. So U_ACBN0, V_ACBN0 and
U_DFPT all come from a single pw.x run.

Four latent bugs were found and fixed; two would have returned a plausible
wrong number rather than crashing. One capability was deliberately traded away
(section 7).

| Quantity | Reference | Our pw.x | Our hp.x |
|---|---:|---:|---:|
| Mn_A 3d on-site `V[1,1]` cell `[0 0 0]` | 5.2344 | **5.2344** | **5.2344** |
| Mn_B 3d on-site `V[2,2]` | 5.2344 | **5.2344** | **5.2344** |
| O 2p on-site `V[3,3]`, `V[4,4]` | 3.4161 | **3.4161** | **3.4161** |
| Mn–O inter-site, cell `[-1 -1 0]` | 2.9031 | **2.9031** | **2.9031** |
| Mn–O inter-site, cells `[-1 0 -1]`, `[-1 0 0]` | 2.9042 | **2.9042** | **2.9042** |
| absolute magnetisation (μB/cell) | 9.36 | **9.36** | — |
| SCF iterations | 33 | 30 | — |

Source: `example/MnO.acbn0.out`, QE v7.3, 128 procs, 5 Jan 2024,
`acbn0_type = 2`.

### Parameter audit

`MnO_ref/stepB.in` differs from `example/MnO.acbn0.in` on **one line**:
`pseudo_dir`. Pseudopotentials verified equivalent rather than assumed —
`psp/Mn.upf` and `psp/O.upf` declare `pseudo_type="NC"`,
`relativistic="scalar"`, `functional="PBE"` (= pseudodojo `nc_sr`) with
`z_valence` 15 and 6, giving the 42 electrons/cell their output echoes.
Identical throughout: `ecutwfc 100`, `ecutrho 400`, PBE, `tetrahedra_opt`,
15³ → 344 k, `celldm 10.289967`, `conv_thr 1e-8`, `diago_full_acc`,
`mixing_beta 0.7`, `startingwfc random`, `ehub_nn_distance 3`,
`ehub_conv_thr 1e-8`, `ehub_mixing 0.85`, `ehub_pot_fix .true.`,
`read_ehub_ns_file .true.`, `ehub_l_choice` MnA-3d / MnB-3d / O-2p.

### The first MnO attempt failed, informatively

Setting `read_ehub_ns_file = .false.` gave a converged AFM state with
Mn_A = 12.1473 and Mn_B = 4.5768 eV. The two Mn are symmetry-equivalent
partners and the reference has them equal — that asymmetry was the tell.

Their input is **step B of a documented two-step recipe** (`Howto2.pdf`):
`ehub_pot_fix` *freezes* the occupations entering the Hubbard energy, and
`read_ehub_ns_file` is what loads step A's DFT-level occupations to freeze.
Turning it off froze whatever `init_ns_v()` guessed — a physically different
calculation, not a different starting point. Preserved in
`MnO_ref/attempt1_from_scratch/`.

---

## 2. Method — one variable per rung

The dangerous failure is not a crash but a plausible wrong number, so no stage
changes two things at once.

| Stage | Driver | Kernel | Bookkeeping | Ground state | Result |
|---|---|---|---|---|---|
| **S0** | pw.x `sum_band` | libpw `new_ns_a0` | ldaV | ACBN0-converged | 9.2501 / 9.2509 / 4.2471 |
| **S1** | ***eACBN0.f90*** | libpw `new_ns_a0` | ldaV | same | identical |
| **S2** | eACBN0.f90 | ***our `hp_acbn0_ns_a0`*** | ldaV | same | byte-identical |
| **S3a** | eACBN0.f90 | our kernel | ***ldaU*** | same | identical |
| **S3b** | eACBN0.f90 | our kernel | ldaU | ***official bare-DFT*** | 10.5455 / 4.2717 |
| **V0** | pw.x `sum_band` | libpw `new_ns_e0_1` | ldaV | eACBN0-converged | 7.7065 / 2.6538 + V |
| **V1** | ***eACBN0.f90*** | libpw `new_ns_e0_1` | ldaV | same | identical, 456 lines |
| **V2** | eACBN0.f90 | ***our `hp_acbn0_ns_e0_1`*** | ldaV | same | identical, 456 lines |
| **V3a** | eACBN0.f90 | our kernel | ***ldaU*** | same | identical, 456 lines |
| **V3b** | eACBN0.f90 | our kernel | ldaU | ***official bare-DFT*** | 10.5455 / 4.2717 + V |

*Italic* = the one thing that changed on that rung.

**S3b's direct reference.** S3b changes the input data, not the code, so its
number legitimately moves and needed its own pw.x check. Blocked at first:
with `lda_plus_v` set, `read_scf` reads `occupation_uv.txt` unconditionally
(`io_rho_xml.f90:644`) and the official save has only the ldaU `occup.txt`.
Solved in two steps — step 1 with `lacbn0 = .false.` (so `ehub_u = ehub_v = 0`
and the Hamiltonian is the official one) writes the file; step 2 does one ACBN0
evaluation on it.

```
                       pw.x one-shot   hp.x S3b
    Ni1 3d               10.5455        10.5455    exact
    Ni2 3d               10.5455        10.5455    exact
    O   2p                4.2718         4.2717    last printed digit
    Tr[ns] atom 1         8.46278        8.4628    exact

    control: step1 E_tot = official ../NiO E_tot = -737.16322501 Ry
```

### Static faithfulness of the transplant

Diffed against `PW/src/new_ns_v.f90` with renames undone and blank/comment
lines stripped. **2472 code lines, five differences** — the three `l_separation`
guards, plus (in `symm_n5`, added at V3) the `nnc == 0` guard.

| Our routine | From | Code lines | Difference |
|---|---|---:|---|
| `hp_acbn0_ns_a0` | `new_ns_a0` 233–633 | 401 | identical |
| `hp_acbn0_orberi` | `orberi` 3780–3818 | 39 | identical |
| `hp_acbn0_ns_e0_1` | `new_ns_e0_1` 636–1808 | 833 | identical |
| `hp_acbn0_ns_v` | `new_ns_v` 10–232 | 159 | identical |
| `hp_acbn0_symm_n4` | `symmetrize_n4` 2964–3165 | 161 | +1 bugfix |
| `hp_acbn0_symm_n4m` | `symmetrize_n4m` 3168–3399 | 179 | +1 bugfix |
| `hp_acbn0_symm_n5` | `symmetrize_n5` 3400–3779 | 307 | +2 bugfixes |

---

## 3. Root cause of every bug

The fork's patch touches 60 files. **None is under `HP/`** — the authors did
DFT+U+V in `pw.x` only. `pw.x` initialises `MODULE ldaV` from *two* paths:

- `scf.in` → `dftU_iosys_v1/_v2` (`input.f90:2083, 2188`)
- the XML → `qexsd_copy` + `init_lda_plus_v` (`read_file_new.f90:325`)

`hp.x` only ever walks the second. Anything set only on the input path arrives
undefined, and the XML does not round-trip it.

A sweep of all 35 `ldaV` module variables against that criterion turned an
ad-hoc hunt into a closed list. Genuine risks, all handled:
`remove_ehub_u`/`_v`, `ehub_mixing`, `ehub_u`, `ehub_v`, `ehub_conv_thr`
(guarded by `stop_acbn0`), `pot_fix` (never read by `new_ns_e0_1`).
Unreachable from the kernels: `ehubene`, `ehubene_old`, `dehubene`,
`d_spin_ldav`.

### Bug 1 — `remove_ehub_u` / `remove_ehub_v` undefined (silent)

No initialiser in `ldaV.f90`; assigned only at `input.f90:2083`. The kernels
multiply the final *U* by them **twice**. Undefined, this does not crash — it
silently scales or zeroes every Hubbard parameter. Set to the documented
default of 1.

Not hypothetical: `example/GaAs.acbn0.in_noncollienar` sets
`remove_ehub_u(1,1) = 0` to kill the s-channel *U*.

### Bug 2 — `iunwfc` is closed at the insertion point (crash)

`read_file` buffers wavefunctions to `iunwfc` then closes it with
`close_buffer(..., 'KEEP')` (`read_file_new.f90:62`). `'KEEP'` preserves the
*files*, but the call also deregisters the unit, so `get_buffer` falls through
to `davcio` on a closed unit. Official `hp.x` never trips this because DFPT reads
through `iuwfc`, opened later in `hp_openfil_q`.

### Bug 3 — `hubbard_get_nnc` returns an undefined index

`PW/src/hubbard_nns.f90:241`. `nnc` is assigned only inside
`IF (MAXVAL(ehub_nn_vect(...)) == 1)`, so a cell carrying no neighbour pair
falls through untouched. `hubbard_symmetry` asks for the index of a *rotated*
cell — not guaranteed to be one of the occupied ones — and `symmetrize_n5`
then indexes `n(nnc)` with whatever was on the stack.

Never fires on the fork's own NiO ground state; fires on the official one, with
an identical neighbour table, identical symmetry (12 ops) and identical
manifold. In an optimised build a garbage index is silent whenever it lands in
a mapped page. The optimised build gave only a UCX backtrace; a rebuild with
`-check bounds,uninit -fpe0` named it at once:

```
forrtl: severe (408): Subscript #1 of the array N has value 1059660741
                      which is greater than the upper bound of 19
  hp_acbn0_symm_n5_   3949  eACBN0.f90
  hp_acbn0_ns_e0_1_   2587  eACBN0.f90
```

Fixed by returning 0 and skipping the pair — the correct contribution, since a
rotated image outside the neighbour set contributes nothing to the symmetrised
occupation matrix.

### Bug 4 — `l_separation` writes one element before the array (latent)

In `symmetrize_n4`, `_n4m` and `_n5` alike, a type with no Hubbard manifold has
`ehub_l(1,nt) = 0` and the next line writes `l_separation(0,nt)`. Dormant when
every type is Hubbard-active — true of NiO, MnO and GaAs — and fires the moment
a spectator species is added. **Live in the fork.**

### Two gaps in hp.x's own guards

`hp_readin` refused to start without a `HUBBARD` card and refused
`tetrahedra`. Both are linear-response preconditions; ACBN0 solves no
Sternheimer equation. But **neither authors' reference input has a `HUBBARD`
card**, and `MnO.acbn0.in` uses `tetrahedra_opt` — between them these blocked
the entire validation. Relaxed for `acbn0 > 0` only; ACBN0-side preconditions
enforced in `hp_acbn0_check` / `hp_acbn0_check_v`. The dispatch also moved
above `hp_init`, which counts ldaU Hubbard atoms and aborts when a ldaV-only
ground state legitimately has none.

---

## 4. What the two estimators say

Both from the *same* NiO ground state (official `../NiO`, U ≈ 0, ortho-atomic,
identical manifold):

| | ACBN0 (this work) | DFPT (official hp.x) |
|---|---:|---:|
| Ni 3d | **10.5455** | 6.7427 |
| O 2p | **4.2717** | 7.3243 |

They disagree in *opposite directions* on the two species. That divergence is
the result.

- **Self-consistency:** one-shot on bare DFT gives Ni 10.55; at the ACBN0 fixed
  point 9.25. O barely moves, 4.27 → 4.25.
- **Turning on V** (3.0 Å cutoff): on-site drops to Ni 7.71, O 2.65, with
  V(Ni–O) = 3.18, V(Ni–Ni) = 1.05, V(O–O) = 0.54 eV. The same screening is now
  shared between U and V.
- **Cost:** ~12 s on 128 cores against an existing ground state.

### With inter-site V, from that same single run

Once V3 put the +V branch on ldaU it could be pointed at the same ground state.
The on-site U comes back at **10.5455** and **4.2717** — identical to the
+U-only branch to all six figures, from a different kernel. Correct, not a
coincidence: the on-site ratio depends only on on-site occupations, so the
cutoff cannot move it in a one-shot evaluation; it enters only through
self-consistency.

| Pair | d (Å) | V (eV) | shell |
|---|---:|---:|---|
| Ni–O | 2.0603 | **3.1859** | octahedral nearest neighbour |
| Ni–Ni | 2.9020 | 1.4359 | second shell, split by the |
| O–O | 2.9020 | 0.8741 | rhombohedral AFM-II distortion |
| Ni–Ni | 2.9254 | 1.4242 | |
| O–O | 2.9254 | 0.8671 | |

Cutoff d_IJ = 3.0 Å, cleanly excluding the next shell at 3.5398 Å.

---

## 5. Corroboration

**The two kernels agree with each other.** `Howto2.pdf` states `acbn0_type = 1`
is *not* the maintained path and that +U should be done as `acbn0_type = 2`
with `ehub_nn_distance = 0.0`. Run that way NiO gives 9.2501 / 9.2509 / 4.2471
— identical to `acbn0_type = 1`, same 37 iterations, zero neighbour pairs
found. `new_ns_a0` and `new_ns_e0_1` are independent kernels (~400 vs ~1170
lines, different occupation-matrix rank), so this is real cross-validation.

**Multi-channel manifolds work.** GaAs with Ga 4s+4p and As 4s+4p is the only
case where `ehub_l(m,nt)` varies with *m*, so it is the first real exercise of
`l_separation`. hp.x reproduces pw.x across 204 lines, with the block
structure `A0%U(l(m1),l(m2))` predicts:

```
Ga on-site [1,1,4,4] cell [0 0 0]
  0.9709  0.5808  0.5808  0.5808     U[s,s] = 0.9709
  0.5808  0.3441  0.3441  0.3441     U[s,p] = 0.5808
  0.5808  0.3441  0.3441  0.3441     U[p,p] = 0.3441
  0.5808  0.3441  0.3441  0.3441
```

---

## 6. Usage

```fortran
 &inputhp
   prefix      = 'NiO',
   outdir      = './outdir',
   acbn0       = 2,     ! 0 = DFPT (default) | 1 = ACBN0 U | 2 = eACBN0 U+V
   nn_distance = 3.0,   ! d_IJ in Angstrom, for acbn0 = 2
 /
```

`acbn0 > 0` evaluates the parameters from the ground state in `outdir`, writes
`eACBN0_hub.dat`, and stops before the q-point loop. The kernels also print the
same `V[atom i, atom j, ...]` tables `pw.x` prints, which is what makes every
comparison here a plain `diff`.

**Both branches now read official QE's `ldaU`**, so both need an ortho-atomic
`HUBBARD` card and one channel per type, and both run on **any** such ground
state — including the one DFPT consumes. Neither needs `lacbn0` any more.

`nn_distance` is the inter-site cutoff d_IJ, required for `acbn0 = 2` unless
the run carries `lda_plus_v` (then `ehub_nn_distance` comes from the XML).
Because hp.x now owns the neighbour table, **the cutoff can be swept without
re-running pw.x** — which the fork cannot do.

Note on the reference inputs: `example/MnO.acbn0.in` declares its manifold
through `ehub_l_choice` and carries no `HUBBARD` card, so under V3 it needs one
added — same manifold, `U = 0.0001` eV, numerically nothing. With that one
addition it still reproduces the reference to all 168 lines
(`MnO_ref_v3/`).

---

## 7. Open items

**Multi-channel manifolds no longer run on the +V branch.** V3's whole point
is that `ldaU` describes the manifold — and `ldaU` cannot express two Hubbard
channels on one atomic type. Its scalar `Hubbard_l(nt)` names one channel and
puts any second in a *separate* `offsetU_back` block, where `ldaV` concatenates
them into a single manifold. So GaAs with Ga 4s+4p, which ran under V2, is now
refused with an explicit error rather than silently projected onto one channel.
The V2 build is kept at `eACBN0.f90.v2bak` for that case; making it work under
`ldaU` means teaching the kernels about background channels.

**GaAs cannot validate against the reference.** Theirs is noncollinear +
spin-orbit (→ `new_ns_e0_nc_1`, a separate ~2300-line kernel set, not
transplanted) with PBEsol `nc-fr` pseudos where we have PBE `nc_sr`. Wrong
functional and wrong relativistic treatment — unreachable by construction. Our
GaAs run is a smoke test comparing our pw.x to our hp.x.

**`remove_ehub_u` has no route into hp.x.** Input-only, not round-tripped
through the XML. We default it to 1, correct whenever the pw.x run did — but
any calculation masking a channel (as their GaAs does) will disagree. Fix: an
`&INPUTHP` tag.

**`nspin = 1` is rejected, correctly.** Not our limitation:
`init_lda_plus_v` refuses it at `ldaV.f90:103`, and the A0/HF assembly
hard-codes *both* spin slots — `p%O(m3,m1,2,na)` on an array dimensioned
`(...,nspin,...)` — so `nspin = 1` reads past the end silently. Use
`nspin = 2` with `starting_magnetization = 0.0`.

---

## 8. Artefacts

- `HP/src/eACBN0.f90` — 3763 lines, 15 subroutines
- Modified: `hpcom.f90`, `hp_readin.f90`, `hp_bcast_input.f90`, `hp_main.f90`,
  `hp_print_clock.f90`, `Makefile`, `make.depend`
- Backups kept: `eACBN0.f90.{s1,s2,s3,v1}bak`, `hp_readin.f90.bak`,
  `hp_main.f90.bak`
- `logs/260901-eacbn0-hpx.txt` — full session log
- Runs under `/gpfs/yjkang/examples/hp/`: `NiO_acbn0`, `NiO_acbn0_v`,
  `NiO_dual`, `NiO_ref1shot`, `NiO_type2_nn0`, `MnO_ref`, `GaAs_smoke`,
  `_final_check`
