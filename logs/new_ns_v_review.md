# `new_ns_v.f90` — Integrated Paper + Code Review

**Target file:** `~/projects/qe_aCBN0-7.3/PW/src/new_ns_v.f90` (6751 lines, 269 kB)
**Paper reviewed:** N. Tancogne-Dejean, M. J. T. Oliveira, A. Rubio,
*"Self-consistent DFT+U method for real-space time-dependent density functional theory calculations"*,
**Phys. Rev. B 96, 245133 (2017)** — `logs/PhysRevB.96.245133.pdf`
**Repo:** `qe_aCBN0-7.3` = QE 7.3 + an out-of-tree aCBN0 / Hubbard-V patch (git: `1711488 Initial commit: QE 7.3 base + aCBN0/Hubbard-V patch`)
**Date:** 2026-08-27

---

## 0. Read this first — one premise in the request needs correcting

> *"new_ns_v.f90 is the script where that paper's method of calculating +V is implemented"*

Half right, and the half that is wrong matters.

**PRB 96, 245133 (2017) contains no intersite V at all.** It is the Octopus implementation of the
*on-site* ACBN0 functional (self-consistent `U` and `J`), plus two genuinely new pieces: a
**noncollinear-spin extension** (Appendices A and B) and a **real-time TDDFT+U** extension. Search
the PDF for "V_{IJ}", "intersite", or "extended Hubbard" and you get nothing. Its Hubbard energy is
Dudarev's on-site form, Eq. (2), with a single `U_eff = Ū − J̄` per `(I,n,l)` manifold.

So what *is* `new_ns_v.f90`? It is a **superset**:

| Layer | Source of the physics | Where in the file |
|---|---|---|
| On-site ACBN0 (`U`, `J`) — collinear | Agapito–Curtarolo–Buongiorno Nardelli, **PRX 5, 011006 (2015)**, restated as Eqs. (6)–(10) of your PDF | `new_ns_a0` |
| On-site ACBN0 — **noncollinear** | **Appendix A + B of your PDF**, Eqs. (A1)–(A5), (B1)–(B6) — transcribed essentially literally | `new_ns_e0_nc_1` |
| **Intersite `V` (extended/"eACBN0")** | *Not in your PDF.* This is the extended-ACBN0 generalisation (Lee & Son / Tancogne-Dejean & Rubio 2020 line of work — see the citation caveat at the end of this section) | `new_ns_e0_1`, `new_ns_e0_2` |
| Bare occupation matrices `n^{IJσ}_{mm'}(R)` | Eq. (3) generalised to two sites + lattice vector | `new_ns_v`, `new_ns_v_nc` |

The strongest direct evidence that the fork authors read *this exact paper* is a comment they left in
the source at `new_ns_v.f90:938-948`:

```fortran
! Our definition
num(ibnd)%T(l1,l2,:,nt1,nt2) = num(ibnd)%T(l1,l1,:,nt1,nt1) + &
                               num(ibnd)%T(l2,l2,:,nt2,nt2)
IF ( upf(nt1)%psd == upf(nt2)%psd ) num(ibnd)%T(l1,l2,:,nt1,nt2) = ... / 2.0_DP
! The definition of Nicolas Tancogne-Dejean & Angel Rubio
! num(ibnd)%T(l1,l2,:,nt1,nt2) = SQRT(num(ibnd)%T(l1,l1,:,nt1,nt1) * &
!                                     num(ibnd)%T(l2,l2,:,nt2,nt2))
```

That is the fork author explicitly choosing an **arithmetic-mean** intersite renormalisation over the
**geometric-mean** one used by the paper's authors in their own `+V` work. This one commented-out
line is the single most important thing in the file if you intend to publish: it is a documented,
uncited methodological divergence.

> **Citation caveat (be honest with yourself here):** I could not verify the exact bibliographic
> details of the `+V` extension papers offline. The `+V` generalisation is standardly attributed to
> S.-H. Lee & Y.-W. Son (*Phys. Rev. Research* **2**, 043410 (2020)) and, independently, to
> Tancogne-Dejean & Rubio (*Phys. Rev. B* **102**, 155117 (2020)). **Verify both before citing.**
> What is *certain* from the source tree is: (i) your PDF is the on-site + noncollinear paper,
> (ii) the code implements a `+V` generalisation the PDF does not contain.

---

## 1. Paper review — the formalism you actually need

### 1.1 The DFT+U skeleton (Eqs. 1–5)

```
E_{DFT+U}[n, {n^{Iσ}_{mm'}}] = E_DFT[n] + E_ee[{n^{Iσ}_{mm'}}] − E_dc[{n^{Iσ}_{mm'}}]        (1)
```

With Dudarev's rotationally invariant approximation this collapses to a single `U_eff = U − J`:

```
E_U = Σ_{I,n,l} (U^eff_{I,n,l}/2) Σ_{m,σ} ( n^{I,n,l,σ}_{mm} − Σ_{m'} n^{I,n,l,σ}_{mm'} n^{I,n,l,σ}_{m'm} )   (2)
```

The occupation matrix for a periodic system (**this is the single most important equation for the
code review**):

```
n^{I,n,l,σ}_{mm'} = Σ_n Σ_k^BZ w_k f^σ_{nk} ⟨ψ^σ_{n,k}|φ_{I,n,l,m}⟩⟨φ_{I,n,l,m'}|ψ^σ_{n,k}⟩          (3)
```

Minimising (1) at fixed occupations gives the nonlocal Hubbard potential:

```
V^σ_U|ψ^σ_{n,k}⟩ = Σ_{I,n,l} Σ_{m,m'} V^{I,n,l,σ}_{m,m'} P^{I,n,l}_{m,m'} |ψ^σ_{n,k}⟩               (5)
V^{I,n,l,σ}_{m,m'} = U^eff_{I,n,l} ( ½δ_{mm'} − n^{I,n,l,σ}_{mm'} ),   P^{I,n,l}_{mm'} = |φ_{I,n,l,m}⟩⟨φ_{I,n,l,m'}|
```

**Practical caveat the paper makes loudly (Sec. II, Fig. 2):** the computed `U_eff` and the resulting
band gap depend *strongly* on how far the localised orbital is allowed to extend. In Octopus this is
a real-space sphere radius; in a PW code it is the choice of projector (`ortho-atomic` vs `atomic`
vs `pseudo`) and the pseudopotential's `PP_CHI`. Two ACBN0 codes agreeing on the functional can
still disagree on `U` by ~1 eV purely through the projector. Table I of the paper shows exactly this:
their O-2p `U_eff` for MnO is 5.18 eV vs 2.68 eV in the original ACBN0 paper. Same functional,
different basis. **Keep this in mind before you claim a discrepancy is a bug.**

### 1.2 ACBN0 (Eqs. 6–10) — the part `new_ns_a0` implements verbatim

Screened (pseudo-hybrid) Hartree–Fock energy on the Hubbard manifold:

```
E_ee = ½ Σ_{m} Σ_{αβ} P̄^α_{mm'} P̄^β_{m''m'''} (mm'|m''m''')
     − ½ Σ_{m} Σ_{α}  P̄^α_{mm'} P̄^α_{m''m'''} (mm'''|m''m')                                    (6)
```

*(The paper explicitly notes it fixes a sign typo in the original ACBN0 paper here.)*

The two ACBN0 objects — memorise these, they are the whole game:

```
P̄^{I,n,l,σ}_{mm'} = Σ_{nk} w_k f_{nk} N̄^{I,n,l,σ}_{ψnk} ⟨ψ^σ_{nk}|φ_{I,n,l,m}⟩⟨φ_{I,n,l,m'}|ψ^σ_{nk}⟩   (7)

N̄^{I,n,l,σ}_{ψnk}  = Σ_{{I}} Σ_m ⟨ψ^σ_{nk}|φ_{I,n,l,m}⟩⟨φ_{I,n,l,m}|ψ^σ_{nk}⟩                          (8)
```

`N̄` is a **per-band, per-(n,l), per-spin scalar in [0,1]** — a Mulliken-like measure of "how much of
Bloch state `ψ_nk` lives in this localised manifold". The `{I}` sum runs over *all atoms of the same
species*. It is the screening knob: `N̄ → 0` recovers plain DFT, `N̄ → 1` recovers bare Hartree–Fock.
`P̄` is just the ordinary occupation matrix (3) with every band reweighted by `N̄`.

Finally,

```
        Σ_{{m}} Σ_{αβ} P̄^α_{mm'} P̄^β_{m'm''} (mm'|m''m''')
Ū =  ────────────────────────────────────────────────────────────────                            (9)
     Σ_{m≠m'} Σ_α N^α_m N^α_{m'}  +  Σ_{{m}} Σ_α N^α_m N^{-α}_{m'}

        Σ_{{m}} Σ_α P̄^α_{mm'} P̄^α_{m'm''} (mm'''|m''m')
J̄ =  ────────────────────────────────────────────────────                                       (10)
              Σ_{m≠m'} Σ_α N^α_m N^α_{m'}

U^eff = Ū − J̄ ,     N^α_m ≡ n^{I,n,l,α}_{mm}
```

**Structural reading:** numerator = screened HF energy (built from `P̄`); denominator = the same
energy expression evaluated at the *bare* mean-field level (built from the diagonal occupations `N`).
`U` is the ratio. **Every single `U`/`V` in this code is an `E0/HF` ratio of exactly this shape** —
once you see that, the 6751 lines become readable.

The paper is candid that this is not variationally consistent: the true ACBN0 potential would need an
OEP treatment (`P̄` depends on the orbitals); they and the code both just use the plain DFT+U
potential of Eq. (5). Worth a sentence in your paper.

### 1.3 Noncollinear ACBN0 (Appendices A–B) — what `*_nc` implements

Occupations become 2×2 matrices in spin space:

```
n^{I,m,l,σσ'}_{mm'} = Σ_n Σ_k^BZ w_k f_{nk} ⟨ψ^σ_{n,k}|φ_{I,n,l,m}⟩⟨φ_{I,n,l,m'}|ψ^{σ'}_{n,k}⟩       (A2)
```

Rotationally invariant Hubbard energy (Dudarev, noncollinear):

```
E_U = (U−J)/2 [ Σ_σ Σ_m n^{σσ}_{mm} − Σ_{mm'} Σ_{σσ'} n^{σσ'}_{mm'} n^{σ'σ}_{m'm} ]               (A5)
```

Renormalised density matrix, **spin-traced `N̄`**:

```
P̄^{I,n,l,σσ'}_{mm'} = Σ_{nk} N̄^{I,n,l}_{Ψnk} w_k f_{nk} ⟨ψ^σ_{nk}|φ_m⟩⟨φ_{m'}|ψ^{σ'}_{nk}⟩          (B3)
N̄^{I,n,l}_{ψnk}    = Σ_{{I}} Σ_m ⟨ψ_{nk}|φ_{I,n,l,m}⟩⟨φ_{I,n,l,m}|ψ_{nk}⟩                          (B4)
```

Note `N̄` carries **no spin index** in (B4) — the spinor overlap is already spin-traced. This is
exactly what the code does (`num(i)%O(4,4,ntyp,ntyp)` — no spin dimension).

```
        Σ_{σσ'} Σ_{{m}} P̄^{σσ}_{mm'} P̄^{σ'σ'}_{m''m'''} (mm'σ|m''m'''σ')
Ū = ──────────────────────────────────────────────────────────────────────────────────────────   (B5)
    Σ_σ Σ_{mm'} N^{σσ}_m N^{-σ-σ}_{m'} + Σ_σ Σ_{m≠m'} N^{σσ}_m N^{σσ}_{m'} − Σ_σ Σ_m N^{σ-σ}_m N^{-σσ}_m

        Σ_{σσ'} Σ_{{m}} P̄^{σσ'}_{mm'} P̄^{σ'σ}_{m''m'''} (mm'''σ|m''m'σ')
J̄ = ──────────────────────────────────────────────────────────────────                          (B6)
        Σ_σ Σ_{m≠m'} [ N^{σσ}_m N^{σσ}_{m'} + N^{σ-σ}_m N^{-σσ}_{m'} ]
```

**These two equations are reproduced term-for-term in `new_ns_e0_nc_1`.** See §5.1 — the mapping is
one-to-one, including the minus sign in the (B5) denominator. If you need one sentence linking the
PDF to this file, it is: *"the noncollinear ACBN0 kernel of `new_ns_e0_nc_1` is a literal
transcription of Eqs. (B5)–(B6) of PRB 96, 245133."*

### 1.4 Hubbard forces (Eqs. 11–13)

```
F^U_{α,i} = − Σ_{I,m,m',σ} (U/2)( δ_{mm'} − 2 n^{I,n,l,σ}_{mm'} ) ∂n^{I,n,l,σ}_{mm'} / ∂R_{α,i}     (11)
```

with the derivative moved onto the *orbitals* rather than the wavefunctions (Eq. 13) to suppress the
egg-box effect — a real-space-grid concern that does not apply to a PW code. In the fork this
corresponds to `PW/src/force_hub_v.f90` / `stres_hub_v.f90`, which are *outside* the scope of this
review but consume exactly the same `ne0` array.

### 1.5 Numerical results worth borrowing

| Material | `U_eff` (TM 3d), this work / orig. ACBN0 | `U_eff` (O 2p), this work / orig. |
|---|---|---|
| MnO | 4.68 / 4.67 | 5.18 / 2.68 |
| NiO | 6.93 / 7.63 | 2.87 / 3.0 |
| TiO₂ | 0.96 / 0.15 | 10.18 / 7.34 |
| ZnO | 13.3 / 12.8 | 5.95 / 5.29 |

Useful sanity targets if you ever run this fork's `acbn0_type=1` on NiO: **you should land near
6.9–7.6 eV on Ni-3d and ~3 eV on O-2p.** If you get 15 eV, your projector or your `pao_3g` fit is
wrong, not the physics.

---

## 2. What `new_ns_v.f90` contains — file map

| Lines | Routine | Role | Driven by |
|---:|---|---|---|
| 10–230 | `new_ns_v(ns)` | Bare `n^{IJσ}_{mm'}(R)`, collinear. **No `U`/`V` update.** | `lda_plus_v .AND. .NOT.lacbn0`, or converged ACBN0 |
| 233–633 | `new_ns_a0(ns)` | On-site ACBN0 → `ehub_u`. Eqs. (7)–(10) | `acbn0_type == 1` |
| 636–1806 | `new_ns_e0_1(ns)` | eACBN0, **`(l,l')`-resolved** `U` and `V` → `ehub_v` | `acbn0_type == 2` |
| 1809–2961 | `new_ns_e0_2(ns)` | eACBN0, **`(m,m')`-resolved** `V` ("orbital-dependent") | `acbn0_type == 3` |
| 2964–3165 | `symmetrize_n4(n)` | Point-group symmetrisation of a rank-4 on-site array | internal |
| 3168–3397 | `symmetrize_n4m(n)` | Same, for the `dm1(ntyp,4)` family | internal |
| 3400–3777 | `symmetrize_n5(n)` | Same, for `(m,m',σ,I,J,R)` intersite arrays | internal |
| 3780–3818 | `orberi(ehub_l, eri_l)` | Stateful `m → ERI lm-label` walker | internal |
| 3820–4085 | `new_ns_v_nc(ns)` | Noncollinear bare `n^{IJσσ'}_{mm'}(R)` | `noncolin .AND. .NOT.lacbn0` |
| 4087–5200 | `new_ns_e0_nc_1(ns)` | **Noncollinear eACBN0** — Eqs. (B3)–(B6) | `noncolin .AND. acbn0_type == 2` |
| 5202–6748 | `symmetrize_n{4,4m,5}_nc` | Noncollinear symmetrisation (spin rotations via `d_spin_ldav`) | internal |

`new_ns_a0_nc` and `new_ns_e0_nc_2` are **referenced but do not exist** — the call sites in
`sum_band.f90:159,163` and `electrons.f90:1193,1197` are commented out. Noncollinear therefore only
supports `acbn0_type = 2`.

Roughly **70 % of the file is boilerplate**: the four `symmetrize_*` routines are 3800 of the 6751
lines and consist almost entirely of a 16-way `IF/ELSE IF` ladder over `(l1,l2) ∈ {1..4}²` selecting
which of `d1`/`d2`/`d3` to sandwich. §7.2 shows how that collapses to ~40 lines.

---

## 3. Symbol dictionary — paper ↔ code

### 3.1 Physics quantities

| Paper | Code | Declared in | Shape / meaning |
|---|---|---|---|
| `\|φ_{I,n,l,m}⟩` | `wfcv(:,:)` | `ldaV.f90:12` | Löwdin-orthogonalised atomic wfcs × S, only the Hubbard-corrected channels. Built by `orthowfcv` (`orthoatwfc_v.f90`), buffered on unit `iunhubv` |
| `⟨φ_{I,n,l,m}\|ψ^σ_{nk}⟩` | `proj%k(ofstv(na)+m, ibnd)` | local | `= calbec(npw, wfcv, evc, proj)`. **Note: `proj%k(x,i)` holds `⟨φ_x\|ψ_i⟩`**, so `⟨ψ\|φ⟩ = CONJG(proj%k)` |
| `w_k f^σ_{nk}` | `wg(ibnd,ik)` | `wvfct` | QE folds k-weight × occupation into one array |
| `n^{I,n,l,σ}_{mm'}` (Eq. 3) | `n(1)%O(m1,m2,is,na,na)` → `ns(...,1)` → `rho%ne0` | local / `scf_mod.f90:68` | on-site block |
| `n^{IJσ}_{mm'}(R)` (**not in paper**) | `n(nnc)%O(m1,m2,is,na1,na2)` → `ns(:,:,:,:,:,nnc)` | ditto | `nnc=1` is `R=0` |
| `N̄^{I,n,l,σ}_{ψnk}` (Eq. 8) | `num(ibnd)%O(l1,l2,is,nt1,nt2)` | local | index by `(l+1)`, **not** by `m`; band-resolved; deallocated each `ik` |
| `P̄^{I,n,l,σ}_{mm'}` (Eq. 7) | `p%O` (`new_ns_a0`); `dm0%O`, `dm1(nt,l)%O`, `p(nnc)%O` (`new_ns_e0_*`) | local | see §3.3 |
| `(mm'\|m''m''')` (Eq. B2) | `compute_eri(lm1,lm2,lm3,lm4, x1..x4, c1..c4)` | `compute_eri.f90:31` | chemists' notation, returns **Ry** (`×2` at line 174) |
| `Ū` numerator | `A0%U(l1,l2,na)` / `E0(i)%V(...)` | local | |
| `Ū` denominator | `HF%U(l1,l2,na)` / `HF(i)%V(...)` | local | |
| `J̄` numerator / denominator | `A0%J` / `HF%J`, `E0(i)%J` / `HF(i)%J` | local | |
| `U^eff = Ū − J̄` | `ehub_u(na,m1,m2)`; `ehub_v(na,na,m1,m2,1)` | `ldaV.f90:37,36` | in **Ry**; printed in eV via `rytoev` |
| `V^{IJ}` (**not in paper**) | `ehub_v(na1,na2,m1,m2,nnc)` | `ldaV.f90:36` | |
| `V^{I,n,l,σ}_{m,m'}` (Eq. 5) | `v%ne0(m1,m2,is,na1,na2,nnc)` | `v_hubbard_v.f90` | applied by `vhpsi_v.f90` |

### 3.2 Bookkeeping variables (all in `PW/src/ldaV.f90`)

| Variable | Meaning |
|---|---|
| `ehub_lm(nt)` | number of Hubbard `m`-channels on type `nt`, `= Σ_{chosen l} (2l+1)`. `0` ⇒ atom is not corrected |
| `ehub_lm_max` | `MAX_nt ehub_lm(nt)` — leading dimension of every occupation array |
| `ehub_l(m,nt)` | **`l+1`** for channel `m` (so `1=s, 2=p, 3=d, 4=f`). Blocks are contiguous in `m` |
| `ofstv(na)` | offset of atom `na` inside `wfcv` / `proj` |
| `nwfcv` | total number of columns in `wfcv` |
| `ehub_l_choice(nt,1:6)` | input mask: which `l` channels get corrected (`1=s,2=p,3=d,4=f,5=semicore d`). Rewritten in-place by `init_lda_plus_v` from a 0/1 mask into a `PP_CHI` index |
| `ehub_csize = 10`, `ics = 11` | cell-index array is `(-10:10)³` stored as `1:21`; `ics` is the origin |
| `ehub_nn_vect(a,b,c,na1,na2)` | 1 if `na2` in cell `(a-11,b-11,c-11)` is within `ehub_nn_distance` of `na1` in the home cell. **Diagonal `(ics,ics,ics,na,na)` is 0** (`hubbard_nns.f90:46` rejects `r² < 1e-8`) |
| `ehub_nn_dim(1:3)` | max `\|n_i\|` over all neighbour pairs — the loop bound |
| `ehub_nnnc` | number of cells hosting ≥1 pair, **+ the home cell**, computed in `input.f90:2142-2180` |
| `nnc` | running cell index, `1 = home`, `2..ehub_nnnc` in lexicographic `(n1,n2,n3)` order |
| `ehub_mixing` (0.7) | linear mixing of `U`/`V` between SCF iterations |
| `ehub_conv_thr` (1e-8), `stop_acbn0` | freeze `U`/`V` once `max\|ΔV\| < sqrt(ehub_conv_thr)` |
| `remove_ehub_u(nt,l)` | 0/1 multiplier to kill `U` on a channel |
| `remove_ehub_v(nt,l)` | **declared, broadcast, never used** — see Finding #4 |
| `acbn0_type` | `1` = on-site `U` only; `2` = `U+V`, `(l,l')`-resolved; `3` = `U+V`, `(m,m')`-resolved |
| `pot_fix` (`ehub_pot_fix`) | keep `U`/`V` fixed, recompute occupations only |

### 3.3 The three flavours of `P̄` in `new_ns_e0_*` — the subtlest thing in the file

| Code object | Renormalisation weight used | Represents |
|---|---|---|
| `dm0%O(m1,m2,σ,na)` | `num(ibnd)%O(l(m1), l(m2), σ, nt, nt)` | `P̄^{II}` for the **on-site `U`** term |
| `dm1(nt2,l2)%O(m1,m2,σ,na1)` | `num(ibnd)%O(l(m1), l2, σ, nt1, nt2)` | `P̄^{I}` renormalised with the **cross-species/cross-channel** `N̄^{IJ}` — the direct (Hartree-like) half of `V^{IJ}` |
| `p(nnc)%O(m1,m2,σ,na1,na2)` | `num(ibnd)%O(l(m1), l(m2), σ, nt1, nt2)` | genuinely **intersite** `P̄^{IJ}(R)` — the exchange half of `V^{IJ}` |

`dm1` is a `(ntyp,4)`-indexed family: `dm1(nt,l)` is "the on-site `P̄` of every atom, screened as if it
were interacting with channel `l` of species `nt`". That is the fork's answer to a question the paper
never asks: *what `N̄` do you use when the two orbitals sit on different atoms?*

---

## 4. Routine-by-routine walkthrough (collinear)

### 4.1 `new_ns_v(ns)` — lines 10–230 — the pure occupation-matrix builder

**Signature:** `REAL(DP), INTENT(OUT) :: ns(ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc)`

Local storage is a derived-type array `TYPE(occmat_5) :: n(ehub_nnnc)` with an allocatable component
`%O(ehub_lm_max, ehub_lm_max, nspin, nat, nat)`. One `%O` per neighbour cell.

| Lines | What happens | Paper |
|---|---|---|
| 58–66 | `start_clock`, allocate `proj(nwfcv,nbnd)` and zero every `n(i)%O` | — |
| 68–80 | k-loop; `get_buffer(evc)`, `get_buffer(wfcv)`, then `calbec(npw, wfcv, evc, proj)` ⇒ `proj = ⟨φ\|S\|ψ⟩` | the projections in Eq. (3) |
| 77 | hard `errore` if `Hubbard_projectors == 'pseudo'` — **dead branch**, `orthowfcv` forcibly sets `'ortho-atomic'` (`orthoatwfc_v.f90:33`) | paper's orthogonality assumption is thereby satisfied by construction |
| 83–122 | **home-cell block.** `DO na1; DO na2 = na1, nat` with `IF (na1/=na2 .AND. ehub_nn_vect(ics,ics,ics,na1,na2)==0) CYCLE`. Accumulates `n(1)%O(m1,m2,σ,na1,na2) += wg · Re[ proj_k(J,m2) · conj(proj_k(I,m1)) ]` | **Eq. (3)** with `I≠J` allowed |
| 124–143 | fills the transposed `(na2,na1)` half. Skips `na1==na2` (already complete) | Hermiticity |
| 145–203 | **neighbour-cell block.** Triple loop over `(n1,n2,n3)` in `±ehub_nn_dim`, skipping `(0,0,0)`. Phase `arg = 2π k·(n1 a1 + n2 a2 + n3 a3)`, `kphase = e^{i·arg}`, accumulate `n(nnc)%O += Re[kphase · wg · proj_k(J,m2) · conj(proj_k(I,m1))]` | Bloch sum of Eq. (3) over `R` |
| 199 | `IF (MAXVAL(ehub_nn_vect(...)) == 1) nnc = nnc + 1` — **the cell counter only advances for cells that actually host a pair**. This idiom appears **17× in this file alone and 30× across 7 files**, and every copy must stay bit-identical with `input.f90:2148-2178`, `vhpsi_v.f90:175`, `v_hubbard_v.f90`, `force_hub_v.f90`, `stres_hub_v.f90`, `io_rho_xml.f90`, `ldaV.f90`, `hubbard_get_nnc` | — |
| 209–212 | `mp_sum(n(i)%O, inter_pool_comm)` — k-point-pool reduction. `calbec` already did the G-vector (`intra_bgrp`) sum | — |
| 214 | `IF (nspin == 1) n(i)%O = 0.5 · n(i)%O` — undoes QE's spin-degeneracy factor of 2 in `wg`. **Dead code:** `init_lda_plus_v` aborts on `nspin==1` (`ldaV.f90:104`) | — |
| 217 | `CALL symmetrize_n5(n)` | crystal-symmetry averaging |
| 220–224 | copy `n(i)%O → ns(:,:,:,:,:,i)`, deallocate | — |

**Taking `DBLE(...)` per k-point** (lines 179 and 187) discards `Im n^{IJ}(R)` before the BZ sum. That is legitimate
only if the k-mesh is time-reversal complete, which QE guarantees. Same convention as QE's own
`new_nsg.f90`. Fine, but worth a comment in the source.

### 4.2 `new_ns_a0(ns)` — lines 233–633 — on-site ACBN0, `acbn0_type = 1`

This is the routine that maps **one-to-one onto Eqs. (7)–(10)**. Output: `ns(...,nat)` (rank 4, no
site pairs) and, as a side effect on the module, `ehub_u(na,m1,m2)`.

| Lines | Code | Paper |
|---|---|---|
| 310–341 | `n%O(m1,m2,σ,na) += wg · Re[proj(m2)·conj(proj(m1))]` | **Eq. (3)** |
| 343–376 | `num(ibnd)%O(l,l,σ,nt) += \|proj(ofstv(na)+m,ibnd)\|²`, summed over **all atoms `na` of type `nt`** and all `m` in the `l` block. **No `wg` factor** | **Eq. (8)**, inner sums |
| 378–395 | `num(i)%T(:,:,:,nt1) += num(i)%O(:,:,:,nt2)` for every `nt2` with `upf(nt2)%psd == upf(nt1)%psd` | **Eq. (8)**, the `{I}` sum over *"atoms of the same type as atom I"* — implemented as *same chemical element*, so `Ni1`/`Ni2` in an AFM cell are correctly merged |
| 399–441 | `p%O(m1,m2,σ,na) += num(ibnd)%O(l1,l2,σ,nt) · wg · Re[proj(m2)·conj(proj(m1))]` | **Eq. (7)** exactly: `P̄ = Σ_nk w_k f_nk N̄_{ψnk} ⟨ψ\|φ_m⟩⟨φ_m'\|ψ⟩` |
| 452–464 | `symmetrize_n4(n)`, `symmetrize_n4(p)`; `ns = n%O` | — |
| 468–532 | **`Ū`/`J̄` numerators.** Quadruple `m1..m4` loop, `m3` restricted to `l(m3)==l(m1)`, `m4` to `l(m4)==l(m2)`. Accumulate into `A0%U(l1,l2,na)` and `A0%J(l1,l2,na)` | **Eqs. (9)/(10) numerators** |
| 536–563 | **denominators** `HF%U`, `HF%J` from the diagonal `ns(m,m,σ,na) = N^σ_m` | **Eqs. (9)/(10) denominators** |
| 570–598 | `ehub_u_ = A0%U/HF%U − A0%J/HF%J`, times `remove_ehub_u` masks | **`U^eff = Ū − J̄`** |
| 596–598 | `ehub_u = (1−ehub_mixing)·ehub_u_old + ehub_mixing·ehub_u_new` | *not in paper* — SCF damping |
| 605–620 | symmetrise `ehub_u(m1,m2) ← ½(U_{m1m2}+U_{m2m1})` | — |

**Line-by-line verification of the spin sums against Eq. (9):**

```fortran
! lines 507-513, the generic (non-s) branch
A0%U(l1,l2,na) += ( p%O(m3,m1,1,na)*p%O(m4,m2,1,na)     ! α=↑ β=↑
                  + p%O(m3,m1,2,na)*p%O(m4,m2,2,na)     ! α=↓ β=↓
                  + p%O(m3,m1,1,na)*p%O(m4,m2,2,na)     ! α=↑ β=↓
                  + p%O(m3,m1,2,na)*p%O(m4,m2,1,na) )   ! α=↓ β=↑
                * compute_eri(l3,l1,l4,l2, τ,τ,τ,τ, c3,c1,c4,c2)
```
All four `(α,β)` ⇒ `Σ_{αβ}` of Eq. (9). ✔
Index identification: `(m,m') = (m3,m1)`, `(m'',m''') = (m4,m2)`, and the ERI argument order
`(l3,l1,l4,l2)` is literally `(m m' | m'' m''')`. ✔

```fortran
! lines 516-522
A0%J(l1,l2,na) += ( p%O(m3,m1,1,na)*p%O(m4,m2,1,na)
                  + p%O(m3,m1,2,na)*p%O(m4,m2,2,na) )   ! Σ_α only, same spin
                * compute_eri(l3,l2,l4,l1, ...)         ! = (m m''' | m'' m')
```
Same-spin only + swapped ERI ⇒ Eq. (10). ✔

```fortran
! lines 545-563 denominators
m1 == m2 : HF%U += N↑_m·N↓_m + N↓_m·N↑_m                          ! the Σ_{{m}} N^α N^{-α} term at m=m'
m1 /= m2 : HF%U += N↑N↑ + N↓N↓ + N↑N↓ + N↓N↑                       ! Σ_{m≠m'} N^α N^α  +  Σ_{{m}} N^α N^{-α}
           HF%J += N↑N↑ + N↓N↓                                     ! Σ_{m≠m'} N^α N^α
```
Sum over both branches reproduces Eq. (9)'s denominator exactly, and the `J̄` denominator is
Eq. (10)'s. ✔ **This routine is a correct implementation of the paper's Eqs. (7)–(10).**

**The `l=0` special case (lines 497–505, 578–584).** When both orbitals are `s` (`ehub_l == 1`):
only the opposite-spin products enter `A0%U`, and `A0%J`/`HF%J` are never touched. Correct and
*necessary*: an `s` manifold has a single `m`, so `Σ_{m≠m'}` is empty and Eq. (10) would be `0/0`.
Physically, same-spin `s`–`s` self-interaction is exactly cancelled by exchange. Undocumented in the
source — **add a comment**.

**Consequence worth stating in your paper:** `A0%*`/`HF%*` are indexed only by `(l1,l2)`, so
`ehub_u(na,m1,m2)` is **block-constant in `m` within each `l` block**. For a single corrected `d`
manifold this is one scalar `U^eff` per atom — precisely the paper's single `U^eff_{I,n,l}`. The
`(m1,m2)` storage is forward compatibility for `acbn0_type=3`, not extra physics.

### 4.3 `new_ns_e0_1(ns)` — lines 636–1806 — eACBN0, the `+V` routine

This is the routine your project actually cares about. Same skeleton as `new_ns_a0`, but every
quantity carries a site pair `(I,J)` and a lattice vector `R`.

**Early exit (lines 715–725):**
```fortran
IF ( stop_acbn0 ) THEN
   IF ( MAXVAL(ehub_v_diff(:,:,:,:,:)) < SQRT(ehub_conv_thr) ) THEN
      CALL new_ns_v(ns)      ! occupations only — U,V frozen
      GOTO 1                 ! jump to write_ehub_uv
   ENDIF
   ehub_v_diff = 0.0_DP
ENDIF
```
Once `V` stops moving, the expensive ERI machinery is bypassed for the rest of the SCF. Good design.
(The `GOTO 1` is ugly but harmless — label `1` is on `write_ehub_uv` at line 1800.)

**Stage 1 — occupations (lines 756–877).** Identical to `new_ns_v`. Produces `n(1..ehub_nnnc)%O`.

**Stage 2 — `N̄` (lines 879–962).** Two sub-stages:

```fortran
! 879-914 : per-band, per-(type,l) projection weight, summed over atoms of that type
num(ibnd)%O(l,l,σ,nt,nt) += |proj(ofstv(na)+m, ibnd)|²

! 916-931 : merge types sharing a chemical element  →  Eq. (8)'s {I} sum
num(ibnd)%T(:,:,:,nt1,nt1) += num(ibnd)%O(:,:,:,nt2,nt2)   ! if psd(nt1)==psd(nt2)

! 936-948 : THE INTERSITE GENERALISATION — not in the paper
num(ibnd)%T(l1,l2,:,nt1,nt2) = num(ibnd)%T(l1,l1,:,nt1,nt1) + num(ibnd)%T(l2,l2,:,nt2,nt2)
IF ( psd(nt1) == psd(nt2) ) &
   num(ibnd)%T(l1,l2,:,nt1,nt2) = num(ibnd)%T(l1,l2,:,nt1,nt2) / 2.0_DP
```

Read that carefully:
* **different elements** ⇒ `N̄^{IJ} = N̄^I + N̄^J` (a *sum*): "how much of band `nk` lives in the
  **union** of the two manifolds". This is the honest extension of Eq. (8)'s Mulliken reading.
* **same element** ⇒ `(N̄^I + N̄^J)/2` (an *average*): the union is the same manifold, so the `/2`
  removes the double count. For `l1==l2` this correctly returns `N̄` itself.
* the commented alternative `SQRT(N̄^I · N̄^J)` is the geometric mean attributed in the source to
  Tancogne-Dejean & Rubio.

**This sequential in-place update is safe but fragile.** The `(l1,l2)` loop writes `T(l1,l2)` while
reading `T(l1,l1)` and `T(l2,l2)`. It happens to be correct because the diagonal is self-preserving
(`T(l,l) = 2·T(l,l)/2`) and is written before any later iteration reads it. Change the loop order and
you silently corrupt `N̄`. **Add a warning comment or use a separate output array.**

**Stage 3 — the three `P̄` families (lines 963–1195).**

| Lines | Object | `N̄` used | Purpose |
|---|---|---|---|
| 967–1006 | `dm0%O(m1,m2,σ,na)` | `num(l1,l2,σ,nt,nt)` | on-site `U` |
| 1008–1049 | `dm1(nt,l)%O(m1,m2,σ,na1)` | `num(ehub_l(m1,nt1), l, σ, nt1, nt)` | direct part of `V^{IJ}` |
| 1052–1119 | `p(1)%O(m1,m2,σ,na1,na2)` | `num(l1,l2,σ,nt1,nt2)` | exchange part, `R = 0` |
| 1121–1191 | `p(nnc)%O(...)` + `kphase` | ditto | exchange part, `R ≠ 0` |

**Stage 4 — reductions + symmetrisation (lines 1197–1243).** `mp_sum` over pools, optional `nspin==1`
halving, then `symmetrize_n5(n)`, `symmetrize_n4(dm0)`, `symmetrize_n4m(dm1)`, `symmetrize_n5(p)`.
Symmetrising `P̄` is legitimate because `N̄` is a symmetry scalar per band.

**Stage 5 — on-site `U` (lines 1245–1385).** Byte-for-byte the `new_ns_a0` kernel with `p%O → dm0%O`
and `A0 → E0(1)`, writing `ehub_v(na,na,m1,m2,1)`. **So the on-site block of `ehub_v` *is* the
Hubbard `U`.** `ehub_u` is not used in this branch at all.

**Stage 6 — intersite `V`, `R = 0` (lines 1387–1562).** The core new physics:

```fortran
E0(1)%V(l1,l2,na1,na2) +=
   (  dm1(nt2,l2)%O(m3,m1,1,na1) * dm1(nt1,l1)%O(m4,m2,1,na2)     ! ↑↑ ┐
    + dm1(nt2,l2)%O(m3,m1,2,na1) * dm1(nt1,l1)%O(m4,m2,2,na2)     ! ↓↓ │ Σ_{αβ} direct
    + dm1(nt2,l2)%O(m3,m1,1,na1) * dm1(nt1,l1)%O(m4,m2,2,na2)     ! ↑↓ │
    + dm1(nt2,l2)%O(m3,m1,2,na1) * dm1(nt1,l1)%O(m4,m2,1,na2)     ! ↓↑ ┘
    -            p(1)%O(m3,m2,1,na1,na2) * p(1)%O(m4,m1,1,na2,na1)  ! ┐ Σ_α exchange
    -            p(1)%O(m3,m2,2,na1,na2) * p(1)%O(m4,m1,2,na2,na1)) ! ┘
 * compute_eri(l3,l1,l4,l2, τ_I,τ_I, τ_J,τ_J, c3,c1,c4,c2)
```

With `(m,m',m'',m''') = (m3,m1,m4,m2)` this reads

```
Numerator^{IJ} = Σ_{mm'm''m'''} [ Σ_{αβ} P̄^{Iα}_{mm'} P̄^{Jβ}_{m''m'''} − Σ_α P̄^{IJα}_{mm'''} P̄^{JIα}_{m''m'} ] (mm'|m''m''')
```

**This is exactly Eq. (6) with the site labels split.** To see the equivalence of the exchange term,
relabel the dummies `m' ↔ m'''` in Eq. (6)'s second term:
`Σ P̄_{mm'} P̄_{m''m'''} (mm'''|m''m')  ≡  Σ P̄_{mm'''} P̄_{m''m'} (mm'|m''m''')`.
Same integral, indices moved from the ERI to the density matrices. The code uses this second form so
that a **single** ERI serves both terms — which is why only `compute_eri(l3,l1,l4,l2)` appears here
while `new_ns_a0` needs two. Clean and deliberate.

Denominator (lines 1477–1505) is the same expression at the bare mean-field level:
```
HF(1)%V(l1,l2,na1,na2) += Σ_{αβ} N^{Iα}_{m1} N^{Jβ}_{m2} − Σ_α n^{IJα}_{m1m2} n^{JIα}_{m2m1}
```

and (line 1550):
```fortran
ehub_v_(na1,na2,m1,m2) = ( E0(1)%V(l1,l2,na1,na2) / HF(1)%V(l1,l2,na1,na2) ) / 2.0_DP
```

**The `/2` has no counterpart in the paper.** It compensates the fact that `v_hubbard_v` sums over
*ordered* pairs `(I,J)` and `(J,I)`, so each physical bond is counted twice. Definitional, not
derivable from Eqs. (6)–(10). **Flag it explicitly if you write this up.**

**Stage 7 — intersite `V`, `R ≠ 0` (lines 1566–1758).** Same expressions with:
* `ntau(:,na2) = n1·a1 + n2·a2 + n3·a3 + tau(:,na2)` — the image atom's position, fed to `compute_eri`;
* the exchange pair becomes `p(nnc)%O(m3,m2,na1,na2) · p(nnc)%O(m1,m4,na1,na2)` — using the
  `m`-transpose of the *same* `(I,J,R)` block instead of `(J,I,−R)`, which is correct because
  `P̄^{JI}_{m4m1}(−R) = [P̄^{IJ}_{m1m4}(R)]* = P̄^{IJ}_{m1m4}(R)` for the real-part-only storage;
* likewise `HF(nnc)%V` uses `ns(m1,m2,...,nnc)²`.

**Stage 8 — Hermiticity + cleanup (lines 1760–1800).** Only the `nnc=1` block is explicitly
symmetrised (`ehub_v(I,J,m1,m2,1) ← ½[V_{IJ,m1m2} + V_{JI,m2m1}]`). For `nnc>1` no symmetrisation is
applied — **and none is needed**: `E0(nnc)%V` and `HF(nnc)%V` are already invariant under
`(I,m1,R) ↔ (J,m2,−R)` because (a) `(12|34) = (34|12)` and the ERI is translation invariant,
(b) the `dm1` product is symmetric under the swap, and (c) the `p(nnc)` product maps onto itself.
So `V^{IJ}(R) = V^{JI}(−R)` holds by construction and `vhpsi_v` builds a Hermitian operator. Worth a
comment in the source; a future reader will otherwise assume it is an omission.

### 4.4 `new_ns_e0_2(ns)` — lines 1809–2961 — `acbn0_type = 3`

A near-clone of `new_ns_e0_1`. Verified differences (`diff` of the two extracted bodies):

1. `E0(i)%V/J`, `HF(i)%V/J` are dimensioned `(ehub_lm_max, ehub_lm_max, nat, nat)` and indexed by
   `(m1,m2)` instead of `(4,4,...)` indexed by `(l1,l2)` ⇒ **`V` becomes genuinely orbital-resolved**
   (`V_{d_{xy}, p_x} ≠ V_{d_{z²}, p_x}`).
2. The whole `stop_acbn0` / `ehub_v_diff` convergence machinery is **absent** — `acbn0_type=3`
   cannot freeze `U`/`V` and will run the full ERI cost every SCF step.
3. `num(...)%T(:,:,current_spin,...)` instead of `(:,:,:,...)` — cosmetic (the other spin slot is
   zero at this point anyway).
4. `P^{IJ}` construction includes `na1==na2` (`IF (na1 /= na2 .AND. ...) CYCLE`), so `p(1)` also holds
   the on-site block.
5. `E0(1)%V = 0.0_DP` before the intersite accumulation was **dropped**. Not a bug — the on-site
   entries live at `(na,na)` and the intersite loop never touches those (`ehub_nn_vect` diagonal is
   0) — but it is a landmine if anyone ever makes the diagonal a neighbour.
6. The on-site `U` block still writes `E0(1)%V(ehub_l(m1),ehub_l(m2),na,na)`, i.e. it indexes an
   `m`-dimensioned array with `l` values. In-bounds always (`ehub_lm(nt) = Σ(2l+1) ≥ l+1`), and
   internally consistent, but it means **`acbn0_type=3` is `(m,m')`-resolved for `V` and still
   `(l,l')`-resolved for on-site `U`.** Undocumented asymmetry.

### 4.5 `symmetrize_n4` / `n4m` / `n5` — lines 2964–3777

All three implement the same textbook operation:

```
n^{sym}_{m1m2}(I)  = (1/N_sym) Σ_S Σ_{m3m4} D^{l1}_{m3m1}(S) · n_{m3m4}(S·I) · D^{l2}_{m4m2}(S)
n^{sym}_{m1m2}(I,J,R) = (1/N_sym) Σ_S Σ_{m3m4} D^{l1}_{m3m1}(S) · n_{m3m4}(S·I, S·J, S·R) · D^{l2}_{m4m2}(S)
```

with `D^1 = d1`, `D^2 = d2`, `D^3 = d3` from `symm_base` (real-spherical-harmonic rotation matrices),
`S·I = irt(isym,na1)`, and for the pair case `(S·I, S·J, S·R)` from `hubbard_symmetry`
(`hubbard_nns.f90:85-238`, a lightly modified `symonpair`).

The `l_separation(l,nt)` array maps a global channel index `m` to its position inside its own `l`
block, so `d1`/`d2`/`d3` get local indices `1..2l+1`. Necessary because `ehub_lm` concatenates
several `l` manifolds into one index range.

Both routines finish with an explicit `n(m1,m2) ← ½[n(m1,m2) + n(m2,m1)]` symmetrisation, and
`symmetrize_n5` additionally zeroes any `(I,J,R=0)` block that is not a listed neighbour.

**See Findings #1 and #7 — there are two real defects in this family.**

### 4.6 `orberi(ehub_l, eri_l)` — lines 3780–3818

A **stateful iterator**, not a pure function. Given `ehub_l = l+1` and the *current* `eri_l`, it
returns the *next* ERI label in that `l` block:

```
l = 0 (s) : always 1
l = 1 (p) : 2 (p_z) → 3 (p_x) → 4 (p_y) → 2 → …
l = 2 (d) : 5 (d_{z²}) → 6 (d_{xz}) → 7 (d_{yz}) → 8 (d_{x²−y²}) → 9 (d_{xy}) → 5 → …
```

Every call site initialises `l1 = 0` (or `l3 = 0`) immediately before the corresponding `m` loop and
relies on the loop visiting the block's channels in order. Fragile — see Finding #6 — but **the
ordering it produces is correct**: it matches QE's `ylmr2` real-harmonic ordering exactly
(`l=1: p_z, p_x, p_y`; `l=2: 3z²−r², xz, yz, x²−y², xy`), which is what `atomic_wfc` puts into
`wfcv`. `f` orbitals (`ehub_l == 4`) fall through every branch and leave `eri_l` **unchanged** — so
`compute_eri` would silently be handed a stale label. `exponents_and_normalization_constants`
(`compute_eri.f90:181`) has no `CASE` beyond 9 either. **`f`-electron systems are not supported and
fail silently, not loudly.**

### 4.7 The Coulomb integrals — `compute_eri.f90` + `pao_3g.f90`

```
compute_eri(lm1,lm2,lm3,lm4, x1,x2,x3,x4, c1,c2,c3,c4)
  = 2 · ∫∫ d³r₁ d³r₂  φ*_{lm1}(r₁−x1) φ_{lm2}(r₁−x2) φ*_{lm3}(r₂−x3) φ_{lm4}(r₂−x4) / |r₁−r₂|
```
(the `2` converts Hartree → Ry, `compute_eri.f90:175`) — i.e. **Eq. (B2) of the paper**, in chemists'
notation.

Each `φ` is a **3-term Cartesian Gaussian contraction**: `pao_3g(psd, l, c)` returns
`c(1:3)` = contraction coefficients, `c(4:6)` = exponents, hard-coded per element (H … ~Bi) and per
`l`. The angular part is expanded into Cartesian monomials by
`exponents_and_normalization_constants`, with `d_{z²}` needing 3 terms (`2z² − x² − y²`) and
`d_{x²−y²}` needing 2. The inner `coulomb_repulsion` is a C routine (`BIND(C)`) — the standard
Gaussian ERI kernel.

**Three things you must say in a paper about this:**
1. The ERI is **not** computed from the actual pseudo-atomic `PP_CHI` radial function; it uses a
   fitted 3-Gaussian STO-3G-like proxy that depends only on `(element, l)`. Two different
   pseudopotentials for Ni give the *same* `pao_3g` radial function. This is the largest uncontrolled
   approximation in the whole `U`/`V` evaluation and it is invisible in the output.
2. `pao_3g` fails hard (`errore`) on elements not in its table — check coverage before running.
3. Some tabulated fits are numerically nasty: e.g. one entry has `c(2) = −309.43`, `c(3) = +309.50`
   with near-degenerate exponents `0.3412` / `0.3410` — a catastrophic-cancellation fit. Any `U`
   involving that element deserves a hand check.

---

## 5. The noncollinear family

### 5.1 `new_ns_e0_nc_1` — the literal implementation of Eqs. (B5)/(B6)

Spin packing: `j = npol*(is1−1) + is2`, i.e. `ns(...,1)=↑↑`, `2=↑↓`, `3=↓↑`, `4=↓↓` (`new_ns_v.f90:4036-4037`).
Projector index arithmetic (`new_ns_v.f90:4002-4004`):

```fortran
proj( ofstv(na) + m + counter_lm(m,nt) + (2*(ehub_l(m,nt)-1)+1)*(is-1), ibnd )
```
`counter_lm(m,nt)` (built in `ldaV.f90:379-390`) skips the up-spin half of every *preceding* `l`
block; the last term selects the spin half of the *current* block. Necessary because
`atomic_wfc_v_nc` lays out each `l` manifold as `[(2l+1) up | (2l+1) down]`.

**`N̄` (lines 4309–4352, active branch at 4340–4346):** summed over `is1` with `proj·CONJG(proj)` at the *same* spinor index ⇒
spin-traced, **no spin dimension on `num`** ⇒ **Eq. (B4)**. ✔ The commented-out alternative just
above it (with the off-diagonal `is1≠is2` terms) carries the note *"question: How to consider
renormalization number in noncollinear case? The definition of Nicolas Tancogne-Dejean & Angel
Rubio"* — again an explicit, deliberate divergence from the paper.

**`Ū` numerator (lines 4653–4670):**
```fortran
E0(1)%V += DBLE( dm0(m3,m1,↑,↑)·dm0(m4,m2,↑,↑) + dm0(m3,m1,↓,↓)·dm0(m4,m2,↓,↓)
                + dm0(m3,m1,↑,↑)·dm0(m4,m2,↓,↓) + dm0(m3,m1,↓,↓)·dm0(m4,m2,↑,↑) ) · (m3 m1|m4 m2)
```
`Σ_{σσ'} P̄^{σσ} P̄^{σ'σ'}` ⇒ **Eq. (B5) numerator**. ✔

**`J̄` numerator (lines 4673–4681):**
```fortran
E0(1)%J += DBLE( dm0(m3,m1,↑,↑)·dm0(m4,m2,↑,↑) + dm0(m3,m1,↓,↓)·dm0(m4,m2,↓,↓)
                + dm0(m3,m1,↓,↑)·dm0(m4,m2,↑,↓) + dm0(m3,m1,↑,↓)·dm0(m4,m2,↓,↑) ) · (m3 m2|m4 m1)
```
`Σ_{σσ'} P̄^{σσ'} P̄^{σ'σ}` with the swapped ERI ⇒ **Eq. (B6) numerator**. ✔

**Denominators (lines 4701–4723):**
```fortran
m1 == m2 : HF%V += N^{↑↑}N^{↓↓} + N^{↓↓}N^{↑↑} − N^{↑↓}N^{↓↑} − N^{↓↑}N^{↑↓}
m1 /= m2 : HF%V += N^{↑↑}N^{↑↑} + N^{↓↓}N^{↓↓} + N^{↑↑}N^{↓↓} + N^{↓↓}N^{↑↑}
           HF%J += N^{↑↑}N^{↑↑} + N^{↓↓}N^{↓↓} + N^{↑↓}N^{↓↑} + N^{↓↑}N^{↑↓}
```
Summing both branches gives
`Σ_σ Σ_{mm'} N^{σσ}_m N^{−σ−σ}_{m'} + Σ_σ Σ_{m≠m'} N^{σσ}_m N^{σσ}_{m'} − Σ_σ Σ_m N^{σ−σ}_m N^{−σσ}_m`
⇒ **Eq. (B5) denominator, including the minus sign**. ✔ And `HF%J` ⇒ **Eq. (B6) denominator**. ✔

> **This is the tightest paper↔code correspondence in the file.** If you cite PRB 96, 245133 anywhere
> in your work, cite it here.

The intersite `V` blocks (lines 4769–5153) then generalise the same expressions to `(I,J,R)` in the
manner of §4.3.

### 5.2 `new_ns_v_nc` — lines 3820–4085

Noncollinear analogue of §4.1. Occupations are `COMPLEX` and *no* `DBLE()` is taken — the full
`n^{IJσσ'}(R)` is kept, correct for spinors. Hermiticity is imposed as
`n(m2,m1,is2,is1,na2,na1) = CONJG(n(m1,m2,is1,is2,na1,na2))` (line 3957). The Hermiticity *check* at
lines 4044–4076 is entirely commented out — see Finding #8.

---

## 6. Workflow — inputs, outputs, and where this sits in the SCF loop

### 6.1 Call graph

```
PW/src/electrons.f90 :: electrons_scf
  │
  ├─► c_bands()                          ! diagonalise H, incl. vhpsi_v (applies V_hub)
  │
  ├─► sum_band()                          ! PW/src/sum_band.f90:154-181
  │     └── IF (lda_plus_v):
  │           noncolin ─ lacbn0 ─ acbn0_type==2 ─► new_ns_e0_nc_1(rho%ne0_nc)
  │           noncolin ─ .NOT.lacbn0      ─────► new_ns_v_nc  (rho%ne0_nc)
  │           collinear─ lacbn0 ─ type==1 ─────► new_ns_a0    (rho%na0)
  │           collinear─ lacbn0 ─ type==2 ─────► new_ns_e0_1  (rho%ne0)     ◄── THE +V ROUTINE
  │           collinear─ lacbn0 ─ type==3 ─────► new_ns_e0_2  (rho%ne0)
  │           collinear─ .NOT.lacbn0      ─────► new_ns_v     (rho%ne0)
  │
  ├─► mix_rho()                           ! rho%ne0 IS part of the Broyden mixing vector
  │
  ├─► v_of_rho()                          ! PW/src/v_of_rho.f90:138-156
  │     └── v_hubbard_v(rho%ne0, v%ne0, ehubene)      ! E_UV and the potential matrix
  │
  └─► (convergence) ─► IF (pot_fix) recompute occupations ─► print_ks_energies
```

Two facts that matter:
1. **`U`/`V` are updated once per SCF iteration**, inside `sum_band`, *before* mixing. The Hubbard
   parameters therefore ride along with the electronic SCF rather than in an outer loop —
   `ehub_mixing = 0.7` is what keeps that stable.
2. **`rho%ne0` participates in density mixing** (`scf_mod.f90:707,753,842,1671-1740`), exactly like
   `rho%ns` in official QE. Its contribution enters the SCF error estimate.

### 6.2 `new_ns_v` I/O contract

**INPUT (all implicit, via `USE` and buffers):**

| Source | Object | Note |
|---|---|---|
| `wavefunctions` / unit `iunwfc` | `evc(npwx*npol, nbnd)` | current KS orbitals at `ik` |
| `ldaV` / unit `iunhubv` | `wfcv(npwx*npol, nwfcv)` | `S·φ^ortho`, written once by `orthowfcv` |
| `wvfct` | `wg(nbnd,nks)` | `w_k · f_{nk}` |
| `klist` | `xk`, `ngk`, `nks` | k-points |
| `cell_base` | `at`, `alat`; `ions_base`: `tau`, `ityp` | geometry |
| `ldaV` | `ehub_lm`, `ehub_l`, `ofstv`, `ehub_nn_vect`, `ehub_nn_dim`, `ehub_nnnc`, `ehub_csize` | Hubbard bookkeeping, set in `init_lda_plus_v` |
| `symm_base` | `d1,d2,d3`, `nsym`, `irt`, `s`, `ft` | symmetrisation |
| `ldaV` | `ehub_v`, `ehub_u` | **read-modify-write** in the ACBN0 variants (mixing) |
| `uspp_param` | `upf(nt)%psd` | element label, drives the `{I}` sum of Eq. (8) |

**OUTPUT:**

| Object | Shape | Consumer |
|---|---|---|
| `ns` → `rho%ne0` | `(ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc)` | `v_hubbard_v`, `mix_rho`, `force_hub_v`, `stres_hub_v`, `io_rho_xml` |
| `ns` → `rho%na0` (type 1) | `(ehub_lm_max, ehub_lm_max, nspin, nat)` | `v_hubbard_a0` |
| `ldaV::ehub_v(na1,na2,m1,m2,nnc)` | Ry | `v_hubbard_v`, `write_ehub_uv`, XML |
| `ldaV::ehub_u(na,m1,m2)` | Ry | `v_hubbard_a0` |
| `ldaV::ehub_v_diff` | — | early-exit test |
| stdout | `U[atom i, #orb, #orb]` / `V[atom i, atom j, ...]` tables in eV | `write_ehub_uv` (`ldaV.f90:617`) |

**Then, downstream:**

```fortran
! v_hubbard_v.f90 — the DFT+U+V energy actually minimised
E = Σ_I Σ_σ Σ_m [ ½U^I_{mm} n^{IIσ}_{mm} − Σ_{m'} ½U^I_{mm'} n^{IIσ}_{mm'} n^{IIσ}_{m'm} ]
  − ½ Σ_{I≠J} Σ_R Σ_σ Σ_{mm'} V^{IJ}_{mm'}(R) n^{IJσ}_{mm'}(R) n^{JIσ}_{m'm}(−R)

! the potential matrix
v%ne0(m1,m1,σ,I,I,1) += ½ U^I_{m1m1}
v%ne0(m1,m2,σ,I,J,R) -=   V^{IJ}_{m1m2}(R) · n^{JIσ}_{m2m1}(−R)

! vhpsi_v.f90 — applied to |ψ⟩
H|ψ⟩ += Σ_{I,J,R} e^{i k·R} |φ_{I,m1}⟩ v%ne0(m1,m2,σ,I,J,R) ⟨φ_{J,m2}|ψ⟩
```

The on-site part is **exactly Eq. (2) / Eq. (5) of the paper** with `U^eff → U_{mm'}`; the intersite
part is the Campo–Cococcioni DFT+U+V term the paper does not have.

### 6.3 Input namelist (`&SYSTEM`)

| Variable | Default | Meaning |
|---|---|---|
| `lda_plus_v` | `.false.` | enable the extended-Hubbard machinery |
| `lacbn0` | `.false.` | compute `U`/`V` self-consistently (forces `lda_plus_v = .true.`, `input.f90:2092`) |
| `acbn0_type` | `1` | `1` = on-site `U`; `2` = `U+V`, `(l,l')`; `3` = `U+V`, `(m,m')` |
| `ehub_l_choice(nt,1:6)` | `0` | which `l` gets corrected: `1=s, 2=p, 3=d, 4=f, 5=semicore d` |
| `ehub_nn_distance` | `3.0` Å | intersite cutoff → `ehub_nn_vect`, `ehub_nn_dim`, `ehub_nnnc` |
| `ehub_mixing` | `0.7` | linear mixing of `U`/`V` across SCF steps |
| `ehub_conv_thr` | `1e-8` | with `stop_acbn0`, freeze when `max\|ΔV\| < sqrt(1e-8)` = `1e-4` Ry |
| `stop_acbn0` | `.false.` | enable that freeze |
| `ehub_pot_fix` | `.false.` | keep `U`/`V` fixed; recompute occupations only |
| `remove_ehub_u(nt,l)` | `1` | set `0` to kill on-site `U` on a channel |
| `remove_ehub_v(nt,l)` | `1` | **inert — never read** (Finding #4) |
| `ehub_u`, `ehub_v` | `0` | starting / fixed values (Ry) |
| `read_ehub_ns_file`, `read_ehub_uv_file` | `.false.` | restart occupations / parameters from file |

Hard constraints enforced at startup (`ldaV.f90:102-131`): `nspin == 1` is rejected; fully
relativistic PPs require `lspinorb = .true.`; `ehub_nn_dim > 10` is rejected; `Hubbard_projectors` is
forced to `'ortho-atomic'`; `gamma_only` is rejected.

---

## 7. Code review — findings

Ordered by how much they should worry you.

### Finding #1 — out-of-bounds write in **all six** `symmetrize_*` routines ★ real bug

`new_ns_v.f90:2995, 3203, 3435, 5237, 5641, 6037`

```fortran
DO nt = 1, ntyp
   temp = ehub_l(1,nt)
   l_separation(ehub_l(1,nt),nt) = 1      ! <-- index 0 when the type has no Hubbard channel
   DO m = 2, ehub_lm(nt)
      ...
```

`l_separation` is declared `INTEGER :: l_separation(5,ntyp)` (lower bound 1). For any type with
`ehub_lm(nt) == 0` — perfectly legal, and guarded against by `IF (ehub_lm(nt) == 0) CYCLE` everywhere
else in the file — `ehub_l(1,nt)` is `0`, so this writes `l_separation(0,nt)`, i.e. **8 bytes outside
the array** when `nt == 1`, or into the unused `l_separation(5,nt-1)` slot otherwise.

Silent on most compilers; instant abort under `-fcheck=bounds` / `-check bounds`. Triggers whenever
you correct only the transition metal and leave, say, H or a spectator species uncorrected — and it
is worst exactly when that species is `ityp == 1`.

**Fix (one line, six places):**
```fortran
DO nt = 1, ntyp
   IF ( ehub_lm(nt) == 0 ) CYCLE          ! <-- add
   temp = ehub_l(1,nt)
   l_separation(ehub_l(1,nt),nt) = 1
```
Also initialise `l_separation = 0` after declaration; it is currently read-before-write for any
`l` value not present on a type (harmless today only because every read is guarded).

### Finding #2 — the ERI table is recomputed every SCF iteration ★ the big performance win

The `m1,m2,m3,m4` quadruple loops in `new_ns_a0`, `new_ns_e0_1`, `new_ns_e0_2`, `new_ns_e0_nc_1` call

```fortran
CALL pao_3g(upf(nt)%psd, ehub_l(m,nt), c)          ! ~100-branch SELECT CASE
... compute_eri(l3,l1,l4,l2, τ,τ,τ,τ, c3,c1,c4,c2)
```

**every time the routine runs, i.e. once per SCF iteration.** But every argument —
`(species, l, m)` of the four orbitals, the four positions, the Gaussian coefficients — is **fixed
for the whole SCF cycle**. Only the `P̄` prefactors change.

Cost, `d`-only, per atom pair per iteration: `5⁴ = 625` `compute_eri` calls; each expands into up to
`3⁴ = 81` primitive quartets × up to `3⁴ = 81` Cartesian-monomial combinations (when all four are
`d_{z²}`) ⇒ up to `~6.5 × 10³` `coulomb_repulsion` calls per ERI, `~4 × 10⁶` per atom pair. Multiply
by `nat × (neighbours + 1)`.

On top of that the assembly loops are **fully serial and replicated on every MPI rank** — no
`inter_pool`/`intra_bgrp` decomposition, no OpenMP.

**Fix:** hoist the ERIs into a table computed once per geometry (once per `scf`, once per ionic step
for relax/md), keyed on `(na1, na2, nnc, l1, l2, l3, l4)`. Also hoist `pao_3g` out of the `m4` loop
(the coefficients depend only on `(species, l)`). Expect a large constant-factor win with zero
physics change. Parallelising the outer `na1` loop over `inter_pool_comm` with an `mp_sum` is the
obvious second step.

### Finding #3 — `rho%ne0` scales as `nat²` and is in the mixing history ★ scalability wall

```fortran
rho%ne0(ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc)
```

**Dense over all atom pairs and all cells**, even though `ehub_nn_vect` is sparse (each atom has
O(10) neighbours regardless of cell size). Concretely, `ehub_lm_max = 5`, `nspin = 2`,
`ehub_nnnc = 27`:

| `nat` | one `ne0` array |
|---:|---:|
| 4 | 173 kB |
| 20 | 4.3 MB |
| 100 | 108 MB |
| 200 | 432 MB |

and there are several live copies (`rho`, `rhoin`, `v`) **plus the Broyden history**
(`mixing_ndim`, default 8) — so multiply by ~12. At `nat = 100` that is **~1.3 GB per MPI rank**,
duplicated on every rank.

Official QE's DFT+U+V avoids this entirely: `ldaU.f90:444` uses
`nsg(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)` — a *neighbour-list* index instead of a
dense `(nat, nat, ncell)` index. Same information, `nat × 20` instead of `nat² × 27`.

**This is the single most important thing to fix before porting to `hp.x`** — see §8.

### Finding #4 — `remove_ehub_v` is dead ★ user-visible

Declared (`input_parameters.f90:466`), defaulted (`read_namelists.f90:249`), broadcast
(`read_namelists.f90:1052`), copied into the module (`input.f90:2083-2084`) — and **never read by any
computation**. Users will set it and get no effect and no warning.

Relatedly, `remove_ehub_u` is applied **only to the on-site blocks** (`new_ns_v.f90:584, 1362, 2531, 4744`), never to intersite `V`. So there is currently *no* way to disable `V` on a given channel.
Either wire `remove_ehub_v` into the `ehub_v_(na1,na2,...)` assignments or delete it and error on use.

### Finding #5 — `stop_acbn0` silently does nothing for `acbn0_type = 1` and `3`

`ldaV.f90:455-462` allocates `ehub_u_diff` when `acbn0_type == 1`, but `new_ns_a0` never writes or
reads it, and has no early-exit block. `new_ns_e0_2` (`type = 3`) has no convergence machinery at
all. Only `new_ns_e0_1` and `new_ns_e0_nc_1` honour the flag. Given Finding #2, this is the
difference between a 3-hour and a 30-hour job. Either implement it for the other three or
`errore('...','stop_acbn0 requires acbn0_type=2',1)`.

### Finding #6 — `orberi` is a stateful iterator masquerading as a subroutine

`INTENT(INOUT) :: eri_l` with a `+1`-mod-cycle body. Correctness depends on four invariants that are
nowhere stated:
1. `eri_l` is reset to `0` immediately before every `m` loop;
2. the loop visits channels in increasing `m`;
3. `l` blocks are contiguous in `m` (true by construction in `init_lda_plus_v`);
4. the `m3`/`m4` loops apply their `IF (ehub_l(m3,nt) /= ehub_l(m1,nt)) CYCLE` **before** calling
   `orberi` (they do — but the ordering is load-bearing and unmarked).

Any refactor that reorders those loops corrupts the ERI labels **silently** and produces plausible,
wrong `U` values. Replace with a pure function of `(l, m_local)`:
```fortran
eri_lm = l*l + m_local          ! l=0→1 ; l=1→2,3,4 ; l=2→5,6,7,8,9
```
That is exactly the mapping the current code produces, in one line, stateless.

**Also:** `f` orbitals (`ehub_l == 4`) fall through every branch and leave `eri_l` untouched, and
`exponents_and_normalization_constants` has no `CASE` above 9. `ehub_l_choice(nt,4)` is accepted by
the input parser. **Add an explicit `errore` for `l = 3`** instead of returning garbage.

### Finding #7 — `hubbard_get_nnc` can return an uninitialised `nnc`

`hubbard_nns.f90:241-282`: if the requested cell hosts no neighbour pair, the `IF/ELSE` never assigns
`nnc`, and the caller gets whatever was on the stack.

Currently harmless — `symmetrize_n5`'s call at line 3595 computes `nnc0` *before* the `na2` loop and
never dereferences it when the `ehub_nn_vect` guard fires — but it is one refactor away from an
out-of-bounds `n(nnc0)%T`. Initialise `nnc = 0` and `errore` if it is still `0` on exit.

### Finding #8 — the noncollinear Hermiticity check is commented out

`new_ns_v.f90:4044-4076` — a full `|ns(m1,m2,j,I,J) − conj(ns(m2,m1,k,J,I))| > 1e-10`
test, entirely disabled, with `!! temp` markers. For a noncollinear intersite implementation this is
exactly the assertion you want live behind `iverbosity > 0`. Restore it (guarded) or delete it.

### Finding #9 — `p(i)%O` is never explicitly deallocated

`new_ns_v.f90:738` (`e0_1`), `:1902` (`e0_2`), `:4198` (`e0_nc_1`) allocate `p(i)%O` and never `DEALLOCATE`
it, while every sibling array (`n`, `dm0`, `dm1`, `E0`, `HF`, `ehub_v_`, `num`) is explicitly freed.

**Not a leak** — `p` is a local automatic array of a derived type, and F2003 auto-deallocates
allocatable components when the variable ceases to exist. But it is the only asymmetry in an
otherwise scrupulous file, which makes a reader stop and check. Add the three `DEALLOCATE`s.

### Finding #10 — ~3800 lines of copy-pasted `IF/ELSE IF` ladders

The four `symmetrize_*` routines contain 16 near-identical branches each, selecting `d1`/`d2`/`d3`
by `(ehub_l(m1), ehub_l(m2))`. This collapses to:

```fortran
FUNCTION dmat(l, i, j, isym) RESULT(d)
   ! l = 1(s), 2(p), 3(d), 4(f)  -> 1.0, d1, d2, d3
   SELECT CASE (l)
      CASE (1) ; d = 1.0_DP
      CASE (2) ; d = d1(i,j,isym)
      CASE (3) ; d = d2(i,j,isym)
      CASE (4) ; d = d3(i,j,isym)
      CASE DEFAULT ; CALL errore('symmetrize','l not implemented',1)
   END SELECT
END FUNCTION

n%T(m1,m2,is,na1) = n%T(m1,m2,is,na1) &
   + dmat(ehub_l(m1,nt1), m3_, m1_, isym) * n%O(m3,m4,is,na2) &
   * dmat(ehub_l(m2,nt1), m4_, m2_, isym) / DBLE(nsym)
```

Same result, ~40 lines instead of ~3800, and `f` support becomes a one-line table entry. Worth doing
before the `hp.x` port — you do not want to carry this into a second codebase.

### Finding #11 — smaller items

| Where | Issue |
|---|---|
| `new_ns_v.f90:214`, `v_hubbard_v.f90:126` | `nspin == 1` branches are dead (`ldaV.f90:104` aborts) |
| `new_ns_v.f90:77, 3891` | `Hubbard_projectors == 'pseudo'` check is dead (`orthoatwfc_v.f90:33` hard-sets `'ortho-atomic'`) |
| `new_ns_v.f90:96-113`, `:172-190` | `gamma_only` branches are dead (`orthoatwfc_v.f90:40-41` aborts) |
| `new_ns_v.f90:1800` | `GOTO 1` / numeric label — replace with a logical guard |
| `new_ns_a0:497`, `new_ns_e0_1:1276` | the `l == 0` special case (no `J̄`, opposite-spin `Ū` only) is correct but **completely uncommented** |
| `new_ns_v.f90:1550` | the `/2.0_DP` in `V = (E0/HF)/2` is a definitional double-counting factor with no comment |
| `new_ns_v.f90:938-948` | the `N̄^{IJ}` sequential in-place update is order-dependent; add a warning or use a scratch array |
| `new_ns_v.f90` throughout | the `IF (MAXVAL(ehub_nn_vect(...)) == 1) nnc = nnc + 1` cell-counting idiom is duplicated **30× across 7 files** (17× in `new_ns_v.f90` alone; also `vhpsi_v.f90`, `v_hubbard_v.f90`, `force_hub_v.f90`, `stres_hub_v.f90`, `ldaV.f90`, `io_rho_xml.f90`, plus the independent re-derivation in `input.f90:2148-2178`). Every copy must agree on the cell ordering or `U`/`V` silently attach to the wrong bond. One `nnc_of_cell(n1,n2,n3)` lookup table built once in `init_lda_plus_v` would remove a whole class of future bugs — **do this before the port** |
| `new_ns_e0_2` | inherits `(l,l')`-resolved on-site `U` while advertising `(m,m')`-resolved `V` — document or unify |

### What is *right* — and it is a lot

* Eqs. (7)–(10) and (B3)–(B6) are implemented **correctly**, including the sign structure and the
  `Σ_{m≠m'}` vs `Σ_{{m}}` distinction that is easy to get wrong.
* The `{I}` sum of Eq. (8) is correctly implemented as *same chemical element* via `upf%psd`, so
  spin-split types (`Ni1`/`Ni2` in AFM NiO) merge properly.
* Time-reversal / point-group symmetrisation is applied to `n`, `P̄`, `dm0`, `dm1` consistently.
* `V^{IJ}(R) = V^{JI}(−R)` and hence Hermiticity of the Hubbard operator hold **by construction**.
* Parallel reductions are correct: `calbec` handles `intra_bgrp`, `mp_sum(inter_pool_comm)` handles
  k-points, and the routines run before `divide(inter_bgrp_comm, ...)` so band groups do not double
  count.
* The `stop_acbn0` early exit is a genuinely good design for an expensive functional.
* The `dm1(nt,l)` construction is a thoughtful, non-obvious answer to a question the paper never
  poses.

---

## 8. What this means for your `hp.x` port

Your stated goal is moving DFT+U+V from this QE fork into `hp.x`. Three concrete observations:

**1. `hp.x` already speaks a different dialect.** QE 7.3's HP works on
`nsg(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)` with `neighood(na)%neigh(viz)` indices into
a supercell atom list `at_sc` (`ldaU.f90:155-196`; used e.g. at `HP/src/hp_ns_trace.f90:84`). This
fork works on `ne0(ehub_lm_max, ehub_lm_max, nspin, nat, nat, ehub_nnnc)` with
`ehub_nn_vect(a,b,c,na1,na2)`. **The two never meet.** The port is, first and foremost, an index
translation:

```
ne0(m1, m2, is, na1, na2, nnc)   ⟷   nsg(m1, m2, viz, na1, is)
      where  viz = position of (na2, R(nnc)) in neighood(na1)%neigh
```

Writing that mapping as a small standalone module — and validating it by round-tripping a converged
`ne0` through `nsg` and back — is the highest-value first step. It also happens to fix Finding #3 for
free.

**2. The physics is orthogonal to HP's.** `hp.x` computes `U`/`V` by **linear response** (the
`χ₀ → χ` route of Cococcioni–de Gironcoli); this file computes them by a **pseudo-hybrid ratio**
(`E0/HF`). They are two independent answers to the same question. Nothing in `new_ns_v.f90` needs to
move into `hp.x`'s response machinery — what `hp.x` needs from this fork is the *occupation matrix
infrastructure* (`new_ns_v` proper, `symmetrize_n5`, `hubbard_symmetry`, `vhpsi_v`, `force_hub_v`),
not the ACBN0 kernels. **`new_ns_v` (lines 10–230) is the routine you port. `new_ns_e0_*` is the
routine you leave behind** — unless you specifically want an ACBN0-vs-linear-response comparison,
which would in fact make a good paper.

**3. Before you build on it, do these four things:**
   * fix Finding #1 (six one-line guards) — it will bite you the moment you run a cell with a
     spectator species;
   * fix Finding #2 (cache the ERIs) — it decides whether your benchmarks are feasible;
   * fix Finding #3 (sparse neighbour indexing) — it decides how large a cell you can publish;
   * decide, and state in writing, whether you keep the arithmetic-mean `N̄^{IJ}` or switch to the
     geometric mean at `new_ns_v.f90:941`. **A referee who knows the Tancogne-Dejean & Rubio `+V`
     paper will ask.**

---

## Appendix A — equation → line index

| Equation | Collinear | Noncollinear |
|---|---|---|
| (3) `n^{Iσ}_{mm'}` | `new_ns_v.f90:83-122` (on-site + in-cell), `:145-203` (`R≠0`) | `:3903-3943`, `:3986-4020` |
| (A2) `n^{Iσσ'}_{mm'}` | — | `:3917-3934` |
| (8)/(B4) `N̄` | `:343-395` (`a0`), `:879-950` (`e0_1`) | `:4309-4400` |
| (7)/(B3) `P̄` | `:399-441` (`a0`); `:967-1006` `dm0`, `:1008-1049` `dm1`, `:1052-1191` `p` | `:4404-4474` `dm0`/`dm1`, `:4475-4596` `p` |
| (9) `Ū` numerator | `:497-515` (`a0`), `:1277-1295` (`e0_1`) | `:4653-4670` |
| (10) `J̄` numerator | `:516-524` (`a0`), `:1297-1304` (`e0_1`) | `:4673-4681` |
| (9)/(B5) denominator | `:536-563` (`a0`), `:1314-1341` (`e0_1`) | `:4701-4715` |
| (10)/(B6) denominator | `:559-562` (`a0`), `:1336-1341` (`e0_1`) | `:4717-4720` |
| `U^eff = Ū − J̄` | `:570-598` (`a0`), `:1348-1384` (`e0_1`) | `:4730-4765` |
| (6) generalised → `V^{IJ}`, `R = 0` | `:1427-1446` | `:4809-4853` |
| (6) generalised → `V^{IJ}`, `R ≠ 0` | `:1618-1635` | `:5006-5045` |
| `V^{IJ}` denominator | `:1494-1505` (`R=0`), `:1691-1697` (`R≠0`) | `:4880-4894`, `:5083-5097` |
| (B2) ERI | `compute_eri.f90:31-179` | same |
| Symmetrisation | `:2964`, `:3168`, `:3400` | `:5202`, `:5602`, `:5998` |
| (2) energy, (5) potential | `v_hubbard_v.f90:1-137` | `v_hubbard_v.f90:200-367` |
| (5) applied to `\|ψ⟩` | `vhpsi_v.f90:1-185` | `vhpsi_v.f90:~300-440` |
| (11)–(13) forces | `force_hub_v.f90` | same |

## Appendix B — reproducing this review

```bash
# paper text (no poppler on this box; ghostscript works)
gs -q -dNOPAUSE -dBATCH -sDEVICE=txtwrite -dTextFormat=3 \
   -sOutputFile=paper.txt logs/PhysRevB.96.245133.pdf
# equations are mangled by the two-column layout -> render the pages instead
gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -r150 -dFirstPage=2 -dLastPage=3 \
   -sOutputFile=pg%d.png logs/PhysRevB.96.245133.pdf     # Eqs. (1)-(10)
gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -r150 -dFirstPage=7 -dLastPage=8 \
   -sOutputFile=pgB%d.png logs/PhysRevB.96.245133.pdf    # Appendices A, B

# the two eACBN0 variants differ only in V resolution -- confirm with:
sed -n '636,1806p'  PW/src/new_ns_v.f90 > /tmp/e1.f90
sed -n '1809,2961p' PW/src/new_ns_v.f90 > /tmp/e2.f90
diff -u /tmp/e1.f90 /tmp/e2.f90
```

---

*Files consulted beyond `new_ns_v.f90`:* `ldaV.f90`, `hubbard_nns.f90`, `compute_eri.f90`,
`pao_3g.f90`, `orthoatwfc_v.f90`, `v_hubbard_v.f90`, `vhpsi_v.f90`, `sum_band.f90`, `electrons.f90`,
`v_of_rho.f90`, `scf_mod.f90`, `init_ns_v.f90`, `input.f90`, `Modules/input_parameters.f90`,
`Modules/read_namelists.f90`, `ldaU.f90`.
