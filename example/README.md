# example/ -- input and output files of dev hp.x (acbn0 = 1 and 2)

Two materials, both run on 2026-10-30 with the pinned binaries
`bin/pw.x.e2801ae` (PW/src identical to official QE 7.3) and
`bin/hp.x.a1ad49f` (dev hp.x, commit a1ad49f), 1 node x 32 MPI ranks,
`-npool 32`.  See `run.sh` in each directory for the exact job.

```
example/
  pseudo/      Ga.upf As.upf Ni.upf   PseudoDojo NC SR PBE standard (md5 in pseudo/README.md)
  GaAs/        zinc-blende GaAs, s+p manifold on both atoms, s on-site U removed
  Ni/          fcc Ni, 3d manifold, one atom per cell
```

Each material directory holds the same nine files:

| file               | what it is                                                                  |
|--------------------|-----------------------------------------------------------------------------|
| `scf.in/.out`      | plain PBE SCF with dev (= official) pw.x. `nspin = 2`, NO HUBBARD card      |
| `hp_U.in/.out`     | dev hp.x, `acbn0 = 1`: one-shot ACBN0 on-site U                             |
| `eACBN0_hub_U.dat` | the file hp.x wrote as `out/eACBN0_hub.dat` in that run (renamed by run.sh) |
| `hp_UV.in/.out`    | dev hp.x, `acbn0 = 2`, `nn_distance = 3.0`: one-shot eACBN0 U + V           |
| `eACBN0_hub_UV.dat`| same, for the U+V run                                                       |
| `run.sh`           | the sbatch script (module, MKL flags, srun lines, renames)                  |

Not shipped: `out/` (the `.save`, the 32 `*.wfc` files, and the forked-format
dumps hp.x writes into `out/<prefix>.save/`).  Re-running `run.sh` regenerates
all of it in well under a minute of compute (GaAs pw.x 20 s, hp.x < 1 s each;
Ni hp.x U+V 9 s).  `run.sh` refuses to overwrite `scf.out`, `hp_*.out` or
`out/`: copy the material directory next to `pseudo/` somewhere under `/gpfs`,
delete those files, and `sbatch run.sh` from inside it.

## Reference numbers (eV unless stated)

GaAs, a = 10.643882 bohr = 5.6325 A (CELL_PARAMETERS components a/2 = 5.321941 bohr), 15x15x15 k, ecutwfc 100 Ry,
Gaussian smearing 0.01 Ry.  Manifold: Ga 4s + 4p, As 4s + 4p; on-site U of
both 4s channels masked by `acbn0_remove_ehub_u(nt,1) = 0`.

| quantity                   | acbn0 = 1 | acbn0 = 2 (d_IJ = 3.0 A) |
|----------------------------|-----------|--------------------------|
| U(Ga 4p)                   | 0.3418    | 0.3418                   |
| U(As 4p)                   | 2.0105    | 2.0105                   |
| V(Ga 4s - As 4s), d = 2.4389 A | -     | 0.8682                   |
| V(Ga 4s - As 4p)           | -         | 1.3418                   |
| V(Ga 4p - As 4s)           | -         | 0.8281                   |
| V(Ga 4p - As 4p)           | -         | 1.6677                   |
| one-shot E(+U) / E(+U+V), Ry | 0.122578611 | -0.041595649           |
| neighbour cells nnnc       | -         | 7 (home + 6)             |

Ni, celldm(1) = 6.48 bohr, 15x15x15 k, ecutwfc 100 Ry, M-V smearing 0.02 Ry,
`starting_magnetization = 0.0` (non-magnetic start).  Manifold: Ni 3d.

| quantity                   | acbn0 = 1 | acbn0 = 2 (d_IJ = 3.0 A) |
|----------------------------|-----------|--------------------------|
| U(Ni 3d)                   | 16.4018   | 16.4018                  |
| V(3d - 3d), d = 2.4247 A, 12 images | - | 2.5890                  |
| one-shot E(+U) / E(+U+V), Ry | 0.616163765 | 0.600010616            |
| neighbour cells nnnc       | -         | 13 (home + 12)           |

Two things these tables show that the main README explains:

- The on-site U of the `acbn0 = 2` route equals the `acbn0 = 1` U for the same
  ground state and manifold.  It does not depend on `nn_distance`.
- These are ONE-SHOT values, U_0 = F[rho_PBE].  They are not the self-consistent
  U and V that Lee and Son's pw.x converges to.  Stock runs with the same
  structure, cutoff, k-mesh, smearing and pseudopotentials, with the manifold
  declared through stock's `&system` flags, are in
  `/gpfs/yjkang/project-01/3-hub_dat_files`; the main README section 9 puts
  their numbers next to these.  Dev hp.x run on the stock converged `.save`
  reproduces the stock raw (unmixed) last-iteration U and V to the last printed
  digit of the dumps, at most 1e-15 Ry (`2-print_wo_mixing`).
- The Ni start is non-magnetic on purpose, copying Lee and Son's
  `example1/Ni.scf.in`, so the two spin blocks of the occupations are equal.

The Ga and As 4s rows print U = 0.0000 because of the mask, not because the
estimator gave zero.  Remove the two `acbn0_remove_ehub_u` lines to see the
unmasked 4s values.
