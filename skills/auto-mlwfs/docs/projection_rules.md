# Projection Rules

This step selects Wannier projection orbitals from a completed VASP SCF calculation.

## Data Source

Read `vasprun.xml` with `pymatgen.io.vasp.outputs.Vasprun(parse_potcar_file=False)`.
Use `complete_dos.get_element_spd_dos(element)` so every element and every available `s/p/d/f` orbital channel is considered automatically.

## Metal vs Semiconductor

Use `Vasprun.eigenvalue_band_properties`:

- Treat the system as metallic when `gap <= metal_gap_tolerance_ev`.
- Otherwise treat it as a semiconductor or insulator.

## Energy Manifold

Use the total DOS to find all contiguous nonzero-DOS segments. Expand each
segment on both sides by `gap/2 + analysis.manifold_padding_ev`, merge
overlapping segments, then select the merged manifold containing `E_F`. If
`E_F` lies in a remaining gap, join the nearest occupied and unoccupied
manifolds around it.

Use `analysis.zero_density_tolerance` to define numerical zero. The default
`0.0` reproduces the reference autoconstruction rule. Fixed Fermi-centered and
band-edge windows are not supported because one energy width does not transfer
consistently between materials.

## Selection

For each `Element:orbital` channel:

1. Integrate the element-resolved PDOS inside the selected DOS manifold.
2. Divide by the number of atoms of that element, matching the reference
   per-atom orbital decision.
3. Keep the channel when this absolute per-atom weight reaches its configured
   threshold. Defaults are `s=0.4`, `p=1.2`, `d=2.0`, and `f=2.8`.

Do not normalize by the total weight of all channels and do not impose a fixed
channel count. If no channel passes, stop and inspect the DOS or thresholds.

The output format for downstream Wannierization is one projection per line:

```text
Fe:d
Mo:d
O:p
```
