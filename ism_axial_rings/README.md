# ism_axial_rings

Can the SPAD detector rings separate membrane dye (in focus) from internalised
dye (above focus), independently of lifetime?

The lifetime route cannot answer this. Beyond roughly 200 nm from the metal the
MIET quenching curve flattens, so a long lifetime means "somewhere above
200 nm" and nothing more. The detector array carries an independent axial
signal over exactly that range: out-of-focus light spreads across more channels
than in-focus light, so the distribution across channels encodes height
regardless of lifetime.

## Contents

- `simulate_ring_weight_vs_height.m` - the ring weight W(ring, z), the photon
  budget needed for an axial call, and the precision on an out-of-focus
  fraction.
- `vendor/SPADarray_AberrationPSF/` - 399 m-files copied from
  `D:\MATLAB\server\Luminosa\Aberration_estimation\SPADarray_AberrationPSF`.
  Package folders (`+membrane_tracking/+curved_miet`, `+fluctuating_miet`,
  `+focused_ism`, `+hop_trap`) are preserved because those functions begin with
  `import membrane_tracking.<pkg>.*` and break if flattened. The vendored copy
  excludes `.git`, a duplicated `.claude/worktrees` tree, and the output/asset
  images.

The forward model comes from that codebase: `spadEffectivePSFArray` returns
`hEff(y,x,z,k)`, the effective excitation-times-detection PSF per detector
channel from a vectorial Richards-Wolf model, and the ring weight is its
lateral sum.

## Optics as configured

| quantity | value | note |
|---|---|---|
| nominal objective NA | 1.45 | oil immersion |
| immersion index | 1.52 | recorded, not propagated |
| sample index | 1.33 | water-like medium |
| **effective in-sample NA** | **1.303** | `0.98 x n_sample` |
| excitation / emission | 0.640 / 0.690 um | |
| detector | honeycomb23, 0.18 um pitch | same array as Luminosa |

The NA cap is not a convenience. `vectorialPSFBessel` forms
`sinTheta = (NA / nMedium) * rho`, which requires `NA <= nMedium`; NA 1.45 in
water would give `sinTheta = 1.09`. Physically, everything above the critical
angle at the glass-water interface is evanescent and only reaches emitters
within about 100 nm of the surface, so it contributes nothing at the 0.2-1 um
heights this study is about. Axial discrimination is therefore slightly weaker
than a naive NA 1.45 calculation would suggest.

## Three limitations to keep in view

1. **The metal film is not modelled.** MIET works because the metal reshapes
   the emission angular distribution, and none of `homogeneous` or
   `airOnGlass` captures that. Near-surface (membrane) emitters are therefore
   modelled less faithfully than out-of-focus (internalised) ones. Since the
   quantity of interest is the contrast between the two, this biases the
   comparison in a knowable direction rather than invalidating it - but the
   membrane ring signature should not be treated as quantitative.
2. **No water-on-glass stratified model exists** in the vendored code;
   `sampleGeometry` accepts only `airOnGlass` and `homogeneous`. Homogeneous
   water is used.
3. **The PSF is symmetric about focus**, so this measures `|z|`, not signed
   `z`. Acceptable when focus sits at the membrane and the dye of interest is
   above it; it cannot separate above from below.

## Status

`simulate_ring_weight_vs_height.m` has **not been executed**. MATLAB `-batch`
aborts on this machine with a graphics assertion whenever an interactive
session holds the display, so not even `checkcode` could be run against the
final version. Two real bugs were found and fixed by reading the forward model
rather than by running it:

- the sweep originally used `beadBottomHeightUm`, which only
  `airSurfaceBeadSlices` reads. The axial coordinate is `sim.z`, used directly
  by `vectorialPSFBessel`. The original would have returned an identical PSF at
  every height - a flat result that looked plausible.
- `defaultParams` builds `sim.x/y/z` at construction, so overriding `fovXY`
  and `nx` afterwards left the default grid in place. The grids are now
  rebuilt.

Run it before trusting anything it prints.
