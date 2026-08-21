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

| quantity | value |
|---|---|
| objective NA | 1.45, full aperture, NOT capped |
| immersion / glass / sample index | 1.52 / 1.518 / 1.33 |
| sample geometry | `airOnGlass` (the stratified model, general in `nSample`) |
| excitation / emission | 0.640 / 0.690 um |
| detector | honeycomb23, 0.18 um pitch |

An earlier version capped the in-sample NA at `0.98 * n_sample = 1.303`,
justified by the claim that supercritical rays are evanescent and so do not
contribute beyond ~100 nm. That was wrong, and the correction matters.

Supercritical-angle fluorescence (SAF) IS collected. A dipole close to the
interface has near-field content with in-plane momentum above `n_sample*k0`,
and at the water-glass boundary that couples into propagating waves in the
higher-index glass - the reciprocal of TIRF. Because SAF decays over roughly
100-200 nm from the surface, its share of the signal is the steepest axial
reporter available in exactly the range the membrane occupies, and it is
nearly absent for internalised dye several hundred nm up. That is closer to a
binary membrane-versus-internalised discriminator than defocus is. Capping the
aperture discarded the best channel.

The stratified model handles it correctly. Where the homogeneous model needs
`sinTheta = (NA/nMedium)*rho` and therefore `NA <= nMedium`, the interface
model uses `q = NA*rho` with `cosSample = sqrt(1-(q/nSample)^2)`, so the pupil
spans the full aperture and components with `q > nSample` acquire a complex
`cosSample`. Its own header states that this retains evanescent decay above
the critical angle. `sampleGeometry` is named `airOnGlass` but the code is
parameterised by `nSample`, `nGlass`, `nImmersion` and `nDesignGlass`, so it
represents water-on-glass by setting `nSample = 1.33`.

The excitation is put through the same stratified model, because a
full-aperture NA 1.45 focus also produces an evanescent excitation
contribution near the surface. That is a second distance-dependent term and it
belongs in the product.

`spadEffectivePSFArray` in the vendored code is hard-wired to the homogeneous
`psfBessel`, so `spadEffectivePSFArrayInterface.m` (new, in this folder) is the
interface counterpart.

## Remaining limitation

**The metal film is still not modelled.** MIET works because the metal
reshapes the emission angular distribution, adding surface-plasmon-coupled
emission into a narrow high-angle cone with its own distance dependence.
Neither `homogeneous` nor `airOnGlass` captures a metal layer. Near-surface
emitters are therefore represented less faithfully than out-of-focus ones, so
the membrane ring signature is not quantitative even though the contrast
between the two pools remains informative.

The PSF is also symmetric about focus in the homogeneous limit, so far from
the surface this measures `|z|` rather than signed `z`. Near the surface the
interface breaks that symmetry, which is a further reason to keep the full
aperture.

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
