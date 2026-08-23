# ism_axial_rings

Can the SPAD detector rings separate membrane dye (in the contact plane) from
internalised dye (above it), independently of lifetime?

**Yes, as a binary decision, from about 0.2 um upwards. No, as a height
measurement above about 0.3 um.** Those two statements are the whole result and
they are measured, not assumed.

The lifetime route alone cannot answer this. Beyond roughly 200 nm from the
metal the MIET quenching curve flattens, so a long lifetime means "somewhere
above 200 nm" and nothing more. The detector array carries an independent axial
signal over exactly that range: out-of-focus light spreads across more channels
than in-focus light, so the distribution across channels encodes height
regardless of lifetime.

## The headline numbers

### The rings decide something the summed decay provably cannot

`study_ring_height_discrimination` is built so the detector-summed decay carries
**exactly zero** information. Two truths share the same three lifetimes; the
third pool sits either at 0.40 um or at the membrane, and its molecule count is
set separately in each case so it contributes the same fraction of detected
photons. The summed decay depends only on lifetimes and detected photon
fractions, so the two truths have an identical expected summed decay - verified
numerically at a relative difference of 1e-16 to 7e-16 in every condition.

The argument is exact, not approximate: a sum of independent Poisson variables
is Poisson in the sum of their means, so the summed data have *identically
distributed* counts under both truths. Any statistic computed from the summed
decay is therefore exactly null, and the summed column below is a true control.

Detection rate at a matched 5% false-positive rate, third pool at 0.40 um,
60 repeats:

| third pool | N=1e3 | N=3e3 | N=1e4 | N=3e4 | summed control |
|---|---|---|---|---|---|
| 2% of photons | 0.10 | 0.28 | 0.45 | 0.88 | 0.02-0.08 |
| 5% of photons | 0.12 | 0.43 | 0.98 | 1.00 | 0.03-0.10 |
| 10% of photons | 0.67 | 0.92 | 1.00 | 1.00 | 0.02-0.15 |

The summed control stays at chance throughout, as it must. Across all 17
conditions it averages 0.054 against a nominal 0.05; the two values at 0.12-0.15
are the expected upper excursions of a 60-repeat binomial over that many
independent tests, not a leak of information.

The test statistic is the deviance gain from letting the components differ in
height (one shared fitted height versus one height per component). It needs no
knowledge of the true heights, so it runs unchanged on real data.

### Minimum detectable displacement

N = 1e4, third pool at 10% of photons, 60 repeats:

| displacement | ring power | summed | normalised pattern gap |
|---|---|---|---|
| 0.08 um | 0.35 | 0.03 | 0.0029 |
| 0.15 um | 0.52 | 0.02 | 0.0084 |
| 0.25 um | **1.00** | 0.05 | 0.0165 |
| 0.40 um | 1.00 | 0.03 | 0.0235 |
| 0.60 um | 1.00 | 0.12 | 0.0259 |

**The threshold sits between 0.15 and 0.25 um.** That is the number that decides
whether this is worth doing, and it is the good case: MIET lifetime is
informative from the surface to roughly 0.2 um, and the rings become reliable
from roughly 0.2 um upwards. The two methods are complementary with little or no
gap between them, which is not something that had to be true.

## What the fitted height is worth

`diagnose_ring_height_identifiability` settles this, with the decisive tests
using no optimiser at all.

An earlier version of `study_ring_resolved_unmixing` reported a recovered
third-pool height of 0.488 um at every photon count from 1e3 to 3e4, against a
truth of 0.40 um. A bias that stable across a 30-fold change in photon count is
the signature of mis-specification, not of statistical error, and the default
height seed happened to be 0.5 um. Three tests separate the possibilities:

1. **The forward model is correct.** On noiseless data the deviance profile over
   the third height has its minimum at exactly 0.400 um, with a residual of
   3e-11. So the stratified PSF, the weight table and the design matrix are all
   right. The problem is information, not implementation.

2. **The normalised ring pattern saturates above about 0.3 um.** Distance from
   the pattern at 0.40 um is 0.039 at z = 0.30 but only **0.037 at z = 1.00** -
   0.4 um and 1.0 um are less distinguishable from each other than 0.4 is from
   0.3. Absolute collection efficiency does keep falling steeply (summed weight
   2.12 at the surface, 0.31 at 0.40 um, 0.05 at 1.0 um), but that is degenerate
   with the freely fitted amplitude, so only the pattern *shape* identifies
   height.

3. **So height has an interquartile range of about 0.6 um at N = 1e4**, on a
   one-dimensional profile with the lifetimes and other heights pinned at truth
   and no fitting involved, and it does not improve at 3e4.

The 0.488 um was the objective being flat, not the parameter being determined.

The displacement sweep refines this usefully, though. The fitted height spread
tracks the truth closely **below** saturation - 0.089 um recovered for a true
0.08, 0.148 for a true 0.15 - and only overshoots above it (0.480 for 0.40,
0.516 for 0.60). So height is quantitative in the 0-0.2 um band and merely
qualitative above 0.3 um, which is the same physics that sets the detection
threshold.

## Which model to use, and which not to

`study_ring_resolved_unmixing` compares three ways of detecting a third pool at
matched false-positive rate and matched total photons:

| model | ring dimension | ring physics |
|---|---|---|
| `summed` | collapsed | none - what the current pipeline does |
| `ring-free` | used as data | none: `nComp * nRing` free amplitudes |
| `ring-w(z)` | used as data | imposed: pattern forced to `a_j * w(:, z_j)` |

The middle model is there to isolate the two effects, and it earns its place by
failing. Adding the ring dimension **without** the physics is actively harmful:
its power reaches only 0.43 at N = 3e4 where the summed decay reaches 0.95, and
its lifetime error does not improve with photon count at all (0.247 to 0.302 ns
from 1e3 to 3e4). With `nRing` free amplitudes per component the model can
absorb a third pool's signature without ever committing to a third lifetime, so
its null deviance-gain distribution is inflated (median 4.8-6.2 against ~0 for
the other two) and its threshold rises accordingly.

The constrained model is the one that helps. Its clearest gain is on lifetime
accuracy - roughly a factor of 2.6 to 3 better than the summed decay at the same
photons - which is a real information gain even though the height itself is not
identified: constraining the ring pattern to a one-parameter physical family
reduces the effective parameter count and tightens the lifetimes.

Free amplitudes are the wrong choice, and the reason generalises: extra
parameters that the artefact can also exploit do not discriminate. Only a
constraint the artefact cannot satisfy does.

## Contents

Forward model and diagnostics:

- `simulate_ring_weight_vs_height.m` - the ring weight W(ring, z), the photon
  budget needed for an axial call, and the precision on an out-of-focus
  fraction.
- `spadEffectivePSFArrayInterface.m` - per-channel effective PSF across the
  interface. Needed because the vendored `spadEffectivePSFArray` is hard-wired
  to the homogeneous `psfBessel`, which cannot represent NA 1.45 into water.
- `diagnose_ring_height_identifiability.m` - is the fitted height real?

Studies:

- `study_ring_height_discrimination.m` - displaced versus coplanar third pool,
  with the summed decay as an exact null control. **This is the one that
  answers the experimental question.**
- `run_ring_discrimination_sweeps.m` - the photon/abundance and displacement
  sweeps above, sharing one forward-model evaluation.
- `study_ring_resolved_unmixing.m` - the three-model comparison.

Shared fitting code, factored out so the diagnostics exercise the same code the
studies measure rather than a copy that can drift:

- `ring_fit_free_amplitudes.m`, `ring_fit_height_constrained.m` - the two
  fitters, by variable projection with an inner non-negative least squares.
  `ring_fit_height_constrained` accepts `heightGroups`, which is what makes the
  shared-height null available.
- `ring_constrained_deviance.m` - the objective in physical units, so it can be
  profiled directly. This is what makes the no-optimiser tests possible.
- `ring_build_constrained_design.m`, `ring_build_patterns.m`,
  `ring_periodic_decay.m`, `ring_interpolate_weights.m`,
  `ring_poisson_deviance.m`, `ring_poisson_sample.m`, `ring_fit_defaults.m`,
  `ring_seed_tau.m`, `ring_search_options.m`.

Vendored:

- `vendor/SPADarray_AberrationPSF/` - 399 m-files copied from
  `D:\MATLAB\server\Luminosa\Aberration_estimation\SPADarray_AberrationPSF`.
  Package folders (`+membrane_tracking/+curved_miet`, `+fluctuating_miet`,
  `+focused_ism`, `+hop_trap`) are preserved because those functions begin with
  `import membrane_tracking.<pkg>.*` and break if flattened. The copy excludes
  `.git`, a duplicated `.claude/worktrees` tree, and the output images.

See `SAF_AUDIT.md` for whether supercritical-angle fluorescence is handled
correctly across the PSF, aberration and tracking code.

## Optics as configured

| quantity | value |
|---|---|
| objective NA | 1.45, full aperture, NOT capped |
| immersion / glass / sample index | 1.52 / 1.518 / 1.33 |
| sample geometry | `airOnGlass` (the stratified model, general in `nSample`) |
| excitation / emission | 0.640 / 0.690 um |
| detector | honeycomb23, 0.18 um pitch |
| rings | 5, edges at 0.09 / 0.246 / 0.336 / 0.418 um, members 1 / 6 / 6 / 6 / 4 |

An earlier version capped the in-sample NA at `0.98 * n_sample = 1.303`,
justified by the claim that supercritical rays are evanescent and so do not
contribute beyond ~100 nm. That was wrong, and the correction matters.

Supercritical-angle fluorescence (SAF) IS collected. A dipole close to the
interface has near-field content with in-plane momentum above `n_sample*k0`, and
at the water-glass boundary that couples into propagating waves in the
higher-index glass - the reciprocal of TIRF. Because SAF decays over roughly
100-200 nm from the surface, its share of the signal is the steepest axial
reporter available in exactly the range the membrane occupies, and it is nearly
absent for internalised dye several hundred nm up. Capping the aperture
discarded the best channel.

The stratified model handles it correctly. Where the homogeneous model needs
`sinTheta = (NA/nMedium)*rho` and therefore `NA <= nMedium`, the interface model
uses `q = NA*rho` with `cosSample = sqrt(1-(q/nSample)^2)`, so the pupil spans
the full aperture and components with `q > nSample` acquire a complex
`cosSample`. `sampleGeometry` is named `airOnGlass` but the code is
parameterised by `nSample`, `nGlass`, `nImmersion` and `nDesignGlass`, so
water-on-glass is just `nSample = 1.33`.

The excitation is put through the same stratified model, because a
full-aperture NA 1.45 focus also produces an evanescent excitation contribution
near the surface. That is a second distance-dependent term and it belongs in the
product.

## Limitations that affect the numbers above

**The metal film is not modelled.** MIET works because the metal reshapes the
emission angular distribution, adding surface-plasmon-coupled emission into a
narrow high-angle cone with its own distance dependence. Neither `homogeneous`
nor `airOnGlass` captures a metal layer. Near-surface emitters are therefore
represented less faithfully than out-of-focus ones, so the ring signature at
z ~ 0 is not quantitative. The *contrast* between a surface pool and a displaced
pool is the robust part; the exact 0.15-0.25 um threshold is the part most
exposed to this.

**Single-height pools.** Each simulated pool sits at one exact height, while
real membrane and internalised dye occupy a continuum. The constrained model
assumes exactly the single-height form `a_j * w(:, z_j)`, so the simulation is
structurally generous to it. A continuum would broaden the pattern and reduce
power by an amount not measured here.

**No background.** The truths contain no dark counts or afterpulsing, while the
production pipeline models a background term.

**There is no data path yet.** The pipeline saves only detector-summed TCSPC, so
none of this can be run on the real acquisitions as they stand. The hook is
small and localised: `shiftIndex` is already computed per photon in
`immune_cell_MIET_reassigned_sliding_tcspc.m`, so a ring lookup plus one extra
stride term in the destination index gives per-ring TCSPC in both the MATLAB and
the MEX branch, at roughly 214 MB for a 137k-window 4x4 stage. Ring membership
should be taken from the DATA rather than from the simulated layout: the
measured per-channel `shiftsToCenter` is the registered lateral displacement of
each channel, so `norm(shiftsToCenter(k,:))` orders the channels radially
without assuming the vendored `honeycomb23` ordering matches the Luminosa
`im_chan` values. Adding it requires re-reading every PTU file.

## Status

All results above were produced by execution, not by reading the source. Two
bugs were found and fixed by reading the forward model before it ever ran (the
sweep originally varied `beadBottomHeightUm`, which the homogeneous PSF never
reads, and `defaultParams` builds `sim.x/y/z` at construction so overriding
`fovXY`/`nx` afterwards had no effect), and one wrong conclusion was found and
withdrawn by execution (the 0.488 um height).

An adversarial audit of `study_ring_resolved_unmixing` raised findings across
four lenses - identifiability, comparison fairness, statistical validity and
simulation fidelity. Nothing invalidated the model ordering. Two points survived
scrutiny and are worth recording:

- `achievedFP` printing 0.0466667 in every row is an algebraic identity, not a
  check: with 150 repeats the quantile lands on order statistic 143 exactly, so
  the strict inequality always counts 7/150 regardless of the null distribution,
  and it is computed from the same draws that set the threshold.
- "matched false-positive rate" means matched order-statistic index, not matched
  realised rate. The realised rate of a threshold at the 143rd of 150 null draws
  is Beta(8,143): mean 5.30%, standard deviation 1.82%.

Neither changes any reported power, because power is computed from the threshold
and `achievedFP` is never read. The confirmatory re-run uses 300 repeats, where
the quantile is interpolated and the nominal rate is exactly 5%.
