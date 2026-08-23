# SAF audit of the vendored PSF, aberration and tracking code

Does supercritical-angle fluorescence (SAF) get treated correctly? Relevant
because the system is NA 1.45 oil into a water-like sample, so the pupil
extends well past the critical angle, and SAF decays over ~100-200 nm from the
interface - making it the steepest axial reporter available in the range the
membrane occupies.

## 1. PSF models - correct where an interface exists

All four models share the same branch selection, and it is right:

```matlab
function root = positiveRoot(value)
    root = sqrt(complex(value, 0));
    root(imag(root) < 0) = conj(root(imag(root) < 0));
end
```

`complex(value,0)` forces complex arithmetic so a negative argument gives a
purely imaginary root rather than NaN or a clamp, and the `conj` step picks the
positive-imaginary (decaying) branch - the physical evanescent solution.

| model | interface | SAF |
|---|---|---|
| `scalarPSFBessel` | no | not applicable - no boundary to couple through |
| `vectorialPSFBessel` | no | not applicable |
| `scalarPSFBesselAirInterface` | yes | correct |
| `vectorialPSFBesselAirInterface` | yes | correct, and polarisation-resolved |

The interface models carry the whole chain: `cosSample` complex above the
critical angle, complex Fresnel coefficients via
`transmissionCoefficients(nGlass, cosGlass, nSample, cosSample)`, the height
decay `exp(... (nSample*cosSample)*emitterHeight)`, and separate `tp`/`ts` so
dipole orientation is respected - which matters, because vertical dipoles
radiate disproportionately into supercritical angles.

The homogeneous models are not wrong, they simply have no interface. They
cannot represent NA > n_sample as a model of this experiment.

Despite the name, `sampleGeometry = 'airOnGlass'` is a general stratified
Gibson-Lanni model parameterised by `nSample`, `nGlass`, `nImmersion` and
`nDesignGlass`. Water-on-glass is just `nSample = 1.33`.

## 2. Aberration estimators - split, and two are exposed

| estimator | forward model |
|---|---|
| `estimateCenterPointISMWavefront` | interface |
| `estimateFullStackISMWavefront` | interface |
| `estimateExperimentalWavefrontFromPSF` | **homogeneous** (`normalizedStackExplicitDetectorZPlanes`) |
| `estimateTwoPlaneISMWavefront` | **homogeneous** (same) |

The two homogeneous ones are a concrete risk on NA-1.45-into-water data. SAF
imposes a radially symmetric apodisation and phase on the outer pupil, and a
homogeneous fit has no term for it - so it will be absorbed into the retrieved
Zernike coefficients, most naturally into **spherical aberration**, which is
also radially symmetric. Retrieved spherical aberration from those two
estimators should therefore be treated as contaminated until checked.

`normalizedStackAirInterfaceZPlanes` already exists as the interface
counterpart, so the fix is a substitution rather than new physics. The check
worth running first: retrieve aberrations from the same data with both stack
builders and see how much spherical aberration moves.

## 3. FCS on the ISM detector - not present

There is no FCS code in this codebase. The `correlat` matches are Gaussian
process correlation kernels and MSD curves inside the tracking packages.

If it is written, the SAF-specific issue is that the **effective detection
volume becomes per-channel and height-dependent**. SAF adds a shallow,
strongly z-weighted collection term that differs between centre and outer
channels, so channel-pair cross-correlation amplitudes - and the `G(0)` to
concentration conversion - acquire a z-dependent bias for surface-bound
species. A single effective volume per channel would not capture it.

## 4. Particle tracking - fixed-sigma 2D Gaussian, no interface

`+curved_miet/ismDetectorRawResponse` (and the `+fluctuating_miet` twin) is:

```matlab
sigmaSquared = opts.psfSigmaUm^2;
samples = exp(-0.5 * (x.^2 + y.^2) / sigmaSquared);
```

No interface, no height dependence, no dipole orientation. Being 2D is not the
issue; two things still bite:

1. For surface-bound emitters SAF reshapes the lateral PSF relative to
   free space, so a fixed sigma is mis-specified and localisation carries a
   systematic component.
2. More seriously, `+fluctuating_miet` **models membrane height fluctuations
   and couples z to lifetime, while leaving the detector response
   z-independent.** Height therefore moves the modelled lifetime but not the
   modelled channel distribution, although physically it moves both. Real
   height fluctuations then have nowhere to go in the detector model and leak
   into apparent lateral displacement or excess localisation noise, biasing MSD
   and diffusion estimates. `+curved_miet` has the same gap: curvature puts a
   range of z inside one spot.

The minimum fix is to make `psfSigmaUm` a function of height, calibrated from
the interface model, so z enters the detector response as well as the lifetime.
A full fix uses the interface PSF per channel directly.

## Bottom line

The physics is implemented correctly where an interface is modelled. The gaps
are all "the interface-aware model exists but this code path does not call
it" - two aberration estimators, and the tracking detector response. Nothing
requires new theory, only routing the existing stratified model into the
places currently using a homogeneous or Gaussian stand-in.

Provenance: every statement above comes from reading the source, not from
running it. The PSF conclusions have since been corroborated by execution -
`simulate_ring_weight_vs_height` and the studies built on it run, and the
noiseless deviance profile in `diagnose_ring_height_identifiability` recovers
the true emitter height exactly, which exercises the interface PSF, the Fresnel
coefficients and the evanescent branch end to end. The aberration-estimator and
tracking findings in sections 2 and 4 remain unverified by execution.
