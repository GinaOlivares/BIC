# BIC Explorer — Phase 1

An analytical/CMT explorer for wall-controlled TE10 radiation zeros, quasi-BICs and critical coupling. This is not a replacement for a full-wave Maxwell solver.

## Run

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pytest -q
.venv/bin/python -m streamlit run app.py
```

## Reference and implemented equations

`Calculation_II.pdf` is authoritative: TE10 dispersion (pp. 19–20), wall parity (pp. 68–70, 100–103), reflection and absorption (pp. 73–75), critical coupling positions (pp. 77–78), quality factors (pp. 79–80), and TEM/pure-harmonic limitations (pp. 103–113).

- μ=μ0, ε=ε0 εr; fc10=1/(2w√(με)).
- β10=√(ω²με−(π/w)²); λg=2π/β10, only above cutoff.
- Cosine family px=(−1)^m; W̃=1−px exp(2iβL).
- γe=γ+|W̃|², equivalently 4γ+ sin²(βL) for even m and 4γ+ cos²(βL) for odd m.
- LB=qλg/2 (even) or (2q+1)λg/4 (odd).
- Near LB, γe≈4γ+β²(L−LB)².
- Qrad=ω0/γe, using the energy-decay linewidth convention.
- r=[iΔ+(γe−γi)/2]/[−iΔ+(γe+γi)/2], R=|r|², A=1−R.
- Positive equal rates yield critical coupling; zero rates yield no absorption.
- Approximate LPA=LB±√(γi/γ+)/(2β). Exact positions around each zero replace √(γi/γ+)/2 with arcsin(√(γi/γ+)/2). Physical candidates require L>b and 0<γi≤4γ+.

The specification's phrase “half-guide-wavelength shift” is inconsistent with its equations: adjacent even/odd families shift by λg/4; each family's period is λg/2.

## Units and assumptions

Physics uses SI exclusively, with explicit GHz/Hz, mm/m and Hz/rad/s conversions. γ is in s^-1, not Hz; γ/(2π) is the corresponding linewidth in Hz. Normalized |W̃|² and γe/γ+ are dimensionless. The default γ+=1 s^-1 is an explicitly arbitrary reference scale, not a calibrated prediction. An optional user-supplied γ+ changes dimensional linewidths and Q. γi is entered in s^-1.

The entered f is the assumed resonance f0. The localized eigenfrequency, its wall-dependent shift, dielectric mode profiles, absolute overlap normalization and resonator material dependence are unresolved; geometry radii affect only the schematic/clearance checks. The guide εr is the homogeneous background, not the cylinder dielectric constant. A narrowband frequency sweep freezes rates and β at f0 and uses detuning directly to avoid subtractive loss of precision for γ+=1 s^-1.

The symmetric reference model assumes W−=px W+, a PEC termination at x=−L, and a pure cosine harmonic. The dielectric outer cylinder extends through 0≤z≤h; the core occupies 0≤ρ≤a, h/2≤z≤h. It is not a disk.

A TE10 cancellation is only a candidate BIC. The resonance eigenvalue condition must also hold and every other open channel must decouple. In particular m=0 generally leaks through the TEM port. For pure m≥1 an ideal axisymmetric TEM port is symmetry-decoupled, but angular mixing can invalidate this. One-port spectra are conditional predictions; they do not represent the full system when extra ports radiate.

## Safeguards and verification

Below cutoff β is imaginary; no wall-cancellation calculation is performed. At cutoff β=0 and λg is undefined. The single-TE10 check uses the next cutoff min(fc20,fc01), sufficient since every other TE/TM cutoff is at least as high. Full channel enumeration is deferred.

A normalized squared coupling ≤1e-12 is marked as a numerical radiation zero; log10 Q is displayed at a cap of 16, never infinity. The exact ideal mathematical Q diverges; a small finite rate below the plotting threshold is not asserted to be exactly zero. Near-BIC diagnostics use a separate 1e-3 threshold. Critical coupling uses relative tolerance 1e-6 and absolute tolerance 1e-12 s^-1.

Tests cover cutoff, evanescence, even/odd radiation zeros, quadratic scaling, exact critical positions, perfect absorption, passive energy bounds, the zero-rate singularity, invalid inputs, tall-guide exclusion and finite Q plotting. Streamlit smoke tests exercise the default and cutoff states.

## Scope

Only Phase 1 is implemented. Parameter maps, expanded angular-family exploration, full radiation-channel analysis, field solutions and COMSOL import are deferred. Future COMSOL integration should compare eigenfrequencies, decay rates and angular content using documented units and normalization after analytical validation.

## Waveguide Mode Explorer

Open **Waveguide Modes** in the Streamlit sidebar, or visit
`http://127.0.0.1:8501/Waveguide_Modes` after starting the app.

This additional educational page implements the analytical TE/TM fields from
Calculation_II.pdf §§5–9 (pp. 17–21). It includes **10 TE and 10 TM presets**, plus
custom rectangular indices up to 10. These indices differ from the cylindrical
angular index used in the BIC model. TE00 and TM modes with a zero index are rejected.

Change frequency, guide dimensions, background permittivity and longitudinal
amplitude. Inspect all six components on y–z, x–y or x–z planes; the slice slider
moves the plane shown in the 3D locator. Magnitude, instantaneous real field at a
chosen time phase, and imaginary-phasor views are available. Normalized plots use
one common scale for the three E components and another for the three H components,
both defined at x=0 over the entire cross-section. This preserves component ratios
and decay when moving the plane. Disable normalization for V/m and A/m values.

The phasor convention is exp(iβx), with physical time dependence exp(−iωt).
Below cutoff β=iα selects exp(−αx) for x≥0; at cutoff β=0 is evaluated without
using an infinite guide wavelength. The arbitrary TE amplitude H0 (A/m) or TM
amplitude E0 (V/m) is not a power normalization. An optional NPZ download contains
complex SI fields, coordinates and model metadata, independent of display phase.

These are analytical solutions for the uniform homogeneous PEC guide, not fields
of the composite cylinder or the terminated BIC structure. This requested addition
does not implement parameter maps or COMSOL integration.

Additional tests check all 20 presets, tangential electric fields on all PEC walls,
Maxwell's curl equations in propagating/evanescent/cutoff regimes, the known TE10
profile, normalization bounds, slice coordinates, invalid indices and UI controls.
