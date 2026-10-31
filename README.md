# ECHO-Rx — responsive biologic production

**Research concept v0.1. No ECHO-Rx biological or clinical results.**

ECHO-Rx is a therapeutic research branch of [ECHO-M](https://github.com/RayPeloquin/ECHO-M-Host-Coupled-Artificial-Chromosome). The proposed first system uses contained cultured mammalian cells and a tissue model: sense one specified inflammatory input, produce one biologically active protein, and measure both its effect and what remains after production is disabled.

The first proposed benchmark output is **interleukin-1 receptor antagonist (IL-1Ra)**. Related engineered-cell approaches already have published preclinical precedents. Those results do not validate this proposed combination or establish a new treatment.

## Read the design

**[Eight-page design brief — architecture, equations, evidence gates, failure analysis and Court disposition](ECHO-Rx-v0.1-Artificial-Medicine.pdf)**

The brief is dated **28 October 2026**; this repository was prepared **29 October 2026**. The contained biological system, carrier implementation and molecular enable/disable mechanism remain unbuilt and unqualified.

## The architecture

| Module | Proposed function | Evidence required |
| --- | --- | --- |
| Sensor | Respond to one declared inflammatory input | Specificity, sensitivity, latency, cross-talk and independent input reference |
| Producer | Reversible production of one active protein | Activity, basal leak, continuing secretion, variable producer count and cell burden |
| External intervention | Enable or disable production through an actual physical or biochemical mechanism | Measured intervention effectiveness; a software label has no molecular force |
| Carrier | Compare a simpler non-HAC cassette with a separately qualified human artificial chromosome (HAC) | Physical location, structure, copy state, function, burden and decision-relevant added value |
| Optional memory | Record a finite exposure history | Recording competence and assay linkage; the memory bit never latches production ON |
| Containment and observation | Account for cells and link each assay to its specimen and observation window | Escape, attrition, compartment activity, uncertainty and independent tissue outcomes |

The HAC branch retains the intention of avoiding host-genome editing. This is a requirement to investigate, not proof of non-integration. Published host-genome-engineered cell systems are functional precedents, not evidence for that proposed boundary.

The first biological decision is deliberately small: can the contained system meet predeclared active-output, leak, response and shutdown criteria across independent preparations? A carrier or memory module advances only after it earns its additional complexity in the declared comparison.

## Production OFF is not clearance

Intracellular stores, a depot, ongoing secretion and cell death can all affect residual active protein. Production, residual activity, physical cell removal, chromosome persistence and tissue recovery are separate endpoints. A falling inflammatory reporter cannot by itself establish disease recovery or safe shutdown.

### A mass-balance hypothesis

Let `N` be the viable producer-cell count and `B` the **total active intracellular protein mass** across the tracked producer population. `A_D`, `A_T` and `A_S` are active mass in depot, target and spillover pools. All `k`, `g` and `d_N` rates have units `1/time`; `r_syn` has units `mass/cell/time`; concentrations are `C_j = A_j / V_j`.

```text
dN/dt   = (g - d_N) N
dB/dt   = N r_syn - (k_sec + k_intra + d_N) B
dA_D/dt = k_sec B + phi d_N B - (k_out + k_deg) A_D
dA_T/dt = k_out A_D - (k_clear + k_spill) A_T
dA_S/dt = k_spill A_T - k_sys A_S
```

Here `0 <= phi <= 1` represents the fraction of active intracellular material released on death. The remaining `(1-phi) d_N B` is assigned to an inactivation/debris sink. Secretion is counted once, through `k_sec B`.

Uniform death content, instantaneous release, homogeneous compartments and linear clearance are assumptions. Competence loss, cell escape, biological processing delays, binding, transport backflow or changing volumes can require more states. No per-cell production cap establishes a total dose cap when producer number and activity can change. This model supplies no clinical exposure limit or automatic conservative bound.

### Reproduce the reduced shutdown example

[shutdown_model.py](shutdown_model.py) uses **Python 3.9+ and the standard library only**:

```bash
python3 shutdown_model.py
```

It checks an analytical two-compartment solution against a separate numerical RK4 integration and prints a JSON record. The expected record is in [shutdown-model-results.json](shutdown-model-results.json).

This script is a reconstruction from the published equations and design brief. Its checks remain active with `python -O` and validate the reduced mathematical illustration.

For an **ideal, instantaneous stop of all further secretion**, fixed coefficients, no leak and no new death-release input:

```text
a = k_out + k_deg; b = k_out; c = k_clear + k_spill

A_T(t) = A_T0 exp(-ct)
       + b A_D0 [exp(-at) - exp(-ct)] / (c-a)

If a=c: A_T(t) = [A_T0 + b A_D0 t] exp(-ct)
```

With the deliberately arbitrary values `A_D0=6`, `A_T0=1`, `a=0.2`, `b=0.15` and `c=0.5`:

```text
A_T(t) = 3 exp(-0.2t) - 2 exp(-0.5t)
Peak = 1.2804815896 at t = 1.7027520792
Initial target slope = +0.4
Integral from shutdown to infinity = 11
```

All values use normalized mass and arbitrary time units. **They are not biological measurements, doses or design targets.** The transiently loaded depot produces a post-stop peak. An exact pre-OFF steady state under these fixed coefficients instead gives a monotonically decreasing target curve. The script also checks that case, the equal-rate limit, mass balance and the integrated tail. This verifies the mathematical illustration, not the full biological model or a therapeutic platform.

## Evidence gates and limits

1. **Define scope:** cell background, carrier, input, output, tissue endpoint, observation window, thresholds, uncertainty and sampling units.
2. **Verify activity:** distinguish active protein and tissue function from fluorescence or total protein abundance.
3. **Measure response:** input dependence, absent-input output, cross-talk, latency, repeated cycles and loss of competence.
4. **Qualify shutdown:** effective disabling intervention, residual synthesis/secretion, stored protein, death release and active-protein decay.
5. **Compare carriers:** test a specific advantage for HAC under the declared boundary requirements; also measure stability and burden.
6. **Requalify coupling:** examine common failures, copy drift, clonal selection, optional memory burden and missing outcomes.

Keep deaths, losses and escaped cells visible in the denominator. Assays on different cells cannot automatically certify one specific living cell. A nominally independent assay can still share biological failure modes.

Inflammatory sensing is not a disease or infection classifier. A future implant would require separate evidence for transport, oxygenation, immune response, manufacturing consistency, retrieval and cell escape. Retrieval or killing cells does not instantly remove protein already released.

This repository contains no laboratory construction protocol, genetic sequence, dosing instruction or patient-use procedure. No ECHO-Rx efficacy, clinical readiness, universal programmability, first-in-field status or replacement of existing treatment is established.

## Primary references

| Reference | What it supports | Boundary |
| --- | --- | --- |
| [Choi et al., Science Advances 7, eabj1414 (2021)](https://doi.org/10.1126/sciadv.abj1414) | Autoregulated anticytokine production by host-genome-engineered constructs in studied mouse inflammatory-arthritis models | No human or HAC-based ECHO-Rx efficacy |
| [Collins et al., Gels 9, 169 (2023)](https://doi.org/10.3390/gels9020169) | Preclinical hydrogel-contained engineered-cell IL-1Ra delivery | No assumed immune isolation, retrievability or escape performance for a new device |
| [Cimino et al., Journal of Biological Engineering 19, 47 (2025)](https://doi.org/10.1186/s13036-025-00519-7) | In-vitro dual inflammatory/circadian response and tissue-model protection | Circadian control is outside the initial ECHO-Rx scope |
| [Birchak et al., Chromosome Research 33, 29 (2025)](https://doi.org/10.1007/s10577-025-09788-w) | Modular HAC construction and formation in studied cells | Does not establish coupled therapeutic function or control in the proposed cells |
| [Birchak et al., bioRxiv preprint (2026)](https://doi.org/10.64898/2026.01.28.702397) | A research direction for protein-only centromere seeding | Preprint; PubMed indexing is not peer review; no quantitative formation claim is imported |

## Credits and review

**Anubis / Ray Peloquin:** concept, questions and design direction.  
**Asteris Solenne / ChatGPT:** AI-assisted synthesis, source checks, mathematics, documentation, code and image generation.  
**Grok, Claude, Gemini and Meta AI:** adversarial review through the supplied Court tabs; feedback was checked and selectively adopted.

The Court is an AI-assisted critique workflow. Its agreement is not experimental validation, independent laboratory replication or institutional/clinical peer review. The brief records the adopted objections, limitations and source corrections. The illustrative shutdown depot value is `A_D0=6` throughout this release.

## Technical discussion

Useful criticism should identify a particular assumption, observable, equation, comparator or failure condition. Especially useful: independent activity measurement, identifiable shutdown dynamics, meaningful HAC advantages and evidence that would falsify the proposed interpretation. Do not post private medical or personal data in public issues.
