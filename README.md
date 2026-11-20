# Finite replay and proof notes for Determinant zeros and method-specific barriers to linear-in-frequency pairing bounds in the Holstein model

Davis Rens  
SuperMatics Inc.  
[ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197).

Paper DOI: [10.5281/zenodo.24384373](https://doi.org/10.5281/zenodo.24384373). Supplement DOI: [10.5281/zenodo.24384385](https://doi.org/10.5281/zenodo.24384385).  
Source repository: https://github.com/SuperMatics-Inc/adiabatic-pairing-barriers

## Scope and evidence

This preprint asks whether an upper bound on electronic pair-ordering temperature in the spinful square-lattice Holstein model can improve to linear order in phonon frequency at fixed hopping and coupling. It neither proves nor disproves such a bound for the full Gibbs state. In its stated parameter range, a gap-conditioned expansion gives $e^{-d/2}$ clustering for specified fixed-support electronic monomials only when every static site field is $g$-gapped; the Gibbs weight near the Fermi level and the required joint buffer estimate remain unresolved. A time-mean-zero six-piece site-field history on an isolated $K_{1,3}$ star makes the single-spin determinant vanish with $\beta\|v\|_\infty<5$ at $\beta t=40,100,1000$. This obstructs the required cutoff only within the specified sup-norm family of the cited two-thirds-power bound, with its uniform retained-site control and free-dynamic-law Peierls probabilities. It is not a no-go theorem for all per-history methods or for the model's actual pair-ordering temperature. The separate Koma--Tasaki static diagonal-twist result gives a power-law correlation bound, not a linear-frequency ceiling.

No material prediction, measurement, or claim about realized pair order is made. This preprint has not been peer reviewed.

## Review status

The manuscript passed an independent review and two delta reviews. This is not peer review.

## Contents and replay

The manuscript source is `paper/main.tex`; its admitted bibliography is unchanged in `paper/refs.bib`. `replay/` contains the complete finite replay and source pins for both byte-identical proof-check scripts. `run_all.py` runs both baselines and every guarded mutant; a successful run rejects all 36 mutants. The replay checks finite calculations and witnesses. It does not formalize the conditional clustering argument, the Gibbs-field buffer estimate, or a full-Gibbs-state pairing ceiling.

`ADMITTED.sha256` records the original source hashes, including the 11-page PDF kept in release custody but not included in this source tree. `MANIFEST.sha256` covers every file in this tree except itself.

```sh
python -m pip install numpy==2.5.5 scipy==1.18.1 mpmath==1.4.1 python-flint==0.9.0
python -B replay/run_all.py
sha256sum -c MANIFEST.sha256
(cd paper && tectonic --keep-logs main.tex)
```

## Licenses and provenance

The manuscript text in `paper/` is distributed under CC-BY-4.0 through its paper record. The replay code and documentation use Apache-2.0. AI agents generated the proofs, replay code, and manuscript drafts under D. Rens' direction. Separate AI agents carried out the independent review. This preprint has not been peer reviewed.
