# SCF Literature Reading List (Phase 4, grounded)

Rule: one strong source per need; deduplicated by recency-subsumption;
every entry resolvable or flagged. Groups map to work packages in
`docs/SCF_Theory_Plan.md`.

## Group 1: spiked covariance deterministic equivalents (TP-1, TP-4, TP-6)

1. Baik, Ben Arous, Peche (2005), "Phase transition in the largest
   eigenvalue of large sample covariance matrices", Annals of Probability
   33(5). BBP locations mu(l) = (1+l)(l+c)/l (real case).
   https://doi.org/10.1214/009117905000000233
   [access note: projecteuclid PDF gated; working copy via library]
2. Benaych-Georges, Nadakuditi (2011), "The eigenvalues and eigenvectors of
   finite, low rank perturbations of large random matrices", Advances in
   Mathematics 227(1): overlaps xi(l,c) = (1-c/l^2)/(1+c/l).
   https://arxiv.org/abs/0910.2120 (working copy lit/bgn2011.pdf);
   https://doi.org/10.1016/j.aim.2011.02.007
3. Knowles, Yin (2017), "Anisotropic local laws for random matrices",
   PTRF 169. THE Route-A backbone for TP-1: anisotropic law
   Sigma^{-1}G(z) - Pi(z)Sigma^{-1} = O_prec(Psi(z)) (their Theorem 3.6),
   outside-spectrum strengthening (their Theorem 3.7), Pi(z) definition
   (their Section 3), outlier-region remark (their Remark 3.8 - the piece
   TP-1 must extend for the z=0-below-bulk evaluation along the spike).
   https://arxiv.org/abs/1410.3516 (working copy
   lit/knowles_yin_anisotropic.pdf); https://doi.org/10.1007/s00440-016-0715-8

## Group 2: resolvent calculus for ridge/min-norm regression (TP-1, TP-4)

4. Dobriban, Wager (2018), "High-dimensional asymptotics of prediction:
   ridge regression and classification", Annals of Statistics 46(1).
   Resolvent fixed-point skeleton; lambda -> 0 c > 1 limit.
   https://arxiv.org/abs/1507.03003 (lit/dobriban_wager2018.pdf);
   https://doi.org/10.1214/17-AOS1549
5. Hastie, Montanari, Rose, Tibshirani (2022), "Surprises in high-
   dimensional ridgeless least squares interpolation", Annals of Statistics
   50(2). Min-norm geometry, E[P_null] behavior at c > 1.
   https://arxiv.org/abs/1903.08560 (lit/hastie2022.pdf);
   https://doi.org/10.1214/21-AOS2133

## Group 3: null laws and power limits for spectral tests (TP-2, TP-3a)

6. Johnstone (2001), "On the distribution of the largest eigenvalue in
   principal components analysis", Annals of Statistics 29(2). TW1 null.
   https://doi.org/10.1214/aos/1009210544 [gated; library copy]
7. Onatski, Moreira, Hallin (2013), "Asymptotic power of sphericity tests
   for high-dimensional data", Annals of Statistics 41(3). THE TP-3(a)
   skeleton: likelihood-ratio asymptotics under the null (their Theorem 7)
   and the contiguity reading that scopes eigenvalue-alarm impossibility;
    template for frontier power.
   https://arxiv.org/abs/1306.4867 (lit/omh2013.pdf);
   https://doi.org/10.1214/12-AOS1050

## Group 4: estimator-class context (TP-4, TP-6)

8. Cevid, Buhlmann, Meinshausen (2020), "Spectral Deconfounding via
   Perturbed Sparse Linear Models", AoS 48(5). Transform class definition.
   https://arxiv.org/abs/1811.05352 ; journal DOI:
   https://doi.org/10.1214/19-AOS1873
9. Nava, Buhlmann, Sigrist (2026), "Spectrally Deconfounded Gradient
   Boosting", arXiv:2607.09371. EB tuning + linear special case read
   backwards for the collapse lemma.
   https://arxiv.org/abs/2607.09371

## Group 5: factor-number selection (used by TP-4 statement)

10. Onatski (2010), "Determining the number of factors from empirical
    distribution of eigenvalues", REStat 92(4). Ratio rule + consistency.
    https://doi.org/10.1162/rest.2010.11253

## Notes on dedup choices

- Bloemendal-Knowles-Yin isolated-spikes local laws are subsumed for our
  purposes by Knowles-Yin (entry 3) unless TP-1.2 requires edge-of-bulk
  refinements; revisit only if the glue resists.
- OMH (2010) "Optimal testing..." variants not needed beyond entry 7.
- The classical MP law needs no separate entry (Bai-Silverstein chapter
  would be the reference if a formal citation is required later).

## Acquisition snippet (open-access items; VERIFIED ids 2026-08-25)

```bash
mkdir -p scf/lit && cd scf/lit
curl -L -o knowles_yin_anisotropic.pdf "https://arxiv.org/pdf/1410.3516v4"
curl -L -o bgn2011.pdf                 "https://arxiv.org/pdf/0910.2120v3"
curl -L -o dobriban_wager2018.pdf      "https://arxiv.org/pdf/1507.03003v2"
curl -L -o hastie2022.pdf              "https://arxiv.org/pdf/1903.08560v5"
curl -L -o omh2013.pdf                 "https://arxiv.org/pdf/1306.4867v1"
# gated, use library for the published versions:
#   BBP 2005 (Ann. Probab. 33(5)), Johnstone 2001 (Ann. Statist. 29(2)),
#   Cevid et al. AoS version (journal DOI above), Onatski 2010 REStat.
```

NOTE: three earlier recalled arXiv ids in this list were WRONG
(1410.6033 / 1903-vs-2111 mixups / 1605.03239) and were caught by opening
the downloads - titles did not match. All ids above were re-verified via
the arXiv API title search on 2026-08-25 and by inspecting the fetched PDFs.
Keep this check whenever adding entries.

## Notes on dedup choices

- Bloemendal-Knowles-Yin isolated-spikes local laws are subsumed for our
  purposes by Knowles-Yin (entry 3) unless TP-1.2 requires edge-of-bulk
  refinements; revisit only if the glue resists.
- OMH (2010) "Optimal testing..." variants not needed beyond entry 7.
- The classical MP law needs no separate entry (Bai-Silverstein chapter
  would be the reference if a formal citation is required later).
