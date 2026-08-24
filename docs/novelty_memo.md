# Novelty Memo (WP 1.2, completed 2026-08-23)

Status of prior-art verification after completing the UNVERIFIED items left by planning (research plan Section 4.1): **G1 verdict stands: GO**, with corrections recorded below and one repositioning obligation already reflected in the contribution list.

## 1. Search log additions

| Source / query | Volume | Result |
|---|---|---|
| Semantic Scholar forward citations of arXiv:1811.05352 | 55 works | screened |
| OpenAlex forward citations of both Cevid records (arXiv W2900832492 + published W3117548431) | 24 unique | screened |
| Union S2+OpenAlex forward citations of arXiv:1811.05352 and arXiv:2211.01903 | 65 unique works | merged file `data/lit/all_forward_citations_merged.json` |
| Forward-citation keyword screen (phase transition, detectab, minimax, lower bound, Le Cam, adaptive tun, empirical Bayes, frontier, Tracy-Widom, BBP/Baik/Benaych/Johnstone, bias of, debias, outlier eigenvalue, spiked, spike model, factor structure, confounding strength) | 12 hits | none a direct hit (Section 2) |
| arXiv API: "deconfounding" AND "phase transition" | 0 | gap confirmed |
| arXiv API: "hidden confounding" AND "random matrix" | 0 | gap confirmed |
| arXiv API: "confounding" AND ("Benaych-Georges" OR "Tracy-Widom") | 0 each | gap confirmed |
| arXiv API: "confounder detection" AND "spike" | 0 | gap confirmed |
| arXiv API: "spectral deconfounding" | 4 | exactly the known family; no additions |
| arXiv API: "deconfounding" AND "bias" | 48 | all applied ML/econ titles, none theory-frontier |
| arXiv API: "detectability frontier" | 11 | unrelated fields (watermarking, asteroseismology, robotics) |

## 2. Screening outcome on the 65 forward citations

The only hits worth recording as adjacent work (added to nearest-neighbor table):

1. Guo-Zevid-Buhlmann, Doubly Debiased Lasso (arXiv:2004.03758; JASA): inference for sparse beta under hidden confounding via decorrelation + de-biasing. Same problem class, different target (confidence intervals, not a bias map or detectability frontier), no spectral phase analysis. Becomes a cited neighbor and a natural referee touchstone for C2 assumptions.
2. Decorrelating/debiasing simultaneous inference line (2022-2023): targets loading significance, not our estimands.
3. High-dimensional GLM hidden confounding (2022), single-index latent factors (2025), longitudinal GLM double-debiasing (2025), mediation with latent confounding (2025): extensions into other models; none contains a bias phase diagram or a detectability/removability frontier.
4. Deconfounding via profiled transfer learning (2025): different problem; minimax in a transfer sense.

No paper combining dense hidden confounding with phase-transition/detectability language exists in any searched source, as of 2026-08-23. Fail-rule trigger absent.

## 3. Corrections to the research plan's reference layer

1. Cevid et al. was published in **JMLR 21 (2020), paper 19-545** (https://jmlr.org/papers/v21/19-545.html), not Annals of Statistics as the plan states. Venue correction recorded in references.bib.
2. Wang-Blei is JASA 114(528), 1574-1596 (plan said 527).
3. Onatski (2010) DOI is 10.1162/rest_a_00043.
4. BBP's Annals of Probability 33(5) 2005 paper carries the complex-sample title ("Phase transition of the largest eigenvalue for nonnull complex sample covariance matrices", doi:10.1214/009117905000000233); the real-case statement we use follows from the same analysis and is stated explicitly for the real case in BGN (2011). Bib entry uses the verified locator with a comment.
5. Bai-Ng (2002) Econometrica 70(1) 191-221, doi:10.1111/1468-0262.00272 (Crossref initially matched a different SSRN version; pinned manually against the journal record).
6. Published versions now pinned for the whole spectral-deconfounding family: Scheidegger et al., ACM/IMS JDS 2(1) 1-52 (doi:10.1145/3711116); Ulmer et al., JCGS 35(2) 758-768 (doi:10.1080/10618600.2025.2569602); Schur-Peters DecoR, JRSS-B 88(3) 799-818 (doi:10.1093/jrsssb/qkaf067); Rendsburg et al., Mathematical Statistics Learning 7(3) 189-220 (doi:10.4171/msl/47). Janzing-Scholkopf resolved: Journal of Causal Inference 6(1), doi:10.1515/jci-2017-0013 (flagged item closed).
7. LAVA (Chernozhukov-Hansen-Liao): the resolvable record found is the SSRN/Springer introduction (doi:10.2139/ssrn.1908964); the precise method-paper locator stays flagged until Phase 2 actually needs LAVA as a baseline.

## 4. OMH (2013) theorem-level read (WP 1.2 action 3)

Source inspected: arXiv:1306.4867 (AoS 41(3), doi:10.1214/13-aos1100). Findings that bind our C2 design:

1. Setting: X_i ~ N(0, sigma^2(I_p + h v v')), test H0: h=0 vs h>0 from sample covariance eigenvalues, direction v an unspecified nuisance. This is structurally the closest formal template to our detection problem under A4a (spiked perturbation, unknown direction).
2. Below the BBP threshold the null and alternative are mutually contiguous (Le Cam sense); power of eigenvalue-based tests remains larger than size for local alternatives and tends to 1 as h approaches sqrt(c) from below. The "impossibility threshold" reading of BBP is explicitly called overly pessimistic: information hides in small deviations of the full empirical spectrum from Marchenko-Pastur, captured by Stieltjes-transform-type linear spectral statistics, NOT by the largest-eigenvalue contrast alone.
3. Their Theorem 7 gives the null limit of the log likelihood-ratio process; asymptotic power follows via Le Cam's third lemma; explicit envelope formula beta(theta1; mu) = 1 - Phi[Phi^{-1}(1-alpha) - sqrt((theta1 - 1 + e^{-theta1})/2)] in their normalized parameterization.

Consequences for us:
- Our T3 lower-bound skeleton (two-point/contiguity) matches OMH's own technique; adaptation gap is the cross-moment geometry (our signal lives in Cov(X,Y), not in Var(X)).
- Honest frontier statements must distinguish "undetectable by max-eigenvalue statistics" (which our s_eff captures) from "information-theoretically undetectable" (false strictly below threshold at any distance, per OMH contiguity). Phase 2's numerical Le Cam probe and WP 2.3's alignment sweep are designed with this distinction in mind, and the paper language must respect it. This strengthens rather than weakens the story: the interesting boundary is where practical max-eigenvalue alarms stop working, which can sit strictly inside the contiguous region.
- S2-family statistics (LSS with data-chosen test functions) are not optional decoration: they are the OMH-informed answer inside the contiguity region.

## 5. Frozen contribution list (post-G1)

1. C1 (load-bearing): exact asymptotic bias phase diagram of OLS/ridge/spectral trims under dense hidden confounding as a function of (l_j, c, theta), with the decoupling of visibility and harmfulness.
2. C2 (load-bearing, ledger-conditional): detectability/removability frontier under A4a, TW-calibrated upper bound plus contiguity-aware lower-bound discussion following OMH.
3. C3 (repositioned): optimality verdict for EB-tuned spectral trimming measured against the derived frontier (never "we introduce EB tuning"; SDBoost benchmarked directly).
4. C4 (application): calibrated omnibus confounding alarm package validated on semi-synthetic benchmarks.
Cut by default: C5 universality theory (empirical section instead), time-series extension, M2 generality beyond one benchmark family.

Standing obligation carried into Phase 2: re-run the forward-citation screen (both APIs, same vocabulary) within one week of any submission.
