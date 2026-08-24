# Model Card (SCF Phase 1)

Project: Spectral Confounding Frontier (SCF). Ledger version of this card: see `assumption_ledger.md` (the pair is hashed together into every simulation config).

Status: ACTIVE artifact of WP 1.1. Every simulation config in later WPs carries a hash of this card plus the ledger. Any change to definitions here forces a ledger-version bump and a re-run flag on downstream results.

## 1. Data-generating models

### M1 (vector structural target)

One observational dataset consists of n independent copies of (f_i, u_i, eps_i), f_i in R^r, u_i in R^p, eps_i in R:

    X = Lambda f + u,        Y = beta' X + gamma' f + eps,

with Lambda in R^{p x r} dense loadings, beta in R^p the structural coefficient vector, gamma in R^r the confounder-response link, and X in R^{n x p} the observed design. The observable data are (Y, X) only; f is latent.

The confounding contribution to the conditional mean is entirely through the p-vector b := Lambda gamma ("bias direction"). Two parameterizations coincide by construction and both appear in code: the factor form above and the joint-Gaussian form (x_i, y_i) with Cov(X, Y) = Sigma_X beta + Lambda gamma. They generate the same second-moment behavior when all errors are Gaussian; simulations use the factor form.

### M2 (scalar treatment target)

As M1 plus an observed treatment D with

    D = pi' X + delta' f + nu,        Y = tau D + beta' X + gamma' f + eps,

with nu independent noise. Target of interest: scalar tau. M2 is used in Phase 3 only; all Phase 1 work concerns M1.

## 2. Normalizations (binding for every simulation and formula)

1. Noise scale: sigma_u^2 = 1 unless stated otherwise. All spike strengths below are defined relative to sigma_u^2, so results are scale-equivariant.
2. Spike strengths: l_1 >= ... >= l_r >= 0 are the eigenvalues of Lambda Lambda' / sigma_u^2, i.e., population eigenvalues of Sigma_X / sigma_u^2 are tau_j = 1 + l_j.
3. Aspect ratio: c = p/n in (0, inf).
4. Signal normalization: ||beta||_2 = 1. The response link is gamma = g * dir(theta), where g = ||gamma||_2 is a grid parameter and dir(theta) = cos(theta) e_1 + sin(theta) e_2 with e_j the coordinate basis of factor space R^r (theta = pi/2 gives gamma parallel to e_2). Factor directions e_j correspond to columns of Q in the loading model of Section 3.
5. Loading geometry: Lambda = Q diag(sigma_u sqrt(l_1), ..., sigma_u sqrt(l_r)) where Q (p x r, orthonormal columns) is drawn Haar-isotropically (QR of a seeded Gaussian matrix). "Dense loadings" means exactly this isotropic row distribution.
6. Errors: u_ij, eps_i standard normal at finite n (A2/A5); non-Gaussian variants only in robustness phases.

## 3. Why this geometry (interpretation notes)

With Q as above, the columns u_j := Q e_j are the population eigenvector directions of Sigma_X with eigenvalues tau_j = 1 + l_j. The alignment angle theta is measured between gamma and e_1, the direction of the strongest factor. Because Lambda gamma = sum_j sigma_u sqrt(l_j) gamma_j u_j, theta controls how much confounding mass rides on strong versus weak factor directions while leaving the design's spectral profile (l_1, ..., l_r), hence its visibility to scree/TW inspection, untouched. This is the knob that separates "visible" from "harmful" and it is why theta is a grid dimension everywhere.

Note recorded during Phase 1 planning review: an exactly zero OLS bias with nonzero spikes requires Lambda gamma = 0, which cannot happen for generic full-rank Lambda; "harmless" therefore always means "bias ratio small", operationalized numerically (<= 0.02 at pilot tolerances), never exactly zero. The plan's phrase "gamma orthogonal to the loading subspace" is implemented as "gamma aligned with negligible-strength factors" (auxiliary cells with secondary spikes l_2 = l_3 = 10^{-4}, whose OLS bias coefficient sqrt(l)/(1+l) ~ 0.01).

## 4. Estimands and primary functionals

1. E_beta (M1 target): the structural vector beta in Y = beta'X + gamma'f + eps. Under any estimator T, the bias functional is the loading-conditional expectation

       Bias(T) = E[T(Y,X) | Lambda] - beta,

   where the conditional expectation averages over f, u, eps and the A4a beta draw but CONDITIONS on the realized loading geometry Lambda (equivalently Q). Rationale discovered in Phase 1: if Q were redrawn Haar-isotropically on every replicate, the unconditional mean bias vector would vanish by spherical symmetry while each dataset still suffers an O(1) bias along its own realized spike directions; the conditional functional is the scientifically meaningful one and is what all simulations estimate (simulator draws Q once per config, flag q_fixed). The scalar summary used in gates is rel_bias(T) = ||Bias(T)||_2 / ||beta||_2.
2. Confounding-attributed bias: rel_bias_conf(T) = ||E[T under gamma | Lambda] - E[T under gamma = 0 | Lambda]||_2 / ||beta||_2, estimated with common seeds per rep. This isolates the gamma-induced shift from high-dimensional fitting artifacts that also exist when gamma = 0 (important for c > 1 where min-norm OLS shrinks beta even without confounding; twin arms run only for c > 1 cells because at c <= 1 the exact identity of F1 makes rel_bias_conf = rel_bias).
3. tau (M2 target): defined analogously; not exercised before Phase 3.
4. Visibility statistic: V(X) = lambda_max(X'X/n)/sigma_u^2 compared with the white-noise TW upper quantile threshold tw_threshold99(n, p, 1). A cell is called OUTLIER-POSITIVE if V exceeds the threshold.
5. Frontier quantities (defined now, measured later):
   - s_detect(l, c): minimal effective detection spike s such that the A4a-relative test of Section "Detection problem" in the ledger has power >= 1/2 at size 0.05.
   - s_remove(l, c): minimal spike strength such that worst-case (over theta) rel_bias of the best spectral estimator <= delta, delta = 0.05 by default.
   - Decoupling claim C1: there exist (l, c, g) regions with V below the TW threshold yet rel_bias >= 0.2 (invisible-yet-harmful), and regions with V clearly above the threshold (outlier-positive, l > 2 sqrt(c)) yet rel_bias <= 0.02 (visible-yet-harmless).

## 5. Units and conventions for formulas

All deterministic-equivalent formulas assume the normalizations of Section 2: noise variance 1 on the Sigma_X/sigma_u^2 scale; BBP outlier location mu(l) = (1 + l)(l + c)/l for l > sqrt(c) (this corrects a transcription slip in the research plan Section 2.3 item 4, which wrote (1 + cl)/l; see de_formula_sheet.md derivation F4 and its unit test); BGN overlap xi(l, c) = (1 - c/l^2)/(1 + c/l) for l > sqrt(c), else 0.

## 6. Glossary

- SCF: Spectral Confounding Frontier (this project).
- Spike: a population eigenvalue component l_j > 0 of Lambda Lambda'/sigma_u^2; subcritical if l_j < sqrt(c), supercritical if l_j > sqrt(c) (BBP regime).
- Effective spike: scalar summary of the confounding signal entering a given test statistic; made precise in the ledger's detection problem section.
- DE: deterministic equivalent (large-n deterministic approximation to a random quantity).
- LSS: linear spectral statistic.
- Harmful / invisible: refers to rel_bias >= threshold / V below TW threshold, respectively (thresholds fixed in Section 4.5).
