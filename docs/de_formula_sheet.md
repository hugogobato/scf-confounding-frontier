# DE Formula Sheet (SCF Phase 1, WP 1.3)

Companion to `code/de_formulas.py`. Every formula lists: statement, derivation sketch, source, and the expected finite-n correction order. Status tags: EXACT (algebraic identity), TRANSCRIBED (classical result, numerically self-checked in `tests/test_identities.py` before use), NEW-DE (our derivation, validated against simulation before being relied on).

Conventions: sigma_u^2 = 1 unless shown; l_j eigenvalues of Lambda Lambda'/sigma_u^2; tau_j = 1 + l_j population eigenvalues of Sigma_X/sigma_u^2; c = p/n; gamma in factor coordinates; u_j = Q e_j population eigenvectors.

## F1. OLS bias (EXACT)

Statement: plim beta_OLS = beta + Sigma_X^{-1} Lambda gamma with Sigma_X = sigma_u^2 I + Lambda Lambda'; since Lambda gamma = sum_j sigma_u sqrt(l_j) gamma_j u_j, the component form is bias_j = [sqrt(l_j/sigma2)/(1 + l_j)] gamma_j on u_j (with sigma2 = sigma_u^2).

Sketch: Cov(X,Y) = Sigma_X beta + Lambda gamma is an exact second-moment identity; OLS is the sample projection, consistent for the population projection.

Exact finite-n strengthening (used by the pilot overlay): write the "error" w = Lambda f + eps = X a + zeta where a = Sigma_X^{-1}Lambda gamma and zeta := w - X a. Under A2/A3 (Gaussianity) zeta is independent of X with E zeta = 0, so

    E[beta_OLS] - beta = E[(X'X)^{-1} X' X] a + E[(X'X)^{-1} X' zeta] = a + 0,

exactly, for every n with p <= n (and even p = n). Consequence: at c <= 1 the simulated mean-bias overlay tests the code, not the asymptotics; deviation beyond MC error indicates an implementation bug.

Correction order: none (exact). Source: linear algebra; the zeta-decomposition argument is standard (Gaussian conditioning).

## F2. Ridge bias (population identity; finite-n correction O(1/n))

Statement: plim beta_ridge = beta + (Sigma_X + lam I)^{-1} Lambda gamma; component form bias_j = [sqrt(l_j sigma2)/(sigma2(1+l_j) + lam)] gamma_j.

Finite-n: E[(Sigmahat+lam)^{-1} Sigmahat] a differs from (Sigma+lam)^{-1}Sigma a at O(1/n) (fluctuation of the shrinkage operator around its limit); the zeta term still has mean exactly zero. Expected overlay deviation at n = 2000: well under 1 percent.

Source: standard ridge population algebra; DE fluctuation order from Dobriban-Wager (2018) style calculus.

## F3. PCA-k trim tradeoff (EXACT at population level)

Statement: with V_k the top-k population eigenvectors,

    beta_trim = P_k beta + V_k (V_k' Sigma_X V_k)^{-1} V_k' Lambda gamma,

so trimming removes exactly the bias carried by dropped directions and leaves retained-direction bias unchanged (bias_j for j <= k equals the OLS value). The price is signal loss 1 - ||P_k beta||^2/||beta||^2, which is severe under dense beta (A4a): this is why estimation claims instantiate A4b (perturbed sparsity) while detection claims use A4a (ledger Section A4).

Correction order: sample-eigenvector overlap fluctuations, vanishing for supercritical spikes, non-vanishing near/below sqrt(c); pilot reports trim results descriptively only.

Source: population algebra; overlap caveat from BGN (2011).

## F4. BBP outlier location (TRANSCRIBED, CORRECTED)

Statement: for l > sqrt(c),

    mu(l) = sigma2 (1 + l)(l + c)/l = sigma2 (1 + l)(1 + c/l),

else the sample eigenvalue sticks to the bulk edge (1+sqrt(c))^2.

CORRECTION to research plan Section 2.3 item 4, which wrote (1+l)(1+cl)/l: that expression equals the correct one iff (l-1)(c-1)=0. The plan's own mitigation anticipated this class of error ("misremembered classical constants"). Unit test `test_bbp_location` simulates n=6000, p=3000, l=2 (mu = 3.75; the plan's formula would predict 3.0) and confirms the corrected form within finite-size tolerance, plus a convergence trend at n=12000.

Source: Baik-Ben Arous-Peche (2005), Thm 2.1; cross-checked against Johnstone (2001) conventions. Correction order: sample-spike fluctuation O(n^{-2/3}) around mu(l).

## F5. BGN eigenvector overlap (TRANSCRIBED)

Statement: for a unit vector a independent of the noise, |<a, v_hat_j>|^2 -> xi(l,c) = (1 - c/l^2)/(1 + c/l) for l > sqrt(c), else 0.

Sanity anchors used in review: xi -> 1 as l -> inf; xi -> 0 as l -> sqrt(c)+; xi(2, 0.5) = 0.875/1.25 = 0.7.

Source: Benaych-Georges-Nadakuditi (2011) Adv. Math. 227(1) (real case); rectangular version (2012) JMVA 111. Correction order: O(n^{-1/2}) fluctuation of the squared overlap; pilot averages 40+ reps at n >= 6000 to hold the SE far below the 2-decimal tolerance.

## F6. Marchenko-Pastur bulk edges (TRANSCRIBED)

Support [(1-sqrt(c))^2, (1+sqrt(c))^2] x sigma2. Used for scree context and as the subcritical spike ceiling. Source: Marchenko-Pastur (1967); Bai-Silverstein (2010).

## F7. Tracy-Widom threshold for lambda_max (TRANSCRIBED)

White null: (lam_max(X'X) - mu_np)/sigma_np -> TW1 with mu_np = (sqrt(n-1)+sqrt(p))^2, sigma_np = (sqrt(n-1)+sqrt(p))(1/sqrt(n-1)+1/sqrt(p))^{1/3}; cov-scale threshold divides by n. Quantiles TW1 q95 = 0.9793, q99 = 2.0234 (published tables). Calibration verified by Monte Carlo: rejection rate 0.05 +/- 0.01 over 5000 reps (`test_tw_calibration`). Source: Johnstone (2001).

## F8. Min-norm OLS bias for c > 1 (NEW-DE, PARTIALLY VALIDATED)

For p > n, OLS means the minimum-norm solution beta_hat = P_row Y. Decompose w = Lambda f + eps = X a + zeta, a = Sigma^{-1}b, zeta independent of X (A2/A3). Then E[beta_hat] - beta = E[P_row] a exactly (the zeta term vanishes; E[P_row] bounded for p >= n+2).

Two Phase-1 findings recorded during pre-validation (small-c=5 simulations):

1. Fit-artifact cancellation under A4a: because beta is redrawn per replicate with E[beta] = 0 and P_row depends only on (f, u), the shrinkage artifact (E[P_row] - I) beta has conditional expectation zero across replicates. The pilot's ||mean(beta_hat - beta)|| therefore equals the gamma-attributed bias at every aspect ratio; the artifact lives in per-rep ERROR instead (compared against RMS(rel_err) as a secondary check). An earlier draft prediction conflated the two functionals.
2. Capture coefficients: writing cap_j := coefficient by which the population OLS bias component j survives, the WP 1.5 pilot MEASURED the effective capture via fixed-Q projections and found a clean one-line law, validated at ~0.5-1 percent across supercritical (l = 3 sqrt(c)) and subcritical (l = 0.5 sqrt(c) and l = 1e-4) components, at both n = 400 and n = 2000, c = 5:

       cap_j = (1 + l_j)/(c + l_j)      (c > 1; PILOT-VALIDATED CONJECTURE),

   with correct boundary anchors (l -> 0 gives the uniform rowspace fraction 1/c; l -> inf gives 1). Example: aux cell (l = (3 sqrt(5), 1e-4, 1e-4), theta = pi/6) predicts ||bias|| = 0.1914 against a simulated 0.1927 (0.7 percent); mixed/super cells predict capture 0.6583 against measured 0.655-0.659. The earlier xi-based guess xi(l,c) + (1 - xi(l,c))/c predicted 0.607 there and is superseded (kept in code as bgn_capture_superseded for audit). STATUS: conjecture with strong multi-cell numerical support; derivation is a T1 work item (it is a statement about E[P_row] on spike directions that BGN overlaps alone do not deliver). All c > 1 overlays in the pilot memo use this law.

Sanity anchors of the validated part: r=0 gives E[beta_hat|beta fixed] = beta/c (isotropic min-norm fact, visible in RMS error, not in bias); l -> inf gives coefficient -> sqrt-scale OLS limit; the zeta-decomposition makes the confounding part exact given the capture limits.

Correction order: xi_j and capture limits carry O(n^{-1/2}); overall deviation target 10 percent at n = 2000, p = 10000 for supercritical cells.

## F9. Effective detection spike (working definition, heuristic frontier)

s_eff = ||P_spike Sigma^{-1/2} b||^2 = sum_{supercritical j} [l_j/(1+l_j)] gamma_j^2 (the sigma2 cancels: (sigma_u sqrt(l_j) gamma_j)^2 / (sigma_u^2 (1+l_j))). Only supercritical-aligned confounding produces a coherent rank-r shift in the whitened cross-moment statistic S1; subcritical mass hides in the bulk (s_eff = 0 at leading order regardless of g). The s_eff-to-power curve is Phase 2 material (OMH power template); Phase 1 uses s_eff only to label cells.

Source: composition of Johnstone (2001) null geometry with BGN (2011) spike structure; power template Onatski-Moreira-Hallin (2013).

## F10. Onatski ratio rule (TRANSCRIBED, approximate constants)

ED test: rhat = max{k: eig_k/eig_{k+1} > crit_k}; critical values (2.19, 2.09, 2.04, 2.01, 1.99, ...) transcribed approximately from Onatski (2010) Table 1. Flagged APPROXIMATE; baseline-selection only in Phase 1; exact table transcription scheduled before Phase 2 baselines.

## F11. Bai-Ng selector (PROVISIONAL approximation)

PC-style penalty khat = argmin_k [tail-sum + k * mean(eig) * log(max(n,p)) * (n+p)/(np)]. Exact Bai-Ng (2002) constants deliberately deferred (not load-bearing for Phase 1); flagged here and in code docstring.

## Self-check map (tests/test_identities.py)

F1: test_population_ols_identity (1e-10) and exact-finite-n via simulator smoke.
F2: test_ridge_population_identity (1e-10).
F3: test_pca_trim_tradeoff (1e-10).
F4: test_bbp_location (simulation, 3 percent + trend).
F5: test_bgn_overlap (simulation, +-0.02).
F6/F7: test_tw_calibration (rejection 0.05 +/- 0.01 over 5000 reps).
F8: test_minnorm_isotropic_limit (r=0 reduces to beta/c) + pilot overlay at c=5.
