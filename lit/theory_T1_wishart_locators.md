# Theory T1 bibliographic locators: spiked-covariance bias reduction to central Wishart moments

Researched: 2026-08-25. Purpose: pin citable locators for four classical results used in the elementary proof
that reduces a spiked-sample-covariance bias calculation to ordinary central Wishart moments.

Verification methods used: doi.org resolution fetches, Crossref REST API, OpenAlex API, publisher front-matter PDFs,
author-hosted PDFs, Wikipedia (secondary, stable URLs), Monte Carlo/quadrature checks run locally on this machine.
Anything not confirmable online is marked honestly below.

---

## Item 1. Inverse-Wishart first moment

**Result.** If W ~ Wishart_p(Sigma, df=m) with Sigma > 0 and m > p+1, then

  E[W^{-1}] = Sigma^{-1}/(m - p - 1).

Equivalently, in inverse-Wishart parametrization X ~ W^{-1}_p(Psi, nu) (density prop. |X|^{-(nu+p+1)/2}
exp(-tr(Psi X^{-1})/2)), i.e. X^{-1} ~ W_p(Psi^{-1}, nu):

  E[X] = Psi/(nu - p - 1)   for nu > p+1.

**Status of formula:** CONFIRMED independently:
- Monte Carlo check (this machine, 2026-08-25): p=4, n=10, random SPD Sigma, 100,000 draws;
  max abs deviation of E[W^{-1}] estimate from Sigma^{-1}/(n-p-1) ~ 2e-4 relative. Script output archived in
  session log.
- Secondary sources stating it: Wikipedia "Inverse-Wishart distribution", Moments section:
  E[X] = Psi/(nu-p-1) (attributing to Mardia-Kent-Bibby p. 91 and Press 1982), and the review
  "Moments and Identities Involving Inverted Wishart Distribution" (WSEAS Trans. Math., 2020),
  https://wseas.com/journals/mathematics/2020/a285106-060.pdf , which lists "(m-p-1) E(A^{-1}) = Sigma^{-1}"
  attributed to Muirhead (1982), though WITHOUT a theorem number.

Candidate sources and locators:

1. **Muirhead, R.J. (1982). Aspects of Multivariate Statistical Theory. Wiley, New York.**
   - Book DOI (Wiley Online Library): 10.1002/9780470316559 ; Chapter 3 DOI: 10.1002/9780470316559.ch3
     ("Samples from a Multivariate Normal Distribution, and the Wishart and Multivariate Beta Distributions").
     Chapter-level locator VERIFIED via Crossref (api.crossref.org).
   - Exact theorem/exercise number for E[W^{-1}]: **locator NOT verified** (the result is commonly attributed to
     Chapter 3 of Muirhead, but I could not access any legitimate full text to pin the theorem number;
     Google Books snippet API returned HTTP 429 throughout, and the archive.org scan could not be located).
     Do NOT cite a specific Muirhead theorem number without checking the physical book.

2. **Mardia, K.V., Kent, J.T., Bibby, J.M. (1979). Multivariate Analysis. Academic Press.**
   - Cited by Wikipedia (Inverse-Wishart article) at **p. 91** for both the inverse-of-Wishart distributional
     fact and the mean E[X] = Psi/(nu-p-1).
   - Status: locator VERIFIED only *via the Wikipedia citation* (https://en.wikipedia.org/wiki/Inverse-Wishart_distribution);
     the page itself was NOT directly inspected. Treat as reliable-but-secondary.

3. **Gupta, A.K., Nagar, D.K. (2000). Matrix Variate Distributions. Chapman & Hall/CRC, Boca Raton. ISBN 1-58488-046-5.**
   - A searchable scan exists at archive.org (identifier `matrixvariatedis0000gupt`, lending-restricted):
     https://archive.org/details/matrixvariatedis0000gupt
   - Section/theorem number for the inverse-Wishart moments: **locator NOT verified** (search-inside and text
     derivatives are access-restricted; both html and curl attempts were refused).

4. **Eaton, M.L. (1983; reprint 2007). Multivariate Statistics: A Vector Space Approach. IMS Lecture Notes -
     Monograph Series 28.**
   - Chapter 8 "The Wishart Distribution" DOI: 10.1214/lnms/1196285114 (VERIFIED via Crossref).
   - Whether Chapter 8 states E[W^{-1}] explicitly, and its numbering: **UNVERIFIED**
     (Project Euclid is behind an Incapsula bot wall; neither webfetch nor curl could retrieve the chapter).

**Recommended citation practice:** cite the formula as "standard, e.g. Mardia-Kent-Bibby (1979, p. 91) or
Muirhead (1982, Ch. 3)" and/or prove it in one line (it follows by differentiating |W|^{t}: E[|W|^t] integrability
argument), avoiding an uncheckable theorem number entirely.

---

## Item 2. Wishart representation and orthogonal (rotational) invariance

**Result.** If Z is p x m with iid N(0,1) entries, then ZZ' ~ W_p(I_p, df=m). For any fixed orthogonal U,
UZZ'U' = (UZ)(UZ)' has the same distribution because UZ again has iid N(0,1) entries (orthogonal matrices preserve
the standard Gaussian measure). Consequently, if q is a unit vector in R^p, P = I_p - qq' is an orthogonal projection
of rank p-1, and Q is any orthonormal basis of {q}^{perp}, then for U (? x p) with iid N(0,1) rows,

  K = U (I_p - qq') U' = (UQ)(UQ)' ~ W_{rows}(I, df=p-1),

which reduces the claim to the definition of a central Wishart matrix.

**Locators:**

1. Definition S = GG' with iid Gaussian columns: Wikipedia "Wishart distribution", Definition section
   (stable URL: https://en.wikipedia.org/wiki/Wishart_distribution ) states exactly this construction. VERIFIED
   as an accessible statement (fetched 2026-08-25).
2. Muirhead (1982), Chapter 3 (same chapter DOI as above, 10.1002/9780470316559.ch3, VERIFIED via Crossref)
   is the standard textbook home of the definition and of the Bartlett/partitioned results. Specific theorem numbers
   (definition number, Theorem 3.2.x properties): **locator NOT verified** (full text inaccessible).
   For calibration only: Muirhead's Theorem 3.2.10 (partitioned Wishart) is cited by Bodnar & Okhrin (JMA 2008),
   confirming that the partitioned/properties material lives in Section 3.2, but I could not verify which numbered
   item carries the representation/invariance statement.
3. Eaton (1983/2007), Chapter 8 "The Wishart Distribution", DOI 10.1214/lnms/1196285114 (VERIFIED via Crossref).
   Eaton develops the Wishart distribution through projections of a Gaussian matrix (the vector-space approach),
   which is precisely the GPG' framework used here; his exact proposition number: **locator NOT verified**
   (Project Euclid unreachable behind bot protection).
4. Numerical sanity check of the projection claim (this machine): with q = (1,1,0,0)/sqrt(2), p=4, n=10,
   100,000 draws gave E[K] = (p-1) I to within 3e-3. Formula CONFIRMED.

Note: since orthogonal invariance of iid N(0,I) rows is immediate from the change-of-variables density, the proof can
legitimately cite only the Wishart definition (item 2.1 above) and derive K ~ W(df=p-1) in two lines; no deep
external theorem is required.

---

## Item 3. Marchenko-Pastur negative moments

**Result.** Let MP_a be the Marčenko-Pastur law with density

  f(x) = sqrt((b-x)(x-alpha)) / (2 pi a x),  alpha=(1-sqrt(a))^2, b=(1+sqrt(a))^2,  0 < a < 1,

(the ESD of (1/n)XX', p/n -> a, normalized so that the mean is 1). Then:

  (i)  integral x^{-1} dMP_a(x) = 1/(1-a)
  (ii) integral x^{-2} dMP_a(x) = 1/(1-a)^3

**Verification: INDEPENDENTLY VERIFIED twice over (do this rather than trusting memory):**

*Analytic derivation* from the standard self-consistent (quadratic) equation. With m(z) = integral (x-z)^{-1} dMP_a(x),
the Silverstein-Choi equation (Item 4, eq (1.3) below) specialized to H = delta_1 reads

  m(z) = 1 / ( -z + a/(1+m(z)) )   <=>   a z m^2 + (z + a - 1) m + 1 = 0.

For a < 1 there is no atom at 0, so m is analytic at 0. Setting z = 0:

  (a-1) m(0) + 1 = 0  =>  m(0) = 1/(1-a)  =>  integral x^{-1} dMP_a = m(0) = 1/(1-a).    [identity (i)]

Differentiating a z m^2 + (z+a-1) m + 1 = 0 in z and evaluating at z = 0:

  a m(0)^2 + m(0) + (a-1) m'(0) = 0
  =>  m'(0) = ( a m(0)^2 + m(0) ) / (1-a)
            = ( a/(1-a)^2 + 1/(1-a) ) / (1-a)
            = 1/(1-a)^3.
Since m'(z) = integral (x-z)^{-2} dMP_a(x), we get integral x^{-2} dMP_a(x) = m'(0) = 1/(1-a)^3.             [identity (ii)]

*Numerical quadrature* (SciPy, this machine, 2026-08-25), density integrated against x^{-1}, x^{-2}:

  a=0.10 : m_{-1}=1.11111111 vs 1.11111111 | m_{-2}=1.37174211 vs 1.37174211
  a=0.25 : m_{-1}=1.33333333 vs 1.33333333 | m_{-2}=2.37037037 vs 2.37037037
  a=0.50 : m_{-1}=2.00000000 vs 2.00000000 | m_{-2}=8.00000000 vs 8.00000000
  a=0.75 : m_{-1}=4.00000000 vs 4.00000000 | m_{-2}=64.00000000 vs 64.00000000
  a=0.90 : m_{-1}=10.00000000 vs 10.00000000 | m_{-2}=1000.00000001 vs 1000.00000000

Exponent 3 in (ii): CONFIRMED (both derivations agree to displayed precision).

**Citable base for the quadratic/Stieltjes equation:**

- Bai, Z.D., Silverstein, J.W. (2010). Spectral Analysis of Large Dimensional Random Matrices, 2nd ed.,
  Springer Series in Statistics. Book DOI: 10.1007/978-1-4419-0661-8 (print ISBN 978-1-4419-0660-1).
  Chapter 3 "Sample Covariance Matrices and the Marčenko-Pastur Law": Section 3.1 "M-P Law for the iid case"
  (3.1.1 "Moments of the M-P Law", p. 40), Section 3.2 "Generalization to the Non-iid Case" (p. 51, Theorem 3.10),
  Section 3.3 "Proof of Theorem 3.10 by the Stieltjes Transform" with **Section 3.3.1 "Stieltjes Transform of the
  M-P Law", pp. 52-53**, where the fixed-point equation for the MP Stieltjes transform appears.
  Section/page locators VERIFIED directly from the publisher front-matter (table of contents) PDF
  (link.springer.com/content/pdf/bfm:978-1-4419-0661-8/1.pdf, downloaded 2026-08-25). The internal display-number
  of the MP equation inside Section 3.3.1 was NOT inspectable (paywalled body), so cite "BS 2010, Section 3.3.1".
- Silverstein, J.W., Choi, S.I. (1995), equations (1.3) and (1.4) - see Item 4. VERIFIED from the authors' copy.

**Published statement of the negative-moment identities themselves:** NOT FOUND during this search
(queries for "negative moments" + Marcenko-Pastur returned nothing usable). Status: **external citation UNVERIFIED**;
recommend citing the derivation above together with BS 2010 Section 3.3.1 / Silverstein-Choi eq. (1.3), which is fully
self-contained.

---

## Item 4. Silverstein equation / deterministic equivalent of the resolvent

**Result.** For sample covariance matrices B_N = (1/N) T_N^{1/2} X X^* T_N^{1/2} with ESD of T_N -> H and p/N -> c,
the limiting spectral distribution F = F_{H,c} has Stieltjes transform m = m_F(z) characterized as the unique
solution in C^+ of

  m = ( -z + c integral lambda/(1 + lambda m) dH(lambda) )^{-1},

with inverse function z_F(m) = -1/m + c integral lambda/(1+lambda m) dH(lambda).

**Locator: VERIFIED with exact equation numbers.** Source: author-hosted preprint copy of

- Silverstein, J.W., Choi, S.I. (1995). "Analysis of the Limiting Spectral Distribution of Large Dimensional Random
  Matrices". Journal of Multivariate Analysis 54(2), 295-309. DOI: 10.1006/jmva.1995.1058 (VERIFIED via Crossref and
  ACM DL listing).
  Copy: https://jack.math.ncsu.edu/den.pdf (Jack Silverstein's NCSU page; downloaded and read 2026-08-25).
  Inside it:
    - equation **(1.1)**: F_0 = (1 - 1/c) 1_{[0,infty)} + (1/c) F (mass-at-zero convention),
    - equation **(1.2)**: Stieltjes inversion formula,
    - equation **(1.3)**: "m = m_F(z) is the unique solution for m in D^+ of m = (-z + c int lambda/(1+lambda m) dH(lambda))^{-1}",
    - equation **(1.4)**: z_F(m) = -1/m + c int lambda/(1+lambda m) dH(lambda), m in m_F(D).
  These are exactly the self-consistent/deterministic-equivalent equations requested.

Related, complementary pins:

- Bai & Silverstein (2010), Chapter 3, Section 3.3.1 "Stieltjes Transform of the M-P Law", pp. 52-53 (specialization
  H = delta_1). VERIFIED at section level from the publisher TOC (see Item 3); display-equation number inside the
  section NOT verified.
- Silverstein, J.W. (1995). "On the Empirical Distribution of Eigenvalues of a Class of Large Dimensional Random
  Matrices". J. Multivariate Anal. 54(2), 175-192. DOI: 10.1006/jmva.1995.1051 (VERIFIED via Crossref/OpenAlex).
  This is the general-population deterministic-equivalent paper often quoted alongside the above. Its internal
  equation number for m(z) = integral 1/(t(1-c-czm)-z) dH(t): **NOT verified** (no accessible full text found).

**Correction to the brief:** the suggested reference "Silverstein (1995), 'The smallest eigenvalue of large
dimensional Wishart matrices', Annals of Probability" appears to be a misremembering. Crossref/OpenAlex show NO such
1995 Ann. Probab. paper; Silverstein's 1995 journal papers are the three JMA articles listed above (DOIs
10.1006/jmva.1995.1051, .1058, .1083). The smallest-eigenvalue paper is:

- Silverstein, J.W. (1985). "The Smallest Eigenvalue of a Large Dimensional Wishart Matrix".
  Annals of Probability 13(4). DOI: 10.1214/aop/1176992819 (VERIFIED via Crossref/OpenAlex). It is about the almost
  sure limit of the minimal eigenvalue, NOT the resolvent fixed-point equation; use it only if the smallest-eigenvalue
  limit is needed.

---

## DOI verification requested (both corrected)

1. **Benaych-Georges & Nadakuditi (2011)**, "The eigenvalues and eigenvectors of finite, low rank perturbations of
   large random matrices", Advances in Mathematics 227(1), 494-521.
   - Expected DOI 10.1016/j.aim.2010.08.010 : **WRONG.** It resolves (doi.org, fetched 2026-08-25) to Comes &
     Ostrik, "On blocks of Deligne's category Rep_St(1)", Adv. Math. 226(2), 1331-1377.
   - **Correct DOI: 10.1016/j.aim.2011.02.007** - resolves correctly (VERIFIED via doi.org; metadata cross-checked on
     Crossref and OpenAlex: vol 227, iss 1, pp. 494-521, May 2011). arXiv version exists (arXiv:0910.2318) but was not
     separately verified in this pass.

2. **Knowles & Yin, "Anisotropic local laws for random matrices"**, Probability Theory and Related Fields.
   - Expected DOI 10.1007/s00440-016-0715-9 : **WRONG (HTTP 404; DOI does not exist).**
   - **Correct DOI: 10.1007/s00440-016-0730-4** - resolves correctly (VERIFIED via doi.org, Crossref and OpenAlex).
     Bibliography: PTRF **169**(1-2), 257-352; online August 2016, issue dated 2017 (so "PTRF 169 (2017)" is right
     for the print volume).

---

## Verification status summary

| # | Claim | Status |
|---|-------|--------|
| 1 | E[W^{-1}] = Sigma^{-1}/(m-p-1) (formula) | VERIFIED (Monte Carlo + secondary sources) |
| 1 | Muirhead 1982 exact theorem number | NOT VERIFIED (cite Ch. 3 level only) |
| 1 | Mardia-Kent-Bibby p. 91 | VERIFIED-via-secondary (Wikipedia citation) |
| 1 | Gupta-Nagar 2000 section | NOT VERIFIED (restricted scan) |
| 1 | Eaton Ch. 8 contains it | UNVERIFIED (bot-walled) |
| 2 | W=ZZ' definition + orthogonal invariance | VERIFIED as standard (Wikipedia def; Muirhead/Eaton chapters pinned at chapter level); deeper theorem numbers NOT VERIFIED |
| 2 | K = U(I-qq')U' ~ W(I, df=p-1) | Formula VERIFIED (derivation from definition + Monte Carlo) |
| 3 | int x^{-1} dMP_a = 1/(1-a) | VERIFIED (analytic + numeric) |
| 3 | int x^{-2} dMP_a = 1/(1-a)^3 | VERIFIED (analytic + numeric); exponent 3 CONFIRMED |
| 3 | Published statement of these two identities | UNVERIFIED (not found; use the derivation + BS Sec. 3.3.1 / SC95 eq. (1.3)) |
| 4 | Silverstein-Choi eqs (1.3)/(1.4) self-consistent equation | VERIFIED from author-hosted PDF |
| 4 | BS 2010 Ch. 3 / Sec. 3.1, 3.2, 3.3.1 pages | VERIFIED from publisher TOC PDF |
| 4 | BG&N correct DOI | 10.1016/j.aim.2011.02.007 (expected one wrong) |
| 4 | Knowles-Yin correct DOI | 10.1007/s00440-016-0730-4 (expected one nonexistent) |

## Key URLs

- BG&N: https://doi.org/10.1016/j.aim.2011.02.007
- Knowles-Yin: https://doi.org/10.1007/s00440-016-0730-4
- Silverstein-Choi author copy: https://jack.math.ncsu.edu/den.pdf | DOI 10.1006/jmva.1995.1058
- Bai-Silverstein book DOI: https://doi.org/10.1007/978-1-4419-0661-8 (TOC PDF:
  https://link.springer.com/content/pdf/bfm:978-1-4419-0661-8/1.pdf )
- Muirhead Ch. 3 DOI: https://doi.org/10.1002/9780470316559.ch3 (book: 10.1002/9780470316559)
- Eaton Ch. 8 DOI: https://doi.org/10.1214/lnms/1196285114
- Gupta-Nagar scan (restricted): https://archive.org/details/matrixvariatedis0000gupt
- Silverstein 1985 smallest eigenvalue: https://doi.org/10.1214/aop/1176992819
- Marchenko-Pastur (1967), Math. USSR-Sb. 1(4), 457-483: https://doi.org/10.1070/sm1967v001n04abeh001994
  (VERIFIED via OpenAlex; original source of the MP law.)
- Wikipedia Inverse-Wishart (mean formula): https://en.wikipedia.org/wiki/Inverse-Wishart_distribution
