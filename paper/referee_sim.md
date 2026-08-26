# Referee Simulation (WP 5.1d)

Two hostile simulated reviews were executed against the manuscript only (referees saw nothing but `paper.tex`, `parts/*.tex`, `references.bib`), per the plan's referee-simulation requirement: Referee A = RMT skeptic (score 38/100, reject), Referee B = causal-inference/applied skeptic (score 46/100, major revision). Full reports: `referee_review_A.md`, `referee_review_B.md`. Disposition of every objection:

| # | Objection | Disposition | Where answered |
|---|-----------|-------------|----------------|
| A1 | Fixed-r capture law rests on an unproved lemma; no complete proofs | Lemma given named proof mechanisms (Wishart concentration, odd-symmetry conditional means, inverse-Wishart calculus) with fixed-r concentration scope stated plainly; r=1 theorem is fully proved elementarily; proof-status honesty retained | Sec 3, Lemma (off-diagonal vanishing) |
| A2 | Ridge interpolation contains an internal sign contradiction | Verified algebraically and FIXED: correct identity is $(\mathrm{cap}(\lambda)-1)\langle\beta,q\rangle$ since $\mathrm{cap}(\lambda)-1=-\{1+(1+l)\bar m\}^{-1}$; reconciliation note added | Sec 3, Proposition (ridge capture) |
| A3 | Impossibility pillar contradicted by the authors' own falsifiers | Intro, abstract, roadmap, and limitations rewritten: proved population detachment boundary; pre-registered falsification of universal spec(M_aug) contiguity reported as instructive negative result incl. surrogate-S1 pipeline erratum; claims scoped to shipped alarm class | Abstract; Sec 1 item 3; Sec 5; Sec 8 |
| A4 | "Born blind at c=1" presented as universal when driven by profile calibration | General crossing rule $c_\star=(\omega_\star/w_0)^2$ now stated; coincidence with MP boundary explained as exact consequence of calibration $w_0=\omega_\star$, mechanism survives any calibration | Sec 5 visibility subsection |
| A5 | Simulation precision does not support "exact laws" rhetoric | Exactness claimed only for what is proved (r=1 identity; fixed-r DE); overlay tiers co-reported with the non-monotone tier caveat | Sec 3 validation |
| B1 | A4a does all identification work; promised beta-alignment size sweep absent | Gamma-angle alignment stress reported with numbers (S2 power 1.000 except 0.022 at theta=pi/2); deliberate beta-spike-correlation size sweep flagged as open axis in limitations rather than claimed | Sec 5 Fig alignment; Sec 8 |
| B2 | Abstract asserts impossibility the body refutes | Same fix as A3 | Abstract |
| B3 | Benchmarks cannot falsify dense-loading assumption; kill is a scoped strawman | Added limitation: injections build in assumed Haar geometry (validate calibration/power, not A1); T6 scope guard (collapse only under dense-signal bound) and W+ exclusion rationale already documented | Sec 8; Sec 6 |
| B4 | Practical delta thin; sensitivity analysis unnamed/uncited | Citations added (\citep{cinelli2020making}, \citep{vanderweele2017evalue}); combined workflow paragraph retained; package ships concrete actions (verdict, certificate, Onatski trim) | Sec 8 Discussion; Sec 7 |
| B5 | Projection-type competitors excluded from dominance class | Already answered in text: excluded because without 1/d division they are not estimators of beta (transmit confounding at O(sqrt l)); caught by first falsifier run, archived | Sec 6 |

Post-revision, the manuscript entered the independent quality-review loop (see `review_log.md`): an independent reviewer agent scores the paper as an external reader; revisions continue until the score exceeds 85/100.
