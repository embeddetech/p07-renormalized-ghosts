# Process history — Paper 7 (renormalized) — "The multiplicity horn survives renormalization")

Ghost–normal 1+1-D QFT vacuum-decay falsification test (Krylov + Schur fold +
dimension ladder + MPS). Repo slot paper7; the program's sixth paper.

## 2026-07-14 — FIRST GATE AUDIT (counter 0 -> Opus x3). FINDINGS NOT YET APPLIED.
Fresh-eyes, no-histories auditors (claims / numerics / bibliography), all Opus
per counter-0 rule. Realign prereq (the 6 placeholder cites) was DONE in the
prior session. **Findings surfaced but UNAPPLIED — session was redirected to
p1/p5/p6 before the p7 apply cycle. Apply next (see HANDOFF NEXT ACTION).**

- **numerics (Opus): 1 MAJOR + 2 MINOR + 3 NIT.**
  - **MAJOR (CONFIRMED by re-run): Table 1 (tab:verdict), row "M=4, cap 4,
    λ=0.04": full/naive = 0.8964 is WRONG.** Auditor ran that config: actual
    full/naive = **0.9434** (Schur 0.9972). 0.8964 is identical to the stage3
    (M=4, cap 6, λ=0.02) cell — copied from the wrong cell. Fix: main.tex
    line 247, 0.8964 -> 0.9434 (verify-before-apply: re-run the config first).
  - MINOR: large-M Table 1 rows (M=10/12/16, M=6 cap6) produced only by
    t5_campaign.py `job val` writing results_val_M*_cap*.txt — **none archived**;
    values untraceable to stored output (campaign.py logic is correct). Rerun +
    archive, and add to Code-availability result list.
  - MINOR: Sec-mechanism "0.8620" elimination (diag-only AND vertex-only each
    reproduce 0.8620) has **no backing script or archived .txt**; campaign.py
    prints only a combined renPT2 substitute for |Q|>40k, not the two separate
    eliminations. Unverifiable as written.
  - NITs: extensivity reproduced 7.4% vs "within 7%" (say ~7.5%); ghost
    24% vs "25% above floor"; "L=32 ... negative per-site" point unarchived
    (only L=16/dt=0.04 negative point on disk).
  - Everything cheap/reproducible ELSE checked out: channel census 10/21/36/55;
    calibration 5.68e-16; ladder P1 +1.51% / P2 0.156 / P3 2.263; cap4/cap6
    Schur ratios; λ-deficit 0.0057/0.0227/0.0883/0.3144; stage4 a/b/2S,
    born-vacuum 471->7, d=2 Schur 1.0000. MPS floors traced to archived .txt.

- **claims (Opus, adversarial): 1 MAJOR + 5 MINOR + 2 NIT.**
  - **MAJOR: attribution-by-elimination gap (Sec mechanism + abstract).** "The
    entire deficit is second-order shell-splitting + Z-factor" rests on
    (i) an additivity/non-interference assumption (individual-block nullity does
    not entail the residual block is solely responsible) and (ii) is asserted at
    M=8 cap6 — the very config where the full fold "exceeded memory," so neither
    the full fold nor a shell-shell-only fold was actually computed there. Fix:
    soften abstract/Sec-mechanism to "consistent with" shell-splitting, OR run a
    shell-shell-only fold at a reachable config.
  - MINORs: λ-scan 31% (λ=0.04) vs Table's only λ=0.04 row = 10% deficit —
    scan's M unstated, looks contradictory (state the M); abstract "exact
    Feshbach/Schur folding" vs body's static-E0 approx (qualify "static-E0");
    d=2 exponent 2.26 in abstract lacks the body's "one shell deep, mechanism
    check not scaling" caveat; a=2.03e-4 said to "bracket" plateau 1.56-1.84e-4
    but sits ~10% above it (reword "sits just above"); "not a one-dimensional
    accident" from a single one-shell d=2 spot check (downgrade).
  - NITs: \date "July 13, 2026 --- drafted; not yet audited" status token ->
    plain date (house rule); floor 5.73e-2 vs 5.7e-2 rounding.
  - HELD strong: central thesis (fold tracks exact => renormalization) properly
    scoped to finite-channel/1+1/static-E0; not-claimed items (volume scaling,
    born-vacuum log growth, infinite-channel limit) consistently un-asserted.

- **bibliography (Opus): CLEAN.** All 8 external refs OK (metadata +
  claim-support): Barandes, GSTZ2021 (journal title "metastable" unhyphenated,
  correct), Deffayet2026/EwasiukProfumo2026 (2604.* ids resolve, not
  hallucinated), CJM2004, Lagnese2021, AbelSpannowskyWilliams2025,
  BurshteinGoldstein2024 (2410.21941 = PRB 111 174303; golden-rule-onset
  topical fit holds despite superconducting-impurity title). Companion*
  self-cites consistent. Counter **0 -> 1**.

Counters after this pass (findings pending application):
biblio **1** (clean); claims **0** (MAJOR+MINORs); numerics **0** (confirmed
MAJOR + MINORs). Publish HELD.

## 2026-07-14 — FINDINGS APPLIED (same session, after user "follow up on all").
- **numerics MAJOR — FIXED + INDEPENDENTLY VERIFIED.** Re-ran run_case(4,4,0.04)
  myself (scratchpad verify_p5_table.py, capped): full/naive = **0.9434**,
  full/Schur = **0.9972**; and confirmed the erroneous 0.8964 is exactly the
  M=4/cap6/λ=0.02 cell (wrong-cell copy, as the auditor diagnosed). Table 1
  line 247: 0.8964 -> 0.9434, "(tracked)" -> 0.9972.
- **claims MAJOR — SOFTENED.** Abstract "the entire deficit is second-order
  splitting ... plus a Z-factor" -> "the deficit is consistent with ..."; Sec
  mechanism "located, by elimination" -> "attributed, by elimination (each
  partial dressing individually null; the full fold at this configuration is
  memory-limited, so the attribution leans on the reachable-M folds of Table 1)";
  "what the full fold uniquely captures" -> "captures at reachable M"; abstract
  "exact Feshbach/Schur folding" -> "exact (static-E0) ...".
- **claims MINORs APPLIED:** λ-scan labelled "at M=3, cap 6"; abstract d=2 2.26
  gets "one shell deep, a mechanism check rather than a scaling measurement";
  "not a one-dimensional accident" -> "the same behavior holds in the one d=2
  shell we can reach"; a=2.03e-4 "bracketing" -> "sitting just above" the plateau.
- **NITs APPLIED:** \date "July 13, 2026 --- drafted; not yet audited" -> plain
  "July 13, 2026"; floor 5.7e-2 -> 5.73e-2 (body self-consistency); extensivity
  "within 7%" -> "within ~7.5%" (x2, matches reproduced 7.4%).
- Recompiled clean (**7pp**, 0 errors/undefined).

**REMAINING p7 NUMERICS FOLLOW-UP (needs script reruns — NOT silently dropped;
carry to p7's next numerics pass):**
  (1) Archive the large-M campaign outputs (results_val_M10/12/16*, M=6 cap6
      from t5_campaign.py `job val`) so Table 1's large-M rows are traceable;
      add them to Code-availability. (2) Back the 0.8620 diagonal-only /
      vertex-only partial-fold eliminations with a committed script (no current
      script computes the two separate eliminations) — the mechanism-§ hedge
      now stands in until then. (3) The L=32 negative-per-site point lacks an
      archived file.
Counters after application: biblio **1**; claims **0** (MAJOR applied);
numerics **0** (MAJOR applied). Publish HELD.

## 2026-07-14 — NUMERICS FOLLOW-UP CLOSED (script reruns done, same session).
Per user "do the script reruns, get everything finalized."
- **Large-M campaign rows REPRODUCED + ARCHIVED** (reran t5_campaign.py `job val`,
  λ=0.02, sparse-LU Schur): results_val_M10_cap4.txt (0.8959/0.9951 ✓),
  M12_cap4 (naive **0.8566** now archived, Schur 0.9941 ✓), M16_cap4
  (**0.7737/0.9924** ✓, |P|=561, |Q|=15,416), M6_cap6 (naive 0.8833, Schur
  **0.9884** ✓, |Q|=22,056). **Every value matches Table 1 exactly.**
- **0.8620 elimination BACKED** by a new committed script
  `scripts/t5_partial_folds.py` (matrix-free evolution + sparse matvecs — the
  full Schur fold OOMs at |Q|~81k, the partial folds do not). At M=8/cap6
  (dim **82,065**, |P|=1001, |Q|=**81,064**): full/naive = **0.8620**,
  full/diag-only = **0.8620**, full/vertex-only = **0.8620** — both partial
  dressings fail to move the naive ratio, exactly as the mechanism-§ asserts.
  Cited in Sec mechanism + Code-availability.
- Code-availability updated (t5_partial_folds.py + results_val_*/results_
  partial_folds_* listed); checkpoint blobs (ckpt_H_*.npz) gitignored; result
  .txt kept. p7 recompiled clean (7pp). major-findings 18/19 follow-ups DONE.
- **Numerics provenance gaps from the first pass are CLOSED** — the lane is now
  ready for a clean numerics certification pass. Only residual: the L=32
  negative-per-site MPS diagnostic point lacks an archived file (a diagnostic
  aside, not load-bearing; note at next MPS touch). Publish HELD.

## 2026-07-14 — CLAIMS PASS 2 (Fable, adversarial). MAJOR (L=16 cert) applied.
Per-lane cadence (claims flips Opus->Fable). **1 MAJOR + 5 MINOR + 2 NIT applied.**
MAJOR (major-findings 21): abstract/Discussion claimed the MPS "certifies at
L=16 against ED to 6e-6" but the ED comparison is the L=4 gate run; L=16 is
production-only -> scoped to "against ED at L=4 (to 6e-6), production at L=16"
(**pass-1 Opus MISSED this**; model-diversity). MINORs: intro elimination
"identified as" -> "consistent with" (pass-1 softened abstract+Sec-mech but not
the intro); Sec-3 body first fold gets the "(static-E0)" qualifier (pass-1 only
did the abstract); Sec-mechanism +explicit "a shell-shell-only fold was not
separately computed"; Discussion "Exact" class -> attribution moved to "inferred,
not exact"; λ-scan "λ^4 detuning entering" -> "sub-×4 growth ... marking a λ^4
correction" (deceleration reads as a correction, not an additive term). NITs:
25%->24%; "ED box test" -> "L=4 ED gate test". Counter 0. Clean (7pp). The L=32
MPS diagnostic residual is being resolved this session (run in progress).

## 2026-07-14 — FULL WAVE: claims pass 3 (Opus) + biblio (Fable) + numerics (Opus, trace).
**claims: 1 MAJOR + 1 MINOR.** MAJOR (ledger 27, INCOMPLETE-FIX): the wave-3 L=16-cert fix
hit only 2 of 4 instances -> fixed the intro + Sec-mps closing ("certified at L=4,
production at L=16"). MINOR: χ40→64 "moved the floor nowhere" -> "did not lower the floor
(it rose slightly)". **biblio (Fable): 2 CLAIM-SUPPORT issues (reset 1->0)** — Deffayet2026
mischaracterized ("integral-of-motion vacua don't exist for φ²χ²" -> rephrased: their
construction is (x²−y²)-built, x²y² only a nonintegrable perturbation with bounded motion);
CJM2004 "low-dimensional log structure" detached (CJM has NO logs; the log rung is
Companion2's). Other 6 external refs CLEAN. **numerics (Opus, trace-only): numbers trace
cleanly, ZERO contradictions; PROVENANCE gap** — small-M Table rows, λ-series, dimension
ladder, rate window, born-vacuum, 5.7e-16 calibration have no committed stdout. DEFERRED
behind the MPS (needs runs): commit those stdouts; the L=32 negative-point residual = the
running MPS. Counters: claims 0, biblio 0, numerics 0. Clean (7pp).

## 2026-07-14 — WAVE 5: claims pass 4 (Fable) + biblio pass 2 (Opus). CONVERGENCE.
**claims (Fable): no MAJOR** — the L=16-cert fix CONFIRMED COMPLETE (all 5 occurrences
uniformly "gate at L=4, production at L=16"; no stray "certifies at L=16" — Stage 16 fix
validated). 3 MINOR + 3 NIT; APPLIED: "overshoots by 31%" -> "misses by 31%" (caption
convention, 2 passes flagged); λ⁴ correction marked "negative (saturating)"; census
10/21/36/55 clarified vs |shell| (coupled off-diagonal pairings vs full shell dimension).
No-action (marginal/hedged): abstract fork-(a) wording; floor "saturation"; shell-splitting
operational clause. **biblio (Opus): CLEAN -> 0->1** — the reworded Deffayet ((x²−y²)-based;
biquadratic is a bounded nonintegrable perturbation, their Fig 5) + CJM (pure power-law, no
logs; log via Companion2) CONFIRMED accurate to PRIMARY sources; all 13 refs OK. **SS5_NUMERICS.md
fixture committed** (cores-free provenance). Counters: claims 0, biblio **1**, numerics 0.
Clean (7pp).

## 2026-07-14 — L=32 MPS diagnostic RESOLVED (confirms the paper).
The long-running L=32/dt=0.01/χ=40 MPS diagnostic completed (~5.7h). dt01_L32.txt:
dE_end_per_site = **−4.12e−3**, rate_per_site = −9.33e−4 — the per-site transfer is
genuinely NEGATIVE, confirming the §mps "an L=32 point whose per-site transfer came out
negative, i.e. non-extensive under truncation." The full-wave numerics auditor's flag (no
L=32 archive on disk) is now CLOSED; archive committed + recorded in SS5_NUMERICS.md. No
paper edit needed.


## SESSION 2026-07-18/19 (numbering unification + claim-support)
- Folder paper5/ -> paper7/ (series 7); cite keys -> codenames; release-gate leak removed; header/slug corrected (was "p05").
- Claim-support (Opus): CLEAN; all 8 external attributions SUPPORTED, incl. the load-bearing Deffayet (x^2-y^2)-interaction claim (confirmed vs full text -> the pure-biquadratic choice is "forced, not preferred") and both 2026 preprints. Cosmetic (noted, not fixed): BurshteinGoldstein cite-key year 2024 (arXiv) vs journal 2025 -- harmless convention.

## OPUS 4-PAPER WAVE 2026-07-19 (out-of-rotation audit; claims + numerics + biblio)
- **Claims: 1 MAJOR + 4 MINOR (all applied).** MAJOR (ledger 58): Sec:model mode-set definition wrong/self-inconsistent -- text said "M modes, k=-M/2..M/2" but the model + census (10/21/36/55 = M(2M+1)) + scripts use 2M modes (k_n=2πn/L, n∈{±1..±M}, zero dropped); results unaffected (scripts correct), text unreconstructable. Fixed prose; sibling-swept. MINOR (applied): GSTZ "suppression hypothesis"->"metastability suggestion" (GSTZ only suggested); census "off-diagonal/additionally diagonal" double-book removed; "Exact" strength-bucket split (static-E0 fold agreement re-filed as Empirical, since the fold is an approximation with documented drift; only the shell degeneracy + analytic channel-table calibration to 5.7e-16 stay under Exact); line-221 "exact (static-E0) fold" -> "(static-E0) fold". No theorem spine (paper is empirical-by-design; strengths bucketed) -- no numerical result dressed as theorem-grade.
- **Numerics: CLEAN.** Every quantitative claim reproduces from scripts/committed outputs (census, Schur folds, λ² identification, Trotter tradeoff floors, MPS L=4/16/32, barrier decision, d=2 exponent 2.26); heavy MPS runs verified vs committed .txt rather than relaunched. SS5-flagged prior copy-error confirmed already fixed.
- **Biblio: 7/8 CLEAN + 1 fix.** GSTZ2021 title "metastable"->"meta-stable" (arXiv 2007.05541-verified, hyphenated -- NOT the earlier unchecked-hyphenation mistake). Deffayet2026/EwasiukProfumo2026 descriptors ("(x²-y²)-type", "Fock diagonalization") plausible, flagged provisional per pdf-reads-provisional.
- Housekeeping: SS5_NUMERICS.md header/footer "paper5"->"paper7" (stale codename label, internal file).
- **DEFERRED (MINOR, next pass):** χ (ghost field vs MPS bond dimension) and L (circle length vs lattice site count) notation collisions -> a dedicated p7 notation pass (rushing a multi-site rename risks regressions). Also un-instantiated \newtheorem{proposition}/{remark} in preamble; "negative-energy Fock vacuum" loose phrasing; "6%" rounding (all harmless NITs).
- Recompiled 7pp, 0 undefined; both gates PASS.

## NOTATION MICRO-PASS 2026-07-19 (deferred collisions resolved)
- chi (MPS bond dimension) -> D at all bond-dim sites (chi=40/64/100, log chi); the ghost field chi untouched (29 sites). L (continuum circle length) -> ell at k_n=2*pi*n/ell + "circle of length ell"; the MPS lattice site count stays L (=4/16/32, standard). Recompiled 7pp, 0 undefined.

## OPUS CLAIMS REVALIDATION 2026-07-20 (Wave B; 2nd Opus pass)
- **Mode-set MAJOR fix (ledger 58) VERIFIED CORRECT**: reads "2M momentum modes (k_n=2πn/ℓ, n∈{±1..±M}, zero dropped)", consistent with t5_channel_ladder.py (2M modes) and the census 10/21/36/55=M(2M+1) (the old M-mode reading gives 3/6/10/15 -- definitively 2M). |shell|=561 at M=16 resolved CONSISTENT: 561 = 1 vacuum + 32 two-particle + 528 four-particle(=M(2M+1)); "shell dimension larger than 528" correct. Empirical spine all SOUND, nothing dressed above grade (λ² horn, D-barrier decision, Trotter floor, strength-tiering disciplined). Applied 1 NIT (another L->ℓ rename residue): channel-coupling g=λ/(V√...) used stray V for the 1+1d quantization volume = circle length -> ℓ (V appeared as a symbol only here). Skipped 1 harmless NIT (line 460 "Certified ... L=4 (production at L=16)" -- barrier section already clarifies L=16 is not ED-certified). Recompiled 7pp, 0 undefined; both gates PASS.

## OPUS DEEP CLAIMS 2026-07-20 (Wave E; Nth adversarial pass) -- SOUND, 1 MINOR applied
- No derivational error; strength-tiering / static-E0 fold / census all hold. Applied MINOR: Table 1 rows M=10/12/16 (cap 4) were unlabeled under the interjected lambda=0.04 row, inviting misattribution of the headline 0.995->0.774 mode-slide; the slide is fixed-lambda=0.02 (confirmed vs result files) -> annotated M=10 row "lambda=0.02". Left NIT ("|shell|=561" is |P|=vacuum+shell; shell alone 560 -- defensible vacuum-inclusive E=0 subspace per Wave B; argument survives).

## OPUS DEEP CLAIMS 2026-07-20 (hunt-what-was-rationalized wave, p5-p8) -- SOUND, 2 MINOR applied
- No CRITICAL/MAJOR; the ledger-58 2M mode-set / census $M(2M+1)$ (10/21/36/55) / |shell|=561, |Q|=15416 at M=16 REPRODUCED EXACTLY from t5_channel_ladder.py + t5_stage3.py (independent of the manuscript); grep confirms no stale "M modes"/"k=-M/2..M/2" sibling; all abstract<->body numerics match.
- Applied: (MINOR, "Relation to prior work") "the suppression hypothesis tested here is that of~\cite{GSTZ2021}" attributed a stronger claim to GSTZ than they made (they raised a metastability SUGGESTION / benign-rate possibility; the "suppression" fork is the paper's own pre-registered fork (b)); breaks the program's suggestion-not-hypothesis discipline (cf. L39/L441 "metastability suggestion") -> "the metastability suggestion". (MINOR, Trotter/truncation) χ/bond-dim collision: "re-compression to $\chi$" literally reads as re-compression to the ghost FIELD χ; the bond dimension is $D$ everywhere -> "$D$". (The companion L collision is absent: circle length ℓ vs sites L cleanly separated.)
- HELD/reported (NIT): abstract two-null-partial-folds result established at one config (M=8) without an abstract signal (downstream hedge at L456-457 "inferred, not exact" is honest); "runs in production at $L=16$" mildly oversells a barrier-limited run (clarified in-context each instance).
- Recompiled 7pp, 0 undefined. Both gates PASS. No MAJOR.

## OPUS NUMERICS 2026-07-20 (Wave 6; reproduce-from-scripts) -- 2 NIT applied
- Ran all fast scripts + verified the multi-hour MPS runs against committed fixtures. Every spine claim reproduces: census 10/21/36/55, |P|=561/|Q|=15416 at M=16, mode-slide 0.995->0.774 (lambda=0.02), calibration 5.68e-16, ladder exponent 2.263, d=2 48-channel shell, born-vacuum 471->7, Trotter floors 0.337/0.057/0.231, D:40->64 floor 0.057->0.062, L=4 ED gate 6e-6, 24%-above-floor, L=32 negative transfer -4.12e-3. ledger-58 2M mode-set RE-CONFIRMED reproducing.
- Applied 2 NITs: extensivity ~7.5% -> ~7.4% (reproduces at 7.44%, 2 sites) (#48); rate-window "(6%)" -> "(~6%)" (actual 5.6%) (#49). Held-defensible NIT (#50): tab:verdict "|shell|=561" label = |P|=vacuum+shell (shell proper 560); value correct, re-confirm of prior Wave-B/E disposition.
- NOTE: the numerics agents' re-runs appended to one committed fixture (results_partial_folds_M8_cap6.txt); restored via git checkout.
- Recompiled 7pp, 0 undefined. Both gates PASS. No MAJOR.

## FABLE CLAIMS 2026-07-20 (cross-model wave p2/p4/p7) -- 6 MINOR + 1 NIT applied (2 Fable-MAJOR downgraded)
- Fable adversarial: 2 claimed MAJORs on the central inference, both DOWNGRADED to MINOR on verify (conclusion survives; Fable agrees). Applied:
  - (F1 #72) central discriminator "if multiplicity erodes, [fold] cannot track" -- false as stated (within-P erosion IS fold-tracked) -> "erodes by a mechanism beyond the static second-order dressing". (Fork operationally self-consistent; conclusion holds via small deficit + lambda^2 scaling.)
  - (F2 #73) the M-growing folded-ratio drift (1.0000->0.9884) attribution was asserted -> hedged "consistent with static-E0 degradation" + scoped "fork (b) excluded only above ~1%". (Bulk 31%->1% absorption is decisive.)
  - (F3 #74, ledger-58 SIBLING MISS) "at most 16 modes" -> "M=16 (32 modes per field)".
  - (F4 #75) abstract fork (a) "self-energy renormalization" (which the elimination nullifies) -> "renormalization (dressing) of the states".
  - (F5 #76) "|shell|=561" -> "|P|=561" -- RESOLVES the twice-held #50 (paper defn P=vacuum+shell => shell=560; Fable's paper-internal argument decisive).
  - (F8 #77, NIT) abstract "shell saturates" -> "saturates with cutoff".
- HELD (#78): F7 undefined symbol "S" in 2S_exact=4.53e-5 (needs a one-clause definition -- FLAG); F6 "exactly" (parenthetical discloses precision); F9 cluster.
- Recompiled 7pp, 0 undefined. Both gates PASS. Ledger stays 61 (both Fable-MAJORs verified MINOR).

## HELD-ITEM CLEARED 2026-07-21 -- undefined symbol S defined
- minor-ledger #78 F7: "$2S_{\rm exact}$" had no definition of S. From t5_stage4.py (S_exact = sum of couplings g over the resonant E~0 channels, b=2*S_exact), added inline "$S_{\rm exact}=\sum_{\rm shell}g$, the summed coupling over the exactly-resonant channels". Recompiled 7pp, 0 undefined, both gates PASS.

## 4-WIDE FABLE CLAIMS AUDIT 2026-07-23 — 1 MAJOR (ledger 72), Opus-verified against the script
- MAJOR (ledger 72): the recently-added (#78 F7) inline S_exact gloss dropped the square. Script t5_stage4.py builds each channel g=abs(colv)**2 (matrix-element modulus SQUARED), so S_exact=∑|g|² and b=2∑|g|² (resonant channels: population g²t² ⇒ slope 2∑g²). Value 2.27e-5/136ch plausible only as ∑|g|² (⟨|g|⟩~4e-4), implausible as ∑g (~1.7e-7). Opus confirmed against the source before applying (not on the auditor's word). Fix: "∑_shell g"→"∑_shell |g|²" (matrix-element modulus squared); number unchanged; no siblings.
- Census/2M lane (prior-MAJOR site, ledger 58) re-verified EXACTLY clean by independent enumeration (10/21/36/55=M(2M+1); |P|=561 M16cap4; |Q|=22056 M6cap6; etc). Symbol collisions resolved (D bond-dim, ℓ length). 
- Recorded MINOR/NIT (minor-ledger #86): "full degenerate shell" prose incomplete (contains 2-/6-particle states); Λ cutoff + w1..w4 undefined; intro erosion sibling least-hedged; t5_stage2.py uncited. Held.

## W21 WIDE-OPEN AUDIT 2026-07-26 (Opus, adversarial; claims+numerics+biblio) — 1 MINOR applied; NO MAJOR
- **NUMERICS CLEAN cross-model:** every quoted signature reproduces AND both prior MAJOR-class defects are disclosed+congruent. #75 `m_chi=1.35` model switch: the only unequal-mass entry is the born-vacuum diagnostic, now explicitly labeled "a deliberately detuned twin"; part (A) rate window + all tab:verdict numbers use `mchi=1.0` (equal); 471->7 reproduces. Gate-cert injection: `t5_mps_gate.py` injects the exact ground state ONLY to isolate Trotter/truncation error (disclosed verbatim); production L=16 ANNEALS (`prepare_vacuum` = annealed 2nd-order imaginary-time TEBD). Census 10/21/36/55=M(2M+1), |P|=561 (M16), Schur 0.7737/0.9924, partial folds 0.8620, MPS floors 0.337/0.057/0.231, ladder 2.263/0.156/+1.51% — all match committed fixtures.
- **MINOR applied (#142):** the P3 d=2 exponent 2.263 is produced by `t5_dimension_ladder.py` P3, which sweeps cutoff `Ms2=(2,3,4,6,8)` and fits the last three (M=4/6/8) — a multi-cutoff sweep, NOT "one shell deep" (that descriptor fits the *separate* t5_stage4 48-channel ED spot check). Reworded 3 sibling sites (abstract L66, body L358, summary L473): "one shell deep" -> "a short-lever cutoff sweep" (body: "(cutoff M=4--8, three-point fit)"). Direction was conservative (understated the lever arm) so no scientific change; accuracy fix. Re-grep "one shell deep": 0 stale.
- **HELD (biblio, machine-provisional, needs human-PDF):** `BurshteinGoldstein2024` in-text role "Fermi-golden-rule onset in discrete-mode systems" vs listed title "Quantum simulation of the microscopic to macroscopic crossover using superconducting quantum impurities" (PRB 111 174303 (2025), arXiv:2410.21941) — possible role/title mismatch; a superconducting impurity in a structured bath CAN show FGR onset so likely DEFENSIBLE, but confirm against the PDF before any edit. Do NOT fabricate. Standing provisional (unchanged): Deffayet2026/EwasiukProfumo2026 2026 preprints.
- **NOTED, not applied (internal script docstrings, paper is correct):** t5_stage3.py print "naive deficit should grow ~lambda^4" is a PRE-REGISTERED a-priori guess (detuning-only) superseded by the lambda^2/Z-factor data — left intact by design (rewriting a pre-registration after data is bad practice; paper correctly reads lambda^2). t5_stage4.py docstring "window OPENS PROGRESSIVELY / open-channel count grows with t" reads backwards vs its own code (|dE|<2pi/t shrinks with t) and the paper's correct "window can close... 471->7" — internal narration only, no paper/repro impact.
- Recompiled 7pp, 0 undefined, 0 overfull >10pt. Both gates PASS. **p7 -> 0/3 (took a MINOR).**

## W22 WIDE-OPEN AUDIT 2026-07-26 (Fable, adversarial; claims+numerics+biblio) — 1 MAJOR (#83) + 3 MINOR + 2 NIT applied
- **MAJOR #83 (undisclosed box-length heterogeneity in Table 1)** — see major-ledger. FIXED by disclosure (§2 ℓ=8/ℓ=16 + m=1; table caption box-split + box-robustness note). Verdict unchanged (Schur band robust to ℓ).
- **MINOR-2 (method mislabel, 3 sites):** "exact Krylov evolution" (abstract L46, §2 L213, Discussion L466) — the code uses scipy `expm_multiply` (Al-Mohy--Higham truncated-Taylor scaling), NOT a Krylov method; evolutions are numerically exact so only the label was wrong -> "exact sparse-propagator evolution" (replace_all).
- **MINOR-3 (reproducibility gap):** Code-availability list omitted `t5_prod_dt01.py` (verified to exist; produced the dt01_*.txt fixtures carrying the dt=0.01 floor 0.231, λ=0.12 ~1.9× separation, λ=0.03 floor-buried, L=32 −4.12e-3 numbers) -> added to the Sec:mps script list.
- **MINOR-4 (edge overclaim):** §mps "agrees to 2--6e-6 at t=1..4---three orders below signal" — at t=1 it's ~49× (under two orders); "three orders" holds only at t=4 -> "up to three orders below signal (by t=4; ~50× below at t=1)".
- **NIT-2 (standing F6):** partial dressings "each reproduce the naive ratio exactly" -> "to the four digits resolved" (0.8620×3 is a 4-digit null).
- **NIT-3:** §3 "naive ratio converges to ~0.956" — cap sequence non-monotonic (0.9950/0.9324/0.9602/0.9560) -> "settles non-monotonically near ~0.956".
- **HELD NIT-1 (script-side):** Table row M=4/cap4/λ=0.04 (0.9434/0.9972) reproduces but no committed script archives its stdout (provenance only in SS5 fixture note) — add an archived run in a script-hygiene pass.
- **CLAIMS effectively CLEAN:** census/|P| bookkeeping re-derived (10/21/36/55=M(2M+1); 561=1+32+528 at M=16); coupling g + 1/√2 pair factor re-derived; b=2∑|g|² gloss survives ledger-72; λ²/Z-factor read sound; both prior MAJOR-class items (#75 mchi=1.35 detuned-twin disclosed; gate exact-eigenvector injection disclosed, production anneals) re-confirmed congruent. **BIBLIO CLEAN** (standing machine-provisionals only; BurshteinGoldstein W21-held for human-PDF).
- Recompiled 7pp, 0 undefined, 0 overfull >10pt. Both gates PASS. **p7 -> 0/3.**

## W23 CROSS-PAPER BIBLIO CORRECTION 2026-07-31 (orchestrator, Crossref-verified) - 1 fix; counter reset
- p7 carried the same GSTZ2021 error as p3: "Is negative kinetic energy **meta-stable**?" vs Crossref (10.1103/PhysRevD.103.115025) "metastable". Third paper affected by the 07-24 biblio wave's arXiv-over-publisher premise. **major-ledger #84.** Found only because the corpus-wide sibling sweep was run rather than fixing the two papers the auditor happened to name. Rebuilt clean (1 overfull, 1.8pt), both gates PASS.
- NOTE: resets p7's counter. The W23 Opus claims/numerics audit of p7 was still running at this entry; its findings will be logged separately.

## W23 WIDE-OPEN AUDIT 2026-07-31 (Opus, adversarial) + FABLE CROSS-VERIFY - 1 MAJOR (#87) + 2 downgraded MINOR + ~10 MINOR/NIT applied; counter -> 0
- The Opus auditor raised **5 MAJOR candidates**. An adversarial Fable verifier (framed to refute by default, and told to MEASURE rather than argue) returned **1 CONFIRMED, 2 REFUTED, 2 DOWNGRADED**. Given p7 is the paper that already produced #75 and #83 -- both model-congruence misses -- applying five MAJORs on one model's read would have been reckless; two of them do not survive contact with the fixtures.
- **MAJOR #87 (CONFIRMED, and it was in the TITLE):** the "entanglement barrier" mechanism for the lambda=0 MPS floor is impossible in the paper's own model. At lambda=0 the prepared state is an exact eigenstate of the generator (verifier: real-time residual 6.1e-15), so exact entanglement is STATIONARY -- measured middle-cut entropy 0.274049 at every t in [0,12] -- and there is no free-field growth to outrun the S~log D capacity. The verifier also corrected the auditor's proposed replacement: it is not pure static representability either (discarded Schmidt weight beyond rank 40 is 2.2e-11, while the D=40 floor is ~70x the D=100 value). Actual mechanism: Trotter-injected excitation amplified by repeated re-truncation -- the SIMULATED state's entropy climbs 0.274 -> 0.768 while the exact state's does not. Reframed at all 7 siblings (title, abstract, intro, section title, body, "barrier-limited", discussion) to "accumulated-truncation floor", and the body now states the measurement instead of the story. **Every operative result survives** (floor values, saturation decision, no-volume-scaling claim). Why it survived so long: the paper's MECHANICAL account was correct and every number reproduced, so numerics lanes had nothing to flag -- and the paper had explicitly said it did not measure S(t), so no lane ever did. **A stated non-measurement is where an unfalsified story hides.**
- **MINOR (downgraded):** the ladder calibration's "every ladder statement" overreached its archive (calibrated at d=1, ell=8, M=4/6, cap 4; rungs run at ell=16 and P3 is d=2). The verifier ran the missing checks and the claim is true in substance (6.01e-16 at ell=16; exact d=2 vacuum column to 2.80e-16), so the sentence was scoped to the archive rather than weakened. **MINOR (downgraded):** the m_chi=1.35 twin detunes H_0 only, keeping m=1 vertex mode functions (|g|^2 off by 0.584 at the lowest modes) -- but the verifier built a self-consistent twin and got IDENTICAL channel counts and an identically-shaped transfer profile rescaled x0.69, and no quoted number depends on it. Disclosed by parenthetical.
- **Also applied:** "forced, not preferred" non sequitur; "what remains IS the shell-shell block" (Z factor + non-vacuum vertex corrections enumerated); "all reachable M/caps/couplings" -> "every configuration where the fold is computable" (M=8/cap 6 is reachable but its fold OOM'd, and the ratio falls monotonically with |Q|); "positive-norm, negative-energy Fock vacuum" self-contradiction; eta=0.1 smoothing width disclosed with its spread; **Garriga-Vilenkin now cited and added to the bibliography** (the scripts named it, the paper never did); EwasiukProfumo claim-support corrected to the source's actual method; box disclosure extended to all three ell=8 items; unused \newtheorem removed; certification line names D=100 vs D=40.
- **REFUTED (see auditor-error-log):** the saturation metric-mix (the auditor's own "21% lower" commits the metric mix it alleges) and the gate-certification claim (both configs disclosed at the certification site). minor-ledger #162.

## W24 SEAM-SWEEP + OPUS ADJUDICATION 2026-08-01 - 2 MAJOR (#91, #92), BOTH implicating W23's own edits; counter -> 0
- The seam-sweep raised two MAJORs and a directly contradicting earlier auditor, so an **Opus adjudicator** (cross-model: both disputants were Fable) settled them **by measurement**, not argument.
- **#91 -- the $Z$-factor attribution is impossible in the paper's own instrument.** Four sites credited the deficit's $\lambda^2$ piece to a vacuum-norm/leakage $Z$-factor. But the static Schur fold is Hermitized and evolved by expm(-i Heff t), hence exactly unitary on P (measured normP_eff = 1.000000 at every lambda) -- it carries $Z \equiv 1$ by construction. Measured leakage is two orders too small (6.28e-3 vs a 0.1380 naive deficit at M=8/cap6). The correct attribution is derived, not fitted: free energies vanish on the shell, so the splitting is FIRST order in lambda and enters as (dE t)^2, giving exactly lambda^2. **W23's orchestrator had propagated this pre-existing error to a fourth site on an auditor's advice** -- it made it worse, it did not originate it. Also corrected t5_stage3.py, whose pre-registered "should grow ~lambda^4" sat next to the right lambda^2 measurement.
- **#92 -- ENTIRELY SELF-INFLICTED, and the fastest catch in the program (zero waves).** The W23 fix for #87 added three diagnostics taken from a single verifier's report without re-derivation. The adjudicator measured the trajectory at three sizes: **L=4 does not climb** (0.274049 -> 0.275188) while **L=16 climbs 0.249 -> 0.889**, passing 0.77 at t~10.8. The published sentence spliced an L=4 starting value onto an L=16 endpoint as "the same window"; the 2e-11 discarded weight is likewise an L=4 quantity placed in an L=16 discussion, and could not be an L=16 number at all. Rewritten with every figure scoped to its size, plus the caveat that the L=16 initial state is itself variational. **Process rule now recorded: a number the orchestrator inserts on a subagent's report must be re-derived and must carry its configuration.**
- Also applied: eta range widened to ~0.01-0.15 (d=2 shifts 0.144 at eta=0.05, outside the old claim); the "free negative-norm-free" garble from W23 repaired; SS5_NUMERICS stale 7.5% -> 7.4%. Held: the Table-1 lambda=0.04 cell reproduces only via the module API, not by running a committed script.
- The auditor explicitly named the six multi-hour runs it skipped and scoped its verdicts accordingly. minor-ledger #169.

## W25 SEAM-SWEEP 2026-08-01 (Fable) - 1 MAJOR (#94, reproducibility, LEFT OPEN); 1 MINOR + 1 NIT applied; the W24 mechanism fixes CONFIRMED; counter -> 0
- **The W24 reframing holds.** The auditor confirmed the accumulated-truncation framing propagated to title, abstract, intro, section heading and discussion with **no stale "entanglement barrier" sibling anywhere**, that every quoted number now carries its measured configuration, and that all cross-references resolve. Both W24 mechanism fixes survive adversarial re-derivation: the lambda=0 lattice vacuum is an exact eigenstate (residual 6e-15, structurally exact via sector decoupling, hence valid at L=16 too), static representability is ruled out, and the fold's Z = 1 is a true statement about the instrument.
- **MAJOR #94, LEFT OPEN:** the eight floor-mechanism and Z-factor diagnostics are emitted by NO committed script and NO archived fixture, which falsifies the paper's two reproducibility statements. The VALUES are all confirmed -- the auditor independently re-measured every one, including re-evolving the L=16 trajectory to six digits at t=0..3 against the drafting log.
- **The orchestrator tried to close it by writing t5_floor_mechanism.py and FAILED TWICE, then withdrew the script rather than commit it.** (1) dense expm on a 6561^2 matrix, hung; (2) SVD'd the left-canonical Q factor -- whose singular values are all 1 -- so it silently returned log(bond dim) = ln(81) = 4.394449 for every state, a routine that runs cleanly and is entirely wrong; (3) corrected version matches the D=100 control (4.69e-5 vs the quoted 5e-5) and lands at 0.768 where the drafting run reported 0.77, but returns **0.642420 for the L=4 exact middle-cut entropy against the paper's 0.2740**, because one auxiliary-Hamiltonian term was reconstructed from a partial read instead of imported. **A wrong diagnostic is worse than a missing one.** Nothing committed.
- **Two lessons banked:** (a) the quoted entropies never say WHICH bond they cut, and different bonds give materially different values at L=4 -- #92's configuration-provenance defect recurring one level deeper, at the definition of the observable itself; (b) the diagnostic must import the gate's own Hamiltonian construction, not re-derive it.
- Applied: "two orders too small" -> "over an order of magnitude" (a factor ~22 is ~1.3 orders; orchestrator wording from #91) and "within 10^-3 (0.2740 -> 0.2752)" -> "to ~10^-3" (self-contradictory at rounding level; orchestrator wording from #92). minor-ledger #175.
