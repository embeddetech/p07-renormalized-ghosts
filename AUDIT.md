# How these papers are audited

Every manuscript in this program is written by an AI (Claude Fable 5,
Anthropic), directed and reviewed by a human author who is an engineer,
not a physicist. That provenance is disclosed on page 2 of every paper.
This document describes the quality gate each paper must pass before it
circulates, so that "audited" is a checkable claim rather than a vibe.
The complete pass-by-pass ledger for each paper — every finding, every
fix, every counter state — is preserved in `histories/` in the program's
main repository.

## The release gate

No paper circulates until it has accumulated **K = 3 consecutive
fully-clean audit passes**. K began at 2; the standing escalation
clause — if any error is ever caught after a pass was declared clean,
by a later pass, a reader, or post-release, K rises for all subsequent
releases — fired on 2026-07-13 (a pass-8 audit found issues present
during a declared-clean pass 7), and the program accepted the
escalation rather than weakening the clause. It remains standing.

**A pass** is a full audit cycle by **three fresh AI agents** working
independently, with no access to any prior audit report (to prevent
anchoring on previously-cleared claims):

1. **Claims/proofs auditor** — reads the manuscript as a hostile
   referee and independently re-derives every theorem, lemma, and
   proposition from scratch, probing specifically for false converses,
   overstated quantifiers, and claims filed above their delivered
   strength (theorem-grade vs computed vs sampled vs asymptotic vs
   conjectural).
2. **Numerics auditor** — re-runs every public script in full and
   traces **every number in the manuscript** (abstract included) to a
   script output line, a re-verified in-text derivation, or a cited
   source, with robustness spot-checks (perturbed tolerances, refined
   grids, changed seeds where seed-robustness is claimed).
3. **Bibliography/prior-art auditor** — verifies every reference
   field-by-field against primary records (for journal-led citations
   the publisher record governs, not the arXiv preprint), checks every
   quotation verbatim against its source, and runs fresh adversarial
   searches attempting to *refute* the paper's novelty claims.
   Crucially, it verifies not only that each reference is *real* but
   that it *supports the specific claim it is attached to*: for every
   citation it reads the cited work and confirms the paper actually
   backs the in-text characterization — catching misattributed
   citations (real authors, wrong or non-existent supporting paper)
   and topical mismatches that a metadata-only check would pass.

**Clean** means zero findings requiring a change to the manuscript or
scripts — any severity, including style nits. Pure no-action
observations do not count. **Any fix resets the counter to zero.**
Factual corrections proposed by an auditor (an author name, a title, a
number attributed to an external source) are independently verified
against the primary source *before* being applied — the corrector can
be as wrong as the original, and this program has caught auditors'
corrections being reversed by stronger primary evidence.

## Cross-model verification

The papers are authored by one model family and audited by agents from
two (Claude Fable and Claude Opus), in varied framings across passes
(checklist verification, fresh-eyes reading, skeptical-referee,
regression). The final confirming pass before release is run at
release grade. In practice the two architectures have independently
re-derived every theorem in the corpus multiple times each.

## The failure-mode checklist

Audits are aimed, not generic. Each pass explicitly probes a checklist
of failure modes grown from errors this program actually caught in its
own work (each entry exists because the generic read-through missed it
at least once):

- false "iff"s and converse claims;
- manuscript claims stronger *or weaker* than what the code proves;
- chimera citations (one entry fusing two real papers), paraphrases
  drifted into quotation marks, preprint titles cited for published
  versions;
- **valid reference, false characterization** — a real, correctly-cited
  source credited with a result/section/assumption it does not contain
  (metadata-passes, content-fails). The biblio lane's claim-support step
  is mandatory and must be *executed*, not assumed: read the cited work
  and confirm it actually says what the in-text sentence attributes to
  it. (A 2026-07 fabrication — "Harrigan--Spekkens Sec. 6 carves out
  nomic wave functions," which they do not — passed a metadata-only
  citation check and was caught only by an argument-level review.)
- **arXiv ID resolving to a different paper than the journal ref** — the
  identifier must resolve to the *same* work as the title/authors/journal;
  a correct journal ref can sit beside a wrong preprint ID that points at
  an unrelated paper;
- claims of exactness where the artifact shows machine epsilon;
- tolerance-masked numerics (ill-conditioned inverses, grid
  coarseness) — key checks re-run under perturbed settings;
- half-domain scans generalized to full-domain claims;
- conclusions written before the data (pre-run-expectation
  contamination) — predictions are pre-registered in the notebook and
  their failures logged;
- incomplete fix propagation — a corrected constant is grepped across
  all scripts, companion papers, and documentation;
- Monte-Carlo "resolved to ±X" claims not backed by a printed
  truth-vs-estimate correlation (a spread alone can be an echo of the
  injected prior);
- convention-relative quantities stated as convention-free;
- prior-art staleness — fast-moving lanes are re-swept.

## What this is and is not

This gate is the program's pre-circulation quality bar. It is not a
substitute for human peer review — it is the reason the papers are fit
to ask for it. Every theorem carries a full proof, every numerical
claim is reproducible from the public scripts, and the audit ledgers
are public precisely so that a human referee can start from the
program's own record of what has been checked, what was found, and
what was fixed. The notebooks also keep an errors-and-misses ledger:
the program records its own failures, including the ones its audits
took multiple rounds to catch.
