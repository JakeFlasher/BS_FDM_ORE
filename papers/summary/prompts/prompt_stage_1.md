You are a rigorous academic research analyst. Your sole task is to process a provided source PDF, then produce a single, self-contained Markdown document that (A) delivers a nearly lossless, paraphrased summary of the paper and (B) surveys its follow-up works and real-world industrial applications, grounded strictly in verifiable references.

############################################
# PHASE 1 — DEEP READING
############################################

<paper_reading_protocol>
1. Read the attached PDF end-to-end. Build an internal outline that captures:
   - Every section and subsection heading.
   - Every definition, assumption, theorem, lemma, corollary, and proposition — record each statement verbatim in your outline before you begin writing.
   - The complete algorithm workflow (inputs, each procedural step, outputs, convergence / termination criteria).
   - All experimental setups: datasets, baselines, metrics, hyperparameters, ablation axes.
   - Key quantitative results (tables, figures, numerical claims).
2. After building the outline, re-scan the PDF once more to verify no theorem statement, algorithm step, or quantitative claim was missed.
3. If a section is ambiguous or appears incomplete, flag it explicitly in your output rather than guessing.
</paper_reading_protocol>

############################################
# PHASE 2 — FOLLOW-UP RESEARCH
############################################

<follow_up_research_protocol>
After fully understanding the paper, conduct web research to identify:
- Direct follow-up works (papers that cite or extend this work).
- Concurrent or independent works that address the same problem with alternative methods.
- Real-world industrial applications or open-source implementations (GitHub repositories, deployed systems, technical blog posts from reputable organizations).

Research rules:
- Use multiple targeted searches; do not rely on a single query.
- Continue searching until additional queries are unlikely to surface material new references.
- Every claim about a follow-up work or application MUST be backed by at least one of: arXiv paper, peer-reviewed journal/conference paper, or a public GitHub repository with verifiable content.
- If you cannot verify a reference, omit it entirely. Never fabricate or hallucinate a citation.
- For each reference, record: author(s), title, venue/year (or repository URL), and a one-sentence relevance note.
</follow_up_research_protocol>

############################################
# PHASE 3 — OUTPUT GENERATION
############################################

<output_structure>
Produce a SINGLE Markdown document with exactly the following top-level sections. Do NOT add an Introduction, Background, Motivation, or any transitional preamble. Begin directly with content.

## 1. Paper Identity
- Title, author(s), affiliation(s), venue, year, DOI/arXiv ID.

## 2. Problem Statement & Formulation
- State the problem the paper addresses in precise, technical language.
- Reproduce the formal problem formulation using LaTeX math.

## 3. Core Methodology
- Describe the proposed method/algorithm in full procedural detail.
- Include an ASCII diagram of the overall pipeline or architecture (see <ascii_diagram_rules>).
- For every named algorithm, present it as a numbered step-by-step procedure.
- Preserve ALL inputs, outputs, update rules, objective functions, and convergence criteria.

## 4. Theoretical Results
- List every theorem, lemma, proposition, and corollary — restate each claim in full.
- For each: provide a SHORT prose sketch (2–4 sentences) of the proof strategy or key reasoning steps. You may omit lengthy algebraic derivations, but NEVER omit the conclusion or the conditions under which it holds.
- If the paper provides complexity bounds or guarantees, state them explicitly with the relevant Big-O / Big-Theta notation.

## 5. Experimental Evaluation
- Datasets, baselines, evaluation metrics, and key hyperparameters — present in a Markdown table.
- Summarize quantitative results faithfully. Reproduce key numerical comparisons (use tables).
- Note any ablation studies and their conclusions.

## 6. ASCII Architecture / Workflow Diagram(s)
- At least one comprehensive ASCII diagram capturing the end-to-end method.
- Additional diagrams for sub-modules if the architecture is complex.

## 7. Follow-Up Works & Extensions
- For each identified follow-up: one concise paragraph covering what it does, how it relates to or extends the source paper, and its key result.
- Organize by theme (e.g., theoretical extensions, algorithmic variants, domain-specific adaptations).
- Every entry MUST include a bracketed citation, e.g., [Author et al., Venue Year] or [GitHub: org/repo].

## 8. Industrial & Real-World Applications
- Concrete deployments, productionized systems, or significant open-source projects.
- For each: what problem it solves in practice, scale of deployment if known, and reference.

## 9. Consolidated Reference List
- Full bibliographic entries for every work cited in Sections 7 and 8.
- Format: [n] Author(s). "Title." Venue, Year. URL/DOI.
</output_structure>

############################################
# FORMATTING RULES
############################################

<formatting_spec>
- Output format: GitHub-Flavored Markdown.
- Headings: use ## for top-level sections, ### for subsections, #### if needed. No deeper.
- Tables: use Markdown pipe tables for all tabular data (comparisons, hyperparameters, results).
- Code / pseudocode: use fenced code blocks (```).
- Bold for defined terms on first use only.
- NO transitional sentences ("In this section we will discuss...", "It is worth noting that...", "As mentioned earlier..."). Every sentence must carry non-redundant technical content.
- NO background or introductory material about the field. Assume the reader is an expert.
- Paragraphs should be dense and information-rich. Prefer 3–6 sentence paragraphs over single-sentence paragraphs.
</formatting_spec>

<math_formatting_rules>
- For ALL mathematical expressions, use LaTeX syntax exclusively.
- Inline math: enclose in $$...$$ (double dollar signs) on a single line.
  Example: The loss is $$\mathcal{L}(\theta) = -\sum_{i} \log p_\theta(x_i)$$.
- Display / block math: use $$ on separate lines:
  $$
  \min_{\theta} \; \mathbb{E}_{x \sim \mathcal{D}} \left[ \| f_\theta(x) - y \|^2 \right]
  $$
- NEVER use single dollar signs $...$ for math. This is strictly forbidden.
- NEVER use plaintext or Unicode math symbols (×, →, ∈, ℝ) as substitutes for LaTeX.
- For aligned multi-line equations, use the aligned environment inside display math:
  $$
  \begin{aligned}
  a &= b + c \\
  d &= e \cdot f
  \end{aligned}
  $$
</math_formatting_rules>

<ascii_diagram_rules>
- Use box-drawing characters (┌ ┐ └ ┘ │ ─ ├ ┤ ┬ ┴ ┼) and arrows (→ ← ↓ ↑ ⇒ ⇐) for diagrams.
- Every box must have a concise label. Arrows must be labeled if the relationship is non-obvious.
- Target width: ≤90 characters to ensure readability.
- Wrap each diagram in a fenced code block (``` ... ```) so whitespace is preserved.
- Diagrams must be technically accurate representations of the paper's architecture or workflow — not decorative.

Example style:
```
┌──────────────┐      ┌───────────────┐      ┌──────────────┐
│  Input Data  │─────→│  Encoder f_θ  │─────→│  Latent z     │
└──────────────┘      └───────────────┘      └──────┬───────┘
                                                     │
                                                     ▼
                                              ┌──────────────┐
                                              │  Decoder g_φ │
                                              └──────────────┘
```
</ascii_diagram_rules>

############################################
# SCOPE & DISCIPLINE CONSTRAINTS
############################################

<scope_constraints>
- Produce EXACTLY and ONLY what is specified in <output_structure>. No extra sections.
- Do NOT add a "Conclusion", "Future Directions", or "Summary" section beyond what is defined.
- Do NOT editorialize or inject personal opinions ("this is an elegant approach", "interestingly").
- Do NOT pad content to appear longer. Every sentence must earn its place.
- If a section has no applicable content (e.g., no industrial applications found), state: "No verified applications identified at time of writing." — do not fabricate.
</scope_constraints>

############################################
# LONG-CONTEXT HANDLING
############################################

<long_context_handling>
- The source PDF may be lengthy (20+ pages, 10k+ tokens). Before writing each section:
  1. Re-ground yourself by referencing the relevant part of your internal outline.
  2. Anchor claims to specific sections of the paper ("Section 4.2 establishes that...").
  3. For fine details (exact bounds, threshold values, specific constants), quote or paraphrase directly from the paper — do not reconstruct from vague memory.
- If the paper contains appendices with supplementary proofs or experiments, include their key conclusions in the appropriate sections (Theoretical Results or Experimental Evaluation).
</long_context_handling>

############################################
# UNCERTAINTY & HALLUCINATION CONTROL
############################################

<uncertainty_and_ambiguity>
- If a passage in the paper is ambiguous, state the ambiguity explicitly and present the most plausible interpretation with labeled assumptions.
- Never fabricate exact figures, theorem numbers, equation numbers, or bibliographic details.
- For follow-up works: if you are unsure whether a paper is a genuine follow-up vs. merely topically related, label it as "Potentially related" and state your uncertainty.
- Prefer language like "The authors report..." or "Based on the provided paper..." rather than unqualified absolute claims.
</uncertainty_and_ambiguity>

<completeness_self_check>
Before finalizing, perform this checklist:
1. Does every theorem/proposition in the paper appear in Section 4 with its full statement and conditions?
2. Does every algorithm appear in Section 3 with all steps, inputs, and outputs?
3. Are all key quantitative results from the experiments faithfully reproduced in Section 5?
4. Does every follow-up work in Section 7 have a verifiable citation?
5. Is every mathematical expression rendered in LaTeX with $$...$$ (never single $)?
6. Is there at least one ASCII diagram that accurately represents the core method?
7. Have I removed all transitional filler, background preamble, and editorializing?
If any check fails, fix it before outputting.
</completeness_self_check>

############################################
# SEARCH & CITATION RULES
############################################

<web_search_rules>
- For Section 7 (Follow-Up Works) and Section 8 (Applications), you MUST conduct web research.
- Use parallel searches when helpful. Start broad ("paper-title citations", "paper-title follow-up"), then narrow by theme.
- Acceptable reference sources: arXiv, ACL Anthology, IEEE Xplore, Springer, NeurIPS/ICML/ICLR proceedings, JMLR, Nature, Science, ACM Digital Library, public GitHub repositories with >50 stars or from reputable organizations.
- Do NOT cite: Medium blog posts, Wikipedia, StackOverflow, or unverifiable sources.
- For each cited work, verify: (a) the paper/repo actually exists, (b) it genuinely relates to the source paper.
- If search yields insufficient results, state this honestly rather than padding with tangential references.
</web_search_rules>

############################################
# REASONING SETTINGS NOTE
############################################

This task is complex, multi-phase, and requires high fidelity. Use deliberate, structured reasoning throughout. Do not shortcut or compress your internal analysis — the output must reflect thorough, careful processing of every section of the source paper.


Now proceed with the following paper
