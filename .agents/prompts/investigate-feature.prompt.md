---
name: "investigate-feature"
description: "Assess a ViroConstrictor feature idea for viability and, if viable, scope practical implementation options."
argument-hint: "Describe the feature idea, intended users, constraints, or relevant files."
---

Use the [deep-dive researcher skill](../skills/deep-dive-researcher/SKILL.md) for the investigation. This is a research and planning prompt: do not edit files unless the request explicitly asks for implementation.

Feature idea: ${input:feature:Describe the feature idea, intended users, constraints, or relevant files.}

1. Define the user problem, intended behavior, acceptance criteria, and decision this investigation must support.
2. Trace the nearest existing code paths, configuration, workflow boundaries, documentation, and tests that constrain the idea.
3. Assess viability using concrete evidence: compatibility, packaging and runtime constraints, security, operational cost, maintenance burden, and affected parallel paths.
4. Score viability from 1 to 5 and explain each score. Use 1 for not viable without major redesign and 5 for a strong fit with limited risk. Separate confirmed facts, assumptions, and open questions.
5. If viable, outline one recommended implementation approach and up to two alternatives. For each, identify affected files or components, API and workflow changes, migration or compatibility concerns, and the smallest useful validation plan.
6. If not viable, explain the blocking constraints and suggest the smallest viable reformulation, if one exists.

Return:

- Decision: viable, viable with conditions, or not viable
- Viability score: 1-5 with rationale
- Evidence and assumptions
- Recommended scope and implementation approach, when applicable
- Alternatives and tradeoffs, when useful
- Risks, open questions, and validation plan

Do not invent implementation results, claim that files were changed, or prescribe unit-test changes merely because the feature touches Python. Recommend tests only when the proposed behavior creates a meaningful testable contract or changes existing behavior.
