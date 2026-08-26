---
name: deep-dive-researcher
description: Investigate repository behavior, issues, pull requests, reviews, and failures to find root causes and produce evidence-based conclusions. Use when a task asks for a deep investigation or implementation plan.
---

# Deep Dive Researcher

Investigate one concrete question deeply enough to reach a decision. Read `AGENTS.md` and the relevant implementation, tests, history, or review artifacts.

## Procedure

1. Define the question, expected behavior, and the decision the investigation must support.
2. Find the nearest code path that computes, mutates, or controls the behavior; do not stop at a wiring layer.
3. Gather the smallest useful evidence set from implementation, callers, tests, configuration, documentation, and version history.
4. Compare the observed behavior with the intended behavior and distinguish confirmed facts from assumptions.
5. Identify edge cases, security implications, compatibility constraints, and affected parallel paths.
6. Choose a cheap discriminating check, such as a focused test, a reproducible command, or a direct call-site comparison.
7. Run the check when execution is requested and use its result to confirm or revise the hypothesis.
8. Produce one conclusion and a concrete next action. Do not repeat an unchanged investigation without new evidence.

## Output

Use this structure when useful:

```markdown
# Investigation

## Conclusion
[Direct answer and confidence]

## Evidence
- [Implementation or test evidence]
- [Relevant constraint or compatibility fact]

## Risks and edge cases
- [Only material risks]

## Recommended action
[Smallest justified fix, follow-up, or reason no change is needed]

## Validation
[Check run, result, or why execution was unavailable]
```

For review or quality findings, classify each finding as valid, invalid, or requiring clarification, explain whether it conflicts with intended design, and state whether to fix, suppress, or document it. Avoid presenting a plan as completed work.
