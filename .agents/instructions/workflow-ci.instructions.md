---
name: "Workflow and CI Practices"
description: "Use for Snakemake workflows, container tooling, and GitHub Actions changes where runtime paths and generated outputs must remain portable."
applyTo: ".github/workflows/**/*,containers/**/*.py,container_manager/**/*.py,ViroConstrictor/workflow/**/*.py"
---

- Trace workflow behavior from the caller through configuration to the runtime command; distinguish developer build-time behavior from end-user execution.
- Keep paths portable for installed packages and working-directory changes. Prefer existing directory helpers or configurable environment inputs over repository-root assumptions and hardcoded scan paths.
- For GitHub Actions changes, use documented runner-provided environment variables and outputs. Validate generated JSON, matrix values, and output-file writes with the same shell/runtime assumptions used by the workflow.
- When a workflow supports multiple environments, containers, schedulers, or targets, inspect and validate every supported branch rather than only the first failing case.
- Prefer a focused local reproduction first, then validate the corresponding CI or container path with mocked or hermetic inputs when the real runtime is unavailable.
