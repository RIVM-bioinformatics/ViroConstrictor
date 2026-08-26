# Agent Customization

This document describes how repository customization is organized. It is separate from `.agents/AGENTS.md`, which is the concise, client-neutral entry point for working on ViroConstrictor.

## Canonical Layout

Keep reusable guidance in `.agents/`:

- `.agents/AGENTS.md` is the canonical project guide for coding work.
- `.agents/skills/<name>/SKILL.md` contains one focused, reusable workflow. Skills use lowercase-hyphenated names and frontmatter with `name` and `description`.
- `.agents/prompts/` contains short, task-specific `*.prompt.md` templates. Prompts are entry points, not duplicate skill bodies.
- `.agents/agents/` contains Antigravity-native custom-agent adapters discovered from the workspace; portable role procedures belong in `.agents/skills/`.
- `.agents/instructions/` contains optional `*.instructions.md` conditional rules, such as rules for a particular language, directory, or file type.
- `.agents/references/` contains detailed material loaded only when a skill needs it.
- `.agents/scripts/` contains deterministic helpers. Audit scripts before using them and prefer treating them as black boxes.

## Choosing A Customization Type

Use the smallest primitive that matches the job:

| Primitive | Purpose | Activation | Put here |
| --- | --- | --- | --- |
| Always-on instructions | Project facts, architecture, coding standards, and policies that apply to every task | Automatic | `.agents/AGENTS.md` |
| File-based instructions | Rules for a specific language, directory, or file type | Automatic when a path matches | `.agents/instructions/<topic>.instructions.md`, with client adapters only where required |
| Skill | A reusable, multi-step procedure with optional references, scripts, or templates | Agent-selected or explicitly invoked | `.agents/skills/<name>/SKILL.md` |
| Prompt | A manually invoked task template that supplies inputs and expected output | User-selected, usually a slash command | `.agents/prompts/<name>.prompt.md` |
| Custom agent | A role with a focused instruction set and client-specific tool/model permissions | User-selected or delegated | `.agents/agents/<name>.agent.md` for Antigravity, with client adapters where required |

Keep the boundaries strict:

- Put stable project rules such as Python style, commit format, testing requirements, and architecture in `AGENTS.md`. A commit rule is an instruction, not a skill or prompt.
- Use a file-based instruction only when a rule should apply conditionally to matching paths. Do not use `applyTo: "**/*"` for a rule that belongs in the always-on guide.
- Use a skill when the agent should decide that a workflow is relevant or when the workflow needs supporting resources. A skill should explain how to do the work, not encode one particular task.
- Use a prompt when the user wants a repeatable manual entry point with task-specific variables, such as an issue number, target path, or requested output format. A prompt should reference a skill or instruction instead of copying its procedure.
- Use a custom agent when the important difference is role, permissions, tools, model, or handoffs. Keep the reusable procedure in `.agents/skills/`; keep each client's agent metadata in its native adapter.
- Do not create both a prompt and a skill for the same workflow unless the prompt is only a short input wrapper around the skill. If two entries expose the same procedure, keep the skill and remove the prompt.

Keep `AGENTS.md` and skills portable: they must not depend on vendor-specific tool names, slash commands, or model names. VS Code prompt and file-instruction files use the documented frontmatter because that metadata controls discovery. Use relative Markdown links to connect prompts to canonical skills and keep any client adapter to metadata and a pointer only.

Tool names are client-specific. Use Antigravity tool names only in `.agents/agents/`, and use the native VS Code or Claude tool names in their respective adapters. Do not copy a tool list between clients without verifying that client's tool schema.

VS Code does not discover `.agents/instructions/` or `.agents/prompts/` by default. To use these canonical files as workspace customizations, configure `chat.instructionsFilesLocations` and `chat.promptFilesLocations` to point to those two directories in the user's VS Code settings. Keep that setting local when `.vscode/` is ignored.

## Prompt Pattern

```markdown
# Investigate a repository question

Use `.agents/skills/deep-dive-researcher/SKILL.md` for the procedure.

Question: [user-provided question]
Scope: [files, issue, pull request, or subsystem]
Required output: [decision, evidence, risks, and validation]
```

The prompt frames the request; the skill owns the investigation method; `AGENTS.md` supplies the project rules.

## Client Adapters

- Root `AGENTS.md`, `CLAUDE.md`, and `.github/copilot-instructions.md` symlink to `.agents/AGENTS.md` so that VS Code Copilot, Claude Code, Antigravity, and other clients receive the canonical project guide directly without content duplication or unexpanded `@` directives.
- Agent adapters use client-identifying prefixes so users know which harness they belong to:
  - `.github/agents/GH-<Name>.agent.md`: GitHub Copilot adapters in VS Code with native Copilot tools (`[execute, read, edit, search, todo]`).
  - `.claude/agents/Claude-<Name>.md`: Claude Code adapters with native Claude tools (`Read, Grep, Glob, Bash, Edit`).
  - `.agents/agents/AGY-<Name>.agent.md`: Antigravity adapters with native Antigravity tools (`view_file, grep_search, replace_file_content, run_command`).
- `.vscode/settings.json` points `chat.instructionsFilesLocations` and `chat.promptFilesLocations` to `.agents/instructions` and `.agents/prompts` for native VS Code Copilot discovery.
- All agent adapters remain thin wrappers pointing to the reusable functional procedure in `.agents/skills/<name>/SKILL.md`.

Prefer progressive disclosure: keep skill files focused, move large references and scripts beside them, and load those resources only when the workflow needs them. Do not assume vendor-specific tool names, slash commands, models, or frontmatter in portable content.

When adding or updating customization, verify that the canonical file exists, frontmatter parses where applicable, and every adapter points to the same source. Review bundled scripts and external references before enabling them.
