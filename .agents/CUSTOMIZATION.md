# Agent Customization

This document describes how repository customization is organized. It is separate from `.agents/AGENTS.md`, which is the concise, client-neutral entry point for working on ViroConstrictor.

## Canonical Layout

Keep reusable guidance in `.agents/`:

- `.agents/AGENTS.md` is the canonical project guide for coding work.
- `.agents/skills/<name>/SKILL.md` contains one focused, reusable workflow. Skills use lowercase-hyphenated names and frontmatter with `name` and `description`.
- `.agents/prompts/` contains short, task-specific `*.prompt.md` templates. Prompts are entry points, not duplicate skill bodies.
- `.agents/agents/` contains reusable role descriptions where the client supports a shared agent format; client-specific tool and model metadata belongs in adapters.
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
| Custom agent | A role with a focused instruction set and client-specific tool/model permissions | User-selected or delegated | `.agents/agents/<name>.agent.md`, with client adapters only where required |

Keep the boundaries strict:

- Put stable project rules such as Python style, commit format, testing requirements, and architecture in `AGENTS.md`. A commit rule is an instruction, not a skill or prompt.
- Use a file-based instruction only when a rule should apply conditionally to matching paths. Do not use `applyTo: "**/*"` for a rule that belongs in the always-on guide.
- Use a skill when the agent should decide that a workflow is relevant or when the workflow needs supporting resources. A skill should explain how to do the work, not encode one particular task.
- Use a prompt when the user wants a repeatable manual entry point with task-specific variables, such as an issue number, target path, or requested output format. A prompt should reference a skill or instruction instead of copying its procedure.
- Use a custom agent when the important difference is role, permissions, tools, model, or handoffs. Keep the reusable role description in `.agents/agents/`; keep tool lists and frontmatter in the client adapter.
- Do not create both a prompt and a skill for the same workflow unless the prompt is only a short input wrapper around the skill. If two entries expose the same procedure, keep the skill and remove the prompt.

Keep `AGENTS.md` and skills portable: they must not depend on vendor-specific tool names, slash commands, or model names. VS Code prompt and file-instruction files use the documented frontmatter because that metadata controls discovery. Use relative Markdown links to connect prompts to canonical skills and keep any client adapter to metadata and a pointer only.

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

Client files should be thin discovery adapters:

- Root `AGENTS.md` and `CLAUDE.md` import or point to `.agents/AGENTS.md` for clients that require those names.
- `.github/copilot-instructions.md`, `.github/agents/`, and `.github/prompts/` should contain only client-specific metadata or entry points.
- `.claude/agents/` may contain Claude frontmatter adapters, but reusable procedure belongs in `.agents/skills/`.
- Do not mirror a skill or prompt in multiple client directories. If a client cannot discover `.agents/` directly, add a pointer or symlink where supported and document the exception.

Prefer progressive disclosure: keep skill files focused, move large references and scripts beside them, and load those resources only when the workflow needs them. Do not assume vendor-specific tool names, slash commands, models, or frontmatter in portable content.

When adding or updating customization, verify that the canonical file exists, frontmatter parses where applicable, and every adapter points to the same source. Review bundled scripts and external references before enabling them.
