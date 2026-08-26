# Portable Prompts

Store reusable, task-specific `*.prompt.md` templates in this directory. Prompts are manual entry points: they frame a request, collect variables, and define the desired output. They do not provide always-on project rules or own a multi-step workflow.

Use a skill when the agent should recognize and execute a reusable procedure. Use a prompt only when the user benefits from an explicit command or input template. If a prompt invokes a skill, keep it short and reference the skill instead of copying its steps.

Recommended structure:

```markdown
# Task name

Use `.agents/skills/<skill-name>/SKILL.md` for the procedure.

Input: [task-specific value]
Scope: [optional files or subsystem]
Output: [expected result or report shape]
```

Client adapters may expose these templates through their native prompt feature. The `.github/prompts/` location is optional and is not currently present in this repository:

- `.github/prompts/` for VS Code or GitHub Copilot prompt discovery.
- `.claude/` or another client-specific prompt location when required by that client.

Adapters should reference the canonical file with a relative path or `@` import where supported. Do not maintain independent copies of prompt content. A client adapter may add only the frontmatter and syntax needed for discovery.
