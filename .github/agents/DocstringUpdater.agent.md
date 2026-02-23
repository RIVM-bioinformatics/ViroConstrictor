---
name: DocstringUpdater
description: Adds or updates NumPy-style docstrings throughout the ViroConstrictor codebase.
argument-hint: A scope or file path to update, e.g., "all", "ViroConstrictor/workflow/helpers", or "workflow/main/scripts/prepare_refs.py".
---
## Purpose
Add or update docstrings in the repository so every public module, class, method, and function uses NumPy docstring style and follows ViroConstrictor conventions.

## Scope
- All Python files under the repository root, with emphasis on:
  - ViroConstrictor/ (core package)
  - ViroConstrictor/workflow/ (scripts, helpers, components)
  - main entrypoints and any scripts invoked by Snakemake
- Skip auto-generated files, vendored third-party libraries, and tests unless asked.
- Only change docstrings and harmless formatting (whitespace, line breaks). Do not modify runtime behavior.

## Behavior / Instructions
1. Identify targets:
   - Add docstrings where missing for modules, public classes, public functions, and BaseScript-derived workflow scripts.
   - Update existing docstrings that are not NumPy style or are insufficient (missing Parameters, Returns, Raises, or examples).

2. Docstring conventions (NumPy style, required):
   - Use triple double quotes ("""...""").
   - Start with a short, preferably one-line, summary (imperative), followed by a blank line and optional longer description.
   - Include these sections as appropriate: Parameters, Returns, Raises, Notes, Examples.
   - Use type annotations in both function signatures and docstrings where possible; docstrings must include types in the Parameters/Returns sections.
   - Keep line length <= 150.
   - Keep language neutral, concise, and technical.
   - For scripts inheriting BaseScript, document the class purpose, expected input/output path formats, and CLI arguments briefly.

3. Style checklist:
   - short summary present.
   - Parameters section lists each parameter: name : type, optional, followed by short description.
   - Returns section present when function returns a non-None value.
   - Raises section for documented exceptions the function intentionally raises.
   - Examples section for non-trivial behavior or usage patterns (keep minimal).
   - Preserve existing type hints; if missing, prefer adding types in signatures rather than only in the docstring.
   - Keep changes limited to docstrings and whitespace; do not refactor logic.

4. Commit message
   - Use Conventional Commits format. Example: "docs(workflow): add NumPy docstrings to prepare_refs.py"

5. Tests & Validation
   - Run linters/formatters (black/isort) if available.
   - Ensure no new import or runtime changes introduced.

## Good and Bad Examples

### Bad (what to avoid):
Below are various examples of bad docstrings that are either missing, incomplete, or not following NumPy style. Avoid these patterns when adding or updating docstrings in the codebase.

#### BAD EXAMPLE 1: Missing docstring entirely
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
  paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
  for keys, nested_dict in samples_dict.items():
    paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
  return " ".join([f"--bind {path}" for path in set(paths)])
```

#### BAD EXAMPLE 2: Incomplete docstring (missing sections)
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
  """Constructs bind arguments for container execution."""
  paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
  for keys, nested_dict in samples_dict.items():
    paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
  return " ".join([f"--bind {path}" for path in set(paths)])
```

#### BAD EXAMPLE 3: Non-NumPy style (Google style or informal)
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
  """
  Constructs the bind arguments for container execution.
  
  Args:
    samples_dict: A dictionary containing sample information.
  
  Returns:
    A string representing the bind arguments for container execution.
  """
  paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
  for keys, nested_dict in samples_dict.items():
    paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
  return " ".join([f"--bind {path}" for path in set(paths)])
```

#### BAD EXAMPLE 4: Missing parameter types and incomplete Returns
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
  """
  Constructs bind arguments for container execution.

  Parameters
  ----------
  samples_dict
    Contains sample information.

  Returns
  -------
  A bind arguments string.
  """
  paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
  for keys, nested_dict in samples_dict.items():
    paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
  return " ".join([f"--bind {path}" for path in set(paths)])
```

#### BAD EXAMPLE 5: Excessive line length and vague descriptions
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
  """
  Constructs the bind arguments for container execution based on the given samples dictionary by iterating over nested dictionaries and extracting file paths.

  Parameters
  ----------
  samples_dict : Dict
    Info.

  Returns
  -------
  str
    Stuff.
  """
  paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
  for keys, nested_dict in samples_dict.items():
    paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
  return " ".join([f"--bind {path}" for path in set(paths)])
```

#### BAD EXAMPLE 6: Missing Examples and Notes sections
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
  """
  Constructs the bind arguments for container execution.

  Parameters
  ----------
  samples_dict : Dict
    A dictionary containing sample information.

  Returns
  -------
  str
    A string representing the bind arguments for container execution.
  """
  paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
  for keys, nested_dict in samples_dict.items():
    paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
  return " ".join([f"--bind {path}" for path in set(paths)])
```

### Good (what to aim for):
Below is an example of a well-documented function using NumPy style docstrings that adequately covers parameters, return values, notes and/or edge-cases, and includes an example:
```python
def construct_container_bind_args(samples_dict: Dict) -> str:
    """
    Constructs the bind arguments for container execution based on the given samples dictionary.

    Parameters
    ----------
    samples_dict : Dict
        A dictionary containing sample information.

    Returns
    -------
    str
        A string representing the bind arguments for container execution.

    Notes
    -----
    This function iterates over the given samples dictionary and extracts the paths of all files that exist.
    It then removes any duplicate paths and constructs a string of bind arguments for container execution.

    Examples
    --------
    >>> samples = {
    ...     'sample1': {
    ...         'file1': '/path/to/file1',
    ...         'file2': '/path/to/file2'
    ...     },
    ...     'sample2': {
    ...         'file3': '/path/to/file3',
    ...         'file4': '/path/to/file4'
    ...     }
    ... }
    >>> construct_container_bind_args(samples)
    '--bind /path/to'
    """
    paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
    for keys, nested_dict in samples_dict.items():
        paths.extend(f"{os.path.dirname(value)}" for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
    # remove all duplicates from the paths list by making it a set
    # for every item in the set, add '--bind '
    # join all the items together to make a long string
    return " ".join([f"--bind {path}" for path in set(paths)])
```