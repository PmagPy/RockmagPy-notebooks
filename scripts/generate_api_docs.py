"""Generate markdown API documentation for pmagpy.rockmag using griffe.

Parses the rockmag.py source directly (no pmagpy import needed) and writes
a MyST-compatible markdown file to book/rockmag_api.md. Private functions
(leading underscore) are skipped.

The pmagpy package is located through the active environment's import
machinery, so this picks up the editable PmagPy installation.

Usage:
    ~/miniforge3/envs/ess-jbook/bin/python scripts/generate_api_docs.py
"""
from importlib.util import find_spec
from pathlib import Path

from griffe import load, parse_numpy

root = Path(__file__).resolve().parent.parent

# Locate the pmagpy package without executing it (rockmag.py is parsed
# statically by griffe, so the runtime dependencies are not needed)
spec = find_spec("pmagpy")
if spec is None or spec.origin is None:
    raise RuntimeError("pmagpy is not installed in this environment")
pmagpy_root = Path(spec.origin).parent.parent

module = load("pmagpy.rockmag", search_paths=[pmagpy_root])

lines = [
    "# `pmagpy.rockmag` API Reference\n",
]

if module.docstring:
    lines.append(f"{module.docstring.value}\n")

public_functions = {name: func for name, func in module.functions.items()
                    if not name.startswith("_")}

for name, func in sorted(public_functions.items()):
    # Build signature
    params = []
    for param in func.parameters:
        if param.name == "self":
            continue
        if param.default is not None:
            params.append(f"{param.name}={param.default}")
        else:
            params.append(param.name)
    sig = ", ".join(params)

    lines.append(f"## `{name}`\n")
    lines.append(f"```python\n{name}({sig})\n```\n")

    if not func.docstring:
        lines.append("---\n")
        continue

    # Parse the NumPy-style docstring into structured sections
    sections = parse_numpy(func.docstring)

    for section in sections:
        kind = section.kind.value

        if kind == "text":
            lines.append(f"{section.value}\n")

        elif kind == "parameters":
            lines.append("**Parameters**\n")
            for param in section.value:
                annotation = f" (`{param.annotation}`)" if param.annotation else ""
                desc = param.description.replace("\n", " ")
                lines.append(f"- **{param.name}**{annotation} — {desc}")
            lines.append("")

        elif kind in ("returns", "yields"):
            lines.append(f"**{kind.capitalize()}**\n")
            for ret in section.value:
                annotation = f"`{ret.annotation}` — " if ret.annotation else ""
                desc = ret.description.replace("\n", " ")
                lines.append(f"- {annotation}{desc}")
            lines.append("")

        elif kind == "raises":
            lines.append("**Raises**\n")
            for exc in section.value:
                annotation = f"`{exc.annotation}` — " if exc.annotation else ""
                desc = exc.description.replace("\n", " ")
                lines.append(f"- {annotation}{desc}")
            lines.append("")

        elif kind == "examples":
            lines.append("**Examples**\n")
            for _, example in section.value:
                lines.append(f"```python\n{example}\n```")
            lines.append("")

        elif kind == "admonition":
            # NumPy "Notes" and similar sections arrive as admonitions
            title = section.title or "Notes"
            lines.append(f"**{title}**\n")
            lines.append(f"{section.value.description}\n")

    lines.append("---\n")

# Remove trailing separator
if lines and lines[-1] == "---\n":
    lines.pop()

out_path = root / "book" / "rockmag_api.md"
out_path.write_text("\n".join(lines))
print(f"API docs for {len(public_functions)} public functions "
      f"written to {out_path}")
