import re
import sys
from typing import Optional, Tuple, List


def decode_simple_target(_l: str) -> Optional[Tuple[str, List[str]]]:
    m = re.match(r"^([a-zA-Z0-9_\-]+)\s*:([^=]*)$", _l)
    if not m:
        return None
    _target = m.group(1)
    _deps = [d for d in m.group(2).split() if d and not d.startswith(".")]
    return _target, _deps


def decode_variable_assignment(_l: str) -> Optional[str]:
    """
    Support variables assigned using =, :=, ?=, +=
    :param _l:
    :return:
    """
    m = re.match(r"^([a-zA-Z0-9_\-]+)\s*(=|:=|\?=|\+=)\s*(.*)$", _l)
    if not m:
        return None
    _var = m.group(1)
    return _var


if __name__ == "__main__":
    outfmt = sys.argv[1]
    lines = []
    with sys.stdin as f:
        buff = ""
        for l in f:
            l = l.rstrip("\n\r")
            # Skip empty lines
            if not l.strip():
                continue
            # Deal with multi-line Makes
            if l.endswith("/"):
                buff += l[:-1] + " "
            else:
                buff += l
                lines.append(buff)
                buff = ""
        if buff:
            lines.append(buff)

    target_docs = []
    variable_docs = []
    i = 0
    while i < len(lines):
        line = lines[i]
        simple = decode_simple_target(line)
        variable = decode_variable_assignment(line)
        if simple:
            # print("Detected simple target:", simple[0])
            target, deps = simple
            help_lines = []
            j = i - 1
            while j >= 0 and lines[j].strip().startswith("#"):
                help_lines.insert(0, lines[j].strip("#").strip())
                j -= 1
            if help_lines:
                target_docs.append({"target": target, "deps": deps, "help": "\n".join(help_lines)})
        if variable:
            # print("Detected variable assignment:", variable)
            var = variable
            help_lines = []
            j = i - 1
            while j >= 0 and lines[j].strip().startswith("#"):
                help_lines.insert(0, lines[j].strip("#").strip())
                j -= 1
            if help_lines:
                variable_docs.append({"variable": var, "help": "\n".join(help_lines)})
        i += 1

    with sys.stdout as out:
        if outfmt == "md":
            out.write("# Toplevel Makefile Documentation\n\n")
            out.write(f"## Variables\n\n")

            for doc in variable_docs:
                out.write(f"### `{doc['variable']}`\n\n")
                out.write(f"**Help:** {doc['help']}\n\n")

            out.write(f"## Targets\n\n")
            for doc in target_docs:
                out.write(f"### `{doc['target']}`\n\n")
                out.write(f"**Depends on:** {', '.join(f'`{d}`' for d in doc['deps']) if doc['deps'] else 'None'}\n\n")
                out.write(f"**Help:** {doc['help']}\n\n")
        else:
            out.write("Variables: \n")
            for doc in variable_docs:
                out.write(f"- {doc['variable']}: {doc['help']}\n")

            out.write("\nTargets: \n")
            for doc in target_docs:
                out.write(f"{doc['target']}: {' '.join(doc['deps']) if doc['deps'] else 'None'}\n")
                for l in doc["help"].split("\n"):
                    out.write(f"  {l}\n")
                out.write(f"\n")
