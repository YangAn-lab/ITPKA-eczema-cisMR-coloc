#!/usr/bin/env python3
"""Renumber references by order of first appearance (Vancouver style).
Dry-run by default; pass --apply to rewrite the manuscript in place."""
import re, sys, shutil

PATH = "/workspace/P2-02_manuscript_Genomics_submission.md"
APPLY = "--apply" in sys.argv

text = open(PATH, encoding="utf-8").read()
lines = text.split("\n")

# --- locate References section ---
ref_start = next(i for i, l in enumerate(lines) if l.strip() == "## References")
ref_end = next(i for i, l in enumerate(lines) if i > ref_start and l.startswith("## "))
ref_lines = lines[ref_start + 1: ref_end]

# --- parse reference entries ---
entries = {}
for l in ref_lines:
    m = re.match(r"^(\d+)\.\s+(.*\S)\s*$", l)
    if m:
        entries[int(m.group(1))] = m.group(2)
assert len(entries) == 61, f"expected 61 refs, got {len(entries)}"
assert sorted(entries) == list(range(1, 62)), "ref list not 1..61"

# --- citation token regex: [3], [10,11], [5-8], [1, 2, 3], mixed [14,15,16] ---
CIT = re.compile(r"\[(\d+(?:\s*[–-]\s*\d+)?(?:\s*,\s*\d+(?:\s*[–-]\s*\d+)?)*)\]")

def expand(group):
    """'5–8' -> [5,6,7,8]; '10,11' -> [10,11]; '14,15,16' -> [14,15,16]"""
    out = []
    for part in group.split(","):
        part = part.strip()
        m = re.match(r"^(\d+)\s*[–-]\s*(\d+)$", part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            out.extend(range(a, b + 1))
        else:
            out.append(int(part))
    return out

def compress(nums):
    """[1,2,3,5] -> '1–3,5' (runs of >=3 become ranges)"""
    nums = sorted(set(nums))
    runs, cur = [], [nums[0]]
    for n in nums[1:]:
        if n == cur[-1] + 1:
            cur.append(n)
        else:
            runs.append(cur); cur = [n]
    runs.append(cur)
    parts = []
    for r in runs:
        parts.append(f"{r[0]}–{r[-1]}" if len(r) >= 3 else ",".join(map(str, r)))
    return ",".join(parts)

# --- scan body (everything except the reference list) in document order ---
body_idx = [i for i in range(len(lines)) if not (ref_start < i < ref_end)]
order, seen = [], set()
groups_found = 0
for i in body_idx:
    for m in CIT.finditer(lines[i]):
        groups_found += 1
        for n in expand(m.group(1)):
            if n not in seen:
                seen.add(n); order.append(n)

print(f"citation groups in body: {groups_found}")
print(f"distinct refs cited: {len(seen)}")
missing = sorted(set(entries) - seen)
dangling = sorted(seen - set(entries))
print(f"refs never cited: {missing}")
print(f"cited but not in list: {dangling}")
assert not missing and not dangling, "citation/list mismatch — fix before renumbering"

# --- mapping old -> new ---
mapping = {old: new for new, old in enumerate(order, 1)}
print("\nfirst 12 citations in new order:", [mapping[o] for o in order[:12]])
print("old->new sample:", {o: mapping[o] for o in sorted(mapping)[:12]})

# --- rewrite body citations ---
def repl(m):
    return "[" + compress([mapping[n] for n in expand(m.group(1))]) + "]"

new_lines = list(lines)
for i in body_idx:
    new_lines[i] = CIT.sub(repl, lines[i])

# --- rebuild reference list sorted by new number ---
new_ref_lines = [f"{mapping[old]}. {entries[old]}" for old in order]
new_lines[ref_start + 1: ref_end] = new_ref_lines + [l for l in ref_lines if l.strip() == ""][:1]

# --- validation on rebuilt text ---
new_text = "\n".join(new_lines)
recheck = set()
for i, l in enumerate(new_lines):
    if ref_start < i < ref_start + 1 + len(new_ref_lines):
        continue
    for m in CIT.finditer(l):
        recheck.update(expand(m.group(1)))
assert recheck == set(range(1, 62)), f"post-check failed: {sorted(set(range(1,62))-recheck)}"
# order check: first appearances must be 1,2,3,...
seq, seen2 = [], set()
for i, l in enumerate(new_lines):
    if ref_start < i < ref_start + 1 + len(new_ref_lines):
        continue
    for m in CIT.finditer(l):
        for n in expand(m.group(1)):
            if n not in seen2:
                seen2.add(n); seq.append(n)
assert seq == list(range(1, 62)), f"order check failed at {[(k+1,v) for k,v in enumerate(seq) if k+1!=v][:5]}"
print("\nvalidation OK: 61 refs, first-appearance order 1..61, no dangling/uncited")

if APPLY:
    shutil.copy(PATH, PATH + ".bak_prerenumber")
    open(PATH, "w", encoding="utf-8").write(new_text)
    print("APPLIED — manuscript rewritten, backup at .bak_prerenumber")
else:
    print("DRY RUN — use --apply to write")
