#!/usr/bin/env python3
"""check_refs.py - every \\ref has a \\label, every \\cite has a \\bibitem, and nothing is orphaned.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Run:  python check_refs.py     (from part_vi/paper/)
"""
import re

TEX = "su7_proton_row.tex"

s = open(TEX, encoding="utf-8").read()
s = re.sub(r"(?<!\\)%[^\n]*", "", s)

labels = set(re.findall(r"\\label\{([^}]*)\}", s))
refs = set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", s))
bib = set(re.findall(r"\\bibitem\{([^}]*)\}", s))
cited = set()
for g in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", s):
    cited |= {x.strip() for x in g.split(",")}

print("%s: labels %d, refs %d, bibitems %d, keys cited %d"
      % (TEX, len(labels), len(refs), len(bib), len(cited)))
bad = []
for name, s_ in (("refs with no label", refs - labels),
                 ("labels never referenced", labels - refs),
                 ("cites with no bibitem", cited - bib),
                 ("bibitems never cited", bib - cited)):
    print("  %-24s : %s" % (name, ", ".join(sorted(s_)) if s_ else "none"))
    if s_ and name != "labels never referenced":
        bad += sorted(s_)
raise SystemExit(1 if bad else 0)
