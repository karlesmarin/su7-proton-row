#!/usr/bin/env python3
"""check_layout.py - no page of the SU(7) paper carries an internal blank band.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

A tall unbreakable box (the ledger, the verification table, a figure) can push a
page break early and leave a hole in the MIDDLE of a page.  Slack at the FOOT is
what \\raggedbottom is for and is fine; slack in the middle is not.  This measures
the largest gap between consecutive pieces of inked content on every page.

Run:  python check_layout.py     (from part_vi/paper/)
"""
import re
import fitz

PDFS = ["su7_proton_row.pdf", "su7_proton_row_es.pdf"]
LIMIT = 120.0          # pt; a band larger than this is worth looking at
# The page number sits in the footer, below the text block (A4 at margin 2.4cm
# ends at ~774pt).  Counting it as content makes the gap between the last line
# and the folio look like a band, which flags the LAST page of every document --
# a guard that always fires measures nothing.  So the footer is not content.
FOOTER = 780.0
# ...but excluding the footer must not blind the guard to slack ABOVE it. On any
# page but the last, a large gap between the last line and the bottom of the text
# block means an unbreakable box (a table, a figure) did not fit and was pushed
# over, usually leaving a stranded section heading behind. That is a defect and
# the first version of this file could not see it.
TEXT_BOTTOM = 774.0        # A4 at margin 2.4cm
FOOT_LIMIT = 150.0
# ...and the footer filter above had a second, worse blind spot, found by Carles reading the
# PDF on 2026-08-20: it drops EVERYTHING below y = FOOTER, so ink that runs off the BOTTOM of
# the page is dropped too. Figure 1's caption was doing exactly that, in both editions -- the
# last line sat at y = 848.5 on a 841.9 pt page, i.e. printed past the paper edge, and the
# gate said "both editions clean". A float page can overrun without LaTeX warning at all.
# The folio sits at ~785, so anything below BLEED is not a folio: it is spill.
BLEED = 35.0               # pt of true bottom margin that must stay free of ink

ALLBAD = []
for PDF in PDFS:
 doc = fitz.open(PDF)
 print("%s: %d pages" % (PDF, doc.page_count))
 worst = []
 for i, page in enumerate(doc, 1):
     h = page.rect.height
     rows = set()
     for b in page.get_text("blocks"):
         y0, y1 = b[1], b[3]
         for y in range(int(y0), int(y1) + 1):
             rows.add(y)
     for d in page.get_drawings():
         r = d["rect"]
         for y in range(int(r.y0), int(r.y1) + 1):
             rows.add(y)
     for b in page.get_images(full=True):
         for r in page.get_image_rects(b[0]):
             for y in range(int(r.y0), int(r.y1) + 1):
                 rows.add(y)
     # --- spill off the bottom of the sheet, measured BEFORE the footer filter
     if rows and max(rows) > h - BLEED:
         print("  p%-3d ink runs to %6.1f on a %6.1f pt page   SPILLS OFF THE BOTTOM by "
               "%6.1f pt" % (i, max(rows), h, max(rows) - (h - BLEED)))
         ALLBAD.append((PDF, i, "spill %.1f pt past the bottom margin" % (max(rows) - (h - BLEED))))

     rows = {y for y in rows if y < FOOTER}
     if rows and i < doc.page_count:
         slack = TEXT_BOTTOM - max(rows)
         if slack > FOOT_LIMIT:
             # ...but slack at the foot is only a defect if something was STRANDED. If the very
             # next page opens with a full-width float, the slack is just that float not fitting
             # in what was left, which is what \raggedbottom is for and reads correctly. This is
             # measured, not assumed: it asks whether the next page carries an image high up.
             # A \includegraphics of a PDF is a Form XObject, so get_images() is blind to it;
             # what is always there is the caption. A caption in the top half of the next page
             # means that page opens with the float.
             nxt = doc.load_page(i)          # 0-based: the page after this one
             floated = any(re.match(r"\s*(Figure|Figura|Table|Cuadro|Tabla)\s+\d+", b[4])
                           and b[1] < TEXT_BOTTOM * 0.55
                           for b in nxt.get_text("blocks"))
             if floated:
                 print("  p%-3d ink %6.1f..%6.1f   foot slack %6.1f pt, explained: p%d opens "
                       "with a float" % (i, min(rows), max(rows), slack, i + 1))
             else:
                 print("  p%-3d ink %6.1f..%6.1f   FOOT SLACK %6.1f pt  <-- something was "
                       "pushed to the next page" % (i, min(rows), max(rows), slack))
                 worst.append((i, slack, 0, min(rows), max(rows)))
             continue
     if not rows:
         print("  p%-3d EMPTY PAGE" % i)
         worst.append((i, h))
         continue
     ys = sorted(rows)
     gap, at = 0.0, 0
     for a, b in zip(ys, ys[1:]):
         if b - a > gap:
             gap, at = b - a, a
     worst.append((i, gap, at, ys[0], ys[-1]))
     flag = "  <-- internal band" if gap > LIMIT else ""
     print("  p%-3d ink %6.1f..%6.1f   largest internal gap %6.1f pt at y=%d%s"
           % (i, ys[0], ys[-1], gap, at, flag))
 bad = [w for w in worst if len(w) > 2 and w[1] > LIMIT]
 print(" %s: pages flagged (internal band over %.0f pt, or foot slack over %.0f pt): %d"
       % (PDF, LIMIT, FOOT_LIMIT, len(bad)))
 ALLBAD += [(PDF, b) for b in bad]
print()
print("both editions clean" if not ALLBAD else "FLAGGED: %s" % ALLBAD)
raise SystemExit(1 if ALLBAD else 0)
