#!/usr/bin/env python3
"""check_scripts.py - todo script que el artículo nombra existe y tiene salida archivada.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

La sección de disponibilidad de datos promete que cada número regenera desde los
scripts anexos y que su salida está archivada al lado. Es una promesa verificable
y nada la verificaba.

Falla de tres maneras, todas silenciosas:

  - el artículo nombra un script que se renombró o se borró (hoy mismo pasó:
    `su7_can_the_adjoint_host.py` quedó superado por
    `su7_adjoint_is_vectorlike.py` y hubo que cambiarlo en dos ediciones);
  - un script existe pero su salida archivada no, así que la compuerta de números
    no puede encontrar sus cifras;
  - una edición nombra un script y la otra no.

Uso:  python check_scripts.py     (desde part_vi/paper/)
"""
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RAIZ = ".."
OUT = os.path.join(RAIZ, "outputs")


def nombres(tex):
    """los .py que el artículo cita en \\texttt{}, con los escapes de LaTeX
    deshechos: \\allowbreak fuera y \\_ -> _."""
    s = open(tex, encoding="utf-8").read()
    fuera = set()
    for m in re.findall(r"\\texttt\{([^}]*?\.py)\}", s):
        m = m.replace("\\allowbreak", "").replace("\\_", "_")
        fuera.add("".join(m.split()))
    return fuera


EN, ES = nombres("su7_proton_row.tex"), nombres("su7_proton_row_es.tex")

print("  scripts citados: EN %d, ES %d" % (len(EN), len(ES)))
malos = []
if EN != ES:
    malos.append("solo en EN: %s | solo en ES: %s" % (sorted(EN - ES), sorted(ES - EN)))
    print("  *** las dos ediciones citan conjuntos DISTINTOS ***")
else:
    print("  las dos ediciones citan el mismo conjunto")
print()

print("  %-44s %-10s %s" % ("script", "existe", "salida archivada"))
for n in sorted(EN | ES):
    hay = os.path.exists(os.path.join(RAIZ, n))
    if n.startswith("make_figures"):
        # dibuja PDFs, no imprime números: su control es check_figures.py
        print("  %-44s %-10s %s" % (n, "sí" if hay else "*** NO ***", "(dibuja figuras)"))
        if not hay:
            malos.append("%s: citado y NO EXISTE" % n)
        continue
    base = n[:-3]
    # las salidas de sage se archivan como *_sage.txt
    cand = [base + ".txt", base.replace(".sage", "") + "_sage.txt"]
    sal = any(os.path.exists(os.path.join(OUT, c)) for c in cand)
    if not hay:
        malos.append("%s: citado y NO EXISTE" % n)
    elif not sal:
        malos.append("%s: existe pero no hay salida archivada" % n)
    print("  %-44s %-10s %s" % (n, "sí" if hay else "*** NO ***",
                                "sí" if sal else "*** NO ***"))

print()
if malos:
    print("PROBLEMAS:")
    for m in malos:
        print("   " + m)
else:
    print("cada script citado existe y tiene su salida archivada, en las dos ediciones")
sys.exit(1 if malos else 0)
