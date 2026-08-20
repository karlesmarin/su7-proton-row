#!/usr/bin/env python3
"""check_parity.py - las dos ediciones dicen estructuralmente lo mismo.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

La convención de la casa es traducción COMPLETA, no resumen: arXiv exige el
artículo entero en ambos idiomas. Pero una edición se edita después de la otra
--- hoy mismo, seis cambios y un reencuadre se aplicaron primero en inglés ---
y nada avisa si uno se queda sin traducir. Un párrafo que existe en una edición
y no en la otra no rompe ningún build.

Cuenta los elementos estructurales de cada una y los compara. No compara el
texto: compara que haya el mismo número de secciones, puentes, recuadros,
figuras, teoremas y entradas de bibliografía, y que las etiquetas sean las
mismas --- que es lo que se rompe cuando se olvida traducir un bloque.

Uso:  python check_parity.py     (desde part_vi/paper/)
"""
import re
import sys

PAR = [("secciones", r"\\section\{"), ("puentes", r"\\bridge\{"),
       ("keyeq", r"\\begin\{keyeq\}"), ("figuras", r"\\begin\{figure\}"),
       ("observation", r"\\begin\{observation\}"),
       ("proposition", r"\\begin\{proposition\}"),
       ("paragraph", r"\\paragraph\{"), ("bibitem", r"\\bibitem\{"),
       ("longtable", r"\\begin\{longtable\}"), ("tabular", r"\\begin\{tabular\}")]


def carga(p):
    s = open(p, encoding="utf-8").read()
    return re.sub(r"(?<!\\)%[^\n]*", "", s)


en, es = carga("su7_proton_row.tex"), carga("su7_proton_row_es.tex")

print("  %-14s %-6s %-6s %s" % ("elemento", "EN", "ES", ""))
malos = []
for nombre, pat in PAR:
    a, b = len(re.findall(pat, en)), len(re.findall(pat, es))
    igual = a == b
    if not igual:
        malos.append("%s: EN %d, ES %d" % (nombre, a, b))
    print("  %-14s %-6d %-6d %s" % (nombre, a, b, "" if igual else "*** DESIGUAL ***"))

for nombre, pat in (("label", r"\\label\{([^}]*)\}"), ("cite", r"\\cite\{([^}]*)\}")):
    A, B = set(re.findall(pat, en)), set(re.findall(pat, es))
    if A != B:
        malos.append("%s solo en EN: %s | solo en ES: %s" % (nombre, sorted(A - B), sorted(B - A)))
    print("  %-14s %-6d %-6d %s" % (nombre + "s", len(A), len(B),
                                    "" if A == B else "*** DISTINTAS ***"))


# --- filas de cada longtable, una a una.
#
# Contar los ENTORNOS longtable no basta, y esto lo aprendimos por las malas: el 10 de agosto de
# 2026 el libro quedó con 37 filas en inglés y 38 en castellano --- una fila roja partida en dos al
# arreglar un hueco de página en una sola edición --- y este verificador dio verde, porque el número
# de longtables seguía siendo dos. Lo vio Carles leyendo, que es justo lo que las compuertas existen
# para no tener que hacer. Una fila que existe en una edición y no en la otra tampoco rompe ningún
# build, exactamente como un párrafo sin traducir.
def filas_longtable(s):
    out, pos = [], 0
    while True:
        i = s.find("\\begin{longtable}", pos)
        if i < 0:
            return out
        j = s.find("\\end{longtable}", i)
        cuerpo = s[i:j]
        out.append(len([r for r in cuerpo.split("\\\\") if "&" in r]))
        pos = j + 1


fa, fb = filas_longtable(en), filas_longtable(es)
print()
print("  %-14s %-6s %-6s %s" % ("longtable", "EN", "ES", "filas de cada una"))
for n, (a, b) in enumerate(zip(fa, fb)):
    igual = a == b
    if not igual:
        malos.append("filas de la longtable %d: EN %d, ES %d" % (n, a, b))
    print("  %-14s %-6d %-6d %s" % ("  tabla %d" % n, a, b, "" if igual else "*** DESIGUAL ***"))

print()
if malos:
    print("PROBLEMAS:")
    for m in malos:
        print("   " + m)
else:
    print("paridad estructural: las dos ediciones llevan los mismos elementos")
sys.exit(1 if malos else 0)
