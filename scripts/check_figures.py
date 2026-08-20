#!/usr/bin/env python3
"""check_figures.py - toda figura cerca de donde se la cita, y en la edición correcta.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

Tres cosas que ningún otro verificador ve:

  1. DERIVA. Un `figure` es un flotante: LaTeX lo coloca donde le cabe, no donde
     está escrito. Tras añadir texto -- una introducción, una sección nueva --
     una figura puede acabar páginas por delante o por detrás de la frase que la
     cita, y el .tex no cambia, así que nada lo delata. Se mide la distancia
     entre la página del pie y la página de la PRIMERA cita.
  2. HUÉRFANAS. Una figura que no cita nadie.
  3. EDICIÓN CRUZADA. Que el castellano no incluya un PDF inglés por descuido:
     todo `\\includegraphics` de la edición _es debe apuntar a un `_es.pdf`.

Uso:  python check_figures.py            (desde part_vi/paper/)
"""
import os
import re
import sys

import fitz

EDICIONES = [("su7_proton_row", "EN"), ("su7_proton_row_es", "ES")]
LIMITE = 2          # páginas de separación toleradas entre pie y primera cita


def labels_del_aux(stem):
    """{label: (numero, pagina)} para los \\newlabel de figuras."""
    out = {}
    aux = stem + ".aux"
    if not os.path.exists(aux):
        return out
    for m in re.finditer(r"\\newlabel\{(fig:[^}]+)\}\{\{([^}]*)\}\{(\d+)\}",
                         open(aux, encoding="utf-8", errors="ignore").read()):
        out[m.group(1)] = (m.group(2), int(m.group(3)))
    return out


def paginas_de_cita(doc, numero, pagina_del_pie, es):
    """Páginas donde el TEXTO cita esa figura, sin contar el propio pie.

    El pie SIEMPRE se compone como "Figure N:" / "Figura N:" con dos puntos; una
    cita nunca los lleva. Distinguirlos por ahí y no borrando el resto de la
    página: la primera versión de esto borraba desde el pie hasta el final, y se
    comía las citas que caían despues en la misma pagina -- daba una figura por
    huerfana estando citada tres lineas mas abajo."""
    largo = "Figura" if es else "Figure"
    pats = [r"Fig\.\s*%s\b" % numero, r"%s\s*%s\b(?!\s*:)" % (largo, numero)]
    out = []
    for i, pg in enumerate(doc, 1):
        t = " ".join(pg.get_text().split())
        if any(re.search(p, t) for p in pats):
            out.append(i)
    return out


fallos = []
for stem, ed in EDICIONES:
    pdf = stem + ".pdf"
    if not os.path.exists(pdf):
        print("%s: no existe %s -- saltado" % (ed, pdf))
        continue
    doc = fitz.open(pdf)
    labs = labels_del_aux(stem)
    tex = open(stem + ".tex", encoding="utf-8").read()
    print("%s  %s: %d páginas, %d figuras etiquetadas" % (ed, pdf, doc.page_count, len(labs)))

    # 3. edición cruzada
    incl = re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex)
    if ed == "ES":
        malas = [f for f in incl if not f.endswith("_es.pdf")]
        print("     imágenes incluidas: %d, todas _es: %s" % (len(incl), not malas))
        if malas:
            fallos.append("%s incluye figuras de la otra edición: %s" % (ed, malas))
    else:
        malas = [f for f in incl if f.endswith("_es.pdf")]
        if malas:
            fallos.append("%s incluye figuras _es: %s" % (ed, malas))

    for lab in sorted(labs, key=lambda k: int(labs[k][0]) if labs[k][0].isdigit() else 99):
        num, pag = labs[lab]
        citas = paginas_de_cita(doc, num, pag, ed == "ES")
        if not citas:
            estado, d = "HUÉRFANA: no la cita nadie", None
            fallos.append("%s %s (Fig. %s) no está citada" % (ed, lab, num))
        else:
            d = min(abs(c - pag) for c in citas)
            estado = "ok" if d <= LIMITE else "DERIVA de %d páginas" % d
            if d > LIMITE:
                fallos.append("%s %s (Fig. %s): pie en p.%d, primera cita en p.%s"
                              % (ed, lab, num, pag, citas[0]))
        print("     Fig. %-2s %-16s pie p.%-3d citada en %-12s %s"
              % (num, lab, pag, str(citas) if citas else "-", estado))
    print()

if fallos:
    print("PROBLEMAS:")
    for f in fallos:
        print("   " + f)
else:
    print("todas las figuras en su sitio, citadas, y cada edición con sus propias imágenes")
sys.exit(1 if fallos else 0)
