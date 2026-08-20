#!/usr/bin/env python3
"""check_letter_map.py - los punteros de página de la carta, contra el PDF de verdad.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

La carta a Maru vale por una sola cosa: es un mapa. Le dice a una persona ocupada
en qué página está cada cosa que le concierne, para que gaste dos minutos y no dos
horas. Un mapa con las páginas corridas es peor que ningún mapa --- le hace perder
exactamente el tiempo que la carta promete ahorrarle.

Y las páginas se corren solas. Cada párrafo que se añade al artículo puede empujar
una sección entera. Hoy mismo el artículo pasó de 25 a 26 páginas y tres de los
seis punteros caducaron sin que nada avisara.

Esto localiza cada afirmación de la carta por su TEXTO en el PDF compilado --- no
en el fuente --- y devuelve la página donde cae de verdad. Lo que entrega, no lo
que uno cree haber escrito.

Uso:  python check_letter_map.py     (desde part_vi/paper/, tras compilar)
"""
import sys

import fitz

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PDF = "su7_proton_row.pdf"

# Una carta ENVIADA viaja con su PDF adjunto: la copia del destinatario queda
# congelada y su mapa hay que comprobarlo contra ESO, no contra el árbol de
# trabajo, que sigue moviéndose. Sin esto la compuerta grita por páginas que se
# desplazaron después del envío --- pasó el mismo día: promover un recuadro a
# Proposición movió el resultado del 48 de la p.10 a la p.9, y la carta ya
# enviada decía 10, que era correcto para el PDF que él tiene.
#
# Poner a None mientras se redacta; poner el commit el día que se envía.
ENVIADA_EN = "ef0818c"          # commit cuyo PDF se adjuntó, 2026-08-05 23:17 local

# (lo que la carta le promete, texto que lo localiza, sección que la carta cita)
#
# EN EL MISMO ORDEN EN QUE LA CARTA LOS ESCRIBE. La comprobación de abajo empareja
# el n-ésimo puntero de la carta con la n-ésima entrada de aquí. Agrupar por
# sección no basta y pasó en falso una vez: la carta decía "§7, p.14" para el
# erratum cuando el erratum ya estaba en la p.15, y coló porque OTRA cosa de la
# §7 sí estaba en la 14.
PUNTOS = [
    ("eq. (67) demostrada, no muestreada", "is proved: summing", "7"),
    ("eq. (68) desde la tabla del adjunto: el 7 = 4 + 3", "7 = 4 + 3", "2"),
    ("eq. (81) forzada por su eq. (54): m_W = a/2R5", "their eq. (54)", "2"),
    ("Tabla 1 internamente consistente por su eq. (82)", "makes the third column", "7"),
    ("su eq. (67) sostiene tambien la seleccion de la s5", "only w = 3/4 lies inside", "5"),
    # La sonda tiene que sobrevivir a las reescrituras del texto que localiza, o
    # el verificador confunde "lo he reescrito" con "se ha movido de página".
    # Ésta decía "vector-like Dirac pair" y murió al promover el recuadro a
    # Proposición; "vector-like" está en las dos redacciones. Elegir siempre la
    # palabra que sobrevive, no la frase bonita.
    ("su 48 no puede dar materia quiral: es real", "vector-like", "4"),
    ("los dos datos que su sect. 3.2 necesita", "conjugate: read as", "2"),
    ("erratum: los dos terminos con subindice 48 de su eq. (76)", "printed with subscript", "7"),
    ("la sexta fila, y sus dos predicciones", "0.0734", "7"),
]

if ENVIADA_EN:
    import subprocess
    tmp = "_sent_snapshot.pdf"
    ruta = "research/smeft_formalization/part_vi/paper/" + PDF
    with open(tmp, "wb") as fh:
        subprocess.run(["git", "show", "%s:%s" % (ENVIADA_EN, ruta)],
                       stdout=fh, cwd="../../../..", check=True)
    doc = fitz.open(tmp)
    print("CARTA YA ENVIADA: se comprueba contra el PDF de %s, que es el que"
          " tiene el destinatario," % ENVIADA_EN)
    print("no contra el árbol de trabajo. Para volver a redactar, poner ENVIADA_EN = None.")
    print()
else:
    doc = fitz.open(PDF)
paginas = [p.get_text().replace("\n", " ") for p in doc]

# dónde empieza cada sección, para poder decir "S7, p.14" y que las dos sean ciertas
sec = {}
for i, t in enumerate(paginas, 1):
    for n, titulo in ((1, "Introduction"), (2, "What is inherited"), (3, "All six channels"),
                      (4, "A ladder inside"), (5, "The bill, in eighths"),
                      (6, "The selection rule is a half-line"), (7, "The instrument, and where"),
                      (8, "Why the discrepancy does not"), (9, "Scope"),
                      (10, "What this leaves")):
        if titulo in t and n not in sec:
            sec[n] = i

print("%s: %d páginas" % (PDF, doc.page_count))
print()
print("  secciones: " + ", ".join("S%d p.%d" % (n, p) for n, p in sorted(sec.items())))
print()
print("  %-52s %-8s %s" % ("lo que la carta promete", "sección", "página"))
malos = []
for etiqueta, sonda, s in PUNTOS:
    hit = [i for i, t in enumerate(paginas, 1) if sonda in t]
    if not hit:
        malos.append("%s -- NO ENCONTRADO ('%s')" % (etiqueta, sonda))
        print("  %-52s S%-7s *** NO ENCONTRADO ***" % (etiqueta, s))
        continue
    # la sección citada tiene que contener la página
    ini = sec.get(int(s))
    # una sección puede compartir página con la siguiente: el límite es inclusivo
    fin = min([q for n, q in sorted(sec.items()) if q > ini] + [doc.page_count])
    dentro = [i for i in hit if ini <= i <= fin]
    if not dentro:
        malos.append("%s -- cae en p.%s, fuera de S%s (p.%d-%d)" % (etiqueta, hit, s, ini, fin))
    print("  %-52s S%-7s p.%-8s %s"
          % (etiqueta, s, ",".join(str(i) for i in (dentro or hit)),
             "" if dentro else "*** FUERA DE LA SECCIÓN CITADA ***"))


# ---------------------------------------------------------------------------
# Y ahora lo que la primera versión de este fichero NO comprobaba, y por eso
# pasó en verde mientras la carta traía páginas corridas: los NÚMEROS QUE LA
# CARTA IMPRIME. Pertenecer a la sección correcta no basta; el lector va a la
# página que lee.
LETTER = "../../correspondence/email_komori_maru_reply2_DRAFT.md"
try:
    txt = open(LETTER, encoding="utf-8").read()
except OSError:
    print()
    print("  (aviso: no encuentro %s -- no se comprueban sus números)" % LETTER)
    txt = None

if txt is not None:
    import re as _re
    print()
    print("  los números que la carta imprime, contra el PDF:")
    escritos = _re.findall(r"\(§(\d+),\s*p\.\s*(\d+)\)", txt)
    if len(escritos) != len(PUNTOS):
        malos.append("la carta trae %d punteros y este verificador conoce %d: el "
                     "emparejamiento por orden no es fiable hasta arreglarlo"
                     % (len(escritos), len(PUNTOS)))
        print("     *** la carta trae %d punteros, aquí hay %d entradas ***"
              % (len(escritos), len(PUNTOS)))
    for (s, pg), (etiqueta, sonda, ssec) in zip(escritos, PUNTOS):
        real = [i for i, t in enumerate(paginas, 1) if sonda in t]
        ok = s == ssec and int(pg) in real
        if not ok:
            malos.append("«%s»: la carta dice §%s p.%s, está en §%s p.%s"
                         % (etiqueta, s, pg, ssec, real))
        print("     §%-3s p.%-4s %-46s %s"
              % (s, pg, etiqueta[:46], "ok" if ok else "*** CADUCADO ***"))
    if not escritos:
        malos.append("no encuentro ningún puntero '(§N, p.M)' en la carta")

print()
if malos:
    print("PROBLEMAS:")
    for m in malos:
        print("   " + m)
else:
    print("cada punto de la carta cae en la sección que la carta dice, y estas son sus páginas")
sys.exit(1 if malos else 0)
