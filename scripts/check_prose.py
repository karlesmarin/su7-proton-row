#!/usr/bin/env python3
"""check_prose.py - decima compuerta: una frase que se quedo huerfana al mover un parrafo.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

La edicion castellana imprimia, en Data availability:

    ... un verificador lo impone y ha saltado dos veces. Los artefactos principales son
    Un instrumento interactivo de esta serie corre en karlesmarin.github.io/ghu-explorer ...

El parrafo del explorador se inserto ENTRE la entradilla y su lista, y la entradilla se quedo
colgando delante de el.  Ocho compuertas verdes y una lectura completa no lo vieron: ninguna
lee prosa.  Lo vio la comparacion linea a linea de los dos fuentes.

La firma es estrecha y por eso sirve de compuerta: un verbo copulativo seguido de una palabra
capitalizada SIN punto en medio.  Sobre las dos ediciones da tres aciertos legitimos, todos
nombres propios (`is Sym2`, `is Weinberg`, `es Sym2`), asi que la lista blanca es corta y
explicita -- si crece, hay que mirar por que, no ampliarla a ciegas.

  python check_prose.py             comprueba las dos ediciones
  python check_prose.py --falsify   vuelve a meter el fallo real y comprueba que salta
"""
import re
import sys
import pathlib

import fitz

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent
PDFS = ["su7_proton_row.pdf", "su7_proton_row_es.pdf"]

COPULA = {"son", "es", "are", "is", "fue", "era", "eran", "were", "was"}
# nombres propios que SI van detras de una copula, uno a uno y con su motivo
ALLOWED = {
    "Sym": "su 28 es Sym2(7bar) y el 84 es Sym3(7) -- el nombre del funtor",
    "Weinberg": "the dimension-6 statement is Weinberg and Wilczek-Zee",
    "Re": "every term of the potential is Re Li5(...) -- la parte real, un operador",
}
# la cola tiene que ser {1,}: el fallo real era "son Un instrumento", y con {2,} la compuerta
# no lo veia.  Se descubrio con --falsify, que es exactamente para lo que esta.
PAT = re.compile(r"\b([a-zA-Zaeiouñüà-ÿ]{2,6}) ([A-ZÀ-Ý][a-zà-ÿ]{1,})")


def scan(text):
    text = re.sub(r"-\n", "", text)
    text = re.sub(r"\s+", " ", text)
    out = []
    for m in PAT.finditer(text):
        if m.group(1).lower() not in COPULA:
            continue
        out.append((m.group(2), text[max(0, m.start() - 60):m.end() + 40]))
    return out


def main():
    falsify = "--falsify" in sys.argv
    bad = 0
    for f in PDFS:
        d = fitz.open(HERE / f)
        text = "\n".join(d.load_page(i).get_text() for i in range(d.page_count))
        hits = scan(text)
        unknown = [(w, c) for w, c in hits if w not in ALLOWED]
        P("  %-24s %d copula+mayuscula, %d en la lista blanca, %d sin explicar"
          % (f, len(hits), len(hits) - len(unknown), len(unknown)))
        for w, c in unknown:
            P("       *** %s ***  ...%s" % (w, c))
        bad += len(unknown)

    if falsify:
        P("")
        P("  --falsify: se reinyecta el fallo real que motivo esta compuerta.")
        broken = ("un verificador lo impone y ha saltado dos veces. Los artefactos "
                  "principales son Un instrumento interactivo de esta serie corre en")
        hits = [(w, c) for w, c in scan(broken) if w not in ALLOWED]
        P("     %s" % ("DETECTADO: %s" % hits[0][0] if hits else "*** NO DETECTADO ***"))
        if not hits:
            return 1

    P("")
    if bad:
        P("PROSA HUERFANA: %d sitio(s) donde una frase entra en otra sin punto." % bad)
        return 1
    P("prosa: ninguna entradilla colgando en ninguna de las dos ediciones")
    return 0


if __name__ == "__main__":
    sys.exit(main())
