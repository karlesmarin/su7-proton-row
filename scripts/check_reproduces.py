#!/usr/bin/env python3
"""check_reproduces.py - novena compuerta: los guiones citados VUELVEN A DAR lo archivado.

  Autor: Carles Marin <karlesmarin@gmail.com>  (con Claude, Anthropic, como asistente)

Las otras ocho leen el .tex, el .log y los PDF.  Ninguna ejecuta nada.  `check_numbers.py`
comprueba que todo numero impreso es greppable en `../outputs/`, y `check_scripts.py` que
todo guion citado existe y tiene salida archivada -- pero las dos creen al archivo.  Si el
guion deja de producirlo, las dos siguen verdes y la frase de Data availability ---*every
displayed number regenerates from the ancillary scripts*--- pasa a ser falsa sin que nada
suene.

Eso es exactamente lo que paso, y es el motivo de esta compuerta.  `su7_vacuum.py` se edito
despues de escribirse cuatro guiones que rehusan sus funciones haciendo `exec` de sus
primeras 102 lineas; a partir de esa edicion la linea 102 caia DENTRO de `minimise()`, que
compilaba sin su `return` y devolvia `None` a todo el mundo.  `su7_sixth_row.py` -- la fila
sexta pre-registrada de la seccion 7 -- y `su7_wedge_direction.py` -- la figura 5 --
reventaban.  Los numeros del articulo estaban bien; lo roto era poder rehacerlos.

Cada guion se ejecuta a un fichero NUEVO y se compara con el archivo; el archivo no se
toca nunca -- es contra lo que `check_numbers.py` greppa, y reescribirlo con la corrida de
hoy convertiria la comprobacion en una tautologia.

  python check_reproduces.py             ejecuta todo (lento: ~10 min)
  python check_reproduces.py --falsify   demuestra que la compuerta PUEDE fallar
  python check_reproduces.py NOMBRE ...  solo esos guiones

Lo que NO cubre, dicho en voz alta: los `.sage` (necesitan el contenedor) y los guiones que
el articulo cita sin salida archivada.  Los dos casos se imprimen, uno a uno, con su motivo.
"""
import difflib
import os
import pathlib
import re
import subprocess
import sys
import tempfile

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent          # part_vi/paper
PART = HERE.parent                                       # part_vi
OUT = PART / "outputs"
TEXS = ["su7_proton_row.tex", "su7_proton_row_es.tex"]


def cited():
    """los .py y .sage que el articulo nombra en \\texttt{} -- la lectura de check_scripts.py."""
    out = set()
    for t in TEXS:
        s = (HERE / t).read_text(encoding="utf-8")
        for m in re.findall(r"\\texttt\{([^}]*?\.(?:py|sage))\}", s):
            n = "".join(m.replace("\\allowbreak", "").replace("\\_", "_").split())
            if "/" not in n and not n.startswith("."):
                out.add(n)
    return sorted(out)


def run(name):
    """devuelve (returncode, stdout) ejecutando el guion desde part_vi/, como esta escrito."""
    p = subprocess.run([sys.executable, str(PART / name)], cwd=str(PART),
                       capture_output=True, timeout=3600)
    return p.returncode, p.stdout.decode("utf-8", "replace").replace("\r\n", "\n")


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    falsify = "--falsify" in sys.argv
    names = argv or cited()

    P("guiones citados por el articulo: %d" % len(cited()))
    P("")

    same, diff, crash, skipped = [], [], [], []
    for n in names:
        src = PART / n
        arch = OUT / (os.path.splitext(n)[0] + ".txt")
        if not src.exists():
            skipped.append((n, "el guion no existe"))
            continue
        if n.endswith(".sage"):
            skipped.append((n, "Sage: necesita el contenedor, fuera del alcance de esta compuerta"))
            continue
        if not arch.exists():
            skipped.append((n, "no tiene salida archivada en outputs/"))
            continue
        try:
            rc, got = run(n)
        except subprocess.TimeoutExpired:
            crash.append((n, "TIMEOUT a los 3600 s"))
            P("  %-34s TIMEOUT" % n)
            continue
        want = arch.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
        if rc:
            crash.append((n, "salio con codigo %d" % rc))
            P("  %-34s REVIENTA (codigo %d)" % (n, rc))
            continue
        if got == want:
            same.append(n)
            P("  %-34s reproduce su archivo" % n)
        else:
            diff.append(n)
            P("  %-34s *** DIFIERE DEL ARCHIVO ***" % n)
            for line in list(difflib.unified_diff(want.splitlines(), got.splitlines(),
                                                  "archivo", "corrida de hoy",
                                                  lineterm="", n=1))[:20]:
                P("        %s" % line)

    P("")
    for n, why in skipped:
        P("  no comprobado  %-30s %s" % (n, why))
    if skipped:
        P("")
    P("  reproducen %d | difieren %d | revientan %d | no comprobados %d"
      % (len(same), len(diff), len(crash), len(skipped)))
    for n, why in crash:
        P("     REVIENTA  %-28s %s" % (n, why))

    if falsify:
        P("")
        P("  --falsify: la compuerta tiene que poder fallar, asi que se le miente.")
        victim = same[0] if same else None
        if victim is None:
            P("     nada verde que estropear; falsificacion no concluyente")
            return 1
        arch = OUT / (os.path.splitext(victim)[0] + ".txt")
        want = arch.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
        rc, got = run(victim)
        mutated = want.replace("0", "9", 1)
        P("     %s contra un archivo con UN digito cambiado: %s"
          % (victim, "DETECTADO" if got != mutated else "*** NO DETECTADO ***"))
        if got == mutated:
            return 1

    if diff or crash:
        P("")
        P("REPRODUCIBILIDAD ROTA: %d guion(es) ya no dan lo que el articulo cita."
          % (len(diff) + len(crash)))
        return 1
    P("")
    P("todo guion citado que puede ejecutarse aqui vuelve a dar su salida archivada, byte a byte")
    return 0


if __name__ == "__main__":
    sys.exit(main())
