"""Baixa do CVRPLIB as instâncias de calibração (A-n33-k5, A-n54-k7, X-n125-k30) para data/raw/.

Se o download automático falhar, baixe manualmente os arquivos .vrp e .sol em
https://galgos.inf.puc-rio.br/cvrplib/ (conjuntos A e X) e salve em data/raw/.
"""
import urllib.request
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
# id de cada instância no CVRPLIB (links "download/instance/<id>" e "download/bks/<id>")
FILES = {"A-n33-k5": 5, "A-n54-k7": 21, "X-n125-k30": None}
PATTERNS = [
    "https://galgos.inf.puc-rio.br/cvrplib/index.php/en/download/{kind}/{id}",
    "https://raw.githubusercontent.com/PyVRP/Instances/main/CVRP/{name}.{ext}",
]


def valid(text: str, ext: str) -> bool:
    if "<html" in text[:500].lower():
        return False
    return ("DIMENSION" in text and "DEMAND_SECTION" in text) if ext == "vrp" else ("Route #1" in text and "Cost" in text)

ok = True
for name, s in FILES.items():
    for ext in ("vrp", "sol"):
        dest = RAW / f"{name}.{ext}"
        if dest.exists() and valid(dest.read_text(errors="replace"), ext):
            continue
        for pat in PATTERNS:
            if "{id}" in pat and s is None:
                continue
            url = pat.format(name=name, ext=ext, id=s, kind="instance" if ext == "vrp" else "bks")
            try:
                data = urllib.request.urlopen(url, timeout=30).read().decode("utf-8", "replace")
            except Exception:
                continue
            if valid(data, ext):
                dest.write_text(data)
                print("baixado:", dest.name, "de", url)
                break
        else:
            ok = False
            print("FALHOU:", dest.name)
print("tudo certo" if ok else "baixe manualmente os arquivos que falharam (ver docstring)")
