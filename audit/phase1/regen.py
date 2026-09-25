"""Regenere simulation_functions.py depuis le notebook, en reproduisant le
mecanisme de Jupyter (linecache) dont inspect.getsource a besoin."""
import json, linecache, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent / "repo"
NB = REPO / "notebooks" / "02_simulation.ipynb"
DEF_CELL_INDICES = [2, 3, 4, 5, 6]
EXPORT_CELL_INDEX = 8

def regenerate(out_path: Path) -> str:
    nb = json.loads(NB.read_text(encoding="utf-8"))
    ns = {}
    exec("import numpy as np", ns); exec("import pandas as pd", ns)
    sys.path.insert(0, str(REPO / "notebooks"))
    exec("from config_00 import SITES, PARAMS, DE_SETTINGS, SCENARIOS", ns)
    for idx in DEF_CELL_INDICES:
        src = "".join(nb["cells"][idx]["source"])
        fname = f"<ipython-input-{idx}>"
        linecache.cache[fname] = (len(src), None, src.splitlines(keepends=True), fname)
        exec(compile(src, fname, "exec"), ns)
    export_src = "".join(nb["cells"][EXPORT_CELL_INDEX]["source"]).replace(
        'out_path = Path.cwd() / "simulation_functions.py"',
        f'out_path = Path(r"{out_path}")')
    exec(export_src, ns)
    return out_path.read_text(encoding="utf-8")
