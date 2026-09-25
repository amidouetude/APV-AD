"""
Verification ciblee : pourquoi Nairobi (-1.29 deg) est le seul site austral ou
la forme plein sud d'avant donnait PLUS que la forme corrigee.

Hypothese a tester : a quasi-equatorial, une inclinaison de 20 deg est fausse
dans LES DEUX sens (l'optimum radiatif y est ~0 deg), et la convention
"tournee vers l'equateur" y est presque arbitraire. Si c'est bien cela,
l'ecart doit s'annuler quand beta -> 0 et croitre avec beta.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "notebooks"))
sys.path.insert(0, str(REPO / "webapp" / "backend"))
sys.path.insert(0, str(HERE))

from validate_southern_real import poa_with, fetch, SITES  # noqa: E402

name, lat, lon, season, _ = SITES[1]          # Nairobi
df = fetch(lat, lon, name, season)

print(f"{name}  lat = {lat:+.2f} deg")
print(f"{'beta':>6} {'POA corrige':>13} {'POA legacy':>12} {'ecart %':>10}")
print("-" * 45)
for b in [0, 1, 2, 5, 10, 15, 20, 30, 40]:
    c, _ = poa_with(df, float(b), lat)
    l, _ = poa_with(df, float(b), lat, legacy=True)
    print(f"{b:6.0f} {c:13.1f} {l:12.1f} {100*(l-c)/c:+10.2f}")

print()
print("Meme controle a Lusaka (-15.42), pour comparaison :")
name2, lat2, lon2, season2, _ = SITES[2]
df2 = fetch(lat2, lon2, name2, season2)
print(f"{'beta':>6} {'POA corrige':>13} {'POA legacy':>12} {'ecart %':>10}")
print("-" * 45)
for b in [0, 5, 10, 20, 30]:
    c, _ = poa_with(df2, float(b), lat2)
    l, _ = poa_with(df2, float(b), lat2, legacy=True)
    print(f"{b:6.0f} {c:13.1f} {l:12.1f} {100*(l-c)/c:+10.2f}")
