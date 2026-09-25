# Validation australe — données PVGIS réelles

Exécutée **sur le poste utilisateur** (le conteneur d'audit n'a pas accès à
PVGIS). Artefacts : `audit/phase1/southern_validation/`.
Script : `audit/phase1/validate_southern_real.py`.

PVGIS répond en 2,4 s ; Python 3.14.5, pandas 3.0.3, numpy 2.4.4.
Les cinq TMY téléchargées sont conservées avec leurs métadonnées
(`radiation_db`, altitude, horodatage) comme artefacts reproductibles.

---

## 1. POA et productible PV, avant / après correction

Inclinaison fixée à 20°, `d_row` 6 m, `H` 2,5 m — seule la latitude et le
climat changent.

| Site | lat | GHI | Orientation | POA corrigé | POA legacy | Δ POA | Δ PV |
|---|---:|---:|---|---:|---:|---:|---:|
| Ouagadougou | +12,37 | 2 292 | sud | 2 342,1 | 2 342,1 | **0,00 %** | **0,00 %** |
| Nairobi | −1,29 | 2 047 | nord | 1 925,8 | 2 010,7 | **+4,41 %** | +4,18 % |
| Lusaka | −15,42 | 2 198 | nord | 2 288,4 | 1 919,2 | −16,13 % | −15,48 % |
| Windhoek | −22,57 | 2 243 | nord | 2 395,0 | 1 892,4 | −20,99 % | −20,20 % |
| Le Cap | −33,92 | 1 967 | nord | 2 131,0 | 1 632,7 | −23,38 % | −22,43 % |

*(POA en kWh/m²/an, PV en MWh/ha/an. Δ = legacy relatif au corrigé.)*

**Contrôle nord : 0,00 % exactement.** La correction est rigoureusement
neutre là où l'ancienne formule était juste — sur données réelles cette
fois, pas seulement sur le transect à climat constant de la phase 0.

**Hémisphère sud : l'ancien code sous-estimait le POA de 16 à 23 %** sur les
trois sites franchement australs. L'ordre de grandeur mesuré en phase 0
(−8 % à −32 % à climat constant) est confirmé sur climats réels.

---

## 2. L'anomalie de Nairobi — et ce qu'elle révèle

Nairobi est le seul site austral où l'ancienne forme donne **plus** que la
corrigée. Ce n'est pas une erreur : c'est un résultat.

| β | POA corrigé | POA legacy | Δ |
|---:|---:|---:|---:|
| 0 | 2 051,1 | 2 051,1 | **+0,00 %** |
| 5 | 2 035,0 | 2 056,7 | +1,07 % |
| 20 | 1 925,8 | 2 010,7 | +4,41 % |
| 40 | 1 652,7 | 1 810,9 | +9,57 % |

Pour comparaison, Lusaka (−15,42°) : 0,00 % à β = 0, puis −16,13 % à β = 20.

À β = 0 l'orientation n'existe pas, et les deux formes coïncident exactement
— ce qui confirme au passage que rien d'autre ne diffère entre elles.

**Lecture.** À −1,29° de latitude, l'optimum radiatif est β ≈ 0 : une
inclinaison de 20° est mauvaise **dans les deux sens**. Et le climat de
Nairobi est saisonnièrement asymétrique — les mois les plus ensoleillés
(décembre–février) sont ceux où la déclinaison est négative, donc où le
soleil est au sud. La déclinaison moyenne pondérée par l'irradiance y est
légèrement négative, et un panneau très légèrement tourné vers le **sud**
capte marginalement plus, bien que le site soit dans l'hémisphère sud.

Aux inclinaisons optimales de chaque convention, l'écart tombe à **0,27 %**
(2 051,1 contre 2 056,7). Les 4,41 % n'apparaissent qu'à inclinaison
**imposée** à 20° — un mauvais choix pour ce site, et précisément celui que
le plancher d'exploitation force.

**Conséquence de fond.** « Se tourner vers l'équateur » est une heuristique,
pas une loi. Près de l'équateur, et partout où la nébulosité saisonnière est
asymétrique, l'azimut optimal peut être l'autre. L'heuristique code en dur un
choix que la donnée contredit sur certains sites — c'est un argument concret,
mesuré, en faveur de `surface_azimuth` explicite dans une phase ultérieure.
Elle reste néanmoins largement préférable au statu quo : au-delà de ~5° de
latitude australe, elle vaut 16 à 23 points de POA.

---

## 3. Optimum d'inclinaison en climat réel

Ce que la phase 0 avait explicitement marqué « non transférable », faute de
saisonnalité australe.

| Site | lat | β opt. radiatif | 0,45·\|φ\| | Borne DE | Borne active |
|---|---:|---:|---:|---:|---|
| Ouagadougou | +12,37 | 15° | 5,6 | 15 | non |
| Nairobi | −1,29 | **0°** | 0,6 | 15 | **oui** (−3,9 % de POA) |
| Lusaka | −15,42 | 19° | 6,9 | 15 | non |
| Windhoek | −22,57 | 24° | 10,2 | 15 | non |
| Le Cap | −33,92 | 27° | 15,3 | 15 | non |

Le plancher de 15° n'est réellement saturant qu'en zone quasi-équatoriale —
la phase 0, avec un climat unique, en surestimait la portée (elle annonçait
−7,5° à +22,5° de latitude). Sur climats réels, seul Nairobi est concerné
dans cet échantillon, pour une perte de 3,9 % de POA.

**Rectification d'une remarque antérieure.** J'avais relevé que l'optimum
radiatif mesuré (19–27°) dépassait largement la règle 0,45·|φ| du manuscrit
(6,9–15,3°), en laissant entendre que la règle sous-estimait. C'était comparer
deux objectifs différents : l'optimum radiatif maximise le POA seul, tandis
que la règle décrit l'optimum **eLER**, qui arbitre entre production PV et
lumière laissée à la culture et veut donc moins d'inclinaison. Le run DE
ci-dessous le confirme. La règle n'est pas invalidée par ces chiffres.

---

## 4. Run DE complet en hémisphère sud

Lusaka (−15,42°), 50 générations, **256 s**.

| Variable | Valeur | Borne basse | Marge |
|---|---:|---:|---:|
| β | 15,39° | 15,0 | +0,39 |
| d_row | 3,03 m | 3,0 | +0,03 |
| H | 1,50 m | 1,5 | +0,00 |
| V_dig | 4,04 m³ | 2,0 | — |
| HRT | 36,2 j | 20,0 | — |
| f_PV_heat | 0,0500 | 0,05 | +0,000 |

eLER = **1,7110** — LER_crop 0,389 · LER_PV 0,691 · LER_biogas 0,631.
PV 1 504,9 MWh/ha (le plus élevé des sites testés), biogaz 3,79 MWh/ha.

L'optimiseur converge normalement en hémisphère sud et **reproduit le schéma
de convergence aux bornes** : β, `d_row`, `H` et `f_PV_heat` atterrissent tous
sur leur plancher, comme aux quatre sites publiés.

Corollaire méthodologique : puisque le plancher d'inclinaison est saturant,
**ce run ne peut pas tester la règle 0,45·|φ|** — elle prédit 6,9°, en dessous
de la borne. Tester la règle exigerait d'abaisser temporairement le plancher,
ce qui est un choix d'expérience, pas un correctif.

Note : la saison culturale australe utilisée est une **hypothèse** plausible
(semis en début de saison des pluies), pas un calendrier FAO/GIEWS. Elle
n'affecte ni POA ni PV, seulement les volets culture et eau — donc `LER_crop`
et `LER_biogas` ci-dessus sont indicatifs.

---

## 5. Réserves sur `contextvars` — traitées par des tests

Cinq tests ajoutés à `tests/test_concurrency.py`, **127 tests passent** :

| Test | Propriété figée |
|---|---|
| `test_background_thread_carries_its_own_context` | un thread pose sa surcouche lui-même (le motif de `api.py`) |
| `test_outer_context_is_not_visible_to_a_child_thread` | un thread n'hérite **pas** du contexte — contrat explicite |
| `test_no_leak_between_sequential_requests` | aucune référence mutable ne survit à une requête |
| `test_concurrent_asyncio_tasks_are_isolated` | même garantie côté boucle d'événements |
| `test_every_lookup_goes_through_the_registry` | aucun module ne détient un `SITES` autre que le registre |

Le deuxième est le plus important pour la suite : il **documente en le
testant** le fait qu'une tâche de fond doit ouvrir son propre scope. Quand
`JOBS` sera remplacé par une file de travailleurs, ce test échouera si
quelqu'un suppose l'héritage du contexte.

---

## 6. Ce qui reste à faire de votre côté

**Le workflow CI.** `.github/workflows/tests.yml` est protégé en écriture via
les outils distants. Le fichier corrigé vous a été livré ; il faut le copier,
committer, pousser, et **observer le job sur GitHub Actions**. Je ne peux pas
vérifier une exécution GitHub depuis ici — la simulation locale (exit 0) ne
remplace pas cette observation, comme vous l'avez souligné.

**Dette scientifique enregistrée séparément.** `matlab/forcing/poa_irradiance.m`
porte la même hypothèse plein sud. Inerte pour les quatre sites publiés, tous
au nord. Hors périmètre, ne bloque pas la plateforme web.
