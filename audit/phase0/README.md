# Phase 0 — Rapport d'audit (lecture seule)

Date : 2026-09-25. **Aucun fichier existant du dépôt n'a été modifié.**
L'audit a tourné sur une copie, et la correction n'a été appliquée qu'en
mémoire (monkeypatch) pour mesurer les écarts.

Porte de décision à l'issue de cette phase : arbitrer entre **correctif
minimal** et **`surface_azimuth` explicite**, sur la base des mesures
ci-dessous et non d'une préférence a priori.

---

## Synthèse

| # | Question | Réponse mesurée |
|---|---|---|
| A | L'hémisphère Sud est-il faux ? | **Oui.** POA annuel −8,1 % à −32,0 % ; productible PV jusqu'à −30,1 % |
| A | L'hémisphère Nord régresse-t-il ? | **Non.** Écart 0,000 % sur toute la plage +2,5° à +37,5° |
| D | Combien d'endroits portent l'hypothèse ? | **Un seul** : `simulation_functions.py` lignes 201-202 |
| A3 | La borne d'inclinaison est-elle contraignante ? | **Oui**, de −7,5° à +22,5° de latitude ; pénalité POA ≤ 3,9 % |
| B | La concurrence est-elle défaillante ? | **Oui**, mais surtout par plantage, pas par corruption silencieuse |
| C | Les références sont-elles reproductibles ? | **Oui.** Écart max 0,062 % au CSV publié (arrondi du CSV) |

Un bug a été trouvé **dans le correctif que j'avais proposé** — voir §5.

---

## 1. Fiabilité de l'oracle (A1.1)

Avant de mesurer quoi que ce soit, il faut un étalon. Mon implémentation de
la forme complète de Duffie & Beckman (éq. 1.6.5, azimut `gamma` explicite) a
été confrontée à `pvlib.irradiance.aoi` sur 3 715 tirages aléatoires
(latitude, inclinaison, azimut, déclinaison, angle horaire), alimentés par la
**même** position solaire analytique pour isoler la géométrie du modèle de
position solaire.

```
écart absolu max     : 3,386e-14
écart absolu médian  : 1,110e-16
```

C'est la précision machine. La colonne « corrigé » de tous les tableaux qui
suivent est donc un étalon fiable, pas une seconde opinion.

Piège traité au passage : Duffie & Beckman compte `gamma = 0` au **sud**,
pvlib compte l'azimut à partir du **nord** dans le sens horaire. La conversion
`azimut_pvlib = gamma_DB + 180` est faite explicitement.

---

## 2. Ampleur de l'erreur d'hémisphère (A1, A2)

### Protocole

PVGIS n'est pas joignable depuis le conteneur d'audit (proxy, HTTP 403). Cette
contrainte a imposé — et permis — un protocole plus propre : **un seul climat
réel** (TMY de Konya, 8 760 h) appliqué à des latitudes variables. Tout écart
observé est donc imputable à la géométrie, et à rien d'autre. Un jeu austral
réel aurait mélangé climat et géométrie.

*Limite assumée* : les valeurs absolues de POA ne sont pas représentatives des
climats australs. Seuls les **écarts relatifs à latitude donnée** le sont.

### Géométrie pure, intégrée sur l'année (A1.3)

| Latitude | Écart |
|---|---|
| +2,5° à +37,5° | **0,00 %** |
| −5° | −0,62 % |
| −15° | −5,55 % |
| −25° | −15,20 % |
| −35° | **−29,00 %** |

### Vrai pipeline du dépôt, TMY réelle (A2)

Design fixe β = 20°, `d_row` = 6 m, `H_m` = 2,5 m ; seule la latitude varie.

| Latitude | POA actuel | POA corrigé | Écart POA | Écart PV |
|---|---|---|---|---|
| +37,5° … +2,5° | — | — | **0,000 %** | 0,00 % |
| −5,0° | 1 749,5 | 1 944,1 | −10,01 % | −9,00 % |
| −15,0° | 1 610,7 | 1 952,8 | −17,52 % | −16,09 % |
| −25,0° | 1 445,9 | 1 922,9 | −24,80 % | −23,06 % |
| −35,0° | 1 261,3 | 1 855,1 | **−32,01 %** | **−30,08 %** |

### Avec la borne d'inclinaison réellement imposée (β ≥ 15°)

C'est le chiffre décisionnel, puisque l'optimiseur ne descend pas sous 15° :

| Latitude | β imposé | Écart POA |
|---|---|---|
| −5° | 15,0° | −7,50 % |
| −10° | 15,0° | −10,40 % |
| −20° | 15,0° | −16,22 % |
| −30° | 15,0° | −22,03 % |
| −35° | 15,8° | −26,08 % |

### Ce qui n'est *pas* affecté

`F_shad` est **identique au bit près** dans les deux hémisphères, à toutes les
latitudes testées. `compute_shading` n'utilise que l'élévation solaire ; la
ligne `delta_gamma = 0.0` est inerte. C'est mesuré, pas supposé.

---

## 3. Carte des hypothèses d'orientation (D)

Balayage de `notebooks/`, `webapp/backend/`, `matlab/`.

| Hypothèse | Occurrences | Fichiers |
|---|---|---|
| Forme plein sud **codée en dur** | 2 lignes | `simulation_functions.py` **uniquement** (l. 201-202) |
| `delta_gamma = 0` (azimut de rangée) | 1 | `simulation_functions.py` l. 131 — **inerte**, cf. §2 |
| Plein sud mentionné en commentaire | 8 | `simulation_functions.py` |
| `min_beta_deg` | 7 | `config_00.py`, `simulation_functions.py`, `api.py`, `model_bridge.py` |
| `surface_azimuth` | **0** | le concept n'existe nulle part dans le dépôt |

**Ceci tranche la porte de décision telle qu'elle a été formulée.** Le critère
posé était : « si l'orientation plein sud n'est utilisée que dans
`compute_POA`, le correctif minimal peut être acceptable ». La mesure dit :
**une seule fonction, deux lignes**. Le second critère — configurations
est-ouest — reste un choix produit, pas une contrainte technique constatée.

---

## 4. Bornes d'inclinaison (A3)

Séparation conceptuelle demandée — optimum radiatif / auto-nettoyage /
structurel / exploitation — ici seul l'**optimum radiatif** est mesuré.

- La borne basse du DE (15°) est **active de −7,5° à +22,5°** de latitude,
  soit 13 des 30 points testés — l'essentiel de la ceinture intertropicale.
- Pénalité de POA imputable à la borne : **−3,89 % au maximum**. Réelle, mais
  modeste : la courbe POA(β) est plate près de son optimum.
- Les 4 sites publiés sont tous **0,5° à 1,4° au-dessus de leur borne
  effective** :

| Site | Latitude | β optimisé | Borne effective | Marge |
|---|---|---|---|---|
| Konya | 37,87 | 20,71 | 20,0 | +0,71 |
| Almería | 36,83 | 20,51 | 20,0 | +0,51 |
| Ouagadougou | 12,37 | 16,37 | 15,0 | +1,37 |
| Freiburg | 47,99 | 23,08 | 22,0 | +1,08 |

La « boundary-convergence finding » du manuscrit décrit donc le fait que **la
borne est saturante**, pas une propriété de l'optimum agrivoltaïque. À
requalifier dans le texte.

**Réserve importante** : les optima radiatifs par latitude sont calculés avec
le climat de Konya, dont la saisonnalité ne correspond pas aux latitudes
australes. Ces valeurs sont **indicatives seulement** et ne doivent pas être
citées. Les trois constats ci-dessus (plage d'activité de la borne, ordre de
grandeur de la pénalité, saturation des 4 sites) sont robustes ; les valeurs
de β optimal ne le sont pas.

---

## 5. Un bug dans le correctif que j'avais proposé

J'avais proposé `phi_eff = lat - np.sign(lat) * beta`. Mesure :

| Latitude | `phi_eff` obtenu | Attendu | |
|---|---|---|---|
| −0,001° | +20,00° | +20,00° | ok |
| **0,000°** | **0,00°** | **−20,00°** | **BUG** |
| +0,001° | −20,00° | −20,00° | ok |

`np.sign(0) = 0`, donc à la latitude exactement nulle l'inclinaison disparaît
du terme direct : le panneau est traité comme horizontal. Le symptôme se voit
aussi dans A2 (écart de −6,98 % à lat 0, entouré de 0,00 % au nord) et dans A3
(optimum radiatif 0° à lat 0).

Correction : `np.where(lat >= 0, 1.0, -1.0)` au lieu de `np.sign(lat)`.

Ce n'est pas un argument décisif contre le correctif minimal — c'est la
démonstration qu'**une ligne appelle un test autant qu'un refactor**.

---

## 6. Concurrence (B)

### Ce que fait le code

`model_bridge.evaluate_point()` écrit dans le dictionnaire global
`config_00.SITES`, puis le retire dans un `finally`. `api.py` forge la clé
`f"_runtime_{lat}_{lon}"` — **elle ne dépend que des coordonnées**, alors que le
`site_dict` dépend aussi de la culture et du site d'emprunt.

### Résultats

| Scénario | Threads | Corrects | Plantages | Contaminés |
|---|---|---|---|---|
| B1 — mêmes coordonnées, cultures différentes | 12 | 0 | 11 `KeyError` | **1** |
| B2 — mêmes coordonnées, paramètres identiques | 10 | 1 | 9 `KeyError` | 0 |
| B3 — coordonnées différentes (témoin) | 8 | **8** | 0 | 0 |

### Lecture

1. **Le mode dominant est le plantage**, pas la corruption : un thread retire
   la clé pendant qu'un autre calcule → `KeyError` → HTTP 500. Je corrige ici
   ma propre formulation antérieure, qui annonçait des « mauvais résultats
   silencieux » comme mode principal.
2. **La corruption silencieuse existe néanmoins** : 1 résultat sur 12 en B1 a
   renvoyé les paramètres d'une autre requête, sans erreur. C'est le cas le
   plus grave, et il est minoritaire mais réel.
3. **Le protocole de test initialement proposé aurait manqué le bug.**
   « 50 requêtes concurrentes sur Bamako, Lusaka, Nairobi, Dakar, Windhoek » →
   cinq coordonnées distinctes → cinq clés distinctes → c'est exactement le
   témoin B3, qui **passe**. Le test doit **répéter les mêmes coordonnées**, ce
   qui est aussi le cas réaliste : deux utilisateurs sur la même ville, ou un
   double-clic.

---

## 7. Références gelées (C)

Le chemin de code actuel a été rejoué pour les 4 sites publiés, scénario S0,
avec leurs vecteurs de conception optimisés.

- **204 champs** figés dans `c_reference_frozen.json` (résultats physiques +
  eLER + 4E économique).
- Écart maximal au `results_summary.csv` publié : **0,0616 %**, entièrement
  attribuable aux arrondis du CSV. eLER, LER_crop, LER_PV, LER_biogas et
  biogaz : **0,0000 %**.
- Ouagadougou reproduit **eLER 1,723**.

**Tolérance de non-régression retenue : identité stricte (1e-9).** Les quatre
sites sont dans l'hémisphère nord, où la correction doit être rigoureusement
neutre — tout écart non nul après correction signalerait une erreur.

---

## 8. Ce que l'audit n'a pas pu établir

- **Aucune donnée australe réelle.** PVGIS bloqué par le proxy du conteneur.
  Les écarts mesurés sont géométriquement exacts mais n'intègrent pas la
  saisonnalité réelle des climats australs. À refaire depuis un poste ayant
  accès à PVGIS avant de publier des chiffres absolus.
- **Optima d'inclinaison australs** : non transférables, cf. §4.
- **Qualité de la RH du TMY** (prérequis Penman-Monteith) : non testée,
  faute d'accès aux données et à des stations de référence.
- **Comportement de l'optimiseur DE en hémisphère Sud** : non testé
  (un run coûte ~200 s ; à faire en phase 1 après correction).

---

## 9. Éléments pour la porte de décision

| Critère | Correctif minimal | `surface_azimuth` explicite |
|---|---|---|
| Surface de code touchée | 2 lignes, 1 fonction (mesuré) | signature + appelants + config |
| Cas de bord | `np.sign(0)` — trouvé, corrigeable | aucun |
| Hypothèse visible | non : cachée dans un `sign()` | oui : paramètre nommé |
| Rangées est-ouest | impossible | possible |
| Comparaison à `pvlib` | indirecte | directe (même signature) |
| Risque de second refactor | réel | écarté |
| Non-régression nord | identité stricte attendue | identité stricte attendue |

Le critère qui avait été posé comme décisif — « l'orientation est-elle dans
plusieurs couches ? » — est tranché : **non, une seule**. Le choix se reporte
donc sur un arbitrage produit : l'outil doit-il, à terme, représenter des
rangées est-ouest ?

---

## 10. Fichiers

| Fichier | Contenu |
|---|---|
| `a1_aoi_analytical.py` | AOI analytique, contrôle pvlib, erreur par latitude |
| `a2_poa_real_code.py` | POA/PV sur le vrai pipeline ; optima et bornes d'inclinaison |
| `b_concurrency.py` | Les trois scénarios de concurrence |
| `c_reference_freeze.py` | Gel des références + carte des hypothèses |
| `c_reference_frozen.json` | **La référence de non-régression (204 champs)** |
| `a1_*.csv`, `a2_*.csv`, `d_assumption_map.csv` | Résultats détaillés |
| `b_concurrency_results.json` | Résultats de concurrence |

Reproduire : `python3 a1_aoi_analytical.py` puis `a2`, `b`, `c`.
Prérequis : `pandas`, `numpy`, `scipy`, `pvlib`. Les scripts attendent une
copie du dépôt dans `./repo`.
