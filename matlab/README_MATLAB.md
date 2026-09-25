# Modélisation dynamique APV-AD — MATLAB
# Dynamic APV-AD modelling — MATLAB

Modèle dynamique couplé (EDO) du système agrivoltaïque + digestion anaérobie,
complémentaire au pipeline statique Python du dépôt.

Coupled ordinary-differential-equation model of the agrivoltaic + anaerobic
digestion system, complementing the repository's static Python pipeline.

---

## 1. Pourquoi un modèle dynamique / Why a dynamic model

**FR.** Le pipeline Python existant est quasi-statique. Le photovoltaïque est
calculé heure par heure, mais le digesteur est un bilan annuel agrégé : les
solides volatils de toute l'année sont convertis en méthane par un unique BMP
corrigé en température, et cette température est elle-même un proxy,
`T_dig = clip(moyenne mensuelle T_amb + 5, 26, 37)`. La culture est un rapport
d'intégrales de PAR sur la saison. Aucune variable d'état n'existe.

Trois questions ne peuvent donc pas être posées :

1. **Le digesteur tient-il la plage mésophile ?** Le proxy statique place
   Konya à 26,7 °C en moyenne annuelle et ne peut pas dire si le digesteur
   tient 37 °C une nuit de janvier avec le chauffage en marche.
2. **La fraction `f_PV_heat` suffit-elle, au bon moment ?** Dans le modèle
   statique, `f_PV_heat` ne fait que réétiqueter une part de l'énergie
   annuelle PV de « vendue » vers « chaleur » — une répartition comptable qui
   ne peut pas échouer. La chaleur PV n'existe pourtant qu'en journée.
3. **Que se passe-t-il à la récolte ?** Les résidus de toute une saison
   arrivent en quelques jours. Un modèle annuel divise ce total par 365 et
   suppose implicitement un stock infini toujours plein.

**EN.** The existing Python pipeline is quasi-static. Photovoltaics are
resolved hourly, but the digester is an annual lumped balance: a whole year of
volatile solids is converted through a single temperature-corrected BMP, and
that temperature is itself a proxy,
`T_dig = clip(monthly mean T_amb + 5, 26, 37)`. The crop is a ratio of
season-integrated PAR. There is no state variable anywhere.

Three questions therefore cannot be asked: whether the digester holds
mesophilic conditions hour by hour; whether the diverted PV fraction arrives
when the digester is actually cold rather than merely summing to enough over a
year; and what happens when a season's residues land in the space of a few
days. This model carries 19 states and answers all three.

---

## 2. Vecteur d'état / State vector

Temps `t` en **jours**, `t = 0` au 1er janvier 00:00. Voir `config/apvad_states.m`.

| # | État | Unité | Description (FR) | Description (EN) |
|---|------|-------|------------------|------------------|
| 1 | `S1` | g DCO/L | substrat organique soluble | soluble organic substrate |
| 2 | `X1` | g MVS/L | biomasse acidogène | acidogenic biomass |
| 3 | `S2` | mmol/L | acides gras volatils | volatile fatty acids |
| 4 | `X2` | g MVS/L | biomasse méthanogène | methanogenic biomass |
| 5 | `Z` | mmol/L | alcalinité totale | total alkalinity |
| 6 | `C` | mmol/L | carbone inorganique total | total inorganic carbon |
| 7 | `T_dig` | °C | température du digesteur | digester temperature |
| 8 | `T_cell` | °C | température du module PV | PV module temperature |
| 9 | `Dr` | mm | déficit hydrique racinaire | root-zone water depletion |
| 10 | `GDD` | °C·j | temps thermique | growing degree days |
| 11 | `B_dm` | kg MS/ha | biomasse aérienne sur pied | standing crop dry biomass |
| 12 | `M_res` | kg MV | stock de substrat disponible | stored feedstock |
| 13 | `E_pv` | kWh/ha | production PV cumulée | cumulative PV generation |
| 14 | `E_ch4` | kWh/ha | énergie méthane cumulée | cumulative methane energy |
| 15 | `E_heat` | kWh/ha | chaleur livrée cumulée | cumulative heat delivered |
| 16 | `E_pv_sold` | kWh/ha | électricité PV exportée | PV electricity exported |
| 17 | `E_ch4_heat` | kWh/ha | méthane brûlé en chaudière | methane burned in boiler |
| 18 | `W_irr` | mm | irrigation cumulée | cumulative irrigation |
| 19 | `M_harv` | kg MS/ha | matière sèche récoltée | harvested dry matter |

Les états 13 à 19 sont de purs intégrateurs. Les porter dans l'EDO fait que les
totaux annuels sortent directement du solveur, sous son propre contrôle
d'erreur, plutôt que d'être ré-intégrés depuis une grille rééchantillonnée.

States 13–19 are pure integrators, so annual totals come straight from the
solver under its own error control rather than from re-integrating a resampled
output grid.

---

## 3. Équations / Equations

### 3.1 Digestion anaérobie — AM2 (Bernard et al. 2001)

Croissance acidogène (Monod) et méthanogène (Haldane, inhibition par les AGV) :

```
mu1 = fT(T) * mu1max * S1 / (KS1 + S1)
mu2 = fT(T) * mu2max * S2 / (KS2 + S2 + S2^2/KI2)
```

Bilans matière du réacteur parfaitement agité (`Dil = Q_in / V_dig`) :

```
dS1/dt = Dil*(S1in - S1) - k1*mu1*X1
dX1/dt = (mu1 - alpha*Dil - kd1)*X1
dS2/dt = Dil*(S2in - S2) + k2*mu1*X1 - k3*mu2*X2
dX2/dt = (mu2 - alpha*Dil - kd2)*X2
dZ/dt  = Dil*(Zin - Z)
dC/dt  = Dil*(Cin - C) - qC + k4*mu1*X1 + k5*mu2*X2
```

Phase gazeuse et équilibre acido-basique :

```
qM    = k6*mu2*X2                       (mmol CH4 / L / j)
CO2   = C + S2 - Z
B     = Z - S2                          (bicarbonate)
phi   = CO2 + KH*Pt + qM/kLa
P_CO2 = ( phi - sqrt(phi^2 - 4*KH*Pt*CO2) ) / (2*KH)
qC    = kLa*(CO2 - KH*P_CO2)
pH    = -log10( Kb * CO2 / B )
```

**Le terme de Haldane est la raison d'être du modèle dynamique.** Le modèle
statique convertit les MV en méthane par un BMP unique : une surcharge produit
simplement proportionnellement plus de gaz. Ici, une surcharge fait monter `S2`,
`mu2` chute dès que `S2` dépasse `sqrt(KS2*KI2)`, moins d'AGV sont consommés,
`S2` monte encore — et le digesteur peut basculer dans l'état acidifié que les
exploitants observent réellement. L'acidification est une trajectoire, pas une
moyenne annuelle.

The Haldane term is why the dynamic model earns its keep: souring is a
trajectory through a bistable system, and an annual average cannot represent
one. `am2_steady_state.m` exposes both branches in closed form — the smaller
root is the healthy, stable operating point, the larger one is unstable.

Le méthane est un **rendement de sortie** : la fraction CH4 du biogaz et le pH
sont calculés, pas imposés. Le modèle statique fixe `CH4_fraction = 0.60`.

### 3.2 Réponse en température / Temperature response

Deux formes, sélectionnées par `D.temp.model` :

- `'arrhenius'` (défaut) : `fT = min(exp(theta*(T - T_ref)), 1)` — exactement
  la forme exponentielle que le modèle statique applique au BMP
  (`theta = 0.069 °C⁻¹`). C'est le défaut **pour que la comparaison
  statique/dynamique porte sur la dynamique et non sur une loi de température
  changée**.
- `'ctmi'` : modèle à températures cardinales avec inflexion (Rosso et al. 1993),
  qui représente aussi l'inhibition au-dessus de l'optimum. À utiliser à
  Ouagadougou, où le digesteur peut dépasser 40 °C — la forme d'Arrhenius y
  récompenserait une surchauffe au lieu de la pénaliser.

### 3.3 Bilan thermique du digesteur / Digester thermal balance

```
C_th * dT_dig/dt = Q_heat + Q_reac - Q_loss - Q_flow
C_th   = V * rho * Cp
Q_loss = U_eff * A_wall * (T_dig - T_amb)
Q_flow = Q_in * rho * Cp * (T_dig - T_in)
Q_reac = Q_reaction_frac * E_CH4
```

`A_wall` provient du même cylindre (hauteur = diamètre) que le modèle statique.
La constante de temps est de l'ordre de **10 jours** pour un digesteur de 4 m³ :
c'est précisément cette inertie qui tamponne une source de chaleur qui n'existe
qu'en journée.

### 3.4 Régulation et dispatch / Control and dispatch

Loi proportionnelle avec saturation, priorité au PV puis chaudière biogaz :

```
Q_demand = clip( Kp*(T_set - T_dig), 0, Q_max )
Q_pv     = min( Q_demand, f_PV_heat * P_PV * eta_elec_heat )
Q_boiler = min( Q_demand - Q_pv, E_CH4 * eta_boiler )
```

**Le dimensionnement du chauffage est calculé, pas fixé.** `build_forcing.m`
calcule la charge de conception (pertes de paroi à la moyenne mensuelle la plus
froide + chauffe de l'influent) et pose `Q_max = 1.5 × Q_design`. Pour Konya
cela donne ≈ 1 kW, pas les 6 kW qu'un chiffre rond suggérerait. Ce point n'est
pas un détail : un chauffage surdimensionné tient la consigne dans tous les
scénarios et détruit la capacité du modèle à répondre à la question même pour
laquelle le couplage existe.

A sized heater is the difference between a model that can fail and one that
cannot. `F.deficit_W` records demand that neither PV nor the boiler covered.

### 3.5 Photovoltaïque / Photovoltaics

```
T_ss = T_amb + G_eff/(U0 + U1*WS)                    (Faiman 2008)
C_areal * dT_cell/dt = G_eff - (U0 + U1*WS)*(T_cell - T_amb)
P = P_STC * N * (1 - eta_loss) * [1 + beta_p*(T_cell - 25)] * G_eff/1000
```

L'état stationnaire de l'EDO **est** exactement l'expression de Faiman, ce que
vérifie `testPVThermalReducesToFaiman`. Avec `D.pv.dynamic = false` (défaut) on
retrouve le modèle statique au bit près, ce qui fait du PV un témoin et non un
facteur confondant dans la comparaison. Voir §6 pour le compromis mesuré.

### 3.6 Bilan hydrique du sol / Soil water balance (FAO-56)

```
dDr/dt   = ETc_adj - (pluie + irrigation) + percolation
TAW      = 1000*(theta_fc - theta_wp)*Zr
Ks       = 1  si Dr <= RAW ;  (TAW - Dr)/(TAW - RAW)  sinon
ETc_adj  = Ks * Kc * ET0 * (1 - alpha_shade*F_shad)
```

Le terme d'ombrage est celui du modèle statique, mais il agit ici sur un
**état**. La différence est le fond du sujet : le modèle statique rapporte l'eau
*économisée*, une grandeur comptable, parce qu'il ne suit jamais l'eau
réellement présente dans le sol. Porter `Dr` comme état fait apparaître le
bénéfice là où il compte pour l'agronome — moins d'heures sous le seuil de
stress, irrigations différées — et le renvoie dans la croissance via `Ks`.

### 3.7 Croissance de la culture / Crop growth (Monteith RUE)

```
dGDD/dt = max(0, min(T_jour, T_cut) - T_base)     en saison
LAI     = LAI_max * logistique(GDD) * sénescence(GDD)
f_i     = 1 - exp(-k_ext*LAI)                     (Beer-Lambert)
dB/dt   = RUE * fT * Ks * f_i * min(PAR_AV, PAR_sat) - récolte
```

`min(PAR_AV, PAR_sat)` est le plafond de saturation lumineuse du modèle
statique, et c'est ce qui rend l'agrivoltaïsme viable : au-dessus de `PAR_sat`
la culture ne peut pas utiliser la lumière supplémentaire, donc l'ombrer ne
coûte rien sur le plan agronomique tout en produisant de l'électricité. Appliquer
ce plafond instantanément plutôt qu'à une intégrale saisonnière change le
résultat, parce qu'il mord à midi et pas à l'aube.

### 3.8 Remise à zéro inter-saisons / Between-season reset

`GDD` et `B_dm` sont des **états**. En simulation pluriannuelle, ils doivent
revenir à zéro avant la saison suivante, sinon la deuxième saison démarre sur
un couvert déjà sénescent et la culture ne produit pratiquement rien. Hors
saison, `GDD` relaxe vers zéro (τ = 2 j) et toute biomasse encore sur pied est
transférée au stock (τ = 5 j).

Ce défaut mérite d'être signalé parce qu'il est **silencieux** : tous les
*rapports* restent plausibles, puisque la simulation de référence en plein champ
accumule exactement le même temps thermique. `LER_crop` survit ; ce sont les
rendements absolus qui s'effondrent, et le digesteur qui se retrouve
discrètement affamé. C'est le diagnostic de périodicité (§6) qui l'a détecté, et
`testThermalTimeResetsBetweenSeasons` qui garde contre son retour.

La récolte est du premier ordre en biomasse sur pied, donc l'intégrale du taux
de prélèvement fixe la fraction récoltée, `1 - exp(-strength)`. Une impulsion
d'aire unitaire laisserait **37 %** de la culture au champ ; `harvest_strength = 6`
n'en laisse que 0,25 %.

---

## 4. Boucles de rétroaction / Feedback loops

Deux boucles n'existent pas du tout dans le modèle statique :

1. **Boucle thermo-biologique.** `T_dig` fixe les vitesses de croissance ; le
   méthane produit chauffe la cuve (réaction exothermique + chaudière d'appoint)
   et constitue la sortie mesurée. Digesteur froid → méthanogenèse ralentie →
   moins de gaz → moins de chaleur d'appoint disponible. Le modèle statique
   coupe cette boucle en imposant la température.

2. **Boucle ombrage–eau–croissance.** L'ombrage réduit l'évapotranspiration,
   donc le déficit `Dr`, donc `Ks` monte, donc la croissance monte, donc les
   résidus montent, donc l'alimentation du digesteur monte. Le modèle statique
   calcule l'économie d'eau mais ne la laisse jamais agir sur le rendement :
   `LER_crop` y est un pur rapport de lumière.

---

## 5. Calibration

Deux constantes sont recalées, chacune **exactement et en un seul passage**,
parce que le modèle en est linéaire.

| Constante | Recalée sur | Fichier | Pourquoi c'est exact |
|-----------|-------------|---------|----------------------|
| `RUE` | rendement en plein champ déclaré du site | `calibrate_crop_rue.m` | le LAI dépend du temps thermique seul, `Ks` de l'état hydrique seul, la récolte est linéaire en `B` → `B` est strictement proportionnelle à `RUE` |
| `k6` | BMP de référence à 37 °C (300 NmL CH4/g MV) | `calibrate_am2_to_bmp.m` | `k6` n'apparaît que dans `qM = k6*mu2*X2`, jamais dans les bilans matière |

**Point de méthode important.** Les constantes AM2 publiées ont été identifiées
sur des vinasses de distillerie vinicole ; leur rendement absolu en méthane n'a
aucune prétention à valoir pour des résidus de tomate co-digérés avec du fumier
bovin. Quelque chose doit fixer l'échelle. La tentation serait de régler le
modèle dynamique jusqu'à retrouver le méthane annuel du modèle statique — ce
serait circulaire, et cela détruirait par construction le seul résultat
intéressant. L'échelle est donc fixée **à la condition où `BMP_ref` est
effectivement défini** : un test en batch à 37 °C. Tout ce qui suit — ce que le
digesteur produit à Konya en février, la part de l'année passée sous la plage
mésophile — est alors une *prédiction* du modèle dynamique, et son accord ou son
désaccord avec le chiffre annuel statique est un résultat, non un paramètre
ajusté.

The same point in English: `k6` is anchored at the condition where `BMP_ref` is
defined — a 37 °C batch assay — not against the annual figure being compared to.
Otherwise the benchmark would be fitted to its own answer.

### Base de charge organique / Organic loading basis

`main_apvad_dynamic.m` exécute le cas agrivoltaïque **en deux passes**. La passe A
découvre la production réelle de résidus sous les panneaux ; la passe B
reconstruit la dilution de l'influent sur ce chiffre réalisé, de sorte que la
charge organique, la dilution et le substrat réellement disponible soient
mutuellement cohérents. Sans cela, le digesteur est dimensionné pour une récolte
de plein champ qu'il ne reçoit jamais et passe une partie de l'année le stock
vide — un artefact de la base de conception, non une propriété du système.

---

## 6. Numérique / Numerics

**Raideur.** Les constantes de temps vont de ≈ 12 minutes (masse thermique du
module PV) à quelques jours (température du digesteur, croissance des
méthanogènes) puis à la saison entière (temps thermique, stock de substrat) —
quatre ordres de grandeur dans un même système. C'est la définition d'un problème
raide : une méthode explicite serait contrainte à la plus petite constante de
temps toute l'année pour des raisons de stabilité et non de précision. `ode15s`
(NDF à ordre variable) ne réduit le pas que pour la précision.

**Commutations lissées.** Chaque décision de commutation — démarrer
l'irrigation, cesser de puiser dans un stock vide, saturer le chauffage — passe
par `smooth_step` / `smooth_min` plutôt que par un `if`. Un solveur à pas
variable estime l'erreur locale à partir d'un développement de Taylor lisse ;
une vraie discontinuité rend cette estimation caduque, le pas s'effondre et
l'intégration soit rampe, soit enjambe silencieusement l'événement. Les largeurs
de lissage **suivent l'échelle de la grandeur** (2 % de la capacité du chauffage,
2 % de `PAR_sat`, une demi-journée de demande d'alimentation) : une largeur fixe
de 1 W sur un chauffage de 1 kW reste un point anguleux, franchi deux fois par
jour tous les jours de l'année.

**Motif de jacobienne.** Rien dans le modèle ne *lit* un accumulateur : ce sont
des intégrateurs en écriture seule. Leurs colonnes de jacobienne sont donc
structurellement nulles, ce que `ix.JPattern` transmet à `ode15s` — environ un
tiers de coût en moins par jacobienne.

**PV thermique.** `D.pv.dynamic = false` par défaut, et c'est un compromis
*mesuré*, non une préférence. Résoudre une constante de temps de 12 minutes
impose au solveur des pas de l'ordre de la minute toute l'année : une année passe
de ≈ 30 s à plusieurs minutes. `experiment_pv_thermal_lag.m` quantifie ce que
cela achète — face à un forçage horaire, le retard modifie le productible annuel
d'une fraction de pour cent, parce que le module passe l'année à retarder
alternativement un soleil montant et un soleil descendant, et que les deux se
compensent presque exactement. Activez-le pour l'infra-horaire, les vitesses de
rampe, ou toute question où la mémoire thermique du module est le sujet.

**Année de mise en régime.** Une année de spin-up est intégrée puis écartée
(`opts.n_years = 2`). Une AMT est par construction une année représentative et
non une trajectoire : la répéter est la bonne façon d'atteindre un régime
périodique. Intégrer une seule année depuis des conditions initiales arbitraires
reviendrait à publier le transitoire de démarrage du digesteur comme s'il
s'agissait d'un résultat annuel.

**Calendrier AMT.** Les fichiers PVGIS assemblent l'année mois par mois à partir
d'années civiles différentes (Konya : 2005, 2008, 2009, 2010, 2012, 2013, 2017,
2018). Leurs horodatages reculent donc à chaque frontière de mois, et les mois
issus d'une année bissextile décalent le jour julien d'une unité — la colonne
`day_of_year` du fichier contient des valeurs jusqu'à 366 mais seulement 364
valeurs distinctes, avec 24 heures dupliquées. `load_climate.m` reconstruit
l'axe temporel en projetant `(mois, jour, heure)` sur une **année de référence
non bissextile**, puis trie. Sans cela `griddedInterpolant` refuse les données.

---

## 7. Utilisation / Usage

```matlab
cd('<racine du projet>/matlab')

% Un site, avec l'optimum DE lu dans outputs/csv/optimal_params.json
out = main_apvad_dynamic('Konya');

% Les quatre sites + tableau de synthèse (env. 3-4 min par site)
S = run_all_sites();

% Modèle à températures cardinales, pertinent en climat chaud
out = main_apvad_dynamic('Ouagadougou', struct('temp_model','ctmi'));

% Point de conception imposé
d = struct('beta_deg',25,'d_row_m',6,'H_m_m',2.5, ...
           'V_dig_m3',10,'HRT_days',30,'f_PV_heat',0.2);
out = main_apvad_dynamic('Konya', struct('design',d));

% Redessiner les figures sans réintégrer
plot_from_saved('Konya');

% Tests unitaires
run_all_tests
```

Sorties écrites dans `outputs/` :

| Fichier | Contenu |
|---------|---------|
| `csv/dynamic_<site>.mat` | résultat complet (états, sorties, agrégats) |
| `csv/dynamic_<site>_summary.csv` | tableau comparatif statique / dynamique |
| `csv/dynamic_all_sites.csv` | synthèse inter-sites (`run_all_sites`) |
| `figures/dyn_annual_<site>.png` | 8 panneaux, moyennes journalières |
| `figures/dyn_detail_<site>.png` | fenêtre de 10 jours au pas horaire |

**Prérequis.** `outputs/csv/hourly_<site>.csv` doit exister — exécuter d'abord
`01_load_data.ipynb`. `outputs/csv/optimal_params.json` est optionnel : sans lui,
le point de conception est le milieu des bornes DE et un avertissement est émis.
Aucune toolbox n'est requise au-delà de MATLAB de base (testé sous R2024a).

---

## 7bis. Résultats de référence — Konya / Benchmark results — Konya

Optimum DE de Konya, année en régime périodique, `temp_model = 'arrhenius'`.

| Grandeur | Statique | Dynamique | Écart |
|----------|---------:|----------:|------:|
| PV total (MWh/ha) | 1295,2 | 1297,4 | **+0,2 %** |
| PV exporté (MWh/ha) | 1228,7 | 1295,4 | +5,4 % |
| Biogaz (MWh/ha) | 3,66 | 8,76 | +139 % |
| Demande de chaleur (MWh/ha) | 2,26 | 3,67 | +62 % |
| BMP effectif (NmL/g MV) | 148,8 | 302,3 | +103 % |
| T digesteur moyenne (°C) | 26,7 | 35,7 | +33 % |
| GCR (%) | 70,32 | 70,32 | **−0,002 %** |
| ET0 saison (mm) | 845,6 | 846,4 | **+0,09 %** |
| Eau économisée (mm) | 159,11 | 159,04 | **−0,05 %** |
| LER_crop (déf. rapport de PAR) | 0,356 | 0,402 | +12,9 % |
| LER_crop (déf. biomasse) | 0,356 | 0,679 | +90,9 % |
| eLER | 1,641 | 2,200 | +34,0 % |

**Ce qui doit concorder concorde.** Le PV total, le GCR, l'ET0 saisonnier et
l'eau économisée s'accordent à 0,2 % près. Ce sont les grandeurs partagées entre
les deux modèles : leur concordance valide le portage, et un écart matériel y
aurait signalé un bug, pas une découverte.

**Ce qui diffère, et pourquoi.**

1. **Température du digesteur (26,7 → 35,7 °C).** Le proxy statique fait flotter
   le digesteur à l'ambiante + 5 °C. Le modèle dynamique constate que la demande
   de chaleur d'une cuve de 4,3 m³ est d'environ **600 W** — soit 0,05 % d'un
   champ PV de 780 kWc — et que le chauffage tient donc la consigne **99,2 % de
   l'année**, dont 47 % fourni par le PV et le reste par la chaudière biogaz. Le
   BMP effectif double en conséquence, et le biogaz avec lui.

   *Conséquence pour le manuscrit :* le modèle statique sous-estime largement le
   biogaz, non par une erreur de cinétique, mais parce que son hypothèse de
   température ignore que le PV co-implanté peut alimenter le chauffage sans
   effort. C'est un argument **en faveur** du couplage APV-AD, pas contre.

2. **LER_crop — un problème de définition, pas de physique.** La formule
   statique est `∑min(PAR_AV, PAR_sat) / ∑PAR_open` : le **numérateur est
   plafonné à la saturation lumineuse, le dénominateur ne l'est pas**. Ce
   rapport est donc inférieur à 1 même **sans aucun ombrage** — il facture au
   système agrivoltaïque la saturation lumineuse propre à la culture, un effet
   qui n'a rien à voir avec les panneaux.

   Recalculée à l'identique sur la trajectoire dynamique, cette définition donne
   0,402 contre 0,356 pour le modèle statique : les deux modèles sont d'accord,
   c'est la **définition** qui est en cause. La définition par biomasse — rapport
   des matières sèches récoltées, avec le même plafond appliqué aux deux termes —
   donne 0,679, et isole effectivement l'effet d'ombrage.

   *Conséquence pour le manuscrit :* cette asymétrie devrait être corrigée ou
   explicitement justifiée. Elle biaise `LER_crop` vers le bas à tous les sites,
   donc `eLER` aussi, et la comparaison au benchmark de 1,94 (Riaz et al. 2022)
   s'en trouve affectée.

3. **LER_PV (+5,4 %).** Le modèle statique retranche `f_PV_heat` de la production
   annuelle sans condition, soit 5,1 % du champ. Le digesteur n'en consomme
   réellement que **0,13 %** ; le reste est exporté. Le modèle statique facture
   donc la dérivation environ **quarante fois** son coût réel — ce qui, entre
   autres, explique pourquoi l'optimiseur pousse `f_PV_heat` vers sa borne
   inférieure aux quatre sites.

4. **Fraction CH4 (73 % calculée contre 60 % imposée).** Sortie du bilan
   acido-basique AM2 plutôt qu'une constante. À traiter avec prudence : les
   constantes de partage du CO2 d'AM2 ont été identifiées sur vinasses, et 73 %
   est au-dessus de la plage usuelle de 55–70 % pour une co-digestion résidus +
   fumier. C'est un candidat prioritaire pour une ré-identification (§9.1).

**Points de fonctionnement.** Au point de conception DE, le digesteur est très
loin de toute instabilité : AGV à 0,6 mmol/L contre un pic de Haldane à
49 mmol/L, pH stable à 7,20, jamais limité par le substrat, 1,1 % de la saison
sous stress hydrique. Le modèle dynamique confirme ici la robustesse du point
optimal — un résultat négatif, mais un résultat.

---

## 8. Métriques que seul le modèle dynamique produit

Ces grandeurs n'ont pas d'équivalent dans le pipeline statique et sont, en
pratique, ce qui explique les écarts du tableau comparatif :

| Métrique | Signification |
|----------|---------------|
| `heat_coverage` | part de la demande de chaleur effectivement couverte |
| `pv_share_of_heat` | part de la chaleur livrée provenant du PV |
| `frac_time_below_T_min` | temps passé sous la température minimale effective |
| `frac_time_mesophilic` | temps passé au-dessus de 35 °C |
| `frac_time_feed_limited` | temps où le stock de substrat limite l'alimentation |
| `frac_time_VFA_inhibited` | temps passé au-delà du pic de Haldane |
| `pH_min` | pH minimal atteint sur l'année |
| `x_CH4_mean` | fraction CH4 calculée (le statique impose 0,60) |
| `frac_season_water_stressed` | part de la saison sous stress hydrique |

---

## 9. Limites / Limitations

1. **Constantes cinétiques AM2 non spécifiques au substrat.** Identifiées sur
   vinasses (Bernard et al. 2001). Seul `k6` est recalé ; `mu1max`, `mu2max`,
   `KS1`, `KS2`, `KI2` devraient être ré-identifiés sur des mesures de
   co-digestion résidus de tomate + fumier bovin. La position du seuil
   d'inhibition en dépend directement.
2. **Pas de précipitations.** Les fichiers TMY PVGIS n'en contiennent pas : le
   bilan hydrique est piloté par la seule irrigation. C'est acceptable en climat
   semi-aride (Konya, Almería, Ouagadougou) mais **matériellement faux à
   Freiburg** (Cfb, océanique). Ajouter une colonne de pluie et renseigner
   `D.water.rain_column` avant d'exploiter les résultats hydriques de Freiburg.
3. **Digesteur parfaitement agité, monophasé.** Pas de stratification, pas de
   rétention de biomasse différenciée au-delà du paramètre `alpha`.
4. **Ciel isotrope pour le POA**, hérité du modèle statique (écart documenté de
   1,5–3,9 % contre Perez).
5. **Une seule culture, un seul cycle par an.**
6. **Aucun stockage thermique dédié.** L'inertie de la cuve est la seule mémoire
   thermique ; un ballon tampon changerait les conclusions sur `f_PV_heat`.
7. **Pas encore de couplage à l'optimiseur.** Le modèle dynamique est ici un
   banc d'évaluation, pas la fonction objectif de la DE. Le coût par évaluation
   (dizaines de secondes) l'interdit en l'état ; une version réduite ou un
   métamodèle serait nécessaire.

---

## 10. Arborescence / File map

```
matlab/
├── main_apvad_dynamic.m        pilote de haut niveau / top-level driver
├── config/
│   ├── apvad_sites.m           paramètres des 4 sites (portés de config_00.py)
│   ├── apvad_params.m          constantes partagées avec le modèle statique
│   ├── apvad_dyn_params.m      paramètres propres au dynamique (AM2, thermique…)
│   └── apvad_states.m          index, unités, motif de jacobienne
├── forcing/
│   ├── load_climate.m          lecture AMT + reconstruction du calendrier
│   ├── build_forcing.m         précalcul + interpolants + dimensionnement
│   ├── solar_angles.m          déclinaison, équation du temps, élévation
│   ├── shading_fraction.m      fraction d'ombrage au sol
│   ├── poa_irradiance.m        transposition plan des modules, bifacial
│   └── et0_hargreaves.m        ET0 Hargreaves-Samani / FAO-56
├── models/
│   ├── apvad_rhs.m             second membre assemblé / assembled RHS
│   ├── am2_kinetics.m          Monod + Haldane
│   ├── am2_gas.m               CH4, CO2, composition, pH
│   ├── am2_steady_state.m      équilibre en forme close (2 branches)
│   ├── growth_temp_factor.m    Arrhenius | CTMI
│   ├── digester_thermal.m      bilan thermique de la cuve
│   ├── pv_module.m             Faiman dynamique + puissance
│   ├── crop_growth.m           temps thermique, LAI, RUE
│   ├── soil_water.m            FAO-56 déficit racinaire
│   ├── feed_scheduler.m        stock et composition de l'influent
│   ├── smooth_step.m           commutation différentiable
│   └── smooth_min.m            minimum différentiable
├── control/
│   └── digester_thermostat.m   régulation proportionnelle, dispatch PV d'abord
├── run/
│   ├── run_dynamic_sim.m       enveloppe ode15s
│   ├── run_all_sites.m         les 4 sites + synthèse
│   ├── apvad_initial_state.m   condition initiale (stock à l'équilibre)
│   ├── postprocess_sim.m       sorties algébriques, agrégats, périodicité
│   ├── calibrate_crop_rue.m    recalage exact du RUE
│   ├── calibrate_am2_to_bmp.m  recalage exact de k6
│   └── compare_with_static.m   comparaison au pipeline Python
├── figures/
│   ├── plot_apvad_dynamics.m   2 figures publication
│   └── plot_from_saved.m       redessin sans réintégration
└── tests/
    ├── APVADTest.m             suite unitaire (25 tests)
    ├── run_all_tests.m         lanceur
    ├── smoke_test.m            vérification structurelle rapide
    └── experiment_pv_thermal_lag.m   compromis PV dynamique / algébrique
```

Les tests s'appuient autant que possible sur des **formes closes et des
identités** plutôt que sur des valeurs attendues stockées : l'angle d'incidence
contre une identité trigonométrique au midi solaire, l'EDO thermique du module
contre l'expression de Faiman vers laquelle elle doit relaxer, le digesteur
intégré contre son propre équilibre analytique. Ils restent donc valides quand
les paramètres changent, et échouent quand la physique casse.

---

## 11. Références / References

- Bernard, O., Hadj-Sadok, Z., Dochain, D., Genovesi, A., Steyer, J.-P. (2001).
  Dynamical model development and parameter identification for an anaerobic
  wastewater treatment process. *Biotechnology and Bioengineering*, 75(4),
  424–438.
- Rosso, L., Lobry, J.R., Flandrois, J.P. (1993). An unexpected correlation
  between cardinal temperatures of microbial growth highlighted by a new model.
  *Journal of Theoretical Biology*, 162(4), 447–463.
- Faiman, D. (2008). Assessing the outdoor operating temperature of photovoltaic
  modules. *Progress in Photovoltaics*, 16(4), 307–315.
- Duffie, J.A., Beckman, W.A. (2013). *Solar Engineering of Thermal Processes*,
  4th ed. Wiley.
- Spencer, J.W. (1971). Fourier series representation of the position of the sun.
  *Search*, 2(5), 172.
- Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). *Crop
  evapotranspiration*. FAO Irrigation and Drainage Paper 56.
- Monteith, J.L. (1977). Climate and the efficiency of crop production in
  Britain. *Philosophical Transactions of the Royal Society B*, 281, 277–294.
- McCree, K.J. (1971). The action spectrum, absorptance and quantum yield of
  photosynthesis in crop plants. *Agricultural Meteorology*, 9, 191–216.
- Shampine, L.F., Reichelt, M.W. (1997). The MATLAB ODE suite. *SIAM Journal on
  Scientific Computing*, 18(1), 1–22.
