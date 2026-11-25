# pet-scan-analysis

Analyse d'un TEP scan injecté (ex. 18F-FDG) à partir d'une série DICOM :

1. lecture de la série PET (`*.dcm`, unités BQML) ;
2. conversion en **SUV** (poids corporel), avec correction de la décroissance entre l'heure d'injection et l'acquisition ;
3. détection des lésions par seuil SUV (défaut 2.5) → **SUVmax, SUVmean, MTV, TLG** par lésion ;
4. option : MIP coronale en PNG.

```
pip install -e .
python petscan.py chemin/vers/serie_dicom --thr 2.5 --mip mip.png
pytest
```

## Interface web

```
python viewer.py chemin/vers/serie_dicom   # puis http://127.0.0.1:8000
python viewer.py --demo                    # volume synthétique, sans DICOM
```

Coupes axiale / coronale / sagittale (sliders), MIP, fenêtre SUV, seuil et tableau des lésions. Fusion TEP+CT si le dossier contient un CT.
Molette = zoom centré sur la souris. Outils : *Naviguer* (un clic positionne les autres coupes), *Règle* (2 clics, distance en mm),
*Mesure auto* (clic sur une structure : région connexe ≥ X % du max local autour du clic ; plus grand diamètre + étendue perpendiculaire, mesurés sur le TEP dans la coupe affichée). Le serveur écoute uniquement en local (données patient). À lancer depuis la racine du repo (`index.html` est lu à côté de `viewer.py`).

### Fonctions de la visionneuse

- **Dossier** : `python viewer.py` sans argument, puis champ de chemin, « Parcourir… » (sélecteur natif, tkinter) et dossiers récents (`~/.pet-scan-analysis.json`, hors repo). Le chargement se fait en arrière-plan ; une erreur laisse l'ancien volume affiché.
- **Vues** : 4 vues de même hauteur, zoom molette centré sur la souris, molette sur tous les curseurs (Maj = ×10), fusion TEP + CT.
- **Outils** : règle en temps réel (barres aux extrémités), mesure auto sur le TEP (région ≥ X % du max local) ou sur le CT (± tolérance en HU), recalculées en direct quand le curseur bouge, loupe ×4 à viseur creux.
- **Série CT** : si le dossier en contient plusieurs (ex. OS MAR, STANDARD MAR, THORAX APNEE), un menu « Série CT » permet de changer ; seul le CT est rechargé (quelques secondes), la coupe courante est conservée. Les séries qui ne couvrent qu'une partie du TEP sont signalées « partiel ».
- **MIP rotatif** : sur sa propre ligne, curseur d'angle (molette = 10° par cran) et bouton Lecture. Les images (tous les 10°, grille réduite ×2) sont calculées en arrière-plan après le chargement ; l'angle 0 est en pleine résolution. Règle et mesure auto fonctionnent sur le MIP à l'angle affiché.
- **Liste des lésions** : un clic sur une ligne (ou Entrée) amène les trois coupes sur le pic de la lésion et pose un repère (anneau à centre vide) sur chaque vue.
- **Comparaison de deux séances** (même patient, dates différentes) : section « Comparaison de séances », dossier de la séance B puis « Comparer ». B est recalée sur A par une **translation entière** estimée par corrélation du squelette au CT (`compare.py`, quelques secondes ; le décalage et un score de corrélation sont affichés). Affichages : A, B, curseur A | B, différence A − B (bleu = baisse, rouge = hausse). Le tableau des lésions gagne « SUVmax B » et « Δ % » (max de B dans un rayon de 12 mm autour du pic de A). **Alignement fin** : « Alignement fin : décaler chaque séance » donne trois curseurs (z, y, x) par séance, indépendants pour A et pour B (en voxels, avec l'équivalent en mm) ; le tableau des lésions (SUVmax B et Δ %), les mesures et le PDF tiennent compte des décalages. Le bouton **« Aligner B sur A (auto) »** cherche, à partir des décalages actuels, le décalage supplémentaire qui superpose le mieux le squelette de B sur celui de A, puis règle les curseurs de B : « Autour de l'affichage » (zone de ±70 mm autour des coupes affichées, au voxel près, environ 1 s) ou « Corps entier » (pas de 3 voxels, environ 2 s). Une corrélation faible (< 0,3) est signalée. Le panneau de droite est redimensionnable (on glisse son bord gauche). Limites : pas de rotation ni de déformation (tête, bras et jambes peuvent rester décalés, surtout en bord de champ) ; les SUV ne sont comparables que si dose, poids, temps d'incorporation et machine sont proches ; le MIP tourné et les organes ne concernent que A.
- **Patient** : cadre replié par défaut, lu dans les en-têtes DICOM (jamais enregistré ni journalisé).
- **Organes** (approximatif) : poumons, cœur, foie, rate, reins, vessie, cerveau, squelette, par règles sur le CT et le TEP (`organs.py`), hors ligne. Les seuils HU ont été réglés sur un CT injecté ; un CT sans contraste ou un patient atypique peut donner de mauvais résultats.
- **Export PDF** : les 3 coupes (avec mesures et, si activée, la surcouche d'organes), les infos patient et la liste des mesures. Le PDF contient l'identité du patient : ne le partagez pas sans précaution.
- Le serveur n'écoute que sur 127.0.0.1 et refuse les requêtes dont l'hôte n'est pas local.

## Limites (à connaître)

- Recherche / apprentissage uniquement, **pas un dispositif médical**.
- Le seuil fixe capte aussi les fixations physiologiques (cerveau, cœur, vessie, reins) : pas de segmentation d'organes.
- Pas de SUL/SUV lean body mass, pas de SUVpeak, pas de DICOM multi-séries, pas de correction si l'heure d'acquisition n'est pas celle du début de scan.
- Ne committez jamais de DICOM (données patient) ; `.gitignore` les exclut.
