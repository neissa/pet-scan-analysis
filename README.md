# pet-scan-analysis

Analyse d'un TEP scan injecté (ex. 18F-FDG) à partir d'une série DICOM :

1. lecture de la série PET (`*.dcm`, unités BQML) ;
2. conversion en **SUV** (poids corporel), avec correction de la décroissance entre l'heure d'injection et l'acquisition ;
3. détection des points de fixation par seuil SUV (défaut 2.5) → **SUVmax, SUVmean, MTV, TLG** par point de fixation ;
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

Coupes axiale / coronale / sagittale (sliders), MIP, fenêtre SUV, seuil et tableau des points de fixation. Fusion TEP+CT si le dossier contient un CT.
Molette sur une image = coupe suivante / précédente (Maj = ×10, angle ±10° sur le MIP) ; Ctrl + molette = zoom centré sur la souris. Outils : *Naviguer* (un clic positionne les autres coupes), *Règle* (2 clics, distance en mm),
*Mesure auto* (clic sur une structure : région connexe ≥ X % du max local autour du clic ; plus grand diamètre + étendue perpendiculaire, mesurés sur le TEP dans la coupe affichée). Le serveur écoute uniquement en local (données patient). À lancer depuis la racine du repo (`index.html` est lu à côté de `viewer.py`).

### Fonctions de la visionneuse

- **Dossier** : `python viewer.py` sans argument, puis champ de chemin, « Parcourir… » (sélecteur natif, tkinter) et dossiers récents (`~/.pet-scan-analysis.json`, hors repo). Le chargement se fait en arrière-plan ; une erreur laisse l'ancien volume affiché.
- **Vues** : une ligne axial / coronal / sagittal par séance (A ; avec une séance B : lignes A, B et « Mix » A | B ou différence A − B), toutes de même hauteur et remplissant la largeur ; coupes, zoom, repères et mesures synchronisés entre les lignes. Curseur vertical de zoom à gauche de chaque ligne (centré sur le point courant), Ctrl + molette = zoom centré sur la souris, molette = coupes, molette sur tous les curseurs (Maj = ×10), fusion TEP + CT. La roue crantée (à côté d'« Exporter PDF ») masque ou affiche le panneau de réglages.
- **Outils** (barre sous les dossiers, avec « Réinitialiser le zoom » ; région % et TEP/CT apparaissent avec la mesure auto) : règle en temps réel, semi-transparente (barres aux extrémités), mesure auto sur le TEP (région ≥ X % du max local) ou sur le CT (± tolérance en HU), recalculées en direct quand le curseur bouge, loupe ×4 à viseur creux.
- **Série CT** : si le dossier en contient plusieurs (ex. OS MAR, STANDARD MAR, THORAX APNEE), un menu « Série CT » permet de changer ; seul le CT est rechargé (quelques secondes), la coupe courante est conservée. Les séries qui ne couvrent qu'une partie du TEP sont signalées « partiel ».
- **MIP rotatif** : sur sa propre ligne (avec B : MIP de A, de B et mix côte à côte, B et mix de face), curseur d'angle (molette = 10° par cran) et bouton Lecture. Les images (tous les 10°, grille réduite ×2) sont calculées en arrière-plan après le chargement ; l'angle 0 est en pleine résolution. Règle et mesure auto fonctionnent sur le MIP à l'angle affiché.
- **Points de fixation** : tableau triable (colonne « Localisation » : organe dans ou près duquel se trouve le point, une fois les organes calculés) ; un clic sur une ligne (ou Entrée) amène les trois coupes sur le pic et pose un repère (anneau à centre vide) sur chaque vue.
- **Comparaison de deux séances** (même patient, dates différentes) : champ « Séance B » en haut de page puis « Comparer ». B est recalée sur A par une **translation entière** estimée par corrélation du squelette au CT (`compare.py`, quelques secondes ; le décalage et un score de corrélation sont affichés). Affichage : lignes A, B et mix (A | B avec une poignée qui déplace la séparation, ou différence A − B : bleu = baisse, rouge = hausse, au choix dans « Affichage »). Le tableau des points de fixation gagne « SUVmax B » et « Δ % » (max de B dans un rayon de 12 mm autour du pic de A). **Alignement fin** : « Alignement fin : décaler chaque séance » donne trois curseurs (z, y, x) par séance, indépendants pour A et pour B (en voxels, avec l'équivalent en mm) ; le tableau des points de fixation (SUVmax B et Δ %), les mesures et le PDF tiennent compte des décalages. Le bouton **« Aligner B sur A (auto) »** cherche, à partir des décalages actuels, le décalage supplémentaire qui superpose le mieux le squelette de B sur celui de A, puis règle les curseurs de B : « Autour de l'affichage » (zone de ±70 mm autour des coupes affichées, au voxel près, environ 1 s) ou « Corps entier » (pas de 3 voxels, environ 2 s). Une corrélation faible (< 0,3) est signalée. Le panneau de droite est redimensionnable (on glisse son bord gauche). Limites : pas de rotation ni de déformation (tête, bras et jambes peuvent rester décalés, surtout en bord de champ) ; les SUV ne sont comparables que si dose, poids, temps d'incorporation et machine sont proches ; le MIP tourné et les organes ne concernent que A.
- **Patient** : cadre replié par défaut, lu dans les en-têtes DICOM (jamais enregistré ni journalisé).
- **Organes** (approximatif) : poumons, cœur, foie, rate, reins, vessie, cerveau, squelette, par règles sur le CT et le TEP (`organs.py`), hors ligne. Les seuils HU ont été réglés sur un CT injecté ; un CT sans contraste ou un patient atypique peut donner de mauvais résultats.
- **Export PDF** : les 3 coupes (avec mesures et, si activée, la surcouche d'organes), les infos patient et la liste des mesures. Le PDF contient l'identité du patient : ne le partagez pas sans précaution.
- Le serveur n'écoute que sur 127.0.0.1 et refuse les requêtes dont l'hôte n'est pas local.

## Limites (à connaître)

- Recherche / apprentissage uniquement, **pas un dispositif médical**.
- Le seuil fixe capte aussi les fixations physiologiques (cerveau, cœur, vessie, reins) : pas de segmentation d'organes.
- Pas de SUL/SUV lean body mass, pas de SUVpeak, pas de DICOM multi-séries, pas de correction si l'heure d'acquisition n'est pas celle du début de scan.
- Ne committez jamais de DICOM (données patient) ; `.gitignore` les exclut.
