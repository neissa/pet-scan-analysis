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

## Limites (à connaître)

- Recherche / apprentissage uniquement, **pas un dispositif médical**.
- Le seuil fixe capte aussi les fixations physiologiques (cerveau, cœur, vessie, reins) : pas de segmentation d'organes.
- Pas de SUL/SUV lean body mass, pas de SUVpeak, pas de DICOM multi-séries, pas de correction si l'heure d'acquisition n'est pas celle du début de scan.
- Ne committez jamais de DICOM (données patient) ; `.gitignore` les exclut.
