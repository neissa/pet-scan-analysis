<p align="center"><img src="assets/netscan-256.png" alt="NETScan" width="160"></p>

# NETScan

**NETScan** est une visionneuse locale pour études TEP/CT (18F-FDG) au format DICOM : SUV, points de fixation, mesures, comparaison de deux séances.

> Recherche / apprentissage uniquement, **pas un dispositif médical**. Les données restent sur votre machine.

## Téléchargement

Les applications prêtes à l'emploi sont dans la section [**Releases**](../../releases) :

- `PetScan-win64.zip` : Windows, sans installer Python (décompresser, lancer `PetScan.exe` ; SmartScreen peut avertir, exécutable non signé).
- `PetScan-android-debug.apk` : tablette Android 8.0+ (APK de debug, à autoriser en « sources inconnues »).

Ou depuis les sources :

```
pip install -e .
python viewer.py "<dossier DICOM>"   # http://127.0.0.1:8000
python viewer.py --demo              # volume synthétique, sans DICOM
python petscan.py <dossier> --thr 2.5 --mip mip.png   # ligne de commande
pytest -q
```

## Fonctions

### Chargement et SUV
- **Lecture DICOM rapide** : ouvre un dossier de séries TEP (et CT) en quelques secondes. Lecture parallèle adaptée au nombre de cœurs, barre de progression, en-têtes mis en cache.
- **SUV poids corporel** : valeurs quantitatives. Correction de la décroissance entre l'injection et l'acquisition.
- **Sélecteur de dossier au clic** : un clic dans le champ A ou B ouvre le sélecteur natif (Windows, Android). A et B sont indépendants ; « Charger » valide un chemin saisi.
- **Série CT au choix** : plusieurs CT dans l'examen ? Un menu permet de changer (les séries partielles sont signalées).
- **Infos patient** : lues dans les en-têtes, cadre replié par défaut, jamais enregistrées.

### Visualisation
- **Trois coupes** : axial, coronal, sagittal synchronisés, avec fenêtre SUV réglable et fusion TEP + CT.
- **Opacité du TEP** : curseur de 0 à 100 %. En mesure, réglé automatiquement (TEP ≥ 70 %, CT ≤ 20 %), modifiable ensuite.
- **MIP rotatif** : projection d'intensité maximale, curseur d'angle et lecture ; trois MIP alignés sur une ligne en comparaison.
- **Zoom et navigation** : molette = coupe, Ctrl + molette = zoom sous la souris, curseur de zoom, boutons −/+ à répétition, réinitialisation.
- **Autocentrage** : recentre l'image sur le point courant à chaque changement.
- **Thème** : clair, sombre ou identique au système (réglages).
- **Dates des séances** : JJ/mm, ou JJ/mm/aaaa si l'année diffère de l'année en cours.

### Points de fixation
- **Détection par seuil SUV** : SUVmax, SUVmean, MTV, TLG et pic par lésion. Seuil réglable par −/+ ou au clavier.
- **Filtre de volume** : masque les petites lésions sous un volume minimal (mL), avec compteur.
- **Tableau triable** : tri par n'importe quelle colonne, y compris la distance au dernier point cliqué.
- **Favoris** : étoile pour garder des points et y revenir.
- **Survol** : un anneau sur les trois vues repère la lésion survolée dans le tableau.
- **Localisation** : organe dans ou près duquel se trouve le point.
- **PERCIST** : SUL (masse maigre, James ; repli SUV poids signalé sans taille ni sexe), SULpeak (sphère de 1 cm³) et SULmax par lésion, référence foie (sphère de 3 cm placée automatiquement, repli marge 0 puis sphère de 2 cm si le foie segmenté est étroit, moyenne ± écart-type, seuil 1,5 × moy. + 2 σ), seuil de détection PERCIST (référence foie ramenée en SUV) ou manuel, voxels de la vessie, des reins, du cerveau et du cœur (+ 1 cm) retirés avant détection (une lésion à leur contact est signalée), case « retenir » par ligne et totaux MTV / TLG des lésions retenues.

### Mesures
- **Règle** : distance en mm sur 2 clics, recalculée en direct.
- **Mesure auto** : région connexe autour du clic, sur le TEP (≥ X % du max local) ou le CT (± tolérance HU). Plus grand diamètre et étendue perpendiculaire.
- **Loupe ×4** : viseur creux pour placer les points précisément ; au tactile, appui long.

### Comparaison de deux séances
- **Recalage automatique** : B recalée sur A par translation (corrélation du squelette au CT), avec score.
- **Alignement fin** : curseurs z/y/x par séance, et bouton « Aligner B sur A (auto) » (zone locale ou corps entier).
- **Mix et différence** : A | B avec poignée de séparation, ou différence A − B (bleu = baisse, rouge = hausse).
- **Évolution par lésion** : SUVmax B et Δ % dans le tableau.

### Suivi multi-dates (carte de référence)
- **Définir comme référence** : les points retenus (✓) deviennent la carte du patient (position DICOM et voxel, date, SUVmax/SULmax, SULpeak, volume à 50 % du max local, organe) ; ajout d'un point au clic, renommage, note, exclusion, retrait.
- **Patient reconnu** à l'ouverture d'une séance ; « Appliquer le suivi » (jamais automatique) recale la séance sur la référence (squelette au CT, raffiné à ±15 mm par lésion), puis mesure chaque lésion à sa position, sans seuil global (voxels de la vessie, des reins, du cœur et du cerveau écartés ; ⚠ « excrétion possible » si max > 20 ou urine à < 20 mm, la lésion restant mesurée) et propose les nouvelles fixations.
- **Tableau de suivi** : une colonne par date (n dates), Δ % du SUVmax vs référence et vs date précédente ; « Enregistrer cette date » ajoute les mesures (nombres seulement) et les conditions de la séance (activité, poids, incorporation, machine, unité SUL/SUV, référence foie). MTV/TLG ne servent à aucun statut.
- **Méthode d'analyse** : « Type de cancer » (sarcome desmoplastique à petites cellules rondes, autre / général) puis « Méthode » (PERCIST 1.0 par défaut, ou RECIST 1.1), choix mémorisé avec le patient pseudonymisé. Registre extensible (`methods.py`, un module `method_*.py` par méthode).
- **Avis de statut PERCIST** (aide à la lecture, pas un avis médical), calculé sur les lésions suivies à chaque date (et sur la séance affichée si le suivi y est appliqué sans être enregistré) : cible = lésion retenue la plus active (SULpeak ; SULmax si elle fait moins de 1 ml) ; RMC (cibles ≤ moyenne du foie), RMP (≥ 30 % et ≥ 0,8 de baisse), PMP (≥ 30 % et ≥ 0,8 de hausse, nouvelle lésion retenue par l'utilisateur, ou volume d'une lésion ≥ +75 % et ≥ +1 ml), sinon MMS. Lésions « excrétion possible » ou au contact du cœur, des reins, de la vessie ou du cerveau : exclues et listées. Avertissements : conditions d'acquisition différentes ou illisibles, SUL de repli, référence foie absente ou variable (> 20 %), cible sous le seuil de mesurabilité, volume « limité à la fenêtre ». Nadir de la cible (plus basse valeur depuis la référence incluse) et % vs nadir par date, dans le texte et le tableau ; une remontée ≥ 30 % et ≥ 0,8 depuis le nadir est signalée (« remontée depuis le nadir », avertissement) sans changer le statut officiel, calculé vs référence.
- **Réponse hétérogène** (PERCIST et RECIST, aide à la lecture) : à chaque date, chaque lésion de référence retenue (mêmes exclusions) est classée en baisse / stable / hausse d'après son SULpeak (SULmax sous 1 ml) vs référence, seuils PERCIST (≥ 30 % et ≥ 0,8) ; « Lésions : x en baisse, y stables, z en hausse » (et non évaluables), ⚠ « réponse mixte » si au moins une baisse et une hausse à la même date.
- **Avis de statut RECIST 1.1** (aide à la lecture, pas un avis médical) : à chaque mesure d'une lésion suivie, son diamètre est mesuré sur le CT de la séance, dans le plan axial, à la position de la lésion (région par tolérance ± 50 HU ; plus grand diamètre et petit axe), et enregistré avec la date (nombres seulement). Cibles choisies à la référence parmi les lésions retenues, plus grandes d'abord : 5 au plus, 2 par organe ; mesurable si ≥ 10 mm (2 × l'épaisseur de coupe si elle dépasse 5 mm) ; ganglion (bouton « G » de la carte) : petit axe ≥ 15 mm (10–15 mm : non cible). Somme des diamètres par date, % vs référence et vs nadir : RC (cibles disparues, ganglions < 10 mm), RP (≤ −30 % vs référence), PE (≥ +20 % et ≥ +5 mm vs nadir, nouvelle lésion retenue, ou progression sans équivoque des non-cibles cochée par le lecteur), sinon MS ; NE si un diamètre manque. Mêmes exclusions que PERCIST (à une date au moins).

### Organes et export
- **Segmentation d'organes** (approximative) : poumons (sans trachée), cœur, foie, rate (sans la maladie périsplénique hyperfixante), reins, vessie, cerveau, squelette, aorte thoracique descendante (référence pool sanguin PERCIST affichée si le foie est atteint) et liquide libre (ascite, approximatif), par règles CT + TEP, hors ligne. Seuls cœur, reins, vessie et cerveau sont exclus des mesures (marge 10 mm, 5 mm autour du cœur).
- **Export PDF** : coupes avec mesures, infos patient, liste des mesures. Le PDF contient l'identité du patient : à partager avec précaution.

### Tactile (tablette)
- **Appui long** : loupe ×4 qui suit le doigt.
- **Pincement** : zoom et déplacement à deux doigts.
- **Glisser-déposer natif désactivé** : pas de gêne au toucher.

### Confidentialité
- **Local uniquement** : serveur sur 127.0.0.1, contrôle de l'hôte, pas de journal des requêtes.
- **Seule exception au « rien sur disque »** : un dossier unique (`~/.netscan`, à côté de `PETSCAN_CONFIG` sur Android) avec les cartes de référence pseudonymisées (empreinte salée de l'identifiant + naissance, jamais le nom, l'identifiant ni la naissance) et les réglages ; dossiers récents désactivés par défaut ; tout est listé et supprimable (aperçu, confirmation) dans le panneau « Données enregistrées », favoris du navigateur compris.

## Construire soi-même

`build/`, `dist/` et `*.spec` sont ignorés par git.

**Windows** (icône `assets/netscan.ico`) :

```
pip install -e . pyinstaller
python -m PyInstaller --noconfirm --clean --onedir --name PetScan --icon assets/netscan.ico --add-data "index.html;." viewer.py
```

Résultat : `dist/PetScan/` (à garder entier). Le build est propre à l'OS ; à refaire après toute modification.

**Android** : projet Gradle dans `android/` (serveur Python embarqué via [Chaquopy](https://chaquo.com/chaquopy/), Python 3.10, WebView, `arm64-v8a` et `x86_64`, minSdk 26). Prérequis : JDK 17 ou 21 (pas 25), SDK Android (platform 36), Python 3.10, et `android/local.properties` (non commité) :

```
sdk.dir=C:/Users/<vous>/AppData/Local/Android/Sdk
```

```
$env:JAVA_HOME = "C:\chemin\vers\jdk-21"
cd android
.\gradlew.bat assembleDebug
adb install -r app\build\outputs\apk\debug\app-debug.apk
```

À la première ouverture, l'appli demande l'accès à tous les fichiers (le serveur lit les DICOM par chemin). Copier l'examen dans un dossier de la galerie, puis toucher le champ de la séance.

## Limites

- **Android** : un examen réel demande beaucoup de RAM ; lecture DICOM séquentielle (plus lente) ; export PDF non testé dans la WebView ; APK de debug non signé ; carte SD non testée.
- **Mesure** : le seuil fixe capte aussi les fixations physiologiques (cerveau, cœur, vessie, reins) ; segmentation d'organes heuristique.
- **Comparaison** : translation seule (pas de rotation ni de déformation) ; SUV comparables seulement si dose, poids, incorporation et machine sont proches. Comparaison A/B : grilles TEP d'espacement voisin (±10 %) exigées. Suivi multi-dates : grilles d'espacements différents acceptées (ex. 1,37 → 1,56 mm) ; le CT de référence est rééchantillonné (linéaire, même origine DICOM) pour le recalage sur la zone commune des deux CT, les positions sont converties en mm et les SUV mesurées sur la grille native de chaque séance (jamais interpolées) ; un avertissement le signale, car SUVmax et volume dépendent de la taille du voxel (volume partiel).
- **Avis de statut** : PERCIST simplifié (référence = première date, pas de nadir ; croissance jugée sur le volume à 50 % du max local, pas sur le TLG) ; les nouvelles fixations ne comptent que si vous les avez ajoutées à la carte ; dates enregistrées avant cette version sans conditions d'acquisition (comparabilité non vérifiée). RECIST 1.1 : mesures sur le CT basse dose sans injection de la TEP-CT (non diagnostiques), diamètre par tolérance HU approximatif (région fuyante, non bornée ou de densité non tissulaire = lésion non mesurable, jamais de valeur inventée ; une lésion disparue fait souvent « fuir » la région, la RC est donc rarement conclue), cible non ganglionnaire considérée disparue sous 5 mm, aucun ganglion segmenté (à marquer à la main), dates enregistrées avant cette version « diamètre non enregistré » (NE).
- **Absent** : DICOM multi-séries.
- Ne committez jamais de DICOM (`.gitignore` les exclut).

## Licence

[MIT](LICENSE). Le logiciel est fourni tel quel, sans garantie ; ce n'est pas un dispositif médical.
