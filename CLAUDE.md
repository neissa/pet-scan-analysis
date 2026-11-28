# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Visionneuse web locale pour études TEP/CT (FDG) au format DICOM. Python seul (numpy, scipy, pydicom, matplotlib, pillow) ; le front est un unique `index.html` en JS natif, sans build. Le README (en français) liste les fonctions ; l'UI et les commentaires sont en français.

## Commandes

```
pip install -e .
python viewer.py "<dossier DICOM>"   # http://127.0.0.1:8000 ; sans argument : choix du dossier dans la page
python viewer.py --demo              # volume synthétique, sans DICOM
python petscan.py <dossier> --thr 2.5 --mip mip.png   # CLI : SUV + lésions
pytest -q                            # tous les tests (à rejouer aussi sous Python 3.10 + numpy 1.23 / scipy 1.8.1 / matplotlib 3.6, versions de l'APK ; interpréteur `buildPython` de android/local.properties) (~10 s, fantômes synthétiques, pas de DICOM)
pytest -q test_viewer.py -k autoalign   # un seul test
python -m PyInstaller --noconfirm --clean --onedir --name PetScan --icon assets/netscan.ico --add-data "index.html;." viewer.py   # app Windows dans dist/PetScan (voir README)
cd android; $env:JAVA_HOME="<JDK 17/21>"; .\gradlew.bat assembleDebug   # APK Android (Chaquopy, Python 3.10) : android/app/build/outputs/apk/debug (voir README)
```

Lancer depuis la racine : `viewer.py` relit `index.html` à chaque requête (recharger la page suffit pour le front), mais tout changement Python demande de relancer le serveur.

## Architecture

- `petscan.py` : lecture DICOM (modalité `PT`, fichiers sans extension, `_scan` ne lit que les en-têtes), SUV poids corporel avec correction de décroissance, CT rééchantillonné sur la grille du TEP (`load_ct`, choix de série), `analyze` (lésions : SUVmax/mean, MTV, TLG, pic), `measure_region` / `measure_tolerance` (mesure auto TEP ou CT), `_describe` (plus grand diamètre + étendue perpendiculaire), `ct_diameter` (diamètres RECIST au CT à une position, `{"why"}` si non fiable).
- `methods.py` : registre des méthodes d'évolution (`method_percist.py`, `method_recist.py`) ; `/methods`, `/avis`. RECIST lit le diamètre CT enregistré par lésion et par date sous `values[id]["ct"]` (absent = « diamètre non enregistré », None = pas de CT) ; drapeaux `lesion["ganglion"]` et `dates[d]["pe_non_cible"]` posés par `/ref/update`.
- `viewer.py` : tout le serveur (`ThreadingHTTPServer`, routes GET dans `do_GET`, seul POST = `/export.pdf`) et le rendu. L'objet `State` porte le volume courant sous forme de tuple remplacé atomiquement `vol=(suv, spacing, ctv, patient, folder)` ; chargements, organes, MIP rotatifs, CT et séance B tournent en threads d'arrière-plan avec leurs caches. `render_study` compose A/B/curseur/différence.
- `compare.py` : recalage de deux séances par **translation entière** (corrélation du squelette au CT) ; `refine_shift` (auto-alignement local ou global) part des décalages manuels. B est rééchantillonnée sur la grille de A. `resample` change l'espacement à origine fixe : le suivi (`_follow`/`track`, `ref_spacing`) y ramène le CT de référence quand les grilles diffèrent (`same_grid`, ±0,5 %), convertit les positions mémorisées (indices de la grille de référence) en mm et mesure le TEP sur sa grille native.
- `methods.py` : **registre des méthodes d'analyse** d'évolution enfichables par type de cancer (`CANCERS`, `Method` : id, nom, données requises, cancers, `compare(reference, mesures_par_date, contexte) -> {statut, details, avertissements}`, défaut par type). Une méthode = un module `method_*.py` (dict `METHOD`) importé explicitement dans `methods.py` (PyInstaller / copie Android des `*.py` de la racine) ; `method_percist.py` (statut sur les lésions suivies, jamais un MTV à seuil global), `method_recist.py` déclaré seulement (`compare=None`, « bientôt »). viewer.py (`/methods`, `/avis`, `POST /ref/update` op `method`) et index.html n'en codent aucune ; le choix est stocké dans la carte du patient (`analyse`).
- `organs.py` : segmentation d'organes heuristique (règles HU + SUV), approximative, réglée sur un CT injecté ; `trace_ureter` estime les uretères (étiquettes 13/14, Dijkstra CT + TEP), que `petscan.ureter_paths` change en polylignes pour le drapeau « excrétion possible » de `measure_at` (droite rein → vessie à défaut).
- `export_pdf.py` : PDF A4 construit en mémoire à partir de l'état de la page.
- `patient.py` : infos patient lues dans les en-têtes.
- `store.py` : **seule écriture sur disque** (exception encadrée) : dossier unique `~/.netscan` (ou à côté de `PETSCAN_CONFIG`), écriture atomique, cartes de lésions de référence par patient pseudonymisé (SHA-256 salé de identifiant + naissance, jamais l'identité en clair), réglages et récents (désactivés par défaut), tout supprimable depuis le panneau « Données enregistrées ». Le suivi (`viewer.track`, `petscan.measure_at`) mesure aux positions mémorisées sans seuil global.

## Conventions à ne pas casser

- Axes des volumes : `(z, y, x)` partout ; espacement `(z, y, x)` en mm.
- Décalages manuels : `out[z,y,x] = vol[z+dz, y+dy, x+dx]` côté serveur ; la valeur affichée dans l'UI est l'opposé (`-d`). Les index de coupe passés en `c=` sont des indices, pas des décalages (`_off(q, key, lim)`).
- Sécurité (données patient) : écoute sur 127.0.0.1 uniquement, contrôle de l'en-tête `Host` (anti DNS-rebinding), `Sec-Fetch-Site` exigé sur les routes qui changent l'état, `Cache-Control: no-store`, pas de journal des requêtes. Une erreur non gérée dans `do_GET` renvoie du HTML (`send_error`), que le front ne sait pas lire comme du JSON : préférer `_json({"error": …})`.
- APK Android : le code Python doit rester compatible Python 3.10, numpy 1.23, scipy 1.8.1, matplotlib 3.6 (versions fournies par Chaquopy).
- Ne jamais committer de DICOM, PNG ni PDF (`.gitignore`), sauf les icônes de la marque (`assets/`, `android/app/src/main/res/`), ni de chemins ou identités réels dans les tests ou messages de commit ; les tests utilisent des fantômes synthétiques.

## Routeur de travail (modèles et agents)

- **Planification** : la session principale planifie avec **Sonnet, effort moyen**. Elle ne code pas elle-même les tâches indépendantes : elle les confie à des agents.
- **Exécution** : chaque agent travaille dans **son propre worktree** (`isolation: "worktree"`), avec **Haiku, effort élevé**. Deux agents ne touchent jamais les mêmes fichiers.
- **Escalade** : si le sujet est **complexe** (plusieurs modules, algorithme délicat, conception) ou après **2 échecs** sur le même objectif, la tâche passe à **Opus, effort moyen** (agent ou session).
- **Boucle de base de tout plan** :
  1. **Analyser** : lire le code concerné et lister tout ce que le changement touche (appelants, tests, front).
  2. **Fixer un objectif testable** : un critère vérifiable (un test qui échoue puis passe, une réponse d'API, un comportement observé dans la page).
  3. **Coder** le minimum nécessaire.
  4. **Tester** (`pytest -q`, et vérification dans le navigateur pour le front).
  5. **Recommencer** si l'objectif n'est pas atteint ; compter les échecs pour déclencher l'escalade.
- **Worktrees** : `isolation: "worktree"` fonctionne depuis un dépôt à chemin court, par ex. `C:\git\pet-scan-analysis` (worktrees créés dans `.claude/worktrees/agent-*`, branche `worktree-agent-*`). Il échouait quand le dépôt était dans un chemin trop long (`'$GIT_DIR' too big`) : garder le dépôt à un chemin court ; à défaut, `git worktree add -b <branche> C:/wt/<nom>` et un brief avec ce chemin absolu. L'agent commite sans pousser.
- **Ménage obligatoire** une fois le travail d'un worktree fini et repris (merge ou cherry-pick) : vérifier le commit de l'agent, puis `git worktree remove <chemin>`, `git branch -D <branche>` (y compris les `worktree-agent-*` orphelines), `git worktree prune`, et contrôler avec `git worktree list` et `git branch`.
