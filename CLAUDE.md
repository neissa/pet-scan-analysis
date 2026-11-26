# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Visionneuse web locale pour études TEP/CT (FDG) au format DICOM. Python seul (numpy, scipy, pydicom, matplotlib, pillow) ; le front est un unique `index.html` en JS natif, sans build. Le README (en français) liste les fonctions ; l'UI et les commentaires sont en français.

## Commandes

```
pip install -e .
python viewer.py "<dossier DICOM>"   # http://127.0.0.1:8000 ; sans argument : choix du dossier dans la page
python viewer.py --demo              # volume synthétique, sans DICOM
python petscan.py <dossier> --thr 2.5 --mip mip.png   # CLI : SUV + lésions
pytest -q                            # tous les tests (~10 s, fantômes synthétiques, pas de DICOM)
pytest -q test_viewer.py -k autoalign   # un seul test
```

Lancer depuis la racine : `viewer.py` relit `index.html` à chaque requête (recharger la page suffit pour le front), mais tout changement Python demande de relancer le serveur.

## Architecture

- `petscan.py` : lecture DICOM (modalité `PT`, fichiers sans extension, `_scan` ne lit que les en-têtes), SUV poids corporel avec correction de décroissance, CT rééchantillonné sur la grille du TEP (`load_ct`, choix de série), `analyze` (lésions : SUVmax/mean, MTV, TLG, pic), `measure_region` / `measure_tolerance` (mesure auto TEP ou CT), `_describe` (plus grand diamètre + étendue perpendiculaire).
- `viewer.py` : tout le serveur (`ThreadingHTTPServer`, routes GET dans `do_GET`, seul POST = `/export.pdf`) et le rendu. L'objet `State` porte le volume courant sous forme de tuple remplacé atomiquement `vol=(suv, spacing, ctv, patient, folder)` ; chargements, organes, MIP rotatifs, CT et séance B tournent en threads d'arrière-plan avec leurs caches. `render_study` compose A/B/curseur/différence.
- `compare.py` : recalage de deux séances par **translation entière** (corrélation du squelette au CT) ; `refine_shift` (auto-alignement local ou global) part des décalages manuels. B est rééchantillonnée sur la grille de A.
- `organs.py` : segmentation d'organes heuristique (règles HU + SUV), approximative, réglée sur un CT injecté.
- `export_pdf.py` : PDF A4 construit en mémoire à partir de l'état de la page.
- `patient.py` : infos patient lues dans les en-têtes.

## Conventions à ne pas casser

- Axes des volumes : `(z, y, x)` partout ; espacement `(z, y, x)` en mm.
- Décalages manuels : `out[z,y,x] = vol[z+dz, y+dy, x+dx]` côté serveur ; la valeur affichée dans l'UI est l'opposé (`-d`). Les index de coupe passés en `c=` sont des indices, pas des décalages (`_off(q, key, lim)`).
- Sécurité (données patient) : écoute sur 127.0.0.1 uniquement, contrôle de l'en-tête `Host` (anti DNS-rebinding), `Sec-Fetch-Site` exigé sur les routes qui changent l'état, `Cache-Control: no-store`, pas de journal des requêtes. Une erreur non gérée dans `do_GET` renvoie du HTML (`send_error`), que le front ne sait pas lire comme du JSON : préférer `_json({"error": …})`.
- Ne jamais committer de DICOM, PNG ni PDF (`.gitignore`) ni de chemins ou identités réels dans les tests ou messages de commit ; les tests utilisent des fantômes synthétiques.

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
- **Worktrees sous Windows** : le chemin de ce dépôt est trop long pour git (`'$GIT_DIR' too big`), donc `isolation: "worktree"` échoue et laisse une branche `worktree-agent-*` vide. Contournement testé : `git worktree add -b <branche> C:/wt/<nom>`, puis l'agent travaille uniquement dans ce dossier (chemin absolu dans son brief) et commite sans pousser.
- **Ménage obligatoire** une fois le travail d'un worktree fini et repris (merge ou cherry-pick) : vérifier le commit de l'agent, puis `git worktree remove C:/wt/<nom>`, `git branch -D <branche>` (y compris les `worktree-agent-*` orphelines), `git worktree prune`, et contrôler avec `git worktree list` et `git branch`.
