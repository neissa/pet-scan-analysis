"""Identification APPROXIMATIVE des principaux organes sur un volume CT + TEP, hors ligne, par règles simples
(seuils HU, position relative, fixations FDG). Ce n'est pas une segmentation validée : recherche uniquement.
Orientation supposée : patient en décubitus dorsal, tête vers les z croissants, droite du patient vers les x faibles."""
import numpy as np
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

ORGANS = {  # id: (nom, couleur)
    1: ("Poumon droit", "#4aa3ff"), 2: ("Poumon gauche", "#4ad7ff"), 3: ("Cœur (approx.)", "#ff4a4a"),
    4: ("Foie", "#c8802a"), 5: ("Rate", "#a05ad8"), 6: ("Rein droit", "#2ecc71"), 7: ("Rein gauche", "#7dff7d"),
    8: ("Vessie", "#ffe14a"), 9: ("Cerveau", "#ff8fd0"), 10: ("Squelette", "#e8e8e8"),
    11: ("Aorte thoracique", "#ff6f00"), 12: ("Liquide libre (approx.)", "#3a6dff"),
    13: ("Uretère droit (estimé)", "#e6c200"), 14: ("Uretère gauche (estimé)", "#fff07a"),
}
URETER = {6: 13, 7: 14}  # rein -> uretère estimé (ni exclu des mesures ni PHYSIO : seulement le drapeau « excrétion possible »)
# fixation FDG physiologique exclue des mesures : cœur, reins, vessie, cerveau. Ni la rate, ni l'aorte, ni le liquide :
# leur fixation est faible et une lésion à leur contact (périsplénique, ganglion para-aortique, péritoine) doit rester mesurée.
PHYSIO = (3, 6, 7, 8, 9)
# marge d'exclusion propre à un organe (mm), plafonnant celle de `analyze` : réduite autour du cœur pour ne pas avaler
# un ganglion cardio-phrénique collé au péricarde (le myocarde est au bord du masque, son débord reste marqué « physio »)
PHYSIO_MARGIN_MM = {3: 5}


def _largest(mask):
    lab, n = ndimage.label(mask)
    return lab == 1 + int(np.argmax(ndimage.sum(mask, lab, range(1, n + 1)))) if n else mask


def _body(ct):
    """Masque du corps : HU > -300, ouverture (retire la table et le bruit), plus grande composante, trous remplis par coupe."""
    m = ndimage.binary_opening(ct > -300, iterations=3)
    body = _largest(m)
    for k in np.flatnonzero(body.any(axis=(1, 2))):
        body[k] = ndimage.binary_fill_holes(body[k])
    return body


def _core(mask, soft, erode, dilate, region=None):
    """Plus grande masse compacte de `mask` : érosion (coupe les ponts fins vers les muscles voisins),
    plus grande composante, puis dilatation bornée aux tissus mous."""
    if region is not None:
        mask = mask & region
    if not mask.any():
        return mask
    core = _largest(ndimage.binary_erosion(mask, iterations=erode))
    out = ndimage.binary_dilation(core, iterations=dilate) & soft
    return out & region if region is not None else out


def segment_organs(suv, ct, spacing, step=2):
    """suv, ct : volumes [z,y,x] de même grille ; spacing (z,y,x) en mm ; step : sous-échantillonnage de calcul.
    Retourne (labels uint8 de la taille de suv, {id: {nom, couleur, volume_ml, suv_max, suv_mean, hu_mean}}).
    Les organes introuvables sont simplement absents."""
    s3 = (slice(None, None, step),) * 3
    p, c = suv[s3], ct[s3]
    sz, sy, sx = (v * step for v in spacing)
    mm = lambda d, s: int(round(d / s))  # mm -> voxels réduits
    vox_ml = sz * sy * sx / 1000
    lab = np.zeros(c.shape, np.uint8)
    nz = c.shape[0]

    def put(mask, oid):
        lab[mask & (lab == 0)] = oid

    def zr(a, b):  # tranche de z bornée au volume
        return slice(max(0, int(a)), max(0, min(nz, int(b))))

    body = _body(c)
    zs = np.flatnonzero(body.any(axis=(1, 2)))
    xs = np.flatnonzero(body.any(axis=(0, 1)))
    if not zs.size:
        return _upsample(lab, suv.shape, step), {}
    zmin, zmax = zs[0], zs[-1]
    xmid = (xs[0] + xs[-1]) / 2
    yb = np.flatnonzero(body.any(axis=(0, 2)))
    ymid = (yb[0] + yb[-1]) / 2  # l'avant du patient est vers les y faibles
    X = np.arange(c.shape[2])[None, None, :]
    Y = np.arange(c.shape[1])[None, :, None]
    # fenêtres HU larges : le CT peut être injecté (vaisseaux, reins et cavités cardiaques jusqu'à ~250-300 HU)
    bone = ndimage.binary_opening(body & (c > 250))
    soft = body & (c > 10) & (c < 280) & ~bone
    smc = ndimage.uniform_filter(c.astype(np.float32), 3)  # CT lissé : le CT ramené sur la grille TEP est bruité

    # Poumons : grosses poches d'air à l'intérieur du corps, séparées droite/gauche par la ligne médiane. Les voies
    # aériennes (pharynx, larynx, trachée) sont retirées avant : petites poches d'air près de la ligne médiane, coupe par coupe.
    air = body & (c < -320)
    L, n = ndimage.label(air, structure=[[[0] * 3] * 3, [[0, 1, 0], [1, 1, 1], [0, 1, 0]], [[0] * 3] * 3])  # 2D par coupe
    if n:
        idx = range(1, n + 1)
        area, cx = ndimage.sum(air, L, idx) * sy * sx, np.array([m[2] for m in ndimage.center_of_mass(air, L, idx)])
        air &= ~np.isin(L, 1 + np.flatnonzero((area < 800) & (np.abs(cx - xmid) * sx < 25)))
    L, n = ndimage.label(air)
    big = [i + 1 for i, v in enumerate(ndimage.sum(air, L, range(1, n + 1)) * vox_ml) if v > 300] if n else []
    lungs = np.isin(L, big) if big else np.zeros(c.shape, bool)
    has_lungs = lungs.any()
    if has_lungs:
        for side, oid in ((X < xmid, 1), (X >= xmid, 2)):
            # fermeture : vaisseaux et petits nodules dans le poumon (un remplissage des trous y mettrait la coupole du foie)
            put(ndimage.binary_closing(lungs & side, iterations=2) & body & ~bone, oid)
        lungs = (lab == 1) | (lab == 2)
        area = lungs.sum(axis=(1, 2))
        lz = np.flatnonzero(area > 0.2 * area.max())  # sans la trachée, qui prolonge le sommet jusque dans le cou
        zl0, zl1 = lz[0], lz[-1]
        Lh = zl1 - zl0

    # Cerveau : tête, partie haute, tissus mous nettement hypermétaboliques en FDG.
    z0b = zmax - mm(150, sz)
    if has_lungs:
        z0b = max(z0b, zl1 + mm(100, sz))
    if z0b < zmax:
        head = np.zeros(c.shape, bool)
        head[zr(z0b, zmax + 1)] = True
        cand = head & body & (c > -10) & (c < 100)
        if cand.any():
            thr = max(3.0, 0.3 * np.percentile(p[cand], 99.5))
            brain = _largest(ndimage.binary_closing(cand & (p >= thr), iterations=2))
            for k in np.flatnonzero(brain.any(axis=(1, 2))):
                brain[k] = ndimage.binary_fill_holes(brain[k])
            if brain.sum() * vox_ml > 100:
                put(brain, 9)

    if has_lungs:
        zone = np.zeros(c.shape, bool)
        # Cœur : entre les poumons, en avant, dans la moitié basse de leur hauteur, tissus mous compacts.
        zone[zr(zl1 - 0.85 * Lh, zl1 - 0.30 * Lh)] = True
        zl = lungs & zone
        if (zl & (X < xmid)).any() and (zl & (X >= xmid)).any():
            xr = np.nonzero((zl & (X < xmid)).any(axis=(0, 1)))[0].max()           # bord interne du poumon droit
            xl = np.nonzero((zl & (X >= xmid)).any(axis=(0, 1)))[0].mean()         # centre du poumon gauche
            yl = np.nonzero(zl.any(axis=(0, 2)))[0]
            # avant du rachis : os médian dans la moitié arrière du thorax (le sternum est devant)
            sp_ = zone & bone & (np.abs(X - xmid) < mm(15, sx)) & (Y > (yl.min() + yl.max()) / 2)
            ysp = np.nonzero(sp_.any(axis=(0, 2)))[0].min() if sp_.any() else yl.mean() + mm(10, sy)
            # Aorte descendante : tube de sang injecté (> 100 HU) à gauche et juste devant le rachis.
            ao = zone & soft & (c > 100) & (X > xmid - mm(5, sx)) & (X < xmid + mm(50, sx)) & \
                (Y > ysp - mm(30, sy)) & (Y < ysp + mm(25, sy))
            ao = _largest(ndimage.binary_opening(ao, iterations=2)) if ao.any() else ao
            za = np.flatnonzero(ao.any(axis=(1, 2)))
            if za.size and (za[-1] - za[0]) * sz > 40 and 10 < ao.sum() * vox_ml < 150:
                put(ao, 11)
            # Cœur : jusqu'à l'avant du rachis (l'aorte, déjà posée, en est retirée)
            box = np.zeros(c.shape, bool)
            box[zr(zl0 + mm(10, sz), zl1 - 0.30 * Lh)] = True  # la pointe du cœur descend parfois sous la zone
            box &= (X > min(xr - mm(10, sx), xmid - mm(45, sx))) & (X < xl + mm(20, sx)) & (Y > yl.min() + mm(5, sy)) & \
                (Y < ysp - mm(5, sy))
            put(_core(soft & ~lungs & (lab == 0), soft & (lab == 0), 3, 3, box), 3)
            # TEP en respiration libre, CT en apnée : le myocarde fixant déborde souvent de 1-2 cm du cœur du CT (vers le
            # poumon ou le foie). Myocarde nettement fixant : on rattache au cœur ce qui, contigu, fixe plus de la moitié
            # de son max (un ganglion cardio-phrénique, bien moins fixant, reste dehors)
            heart = lab == 3
            if heart.any() and p[heart].max() >= 6:
                hot = (p >= 0.5 * p[heart].max()) & body & ~bone & (lab <= 2)
                lab[ndimage.binary_dilation(heart, iterations=mm(25, sx), mask=heart | hot)] = 3
        # Foie : sous la coupole droite (qui monte bien au-dessus du bas des poumons) jusqu'au lobe gauche, au-delà de la
        # ligne médiane. Noyau = plus grosse masse compacte (CT lissé : le foie non injecté, ~55 HU, est bruité), puis
        # croissance limitée (30 mm) dans les voxels proches de sa densité médiane : rend les bords que l'érosion a coupés.
        liver_box = np.zeros(c.shape, bool)
        liver_box[zr(zl0 - mm(220, sz), zl0 + mm(50, sz))] = True
        liver_box &= (X < xmid + mm(30, sx)) & (Y < ymid + mm(80, sy))  # le lobe droit descend jusqu'à la paroi postérieure
        free = soft & (lab == 0) & liver_box
        liver = _core((smc > 45) & (smc < 170) & free, soft, 3, 3, liver_box)
        if liver.any():
            near = free & (np.abs(c - np.median(c[liver])) < 35)  # CT brut : le lissage rognerait le bord
            liver |= ndimage.binary_dilation(liver, iterations=mm(30, sx), mask=near)
            for k in np.flatnonzero(liver.any(axis=(1, 2))):
                liver[k] = ndimage.binary_fill_holes(liver[k])  # vaisseaux du foie (plus denses une fois injectés)
        put(liver, 4)
        # Rate : sous le dôme gauche, en arrière et en dehors.
        spl_box = np.zeros(c.shape, bool)
        spl_box[zr(zl0 - mm(130, sz), zl0 + mm(10, sz))] = True
        spl_box &= (X > xmid + mm(40, sx)) & (Y > ymid - mm(20, sy))
        spleen = _core((c > 40) & (c < 170) & soft & (lab == 0), soft, 2, 2, spl_box)
        if spleen.any():
            # rate injectée homogène : on retire ce qui n'a pas sa densité (estomac, côlon) ou fixe nettement plus qu'elle
            # (maladie périsplénique, à garder mesurable et hors de la rate)
            ok = (np.abs(c - np.median(c[spleen])) < 40) & (p < 1.4 * np.median(p[spleen]))
            spleen = _largest(ndimage.binary_opening(spleen & ok))
        if spleen.sum() * vox_ml > 30:
            put(spleen, 5)
        # Reins : tissus mous très captants (excrétion du FDG), en arrière, de part et d'autre du rachis.
        kz = np.zeros(c.shape, bool)
        kz[zr(zl0 - mm(280, sz), zl0)] = True
        kc = kz & body & ~bone & (lab == 0) & (p >= 3.0) & (c > 5) & (c < 300) & (Y > ymid - mm(10, sy)) & (np.abs(X - xmid) > mm(25, sx)) & (np.abs(X - xmid) < mm(130, sx))
        def kidney_of(side, band=None):
            cand = kc & side if band is None else kc & side & band
            k = _largest(ndimage.binary_opening(cand)) if cand.any() else cand
            if k.sum() * vox_ml <= 5:
                return None
            # le TEP ne donne que la zone la plus captante (parfois les seules cavités, dilatées) : on prolonge de proche en
            # proche dans le parenchyme rehaussé (> 100 HU, au-dessus du muscle) ; sans assez de parenchyme rehaussé,
            # ce n'est pas un rein (masse hyperfixante rétropéritonéale)
            k = ndimage.binary_dilation(k, iterations=mm(20, sx), mask=k | ((c > 100) & (c < 300) & ~bone & kz & side & (lab == 0)))
            return k if (k & (c > 150)).sum() * vox_ml > 30 else None

        right, left = kidney_of(X < xmid), kidney_of(X >= xmid)
        zc = lambda m: np.flatnonzero(m.any(axis=(1, 2))).mean()
        if right is not None and left is not None and abs(zc(right) - zc(left)) > mm(60, sz):
            # les deux reins sont à peu près à la même hauteur : le candidat le moins captant est refait dans la bande de l'autre
            big, small = (right, left) if p[right].max() >= p[left].max() else (left, right)
            band = np.zeros(c.shape, bool)
            band[zr(zc(big) - mm(60, sz), zc(big) + mm(60, sz))] = True
            redo = kidney_of(X >= xmid if small is left else X < xmid, band)
            right, left = (big, redo) if big is right else (redo, big)
        for k, oid in ((right, 6), (left, 7)):
            if k is not None:
                put(k, oid)

    # Vessie : dans le bas du volume, au centre, liquide (HU faible) très captant (urine marquée).
    bz = np.zeros(c.shape, bool)
    bz[zr(zmin, zmin + 0.32 * (zmax - zmin))] = True
    # L'urine a ~0-40 HU mais le CT ramené sur la grille TEP est bruité (±20 HU) : la densité ne la sépare pas d'une masse
    # pelvienne ; c'est la fixation, bien plus forte dans l'urine (SUV souvent > 30) que dans une masse, qui la sépare.
    bc = bz & body & (np.abs(X - xmid) < mm(70, sx)) & (p >= 8) & (smc > -30) & (smc < 70) & ~bone
    if bc.any():
        bl = _largest(bc)
        bl = _largest(ndimage.binary_closing(bc & (p >= max(8.0, 0.25 * p[bl].max())), iterations=1))
        if bl.sum() * vox_ml > 3:
            put(bl, 8)

    # Uretères (estimés, approximatifs) : chemin de coût minimal rein -> vessie, tube d'environ 6 mm de rayon
    zone = body & ~bone & (c > -200) & (c < 300)
    for kid, oid in URETER.items():
        path = trace_ureter(p, lab, (sz, sy, sx), kid, xmid, zone)
        if path is not None:
            tube = np.zeros(c.shape, bool)
            tube[tuple(path.T)] = True
            put(ndimage.binary_dilation(tube, iterations=int(6 // sx)) if 6 >= sx else tube, oid)

    # Liquide libre (ascite) : poches de densité liquide peu fixantes, sous les poumons, à plus de 25 mm de la peau (hors
    # paroi), hors organes. L'ouverture retire les interfaces graisse/muscle, dont la moyenne tombe aussi dans -10..25 HU.
    # Une poche de moins de 100 mL au contact du foie est la vésicule biliaire. Approximatif : une anse digestive
    # remplie de liquide peut s'y confondre ; une lame de moins d'environ 1 cm d'épaisseur n'est pas vue.
    zf = zl0 if has_lungs else nz
    free = np.zeros(c.shape, bool)
    free[:zf] = (p[:zf] < 2) & (lab[:zf] == 0) & (ndimage.distance_transform_edt(body[:zf], sampling=(sz, sy, sx)) > 25)
    L, n = ndimage.label(ndimage.binary_opening(free & (smc > -10) & (smc < 25), iterations=2))
    if n:
        vol = ndimage.sum(L > 0, L, range(1, n + 1)) * vox_ml
        gb = set(np.unique(L[ndimage.binary_dilation(lab == 4, iterations=2)]))
        keep = np.isin(L, [i + 1 for i, v in enumerate(vol) if v > 20 and not (v < 100 and i + 1 in gb)])
        # le lissage et l'ouverture rognent le bord : on le rend dans le CT brut
        put(ndimage.binary_dilation(keep, iterations=2, mask=keep | (free & (c > -30) & (c < 45))), 12)

    put(bone, 10)  # le squelette passe en dernier : il ne prend jamais la place d'un autre organe

    info = {}
    for oid, (name, color) in ORGANS.items():
        m = lab == oid
        if m.any():
            info[oid] = {"nom": name, "couleur": color, "volume_ml": float(m.sum() * vox_ml), "suv_max": float(p[m].max()),
                         "suv_mean": float(p[m].mean()), "hu_mean": float(c[m].mean())}
    return _upsample(lab, suv.shape, step), info


def trace_ureter(p, lab, spacing, kid, xmid, zone):
    """Trajet APPROXIMATIF de l'uretère du rein `kid` (6 droit, 7 gauche) à la vessie (8) : chemin de coût minimal (Dijkstra,
    26 voisins) du centre du rein au centre de la vessie, dans `zone` (corps hors os et air), hors autres organes, du côté
    du rein (la vessie exceptée). Une fixation TEP tubulaire (composante hors organes
    > 2 × le fond, longue d'au moins 40 mm, section moyenne < 3 cm² : urine dans l'uretère) coûte jusqu'à 10 fois moins ;
    sans elle, le chemin suit à peu près la droite. Retourne les voxels du chemin (N, 3) ou None (rein ou vessie absents,
    pas de chemin, chemin incohérent : plus de 2 fois la ligne droite) ; l'appelant garde alors la droite.
    shortcut: une chaîne ganglionnaire allongée et fixante passe pour un uretère ; ni psoas ni bassinet repérés."""
    if not (lab == kid).any() or not (lab == 8).any():
        return None
    sp = np.asarray(spacing, float)
    a, b = (np.round(ndimage.center_of_mass(lab == k)).astype(int) for k in (kid, 8))
    z0, z1 = max(0, min(a[0], b[0]) - 2), max(a[0], b[0]) + 3
    X = np.arange(p.shape[2])[None, None, :]
    side = (X <= xmid) if kid == 6 else (X >= xmid)  # le rachis lui-même est hors de `zone` (os)
    ok = (zone & np.isin(lab, (0, 12, kid, 8)) & (side | (lab == 8)))[z0:z1]
    a, b = a - (z0, 0, 0), b - (z0, 0, 0)
    ok[tuple(a)] = ok[tuple(b)] = True
    pb, lb = p[z0:z1], lab[z0:z1]
    bg = float(np.median(pb[ok]))
    hot = ok & np.isin(lb, (0, 12)) & (pb > 2 * bg)
    L, n = ndimage.label(hot)
    tub = np.zeros(hot.shape, bool)
    for i, sl in enumerate(ndimage.find_objects(L), 1):
        ext = max((s.stop - s.start) * d for s, d in zip(sl, sp))
        m = L[sl] == i
        if ext >= 40 and m.sum() * sp.prod() / ext < 300:
            tub[sl] |= m
    cost = np.where(tub, 1 / np.clip(pb / max(bg, 1e-3), 1, 10), 1.0)
    idx = np.full(ok.shape, -1)
    idx[ok] = np.arange(int(ok.sum()))
    rows, cols, w = [], [], []
    for d in np.ndindex(3, 3, 3):
        d = np.array(d) - 1
        if tuple(d) <= (0, 0, 0):  # chaque arête une fois (graphe non orienté)
            continue
        src = tuple(slice(max(0, -k), n - max(0, k)) for k, n in zip(d, ok.shape))
        dst = tuple(slice(max(0, k), n - max(0, -k)) for k, n in zip(d, ok.shape))
        A, B = idx[src], idx[dst]
        m = (A >= 0) & (B >= 0)
        rows.append(A[m])
        cols.append(B[m])
        w.append((cost[src][m] + cost[dst][m]) / 2 * float(np.linalg.norm(d * sp)))
    N = int(ok.sum())
    G = csr_matrix((np.concatenate(w), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    ia, ib = int(idx[tuple(a)]), int(idx[tuple(b)])
    dist, pred = dijkstra(G, directed=False, indices=ia, return_predecessors=True)
    if not np.isfinite(dist[ib]):
        return None
    nodes = [ib]
    while nodes[-1] != ia:
        nodes.append(int(pred[nodes[-1]]))
    path = np.argwhere(ok)[nodes]
    if np.linalg.norm(np.diff(path, axis=0) * sp, axis=1).sum() > 2 * np.linalg.norm((b - a) * sp):
        return None
    return path + (z0, 0, 0)


def _upsample(lab, shape, step):
    full = np.repeat(np.repeat(np.repeat(lab, step, 0), step, 1), step, 2)
    return full[:shape[0], :shape[1], :shape[2]]
