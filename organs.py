"""Identification APPROXIMATIVE des principaux organes sur un volume CT + TEP, hors ligne, par règles simples
(seuils HU, position relative, fixations FDG). Ce n'est pas une segmentation validée : recherche uniquement.
Orientation supposée : patient en décubitus dorsal, tête vers les z croissants, droite du patient vers les x faibles."""
import numpy as np
from scipy import ndimage

ORGANS = {  # id: (nom, couleur)
    1: ("Poumon droit", "#4aa3ff"), 2: ("Poumon gauche", "#4ad7ff"), 3: ("Cœur (approx.)", "#ff4a4a"),
    4: ("Foie", "#c8802a"), 5: ("Rate", "#a05ad8"), 6: ("Rein droit", "#2ecc71"), 7: ("Rein gauche", "#7dff7d"),
    8: ("Vessie", "#ffe14a"), 9: ("Cerveau", "#ff8fd0"), 10: ("Squelette", "#e8e8e8"),
}
PHYSIO = (3, 6, 7, 8, 9)  # fixation FDG physiologique : cœur, reins, vessie, cerveau


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

    # Poumons : grosses poches d'air à l'intérieur du corps (la trachée y est rattachée), séparées droite/gauche par la ligne médiane.
    L, n = ndimage.label(body & (c < -320))
    big = [i + 1 for i, v in enumerate(ndimage.sum(body & (c < -320), L, range(1, n + 1)) * vox_ml) if v > 300] if n else []
    lungs = np.isin(L, big) if big else np.zeros(c.shape, bool)
    has_lungs = lungs.any()
    if has_lungs:
        put(lungs & (X < xmid), 1)
        put(lungs & (X >= xmid), 2)
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
            box = zone & (X > xr - mm(10, sx)) & (X < xl + mm(20, sx)) & (Y > yl.min() + mm(15, sy)) & (Y < yl.mean() + mm(10, sy))
            put(_core(soft & ~lungs, soft, 2, 2, box), 3)
        # Foie : sous la coupole droite (qui monte bien au-dessus du bas des poumons) jusqu'au lobe gauche, au-delà de la
        # ligne médiane. Noyau = plus grosse masse compacte (CT lissé : le foie non injecté, ~55 HU, est bruité), puis
        # croissance limitée (30 mm) dans les voxels proches de sa densité médiane : rend les bords que l'érosion a coupés.
        liver_box = np.zeros(c.shape, bool)
        liver_box[zr(zl0 - mm(220, sz), zl0 + mm(50, sz))] = True
        liver_box &= (X < xmid + mm(30, sx)) & (Y < ymid + mm(40, sy))
        smc = ndimage.uniform_filter(c.astype(np.float32), 3)
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
        if spleen.sum() * vox_ml > 30:
            put(spleen, 5)
        # Reins : tissus mous très captants (excrétion du FDG), en arrière, de part et d'autre du rachis.
        kz = np.zeros(c.shape, bool)
        kz[zr(zl0 - mm(280, sz), zl0)] = True
        kc = kz & body & ~bone & (lab == 0) & (p >= 3.0) & (c > 5) & (c < 300) & (Y > ymid - mm(10, sy)) & (np.abs(X - xmid) > mm(25, sx)) & (np.abs(X - xmid) < mm(130, sx))
        def kidney_of(side, band=None):
            cand = kc & side if band is None else kc & side & band
            k = _largest(ndimage.binary_opening(cand)) if cand.any() else cand
            if k.sum() * vox_ml <= 40:
                return None
            # le TEP ne donne que la zone la plus captante : on prolonge dans le parenchyme rehaussé voisin (~12 mm)
            return k | (ndimage.binary_dilation(k, iterations=mm(12, sx)) & (c > 60) & (c < 300) & ~bone & kz & side)

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
    bc = bz & body & (np.abs(X - xmid) < mm(70, sx)) & (p >= 8) & (c > -20) & (c < 40)
    if bc.any():
        bl = _largest(ndimage.binary_closing(bc, iterations=1))
        if bl.sum() * vox_ml > 8:
            put(bl, 8)

    put(bone, 10)  # le squelette passe en dernier : il ne prend jamais la place d'un autre organe

    info = {}
    for oid, (name, color) in ORGANS.items():
        m = lab == oid
        if m.any():
            info[oid] = {"nom": name, "couleur": color, "volume_ml": float(m.sum() * vox_ml), "suv_max": float(p[m].max()),
                         "suv_mean": float(p[m].mean()), "hu_mean": float(c[m].mean())}
    return _upsample(lab, suv.shape, step), info


def _upsample(lab, shape, step):
    full = np.repeat(np.repeat(np.repeat(lab, step, 0), step, 1), step, 2)
    return full[:shape[0], :shape[1], :shape[2]]
