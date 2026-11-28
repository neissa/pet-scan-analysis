"""RECIST 1.1 : déclaré seulement (non implémenté, non sélectionnable). Cibles prévues : lésions TEP, diamètre mesuré sur le CT."""
METHOD = dict(id="recist11", nom="RECIST 1.1", donnees=("TEP", "CT"), cancers=("dsrct", "general"),
              note="bientôt : diamètres au CT (5 cibles, 2 par organe, petit axe des ganglions)")
