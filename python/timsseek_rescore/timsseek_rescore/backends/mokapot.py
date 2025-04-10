

import mokapot
import numpy as np
from rich.pretty import pprint

from ..feateng import to_mokapot


def mokapot_stuff(data, outdir):
    ds = to_mokapot(data)
    pprint("Brewing Mokapot")
    models, scores = mokapot.brew([ds])
    qvals = mokapot.qvalues.qvalues_from_scores(scores[0], ds.targets)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = np.sum(ds.targets[qvals < ct])
        ssc = scores[0][qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(
            f"Mokapot Number of targets at {ct}: {num_at_thresh};"
            f" Score: {score_at_thresh}"
        )
