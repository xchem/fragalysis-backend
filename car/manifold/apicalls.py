import requests
from ratelimit import limits, sleep_and_retry
import os

api_key = os.environ["MANIFOLD_API_KEY"]


@sleep_and_retry
@limits(calls=100, period=60)
def getManifoldRetrosynthesis(smiles: str):
    """Call Manifold API to search for a retrosynthesis for a given smiles

    Parameters
    ----------
    smiles: str
        SMILES for Manifold retrosynthesis search
    """

    data = {
        "smiles": smiles,
        "maxLeadTimeWeeks": 12,
        "maxSearchDepth": 3,
        "maxNumRoutesToReturn": 10,
    }

    response = requests.post(
        url="https://api.postera.ai/api/v1/retrosynthesis/",
        headers={
            "X-API-KEY": api_key,
        },
        json=data,
    )
    return response.json()


def getexactsearch(smiles):
    response = requests.post(
        "https://api.postera.ai/api/v1/exact/",
        headers={
            "X-API-KEY": api_key,
        },
        json={"smiles": smiles},
    )

    return response.json()
