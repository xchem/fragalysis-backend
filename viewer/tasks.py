from __future__ import absolute_import, unicode_literals
from celery import shared_task
from frag.network.query import get_full_graph


@shared_task
def get_my_graph(smiles):
    return get_full_graph(smiles)
