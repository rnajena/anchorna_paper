# (C) 2025, Tom Eulenfeld, MIT license
from sugar import read


def get_24_ids():
    with open('comparison/ids_24.txt') as f:
        return f.read().split()


def get_50_ids():
    with open('comparison/ids_50.txt') as f:
        return f.read().split()


def get_55_ids():
    return read('comparison/pesti_example.gff').ids


def id_order_key():
    order = get_55_ids()
    return lambda s: order.index(s.id)
