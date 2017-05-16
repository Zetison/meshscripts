import click
import numpy as np

from splipy import curve_factory as cf, surface_factory as sf
from splipy.IO import G2


@click.command()
@click.option('--elements', nargs=2, default=(3, 20))
@click.option('--radius', default=4.0)
@click.option('--out', default='out')
def thingy(radius, elements, out):
    right = cf.circle_segment(np.pi/2)
    right.rotate(-np.pi/4).translate((radius-1, 0, 0))

    left = right.clone().rotate(np.pi)
    right.reverse()

    thingy = sf.edge_curves(left, right)
    thingy.raise_order(0, 1)
    thingy.refine(*[e - 1 for e in elements])

    with G2(out + '.g2') as f:
        f.write([thingy])


if __name__ == '__main__':
    thingy()
