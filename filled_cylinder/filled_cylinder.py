import click
import numpy as np
from lxml import etree

from splipy import curve_factory as cf, surface_factory as sf, volume_factory as vf
from splipy.IO import G2


@click.command()
@click.option('--radius', default=1.0)
@click.option('--length', default=2.0)
@click.option('--elements-rad', default=10)
@click.option('--elements-len', default=15)
@click.option('--out', default='out')
def cylinder(radius, length, elements_rad, elements_len, out):
    square = sf.square(size=2*radius/3, lower_left=(-radius/3, -radius/3))
    square.set_dimension(3)
    square.raise_order(2, 2)
    square.refine(elements_rad-1, elements_rad-1)

    pts = np.zeros((elements_rad+1, 3))
    angles = np.linspace(-np.pi/4, np.pi/4, elements_rad+1)
    pts[:,0] = radius * np.cos(angles)
    pts[:,1] = radius * np.sin(angles)
    curve = cf.cubic_curve(pts, t=square.knots('v'))

    sector = sf.edge_curves(square.section(u=-1), curve)
    sector.raise_order(0, 2)
    sector.refine(0, elements_rad-1)
    sector.swap()

    sectors = [sector.clone().rotate(angle) for angle in [0, np.pi/2, np.pi, np.pi*3/2]]
    patches = [vf.extrude(patch, (0, 0, length)) for patch in [square] + sectors]
    for patch in patches:
        patch.raise_order(0, 0, 2)
        patch.refine(0, 0, elements_len)

    with G2(out + '.g2') as f:
        f.write(patches)

    root = etree.Element('geometry')
    etree.SubElement(root, 'patchfile').text = out + '.g2'
    topology = etree.SubElement(root, 'topology')
    for mid, sid, midx, sidx, rev in [(1, 2, 2, 1, False), (1, 3, 4, 1, True),
                                      (1, 4, 1, 1, True),  (1, 5, 3, 1, False),
                                      (2, 3, 4, 3, False), (2, 5, 3, 4, False),
                                      (3, 4, 4, 3, False), (4, 5, 4, 3, False)]:
        etree.SubElement(topology, 'connection').attrib.update({
            'master': str(mid), 'slave': str(sid),
            'midx': str(midx), 'sidx': str(sidx),
            'reverse': 'true' if rev else 'false',
        })

    topsets = etree.SubElement(root, 'topologysets')

    for name, start, idx in [('wall', 2, 2), ('inflow', 1, 5), ('outflow', 1, 6)]:
        topset = etree.SubElement(topsets, 'set')
        topset.attrib.update({'name': name, 'type': 'face'})
        for i in range(start, 6):
            item = etree.SubElement(topset, 'item')
            item.attrib['patch'] = str(i)
            item.text = str(idx)

    with open(out + '.xinp', 'wb') as f:
        f.write(etree.tostring(
            root, pretty_print=True, encoding='utf-8', xml_declaration=True, standalone=False
        ))


if __name__ == '__main__':
    cylinder()
