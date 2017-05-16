import click
import numpy as np
from lxml import etree

import splipy.curve_factory as cf
import splipy.surface_factory as sf
from splipy.IO import G2


@click.command()
@click.option('--width', default=1.0)
@click.option('--height', default=1.0)
@click.option('--radius', default=0.2)
@click.option('--inner-radius', default=0.4)
@click.option('--nel-ang', default=14)
@click.option('--order', default=4)
@click.option('--out', default='out')
def cut_square(width, height, radius, inner_radius, nel_ang, order, out):

    # Compute number of elements along each part
    rest1 = height - radius - inner_radius
    rest2 = width - radius - inner_radius
    nel_cyl = int(np.ceil(4/np.pi * (inner_radius - radius) / radius * nel_ang))
    nel_rest1 = int(np.ceil(4/np.pi * rest1 / radius * nel_ang))
    nel_rest2 = int(np.ceil(4/np.pi * rest2 / radius * nel_ang))

    # Create quarter circles
    theta = np.linspace(0, np.pi/2, 2*nel_ang+1)[::-1]
    pts = np.array([radius * np.cos(theta), radius * np.sin(theta)]).T
    circle = cf.cubic_curve(pts, boundary=cf.Boundary.NATURAL).set_dimension(3)
    knots = circle.knots('u')
    inner1, inner2 = circle.split(knots[len(knots) // 2])

    # Fill the cylinder patches
    factor = inner_radius / radius
    outer1, outer2 = inner1 * factor, inner2 * factor
    cyl1 = sf.edge_curves(inner1, outer1).set_order(4,4).refine(0, nel_cyl-1)
    cyl2 = sf.edge_curves(inner2, outer2).set_order(4,4).refine(0, nel_cyl-1)

    # Create the "curved rectangles"
    dist = np.sqrt(2) * radius
    edge1 = cf.line((0, height), (dist, height)).set_order(4).set_dimension(3).refine(nel_ang-1)
    rect1 = sf.edge_curves(outer1, edge1).set_order(4,4).refine(0, nel_rest1-1)
    edge2 = cf.line((width, dist), (width, 0)).set_order(4).set_dimension(3).refine(nel_ang-1)
    rect2 = sf.edge_curves(outer2, edge2).set_order(4,4).refine(0, nel_rest2-1)

    # Final square
    edge1 = rect2.section(u=0)
    edge2 = edge1 + (0, height - dist, 0)
    rect = sf.edge_curves(edge1, edge2).set_order(4,4).refine(0, nel_rest1-1)

    diff = 4 - order
    patches = [patch.lower_order(diff, diff) for patch in [cyl1, cyl2, rect1, rect2, rect]]

    with G2(out + '.g2') as f:
        f.write(patches)

    root = etree.Element('geometry')
    etree.SubElement(root, 'patchfile').text = out + '.g2'
    topology = etree.SubElement(root, 'topology')
    for mid, sid, midx, sidx, rev in [(1,2,2,1,False), (1,3,4,3,False), (2,4,4,3,False),
                                      (3,5,2,1,False), (4,5,1,3,False)]:
        etree.SubElement(topology, 'connection').attrib.update({
            'master': str(mid), 'slave': str(sid),
            'midx': str(midx), 'sidx': str(sidx),
            'reverse': 'true' if rev else 'false',
        })

    with open(out + '.xinp', 'wb') as f:
        f.write(etree.tostring(
            root, pretty_print=True, encoding='utf-8', xml_declaration=True, standalone=False
        ))

if __name__ == '__main__':
    cut_square()
