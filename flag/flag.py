import click
import numpy as np
from math import ceil, log, sqrt, pi
import sys
from lxml import etree

from splipy import curve_factory as cf, surface_factory as sf
from splipy.utils.refinement import geometric_refine
from splipy.IO import G2


def graded_space(start, step, factor, N):
    def gen_graded_space(start, step, factor, N):
        for _ in range(N):
            yield start
            start += step
            step *= factor
    return list(gen_graded_space(float(start), float(step), float(factor), N))


def find_factor(initial, total, N, tol=1e-7):
    # Solve initial * (1 - a^N) / (1 - a) = total
    overshoot = lambda a: initial * (1 - (1 + a)**N) / (1 - (1 + a)) - total
    l, u = 0.5, 0.5
    while overshoot(l) > 0:
        l /= 2
    while overshoot(u) < 0:
        u *= 2
    while True:
        m = (u + l) / 2
        s = overshoot(m)
        if abs(s) < tol:
            break
        if s > 0:
            u = m
        else:
            l = m
    return 1 + m


@click.command()
@click.option('--diam', default=1.0)
@click.option('--flag-width', default=0.1)
@click.option('--flag-length', default=5.0)
@click.option('--width', default=20.0)
@click.option('--back', default=40.0)
@click.option('--grad', type=float, default=1.08)
@click.option('--flag-grad', type=float, default=1.05)
@click.option('--nel-rad', type=int, default=40)
@click.option('--nel-circ', type=int, default=120)
@click.option('--nel-flag', type=int, default=40)
@click.option('--order', default=4)
@click.option('--out', default='out')
def flag(diam, flag_width, flag_length, width, back,
         flag_grad, grad, nel_rad, nel_circ, nel_flag, order, out):
    assert(back > width)

    rad_cyl = diam / 2
    width *= rad_cyl
    back = back * rad_cyl - width

    # Create a circle
    angle = 2 * np.arcsin(flag_width / rad_cyl / 2)
    pts = rad_cyl * np.array([
        (np.cos(a), np.sin(a)) for a in np.linspace(angle, 2*np.pi - angle, nel_circ + 1)
    ])
    circle = cf.cubic_curve(pts, boundary=cf.Boundary.NATURAL)
    circle.set_dimension(3)

    # Subdivide it
    nels_side = int(round(nel_circ // 2 * 3 * np.pi / 4 / (np.pi - angle)))
    nels_front = (nel_circ // 2 - nels_side) * 2
    S = (- nel_flag * width + nels_side * (width + back)) / (nel_flag + nels_side)

    kts = circle.knots('u')
    kts = kts[nels_side], kts[nels_side + nels_front]
    circ_up, circ_front, circ_down = circle.split(kts)

    # Extend to boundary
    front = cf.line((-width, width), (-width, -width)).set_order(4).refine(nels_front - 1)
    front = sf.edge_curves(front, circ_front).raise_order(0, 2)
    geometric_refine(front, grad, nel_rad - 1, direction='v', reverse=True)

    up = cf.line((S, width), (-width, width)).set_order(4).refine(nels_side - 1)
    up = sf.edge_curves(up, circ_up).raise_order(0, 2)
    geometric_refine(up, grad, nel_rad - 1, direction='v', reverse=True)

    down = cf.line((-width, -width), (S, -width)).set_order(4).refine(nels_side - 1)
    down = sf.edge_curves(down, circ_down).raise_order(0, 2)
    geometric_refine(down, grad, nel_rad - 1, direction='v', reverse=True)

    # Create the flag
    upt = circle(circle.start('u'))
    fl_up = cf.line((flag_length + rad_cyl, upt[1], 0), upt).raise_order(2)
    geometric_refine(fl_up, flag_grad, nel_flag - 1, direction='u', reverse=True)
    ln_up = cf.cubic_curve(np.array([
        ((1-i)*(width+back) + i*S, width) for i in np.linspace(0, 1, nel_flag + 1)
    ]), boundary=cf.Boundary.NATURAL, t=fl_up.knots('u'))
    fl_up = sf.edge_curves(ln_up, fl_up).raise_order(0, 2)
    geometric_refine(fl_up, grad, nel_rad - 1, direction='v', reverse=True)

    dpt = circle(circle.end('u'))
    fl_down = cf.line(dpt, (flag_length + rad_cyl, dpt[1], 0))
    geometric_refine(fl_down, flag_grad, nel_flag - 1, direction='u')
    ln_down = cf.cubic_curve(np.array([
        ((1-i)*S + i*(width+back), -width) for i in np.linspace(0, 1, nel_flag + 1)
    ]), boundary=cf.Boundary.NATURAL, t=fl_down.knots('u'))
    fl_down = sf.edge_curves(ln_down, fl_down).raise_order(0, 2)
    geometric_refine(fl_down, grad, nel_rad - 1, direction='v', reverse=True)

    fl_back = cf.line((flag_length + rad_cyl, dpt[1], 0), (flag_length + rad_cyl, upt[1], 0))
    ln_back = cf.line((width + back, -width, 0), (width + back, width, 0))
    fl_back = sf.edge_curves(ln_back, fl_back).raise_order(2, 2).refine(40, direction='u')
    geometric_refine(fl_back, grad, nel_rad - 1, direction='v', reverse=True)

    with G2(out + '.g2') as f:
        f.write([up, front, down, fl_up, fl_down, fl_back])

    root = etree.Element('geometry')
    etree.SubElement(root, 'patchfile').text = out + '.g2'
    topology = etree.SubElement(root, 'topology')
    for mid, sid, midx, sidx, rev in [(1,2,2,1,False), (1,4,1,2,False), (2,3,2,1,False),
                                      (3,5,2,1,False), (4,6,1,2,False), (5,6,2,1,False)]:
        etree.SubElement(topology, 'connection').attrib.update({
            'master': str(mid), 'slave': str(sid),
            'midx': str(midx), 'sidx': str(sidx),
            'reverse': 'true' if rev else 'false',
        })

    topsets = etree.SubElement(root, 'topologysets')
    for name, index, entries in [('inflow', 3, [2]), ('outflow', 3, [6]),
                                 ('top', 3, [1,4]), ('bottom', 3, [3,5]),
                                 ('cylinder', 4, [1,2,3]), ('flag', 4, [4,5,6])]:
        topset = etree.SubElement(topsets, 'set')
        topset.attrib.update({'name': name, 'type': 'edge'})
        for pid in entries:
            item = etree.SubElement(topset, 'item')
            item.attrib['patch'] = str(pid)
            item.text = str(index)
    topset = etree.SubElement(topsets, 'set')
    topset.attrib.update({'name': 'inflow', 'type': 'vertex'})
    item = etree.SubElement(topset, 'item')
    item.attrib['patch'] = '2'
    item.text = '1 2'

    with open(out + '.xinp', 'wb') as f:
        f.write(etree.tostring(
            root, pretty_print=True, encoding='utf-8', xml_declaration=True, standalone=False
        ))

if __name__ == '__main__':
    flag()
