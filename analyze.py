import sys
from math import sqrt, sin, cos, atan2, inf

# Filament mass density constant.
FILAMENT_MASS_KG_PER_E = 1.241238 / 1000000.0
NOZZLE_SIZE_MM = 0.4
INITIAL_LAYER_HEIGHT_MM = 0.3
LAYER_HEIGHT_MM = 0.1

def read_gcode_file(filename):
    """Reads elements from a RepRap .gcode file into a list."""
    elements = []
    x0, y0, z0, e0 = 0.0, 0.0, 0.0, 0.0

    ignore = False

    for line in open(filename):
        line = line.rstrip()

        # skip empty lines
        if not line:
            continue

        # eliminate sections that are not part of the finished product
        if line[:7] == ";LAYER:":
            ignore = False
        if line[:6] == ";TYPE:":
            ignore = line[6:] in ("SKIRT", "SUPPORT")

        # skip comments
        if line[0] == ";":
            continue

        # parse code, skip everything but G0, G1 and G92
        code, *args = line.split(" ")
        if code not in ("G0", "G1", "G92"):
            continue

        # default values for a coordinate, if not specified, is to use the
        # current coordinate
        x1, y1, z1, e1 = x0, y0, z0, e0

        # parse arguments from code
        for arg in args:
            if arg[0] == "X":
                x1 = float(arg[1:])
            if arg[0] == "Y":
                y1 = float(arg[1:])
            if arg[0] == "Z":
                z1 = float(arg[1:])
            if arg[0] == "E":
                e1 = float(arg[1:])

        # we only want extrusion moves from non-ignored sections
        if not ignore and e1 > e0:
            elements.append((x0,y0,x1,y1,z0,e1-e0))

        # remember the current position
        x0, y0, z0, e0 = x1, y1, z1, e1

    return elements

def compute_bounding_box(elements):
    xmin, ymin, zmin = inf, inf, inf
    xmax, ymax, zmax = -inf, -inf, -inf

    r = NOZZLE_SIZE_MM / 2.0

    for x0, y0, x1, y1, z, e in elements:
        # layer thickness may be different for the first layer
        if z <= INITIAL_LAYER_HEIGHT_MM:
            h = INITIAL_LAYER_HEIGHT_MM
        else:
            h = LAYER_HEIGHT_MM
 
        xmin = min(x0-r, x1-r, xmin)
        ymin = min(y0-r, y1-r, ymin)
        xmax = max(x0+r, x1+r, xmax)
        ymax = max(y0+r, y1+r, ymax)
        zmin = min(z-h, zmin)
        zmax = max(z, zmax)

    return ((xmin, ymin, zmin), (xmax, ymax, zmax))

def compute_center_of_mass(elements):
    """Compute center of mass for the input elements."""
    cx, cy, cz = 0.0, 0.0, 0.0
    mass = 0.0

    for x0, y0, x1, y1, z, e in elements:
        # mass of filament
        m = FILAMENT_MASS_KG_PER_E * e

        # height of the extruded filament on the current layer
        if z <= INITIAL_LAYER_HEIGHT_MM:
            h = INITIAL_LAYER_HEIGHT_MM
        else:
            h = LAYER_HEIGHT_MM

        # add the contribution of the filament to center of mass
        cx += m * (x0 + x1) / 2
        cy += m * (y0 + y1) / 2
        cz += m * (z - h / 2)

        # accumulate total mass
        mass += m

    return (cx/mass, cy/mass, cz/mass, mass)

def compute_moments_of_inertia(elements, center):
    """Compute moments of inertia about the X/Y/Z axes through a given point."""
    ixx, iyy, izz = 0.0, 0.0, 0.0
    cx, cy, cz = center

    # width of extruded filament
    w = NOZZLE_SIZE_MM

    for x0, y0, x1, y1, z, e in elements:
        # mass of extruded filament
        m = FILAMENT_MASS_KG_PER_E * e
        # length of extruded filament
        l = sqrt((x1 - x0)**2 + (y1 - y0)**2)
        # angle of extrusion with respect to the x axis
        a = atan2(y1 - y0, x1 - x0)

        # height of the extruded filament on the current layer
        if z <= INITIAL_LAYER_HEIGHT_MM:
            h = INITIAL_LAYER_HEIGHT_MM
        else:
            h = LAYER_HEIGHT_MM

        # filament moment of inertia subterms
        i1 = (1/12) * m * h * h
        i2 = (1/12) * m * l * l
        i3 = (1/16) * m * l * w
        i4 = (1/12) * m * w * w

        # filament moment of inertia about principal axes
        ix = i1 + i2*sin(a)**2 + i3*sin(a)*cos(a) + i4*cos(a)**2
        iy = i1 + i2*cos(a)**2 + i3*cos(a)*sin(a) + i4*sin(a)**2
        iz = i2 + i4
 
        # squared distances along axes from center to filament midpoint
        fxfx = ((x0 + x1) / 2 - cx) ** 2
        fyfy = ((y0 + y1) / 2 - cy) ** 2
        fzfz = (z - h/2 - cz) ** 2

        # add the contribution of the filament using the parallel axis theorem
        ixx += ix + m * (fyfy + fzfz)
        iyy += iy + m * (fxfx + fzfz)
        izz += iz + m * (fxfx + fyfy)

    # millimeters squared to meters squared
    ixx /= 1000000.0
    iyy /= 1000000.0
    izz /= 1000000.0

    return (ixx, iyy, izz)

def usage(program):
    print("usage: %s <gcode>" % program)
    return 255

def main(filename):
    elements = read_gcode_file(filename)
    print("Read %s, %d elements parsed." % (filename, len(elements)))

    (xmin, ymin, zmin), (xmax, ymax, zmax) = compute_bounding_box(elements)

    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2
    zmid = (zmin + zmax) / 2

    cx, cy, cz, m = compute_center_of_mass(elements)

    ixx, iyy, izz = compute_moments_of_inertia(elements, (cx, cy, cz))

    print()
    print("Total mass (g):")
    print("\tM = %.4f" % (m * 1000.0))

    print()
    print("Bounding box dimensions and middle point on the platform (mm):")
    print("\tX = %-15.2f Xm = %-15.2f" % (xmax - xmin, xmid))
    print("\tY = %-15.2f Ym = %-15.2f" % (ymax - ymin, ymid))
    print("\tZ = %-15.2f Zm = %-15.2f" % (zmax - zmin, zmid))

    print()
    print("Center of mass on platform and offset from center of bounding box:")
    print("\tCx = %-15.2f dCx = %+-15.2f" % (cx, cx - xmid))
    print("\tCy = %-15.2f dCy = %+-15.2f" % (cy, cy - ymid))
    print("\tCz = %-15.2f dCz = %+-15.2f" % (cz, cz - zmid))

    print()
    print("Moments of inertia for principal axes through the center of mass:")
    print("\tIxx = %-15.6e Ixx/m = %-15.6e" % (ixx, ixx/m))
    print("\tIyy = %-15.6e Iyy/m = %-15.6e" % (iyy, iyy/m))
    print("\tIzz = %-15.6e Izz/m = %-15.6e" % (izz, izz/m))

if __name__ == '__main__':
    try:
        rc = main(*sys.argv[1:]) or 0
    except TypeError:
        rc = usage(sys.argv[0])
    except KeyboardInterrupt:
        rc = 0
    finally:
        sys.exit(rc)
