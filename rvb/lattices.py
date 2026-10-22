"""Small lattice graphs with explicit singlet orientations.

Changing individual edge orientations changes the relative signs of RVB coverings.
Square bonds point from the even to the odd sublattice. Kagome bonds circulate
counterclockwise around both up and down triangles.
"""

from dataclasses import dataclass
from math import sqrt


@dataclass(frozen=True)
class Lattice:
    name: str
    positions: tuple[tuple[float, float], ...]
    edges: tuple[tuple[int, int], ...]
    labels: tuple[tuple[int, ...], ...]

    @property
    def nsites(self):
        return len(self.positions)


def square(nx=2, ny=2, periodic_x=False, periodic_y=False):
    if min(nx, ny) < 1:
        raise ValueError("Lattice dimensions must be positive")
    if any(p and (n < 4 or n % 2) for n, p in ((nx, periodic_x), (ny, periodic_y))):
        raise ValueError("Periodic square dimensions must be even and at least 4")
    labels = tuple((x, y) for x in range(nx) for y in range(ny))
    edges = []
    for x, y in labels:
        for dx, dy, periodic in ((1, 0, periodic_x), (0, 1, periodic_y)):
            xx, yy = x + dx, y + dy
            if (xx >= nx or yy >= ny) and not periodic:
                continue
            xx, yy = xx % nx, yy % ny
            a, b = x * ny + y, xx * ny + yy
            edges.append((a, b) if (x + y) % 2 == 0 else (b, a))
    return Lattice("square", tuple((float(x), float(y)) for x, y in labels), tuple(edges), labels)


def kagome(nx=2, ny=2, periodic_x=False, periodic_y=False):
    if min(nx, ny) < 1:
        raise ValueError("Lattice dimensions must be positive")
    if (periodic_x and nx < 2) or (periodic_y and ny < 2):
        raise ValueError("Periodic kagome dimensions must be at least 2")
    labels = tuple((x, y, s) for x in range(nx) for y in range(ny) for s in range(3))
    offsets = ((0., 0.), (1., 0.), (.5, sqrt(3) / 2))
    positions = tuple((2*x+y+offsets[s][0], sqrt(3)*y+offsets[s][1]) for x, y, s in labels)

    def site(x, y, s):
        if not periodic_x and not 0 <= x < nx:
            return None
        if not periodic_y and not 0 <= y < ny:
            return None
        return 3 * ((x % nx) * ny + y % ny) + s

    edges = []
    for x in range(nx):
        for y in range(ny):
            triangles = ((site(x, y, 0), site(x, y, 1), site(x, y, 2)),
                         (site(x, y, 0), site(x-1, y, 1), site(x, y-1, 2)))
            for triangle in triangles:
                for k in range(3):
                    a, b = triangle[k], triangle[(k+1) % 3]
                    if a is not None and b is not None:
                        edges.append((a, b))
    if len({frozenset(e) for e in edges}) != len(edges):
        raise ValueError("Geometry produced repeated edges")
    return Lattice("kagome", positions, tuple(edges), labels)
