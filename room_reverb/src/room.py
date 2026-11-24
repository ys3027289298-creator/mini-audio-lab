"""Room geometry, material zones and source/receiver layout.

The room is a rectangular box.  Every face belongs to a material zone, and a
zone may cover several faces.  Absorption coefficients are frequency dependent
(one value per octave band, 125 Hz ... 4 kHz), which is what makes a hard
corner sound darker/muffler than the centre: low-frequency energy survives
while high-frequency energy is absorbed by the foam wall.
"""

from dataclasses import dataclass

import numpy as np

# Octave bands analysed in the experiment, Hz.
BAND_FREQS = np.array([125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0])

# Frequency-dependent random-incidence absorption coefficients.
# Values follow the order of BAND_FREQS.
MATERIALS = {
    # Bare concrete / brick, hard and reflective.
    "concrete": {"alpha": np.array([0.02, 0.03, 0.04, 0.05, 0.06, 0.08])},
    # Wooden wall, slightly more absorbent than concrete.
    "wood": {"alpha": np.array([0.10, 0.12, 0.14, 0.16, 0.18, 0.22])},
    # Acoustic foam: weak at low frequencies, very absorbent above 1 kHz.
    "foam": {"alpha": np.array([0.15, 0.25, 0.45, 0.70, 0.90, 0.95])},
}

SPEED_OF_SOUND = 343.0  # m/s


@dataclass
class Face:
    """One rectangular room face.

    zone:    material zone name
    center:  face centre position, m
    area:    face area, m^2
    normal:  outward unit normal (informational)
    axis:    which room axis the face is perpendicular to (0=x, 1=y, 2=z)
    sign:    -1 for the low-coordinate face, +1 for the high-coordinate face
    """

    zone: str
    center: np.ndarray
    area: float
    normal: np.ndarray
    axis: int
    sign: int


@dataclass
class Room:
    size: np.ndarray  # (Lx, Ly, Lz), m
    faces: list

    @property
    def volume(self) -> float:
        return float(np.prod(self.size))

    @property
    def surface_area(self) -> float:
        return float(sum(f.area for f in self.faces))

    def zones(self):
        """Return {zone: {"area": ..., "center": area-weighted centroid}}."""
        out = {}
        for f in self.faces:
            z = out.setdefault(f.zone, {"area": 0.0, "moment": np.zeros(3)})
            z["area"] += f.area
            z["moment"] += f.area * f.center
        for z in out.values():
            z["center"] = z["moment"] / z["area"]
            del z["moment"]
        return out


def build_room() -> Room:
    """8 x 6 x 3 m room with three differently treated wall groups.

    Zone layout (three distinct materials on the enclosing surfaces):
      foam      : the west wall (x = 0) alone
      wood      : the two side walls (y = 0 and y = Ly)
      concrete  : east wall, floor and ceiling
    """
    Lx, Ly, Lz = 8.0, 6.0, 3.0
    size = np.array([Lx, Ly, Lz])
    faces = []

    def add(zone, axis, sign, center, area, normal):
        faces.append(Face(zone, np.array(center, dtype=float), area,
                          np.array(normal, dtype=float), axis, sign))

    # x faces: west = foam, east = concrete
    add("foam", 0, -1, (0.0, Ly / 2, Lz / 2), Ly * Lz, (-1, 0, 0))
    add("concrete", 0, +1, (Lx, Ly / 2, Lz / 2), Ly * Lz, (+1, 0, 0))
    # y faces: both wood
    for y, n in ((0.0, (0, -1, 0)), (Ly, (0, +1, 0))):
        add("wood", 1, -1 if y == 0.0 else 1,
            (Lx / 2, y, Lz / 2), Lx * Lz, n)
    # z faces: floor and ceiling, both concrete
    for z, n in ((0.0, (0, 0, -1)), (Lz, (0, 0, +1))):
        add("concrete", 2, -1 if z == 0.0 else 1,
            (Lx / 2, Ly / 2, z), Lx * Ly, n)

    return Room(size, faces)


def build_layout(room: Room):
    """Source and the three requested receiver positions.

    Source sits near the foam wall so source-to-zone distances are visibly
    non-uniform.  Receivers: room centre, a concrete/wood hard corner, and
    right in front of the foam-covered wall.
    """
    Lx, Ly, Lz = room.size
    source = np.array([0.7, Ly / 2, Lz / 2])

    receivers = {
        "center": np.array([Lx / 2, Ly / 2, Lz / 2]),
        # East wall (concrete) + north wall (wood) + ceiling (concrete) corner.
        "hard_corner": np.array([Lx - 0.3, Ly - 0.3, Lz - 0.3]),
        # 0.3 m in front of the foam wall, mid-span.
        "foam_wall": np.array([0.3, Ly / 2, Lz / 2]),
    }
    return source, receivers
