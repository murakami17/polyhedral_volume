# coding: utf-8
# Copyright (c) 2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the MIT License.

import logging
import numpy
import pymatgen

"""
Calculate volume of polyhedral clusters from atomic coordinates.
"""

logger = logging.getLogger(__name__)


class PolyhedralVolume(object):
    """
    Calculate volume of polyhedral clusters from atomic coodinates.
    
    Parameters
    ----------
    struct: pymatgen.Structure
        Crystal structure which contains polyhedral clusters.
    """
    
    def __init__(self, struct):
        """
        Arguments
        ---------
        struct: pymatgen.Structure
            Crystal structure which contains polyhedral clusters.
        """
        self.struct = struct
        
    def tetrahedron(self, coords):
        """
        Calculates volume of tetrahedral cluster from atomic coordinates.
        
        Arguments
        ---------
        coords: list of numpy.array
            Atomic coordinates of atoms which form octahedral cluster.
        
        Parameters
        ----------
        ab, ac, ad: numpy.array
            Position vector of coords[3 - 1] atoms from coords[0] atom.
        array: numpy.array
            Array of position vectors ab-ad.
        
        Returns
        -------
        float: Volume of the tetrahedral cluster.
        """
        if len(sites) is not 4:
            raise ValueError("Tetrahedral cluster consists of only four atoms.")
        
        ab = coords[1] - coords[0]
        ac = coords[2] - coords[0]
        ad = coords[3] - coords[0]
        
        array = numpy.array([ab, ac, ad])
        
        return abs(numpy.linalg.det(array)) / 6
    
    def octahedron(self, sites):
        """
        Calculates volume of octahedral cluster from atomic coordinates.
        
        Arguments
        ---------
        sites: tuple of int
            Atomic sites which consists of octahedral cluster.
        sites_position: tuple of tuple of int (0-7)
            Tuple of atomic sites which form triangle cluster in octahedron.
        
        Parameters
        ----------
        coords: list of numpy.array
            Atomic coordinates of atoms which form octahedral cluster.
        volume: float
            Calculated volume of octahedron.
        
        Returns
        -------
        volume: float
            Calculated volume of the octahedral cluster.
        """
        if len(sites) is not 6:
            raise ValueError("Octahedral cluster consists of only six atoms.")
        
        coords = []
        for site in sites:
            coords.append(self.struct.cart_coords[site])
        
        gravity = self._calc_center_of_gravity(coords)
        
        volume = 0
        for triangle in sites_position:
            coords_list = [gravity]
            for site in triangle:
                coords_list.append(coords[site])
            volume += tetrahedron(coords_list)
        
        return volume
    
    def cubic(self, sites):
        """
        Calculates volume of cubic cluster from atomic coordinates.
        
        Arguments
        ---------
        sites: tuple of int
            Atomic sites which consists of octahedral cluster.
        
        Parameters
        ----------
        a, b, c, d, e, f, g, h: int
            Atomic sites which consists of octahedral cluster.
        tetra_a, tetra_b, tetra_c, tetra_d, tetra_e, tetra_f: float
            Volume of the tetrahedral cluster consists of cubic cluster.
        
        Returns
        -------
        float: Volume of the cubic cluster.
        """
        if len(sites) is not 8:
            raise ValueError("Cubic cluster consists of only eight atoms.")
        a = sites[0]
        b = sites[1]
        c = sites[2]
        d = sites[3]
        e = sites[4]
        f = sites[5]
        g = sites[6]
        h = sites[7]
        """Under construction
        tetra_a = self.tetrahedron((a, b, c, d))
        tetra_b = self.tetrahedron((d, e, f, g))
        tetra_c = self.tetrahedron((g, h, a, b))
        tetra_d = self.tetrahedron((b, c, d, e))
        tetra_e = self.tetrahedron((e, h, g, f))
        tetra_f = self.tetrahedron((f, a, b, c))
        return tetra_a + tetra_b + tetra_c + tetra_d + tetra_e + tetra_f
        """
        pass
    
    def _calc_center_of_gravity(self, coords):
        """
        Calculates center of gravity.
        
        Arguments
        ---------
        coords: list of numpy.array
            Atomic coordinates which form cluster.
        
        Parameters
        ----------
        coords_sum: numpy.array
            Summention of atomic coordinates which form cluster.
        
        Returns
        -------
        numpy.array: Coordinates of centern of gravity of cluster.
        """
        coords_sum = numpy.zeros(3)
        for coord in coords:
            coords_sum += coords
        return coords_sum / len(coords)
    
if __name__ == "__main__":
    pass
