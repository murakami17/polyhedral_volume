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
        
    def tetrahedron(self, sites):
        """
        Calculates volume of tetrahedral cluster from atomic coordinates.
        
        Arguments
        ---------
        sites: tuple of int
            Atomic sites which consists tetrahedral cluster.
        
        Parameters
        ----------
        a, b, c, d: numpy.array
            Atomic coordinates which consists tetrahedral cluster.
        ab, ac, ad: numpy.array
            Position vector of b-d atoms from a atom.
        array: numpy.array
            Array of position vectors ab-ad.
        
        Returns
        -------
        float: Volume of the tetrahedral cluster.
        """
        if len(sites) is not 4:
            raise ValueError("Tetrahedral cluster consists of only four atoms.")
        a = self.struct.cart_coords[sites[0]]
        b = self.struct.cart_coords[sites[1]]
        c = self.struct.cart_coords[sites[2]]
        d = self.struct.cart_coords[sites[3]]
        ab = b - a
        ac = c - a
        ad = d - a
        array = numpy.array([ab, ac, ad])
        return abs(numpy.linalg.det(array)) / 6
    
    def octahedron(self, sites):
        """
        Calculates volume of octahedral cluster from atomic coordinates.
        
        Arguments
        ---------
        sites: tuple of int
            Atomic sites which consists of octahedral cluster.
        
        Parameters
        ----------
        a, b, c, d, e, f: int
            Atomic sites which consists of octahedral cluster.
        tetra_a, tetra_b, tetra_c, tetra_d: float
            Volume of the tetrahedral cluster consists of octahedral cluster.
        
        Returns
        -------
        float: Volume of the octahedral cluster.
        """
        if len(sites) is not 6:
            raise ValueError("Octahedral cluster consists of only six atoms.")
        a = sites[0]
        b = sites[1]
        c = sites[2]
        d = sites[3]
        e = sites[4]
        f = sites[5]
        tetra_a = self.tetrahedron((a, b, c, d))
        tetra_b = self.tetrahedron((a, b, d, e))
        tetra_c = self.tetrahedron((f, b, c, d))
        tetra_d = self.tetrahedron((f, b, d, e))
        return tetra_a + tetra_b + tetra_c + tetra_d
    
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
        
if __name__ == "__main__":
    pass
