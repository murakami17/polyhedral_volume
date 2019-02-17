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
        if len(coords) is not 4:
            raise ValueError("Tetrahedral cluster consists of only four atoms.")
        
        ab = coords[1] - coords[0]
        ac = coords[2] - coords[0]
        ad = coords[3] - coords[0]
        
        array = numpy.array([ab, ac, ad])
        
        return abs(numpy.linalg.det(array)) / 6
    
    def octahedron(self, sites, sites_position):
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
            volume += self.tetrahedron(coords_list)
        
        return volume
    
    def cubic(self, sites, sites_position):
        """
        Calculates volume of cubic cluster from atomic coordinates.
        
        Arguments
        ---------
        sites: tuple of int
            Atomic sites which consists of cubic cluster.
        sites_position: tuple of tuple of int (0-6)
            Tuple of atomic sites which form triangle cluster in cubic.
        
        Parameters
        ----------
        coords: list of numpy.array
            Atomic coordinates of atoms which form cubic cluster.
        volume: float
            Calculated volume of cubic.
        
        Returns
        -------
        volume: float
            Calculated volume of the cubic cluster.
        """
        if len(sites) is not 8:
            raise ValueError("Cubic cluster consists of only eight atoms.")
        
        coords = []
        for site in sites:
            coords.append(self.struct.cart_coords[site])
        
        gravity = self._calc_center_of_gravity(coords)
        
        volume = 0
        for triangle in sites_position:
            coords_list = [gravity]
            for site in triangle:
                coords_list.append(coords[site])
            volume += self.tetrahedron(coords_list)
        
        return volume
    
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
            coords_sum += coord
        return coords_sum / len(coords)
    
    def _trans_periodic_coords(self, sites):
        """
        Translate periodical atomic coordinates in range -5 < r <=5.
        
        Arguments
        ---------
        
        Parameters
        ----------
        
        """
        for site in sites:
            site_dict = pymatgen.sites.PeriodicSite.as_dict(self.struct.sites[site])
            for i in range(3):
                if site_dict["abc"][i] > 0.5:
                    site_dict["abc"][i] -= 1
                elif site_dict["abc"][i] < -0.5:
                    site_dict["abc"][i] += 1
                else
                    pass
            self.struct.sites[site] = pymatgen.sites.PeriodicSite.from_dict(site_dict)
        
    
if __name__ == "__main__":
    pass
