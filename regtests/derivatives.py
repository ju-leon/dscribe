# -*- coding: utf-8 -*-
"""Copyright 2019 DScribe developers

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import math
import unittest

import numpy as np

import scipy
import scipy.sparse
from scipy.integrate import tplquad
from scipy.linalg import sqrtm

from dscribe.descriptors import SOAP
from testbaseclass import TestBaseClass

from ase import Atoms
from ase.build import molecule


H2O = Atoms(
    cell=[
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ],
    positions=[
        [0, 0, 0],
        [0.95, 0, 0],
        [0.95*(1+math.cos(76/180*math.pi)), 0.95*math.sin(76/180*math.pi), 0.0]
    ],
    symbols=["H", "O", "H"],
)

H = Atoms(
    cell=[
        [15.0, 0.0, 0.0],
        [0.0, 15.0, 0.0],
        [0.0, 0.0, 15.0]
    ],
    positions=[
        [0, 0, 0],

    ],
    symbols=["H"],
)

H2 = Atoms(
    cell=[
        [15.0, 0.0, 0.0],
        [0.0, 15.0, 0.0],
        [0.0, 0.0, 15.0]
    ],
    positions=[
        [-0.5, 0, 0],
        [0.5, 0, 0],

    ],
    symbols=["H", "H"],
)


class SoapDerivativesTests(TestBaseClass, unittest.TestCase):

    def test_constructor(self):
        return
    def test_number_of_features(self):
        return
    def test_sparse(self):
        return
    def test_symmetries(self):
        return
    def test_analytical(self):
        """Tests if the analytical soap derivatives run
        """
        a = SOAP(
            species=[1],
            rcut=3,
            nmax=2,
            lmax=0,
            sparse=False,
        )
        atoms = a.create(H2)

#        derivatives = a.derivatives(H2, positions =[[0.0, 0.0, 0.0], [-0.5, 0, 0], [0.5, 0, 0], ] , method = "analytical", include=None, exclude=None)
        derivatives = a.derivatives(H2, positions =[[0.0, 0.0, 0.0], ] , method = "analytical", include=None, exclude=None)

        print(derivatives)
        print(derivatives.shape)
        print(a._rcut)


if __name__ == '__main__':
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(SoapDerivativesTests))
    alltests = unittest.TestSuite(suites)
    result = unittest.TextTestRunner(verbosity=0).run(alltests)
