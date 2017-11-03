#! /usr/bin/env python

import unittest
import os
import sys
import logging
import random

from fragcoalsim import frag
from fragcoalsim import stats
from fragcoalsim.test import TestLevel
from fragcoalsim.test.support.fragcoalsim_test_case import FragcoalsimTestCase

_LOG = logging.getLogger(__name__)
GLOBAL_RNG = random.Random()


class FragmentationModelTestCase(unittest.TestCase):

    def test_expected_divergence(self):
        gens_since_frag = 20000.0
        pop_size = 10000
        mutation_rate = 1e-6
        expected_div = frag.get_expected_divergence_between_fragments(
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate)
        fm = frag.FragmentationModel(
                seed = 12345,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 1,
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_fragment = 1.0, # this shouldn't matter
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate,
                migration_rate = 0.0)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi(200000))
        self.assertAlmostEqual(expected_div, pi_summary.mean, places = 3)

if __name__ == '__main__':
    unittest.main()

