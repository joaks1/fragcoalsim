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

    def test_ms(self):
        fm = frag.FragmentationModel(
                seed = 123456,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 1000,
                effective_pop_size_of_fragment = 1000,
                effective_pop_size_of_ancestor = 10000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi, pi_a, pi_w = fm.ms_simulate(locus_length = 100, number_of_replicates = 10000)
        print(sum(pi) / len(pi))
        print(fm.expected_divergence)
        p = list(fm.sample_pi(10000))
        print(sum(p) / len(p))
        self.assertEqual(len(pi), 10000)

    def test_expected_divergence_with_one_sample(self):
        gens_since_frag = 2000.0
        pop_size = 10000
        mutation_rate = 1e-6
        expected_div, d_b, d_w = frag.get_expected_divergence(
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 1,
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_ancestor = pop_size,
                effective_pop_size_of_fragment = 1000,
                mutation_rate = mutation_rate)
        fm = frag.FragmentationModel(
                seed = 12345,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 1,
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_fragment = 1000,
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi(200000))
        print(expected_div)
        print(d_b)
        print(d_w)
        self.assertAlmostEqual(expected_div, pi_summary.mean, places = 3)
        self.assertAlmostEqual(expected_div, d_b, places = 6)
        self.assertIsTrue(d_w == float('nan'))

    def test_expected_divergence_no_frag(self):
        fm = frag.FragmentationModel(
                seed = 1111,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 0,
                effective_pop_size_of_fragment = 10000,
                effective_pop_size_of_ancestor = 100000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi(100000))
        print(fm.expected_divergence_within_fragments)
        print(fm.expected_divergence_between_fragments)
        print(fm.expected_divergence)
        print(pi_summary.mean)
        self.assertAlmostEqual(fm.expected_divergence, fm.expected_divergence_within_fragments, places = 6)
        self.assertAlmostEqual(fm.expected_divergence, fm.expected_divergence_between_fragments, places = 6)
        self.assertAlmostEqual(fm.expected_divergence, pi_summary.mean, places = 3)

    def test_expected_divergence(self):
        fm = frag.FragmentationModel(
                seed = 1111,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 10000,
                effective_pop_size_of_fragment = 10000,
                effective_pop_size_of_ancestor = 100000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi(100000))
        print(fm.expected_divergence_within_fragments)
        print(fm.expected_divergence_between_fragments)
        print(fm.expected_divergence)
        print(pi_summary.mean)
        self.assertAlmostEqual(fm.expected_divergence, pi_summary.mean, places = 3)

    def test_expected_divergence_within(self):
        fm = frag.FragmentationModel(
                seed = 1111,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 10000,
                effective_pop_size_of_fragment = 10000,
                effective_pop_size_of_ancestor = 100000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi_within(100000))
        print(fm.expected_divergence_within_fragments)
        print(pi_summary.mean)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, pi_summary.mean, places = 3)

    def test_expected_divergence_between(self):
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
        self.assertEqual(fm.number_of_fragments, 2)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi(200000))
        self.assertAlmostEqual(expected_div, pi_summary.mean, places = 3)

    def test_number_of_fragments(self):
        gens_since_frag = 20000.0
        pop_size = 10000
        mutation_rate = 1e-6
        expected_div = frag.get_expected_divergence_between_fragments(
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate)
        fm = frag.FragmentationModel(
                seed = 12345,
                number_of_fragments = 10,
                number_of_genomes_per_fragment = 1,
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_fragment = 1.0, # this shouldn't matter
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 10)
        pi_summary = stats.SampleSummarizer(
                fm.sample_pi(400000))
        self.assertAlmostEqual(expected_div, pi_summary.mean, places = 3)

if __name__ == '__main__':
    unittest.main()

