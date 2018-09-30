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

    def test_virus(self):
        nreps = 20000
        fm = frag.FragmentationModel(
                seed = 123456,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 20000,
                effective_pop_size_of_fragment = 5000,
                effective_pop_size_of_ancestor = 5000,
                mutation_rate = 1e-5,
                migration_rate = 0.0)
        pi, pi_a, pi_w = frag.run_mspi_simulations(fm,
                number_of_processes = None,
                number_of_replicates = nreps,
                batch_size = 1000,
                locus_length = 1)
        self.assertEqual(len(pi), nreps)
        self.assertEqual(len(pi_a), nreps)
        self.assertEqual(len(pi_w), nreps)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi_w = sum(pi_w) / len(pi_w)
        print("")
        print(mean_pi)
        print(fm.expected_divergence)
        print(mean_pi_a)
        print(fm.expected_divergence_between_fragments)
        print(mean_pi_w)
        print(fm.expected_divergence_within_fragments)
        self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)

    def test_mspi(self):
        fm = frag.FragmentationModel(
                seed = 123456,
                number_of_fragments = 5,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 1000,
                effective_pop_size_of_fragment = 1000,
                effective_pop_size_of_ancestor = 10000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 10, number_of_replicates = 20000)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi_w = sum(pi_w) / len(pi_w)
        print("")
        print(mean_pi)
        print(fm.expected_divergence)
        print(mean_pi_a)
        print(fm.expected_divergence_between_fragments)
        print(mean_pi_w)
        print(fm.expected_divergence_within_fragments)
        self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 3)
        self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 3)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 3)

    def test_pool(self):
        nreps = 20000
        fm = frag.FragmentationModel(
                seed = 1234,
                number_of_fragments = 5,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 1000,
                effective_pop_size_of_fragment = 1000,
                effective_pop_size_of_ancestor = 10000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        pi, pi_a, pi_w = frag.run_mspi_simulations(fm,
                number_of_processes = None,
                number_of_replicates = nreps,
                batch_size = 100,
                locus_length = 10)
        self.assertEqual(len(pi), nreps)
        self.assertEqual(len(pi_a), nreps)
        self.assertEqual(len(pi_w), nreps)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi_w = sum(pi_w) / len(pi_w)
        print("")
        print(mean_pi)
        print(fm.expected_divergence)
        print(mean_pi_a)
        print(fm.expected_divergence_between_fragments)
        print(mean_pi_w)
        print(fm.expected_divergence_within_fragments)
        self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 3)
        self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 3)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 3)

    # def test_ms(self):
    #     fm = frag.FragmentationModel(
    #             seed = 123456,
    #             number_of_fragments = 5,
    #             number_of_genomes_per_fragment = 10,
    #             generations_since_fragmentation = 1000,
    #             effective_pop_size_of_fragment = 1000,
    #             effective_pop_size_of_ancestor = 10000,
    #             mutation_rate = 1e-6,
    #             migration_rate = 0.0)
    #     pi, pi_a, pi_w = fm.ms_simulate(locus_length = 10, number_of_replicates = 20000)
    #     mean_pi = sum(pi) / len(pi)
    #     mean_pi_a = sum(pi_a) / len(pi_a)
    #     mean_pi_w = sum(pi_w) / len(pi_w)
    #     print("")
    #     print(mean_pi)
    #     print(fm.expected_divergence)
    #     print(mean_pi_a)
    #     print(fm.expected_divergence_between_fragments)
    #     print(mean_pi_w)
    #     print(fm.expected_divergence_within_fragments)
    #     self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 3)
    #     self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 3)
    #     self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 3)

    # def test_expected_divergence_with_one_sample(self):
    #     gens_since_frag = 2000.0
    #     pop_size = 10000
    #     mutation_rate = 1e-6
    #     expected_div, d_b, d_w = frag.get_expected_divergence(
    #             number_of_fragments = 2,
    #             number_of_genomes_per_fragment = 1,
    #             generations_since_fragmentation = gens_since_frag,
    #             effective_pop_size_of_ancestor = pop_size,
    #             effective_pop_size_of_fragment = 1000,
    #             mutation_rate = mutation_rate)
    #     fm = frag.FragmentationModel(
    #             seed = 12345,
    #             number_of_fragments = 2,
    #             number_of_genomes_per_fragment = 1,
    #             generations_since_fragmentation = gens_since_frag,
    #             effective_pop_size_of_fragment = 1000,
    #             effective_pop_size_of_ancestor = pop_size,
    #             mutation_rate = mutation_rate,
    #             migration_rate = 0.0)
    #     self.assertEqual(fm.number_of_fragments, 2)
    #     # pi_summary = stats.SampleSummarizer(
    #     #         fm.sample_pi(200000))
    #     pi, pi_a, pi_w = fm.ms_simulate(locus_length = 1, number_of_replicates = 50000)
    #     mean_pi = sum(pi) / len(pi)
    #     mean_pi_a = sum(pi_a) / len(pi_a)
    #     print("")
    #     print(mean_pi)
    #     print(expected_div)
    #     print(mean_pi_a)
    #     print(d_b)
    #     self.assertAlmostEqual(expected_div, mean_pi, places = 3)
    #     self.assertAlmostEqual(expected_div, d_b, places = 3)
    #     self.assertTrue(d_w is None)

    def test_mspi_expected_divergence_with_one_sample(self):
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
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 1, number_of_replicates = 50000)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        print("")
        print(mean_pi)
        print(expected_div)
        print(mean_pi_a)
        print(d_b)
        self.assertAlmostEqual(expected_div, mean_pi, places = 3)
        self.assertAlmostEqual(expected_div, d_b, places = 3)
        self.assertTrue(d_w is None)

    # def test_expected_divergence_no_frag(self):
    #     fm = frag.FragmentationModel(
    #             seed = 1111,
    #             number_of_fragments = 2,
    #             number_of_genomes_per_fragment = 5,
    #             generations_since_fragmentation = 0,
    #             effective_pop_size_of_fragment = 10000,
    #             effective_pop_size_of_ancestor = 100000,
    #             mutation_rate = 1e-6,
    #             migration_rate = 0.0)
    #     self.assertEqual(fm.number_of_fragments, 2)
    #     # pi_summary = stats.SampleSummarizer(
    #     #         fm.sample_pi(100000))
    #     pi, pi_a, pi_w = fm.ms_simulate(locus_length = 10, number_of_replicates = 50000)
    #     mean_pi = sum(pi) / len(pi)
    #     mean_pi_a = sum(pi_a) / len(pi_a)
    #     mean_pi_w = sum(pi_w) / len(pi_w)
    #     print("")
    #     print(mean_pi)
    #     print(fm.expected_divergence)
    #     print(mean_pi_a)
    #     print(fm.expected_divergence_between_fragments)
    #     print(mean_pi_w)
    #     print(fm.expected_divergence_within_fragments)
    #     self.assertAlmostEqual(fm.expected_divergence, fm.expected_divergence_within_fragments, places = 6)
    #     self.assertAlmostEqual(fm.expected_divergence, fm.expected_divergence_between_fragments, places = 6)
    #     self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
    #     self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)

    def test_mspi_expected_divergence_no_frag(self):
        fm = frag.FragmentationModel(
                seed = 1111,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 5,
                generations_since_fragmentation = 0,
                effective_pop_size_of_fragment = 10000,
                effective_pop_size_of_ancestor = 100000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 10, number_of_replicates = 50000)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi_w = sum(pi_w) / len(pi_w)
        print("")
        print(mean_pi)
        print(fm.expected_divergence)
        print(mean_pi_a)
        print(fm.expected_divergence_between_fragments)
        print(mean_pi_w)
        print(fm.expected_divergence_within_fragments)
        self.assertAlmostEqual(fm.expected_divergence, fm.expected_divergence_within_fragments, places = 6)
        self.assertAlmostEqual(fm.expected_divergence, fm.expected_divergence_between_fragments, places = 6)
        self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)

    # def test_expected_divergence(self):
    #     fm = frag.FragmentationModel(
    #             seed = 1111111,
    #             number_of_fragments = 2,
    #             number_of_genomes_per_fragment = 5,
    #             generations_since_fragmentation = 10000,
    #             effective_pop_size_of_fragment = 100000,
    #             effective_pop_size_of_ancestor = 1000,
    #             mutation_rate = 1e-5,
    #             migration_rate = 0.0)
    #     self.assertEqual(fm.number_of_fragments, 2)
    #     # pi_summary = stats.SampleSummarizer(
    #     #         fm.sample_pi(100000))
    #     pi, pi_a, pi_w = fm.ms_simulate(locus_length = 1, number_of_replicates = 100000)
    #     mean_pi = sum(pi) / len(pi)
    #     mean_pi_a = sum(pi_a) / len(pi_a)
    #     mean_pi_w = sum(pi_w) / len(pi_w)
    #     print("")
    #     print(mean_pi)
    #     print(fm.expected_divergence)
    #     print(mean_pi_a)
    #     print(fm.expected_divergence_between_fragments)
    #     print(mean_pi_w)
    #     print(fm.expected_divergence_within_fragments)
    #     self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
    #     self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 2)
    #     # self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)
    #     self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 1)

    def test_mspi_expected_divergence(self):
        fm = frag.FragmentationModel(
                seed = 1111111,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 5,
                generations_since_fragmentation = 10000,
                effective_pop_size_of_fragment = 100000,
                effective_pop_size_of_ancestor = 1000,
                mutation_rate = 1e-5,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 1, number_of_replicates = 100000)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi_w = sum(pi_w) / len(pi_w)
        print("")
        print(mean_pi)
        print(fm.expected_divergence)
        print(mean_pi_a)
        print(fm.expected_divergence_between_fragments)
        print(mean_pi_w)
        print(fm.expected_divergence_within_fragments)
        self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 2)
        # self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 1)

    # def test_expected_divergence_within(self):
    #     fm = frag.FragmentationModel(
    #             seed = 12345678,
    #             number_of_fragments = 2,
    #             number_of_genomes_per_fragment = 10,
    #             generations_since_fragmentation = 10000,
    #             effective_pop_size_of_fragment = 10000,
    #             effective_pop_size_of_ancestor = 100000,
    #             mutation_rate = 1e-6,
    #             migration_rate = 0.0)
    #     self.assertEqual(fm.number_of_fragments, 2)
    #     # pi_summary = stats.SampleSummarizer(
    #     #         fm.sample_pi_within(100000))
    #     pi, pi_a, pi_w = fm.ms_simulate(locus_length = 1, number_of_replicates = 100000)
    #     mean_pi = sum(pi) / len(pi)
    #     mean_pi_a = sum(pi_a) / len(pi_a)
    #     mean_pi_w = sum(pi_w) / len(pi_w)
    #     print("")
    #     print(mean_pi)
    #     print(fm.expected_divergence)
    #     print(mean_pi_a)
    #     print(fm.expected_divergence_between_fragments)
    #     print(mean_pi_w)
    #     print(fm.expected_divergence_within_fragments)
    #     self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
    #     self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 2)
    #     self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)

    def test_mspi_expected_divergence_within(self):
        fm = frag.FragmentationModel(
                seed = 12345678,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 10,
                generations_since_fragmentation = 10000,
                effective_pop_size_of_fragment = 10000,
                effective_pop_size_of_ancestor = 100000,
                mutation_rate = 1e-6,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 1, number_of_replicates = 100000)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi_w = sum(pi_w) / len(pi_w)
        print("")
        print(mean_pi)
        print(fm.expected_divergence)
        print(mean_pi_a)
        print(fm.expected_divergence_between_fragments)
        print(mean_pi_w)
        print(fm.expected_divergence_within_fragments)
        self.assertAlmostEqual(fm.expected_divergence, mean_pi, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_between_fragments, mean_pi_a, places = 2)
        self.assertAlmostEqual(fm.expected_divergence_within_fragments, mean_pi_w, places = 2)

    # def test_expected_divergence_between(self):
    #     gens_since_frag = 20000.0
    #     pop_size = 10000
    #     mutation_rate = 1e-6
    #     expected_div = frag.get_expected_divergence_between_fragments(
    #             generations_since_fragmentation = gens_since_frag,
    #             effective_pop_size_of_ancestor = pop_size,
    #             mutation_rate = mutation_rate)
    #     fm = frag.FragmentationModel(
    #             seed = 123456,
    #             number_of_fragments = 2,
    #             number_of_genomes_per_fragment = 1,
    #             generations_since_fragmentation = gens_since_frag,
    #             effective_pop_size_of_fragment = 1.0, # this shouldn't matter
    #             effective_pop_size_of_ancestor = pop_size,
    #             mutation_rate = mutation_rate,
    #             migration_rate = 0.0)
    #     self.assertEqual(fm.number_of_fragments, 2)
    #     # pi_summary = stats.SampleSummarizer(
    #     #         fm.sample_pi(200000))
    #     pi, pi_a, pi_w = fm.ms_simulate(locus_length = 1, number_of_replicates = 200000)
    #     mean_pi = sum(pi) / len(pi)
    #     mean_pi_a = sum(pi_a) / len(pi_a)
    #     self.assertAlmostEqual(expected_div, mean_pi_a, places = 3)
    #     self.assertAlmostEqual(mean_pi, mean_pi_a, places = 3)

    def test_mspi_expected_divergence_between(self):
        gens_since_frag = 20000.0
        pop_size = 10000
        mutation_rate = 1e-6
        expected_div = frag.get_expected_divergence_between_fragments(
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate)
        fm = frag.FragmentationModel(
                seed = 123456,
                number_of_fragments = 2,
                number_of_genomes_per_fragment = 1,
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_fragment = 1.0, # this shouldn't matter
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 2)
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 1, number_of_replicates = 200000)
        mean_pi = sum(pi) / len(pi)
        mean_pi_a = sum(pi_a) / len(pi_a)
        self.assertAlmostEqual(expected_div, mean_pi_a, places = 3)
        self.assertAlmostEqual(mean_pi, mean_pi_a, places = 3)

    def test_number_of_fragments(self):
        gens_since_frag = 20000.0
        pop_size = 10000
        mutation_rate = 1e-6
        expected_div = frag.get_expected_divergence_between_fragments(
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate)
        fm = frag.FragmentationModel(
                seed = 12345678,
                number_of_fragments = 10,
                number_of_genomes_per_fragment = 1,
                generations_since_fragmentation = gens_since_frag,
                effective_pop_size_of_fragment = 1.0, # this shouldn't matter
                effective_pop_size_of_ancestor = pop_size,
                mutation_rate = mutation_rate,
                migration_rate = 0.0)
        self.assertEqual(fm.number_of_fragments, 10)
        # pi_summary = stats.SampleSummarizer(
        #         fm.sample_pi(400000))
        pi, pi_a, pi_w = fm.mspi_simulate(locus_length = 1, number_of_replicates = 400000)
        mean_pi_a = sum(pi_a) / len(pi_a)
        mean_pi = sum(pi) / len(pi)
        self.assertAlmostEqual(expected_div, mean_pi_a, places = 3)
        self.assertAlmostEqual(mean_pi, mean_pi_a, places = 3)

if __name__ == '__main__':
    unittest.main()

