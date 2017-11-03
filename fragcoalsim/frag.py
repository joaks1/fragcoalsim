#! /usr/bin/env python

import sys
import os
import logging
import random

import msprime

from fragcoalsim import stats

_LOG = logging.getLogger(__name__)


def get_expected_divergence_between_fragments(
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000,
        mutation_rate = 1e-8):
    expected_coal_time = 2.0 * effective_pop_size_of_ancestor
    expected_branch_length = expected_coal_time + generations_since_fragmentation
    gens_diverging = expected_branch_length * 2.0
    expected_divergence = gens_diverging * mutation_rate
    return expected_divergence


class FragmentationModel(object):
    def __init__(self,
            seed,
            number_of_fragments = 5,
            number_of_genomes_per_fragment = 1,
            generations_since_fragmentation = 10.0,
            effective_pop_size_of_fragment = 100,
            effective_pop_size_of_ancestor = 1000,
            mutation_rate = 1e-8,
            migration_rate = 0.0):
        self._number_of_fragments = number_of_fragments
        self._sample_size = number_of_genomes_per_fragment
        self._generations_since_fragmentation = generations_since_fragmentation
        self._fragment_population_size = effective_pop_size_of_fragment
        self._ancestral_population_size = effective_pop_size_of_ancestor
        self._mutation_rate = mutation_rate
        self._migration_rate = migration_rate

        self._seed = seed
        self._rng = random.Random()
        self._rng.seed(self._seed)

        self._population_configs = []
        self._population_splits = []
        self._fragmentation_size_change = None
        self._migration_rate_matrix = []

        for i in range(self._number_of_fragments):
            self._population_configs.append(
                    msprime.PopulationConfiguration(
                        sample_size = self._sample_size,
                        initial_size = self._fragment_population_size,
                        growth_rate = 0.0))
            if i > 0:
                self._population_splits.append(
                        msprime.MassMigration(
                            time = self._generations_since_fragmentation,
                            source = i,
                            destination = 0,
                            proportion = 1.0))
        self._fragmentation_size_change = msprime.PopulationParametersChange(
                time = self._generations_since_fragmentation,
                initial_size = self._ancestral_population_size,
                growth_rate = 0.0,
                population_id = 0)
        self._migration_rate_matrix = [
                [0.0 for i in range(self._number_of_fragments)
                        ] for j in range(self._number_of_fragments)
                ]
        for i in range(self._number_of_fragments):
            for j in range(self._number_of_fragments):
                if i != j:
                    self._migration_rate_matrix[i][j] = self._migration_rate


    def _get_sim_seed(self):
        return self._rng.randint(1, 999999999)

    def simulate(self, number_of_replicates):
        return msprime.simulate(
                sample_size = None, # will be determined from pop configs
                Ne = 1.0, # reference effective pop size,
                length = 1,
                recombination_rate = 0.0,
                mutation_rate = self._mutation_rate,
                population_configurations = self._population_configs,
                migration_matrix = self._migration_rate_matrix,
                demographic_events = self._population_splits + [self._fragmentation_size_change],
                random_seed = self._get_sim_seed(),
                num_replicates = number_of_replicates)

    def sample_pi(self, number_of_replicates):
        return (x.get_pairwise_diversity() for x in self.simulate(
                number_of_replicates))

    def _get_seed(self):
        return self._seed
    seed = property(_get_seed)

    def _get_number_of_fragments(self):
        return self._number_of_fragments

    number_of_fragments = property(_get_number_of_fragments)

    def _get_sample_size(self):
        return self._sample_size

    sample_size = property(_get_sample_size)

    def _get_generations_since_fragmentation(self):
        return self._generations_since_fragmentation

    generations_since_fragmentation = property(_get_generations_since_fragmentation)

    def _get_fragment_population_size(self):
        return self._fragment_population_size

    fragment_population_size = property(_get_fragment_population_size)

    def _get_ancestral_population_size(self):
        return self._ancestral_population_size

    ancestral_population_size = property(_get_ancestral_population_size)

    def _get_mutation_rate(self):
        return self._mutation_rate

    mutation_rate = property(_get_mutation_rate)

    def _get_migration_rate(self):
        return self._migration_rate

    migration_rate = property(_get_migration_rate)

    def _get_population_configs(self):
        return self._population_configs

    population_configs = property(_get_population_configs)

    def _get_population_splits(self):
        return self._population_splits

    population_splits = property(_get_population_splits)

    def _get_fragmentation_size_change(self):
        return self._fragmentation_size_change

    fragmentation_size_change = property(_get_fragmentation_size_change)

    def _get_migration_rate_matrix(self):
        return self._migration_rate_matrix

    migration_rate_matrix = property(_get_migration_rate_matrix)


class FragmentationDiversityTracker(object):
    def __init__(self,
            seed,
            years_since_fragmentation_to_sample = [
                    0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
            generation_time = 1.0,
            number_of_simulations_per_sample = 1000,
            number_of_fragments = 5,
            number_of_genomes_per_fragment = 1,
            effective_pop_size_of_fragment = 100,
            effective_pop_size_of_ancestor = 1000,
            mutation_rate = 1e-8,
            migration_rate = 0.0):
        self._seed = seed
        self._rng = random.Random()
        self._rng.seed(self._seed)
        self._years_to_sample = [float(y) for y in years_since_fragmentation_to_sample]
        self._generation_time = float(generation_time)
        self._number_of_simulations = number_of_simulations_per_sample
        assert len(self._years_to_sample) > 0
        self.fragmentation_models = []
        self.pi_samples = []
        for year in self._years_to_sample:
            frag_model = FragmentationModel(
                    seed = self._get_model_seed(),
                    number_of_fragments = number_of_fragments,
                    number_of_genomes_per_fragment = number_of_genomes_per_fragment,
                    generations_since_fragmentation = year / self._generation_time,
                    effective_pop_size_of_fragment = effective_pop_size_of_fragment,
                    effective_pop_size_of_ancestor = effective_pop_size_of_ancestor,
                    mutation_rate = mutation_rate,
                    migration_rate = migration_rate)
            self.fragmentation_models.append(frag_model)
            pi_summary = stats.SampleSummarizer(
                    frag_model.sample_pi(self._number_of_simulations))
            self.pi_samples.append(pi_summary)

    def _get_model_seed(self):
        return self._rng.randint(1, 999999999)

    def _get_seed(self):
        return self._seed
    seed = property(_get_seed)

    def _get_years_to_sample(self):
        return self._years_to_sample
    years = property(_get_years_to_sample)

    def _get_years_to_sample(self):
        return self._years_to_sample
    years = property(_get_years_to_sample)

    def _get_number_of_simulations(self):
        return self._number_of_simulations
    number_of_simulations = property(_get_number_of_simulations)

    def _get_generation_time(self):
        return self._generation_time
    generation_time = property(_get_generation_time)
