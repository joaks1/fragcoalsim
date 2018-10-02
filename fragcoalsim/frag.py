#! /usr/bin/env python

import sys
import os
import math
import logging
import random
import re
import subprocess
import multiprocessing
import itertools
from contextlib import contextmanager
import json

# import msprime

from fragcoalsim import stats, check_external_tool
from fragcoalsim import mspi

_LOG = logging.getLogger(__name__)


def get_expected_coal_time_between_fragments(
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000):
    expected_anc_coal_time = 2.0 * effective_pop_size_of_ancestor
    return expected_anc_coal_time + generations_since_fragmentation

def get_expected_divergence_between_fragments(
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000,
        mutation_rate = 1e-8):
    expected_coal_time = get_expected_coal_time_between_fragments(
            generations_since_fragmentation = generations_since_fragmentation,
            effective_pop_size_of_ancestor = effective_pop_size_of_ancestor)
    gens_diverging = expected_coal_time * 2.0
    expected_divergence = gens_diverging * mutation_rate
    return expected_divergence

def get_expected_coal_time_within_fragments(
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000,
        effective_pop_size_of_fragment = 100):
    # Get expected coalescence time of two gene copies conditional on them
    # coalescing within the fragment population

    n2 = 2.0 * effective_pop_size_of_fragment
    expected_coal_time_within_frag = n2
    weight_within = True
    # Above doesn't work, because coalescences that occur in ancestor skew this
    # E.g., the expectation can be greater than time since fragmentation
    # Instead, since the number of generations of the fragment branch is capped
    # by the generations_since_fragmentation, we can take a weighted average
    # over every generation within the fragment population (every generation is
    # simply weighted by the probability of coalescence at that generation)
    if (10 * expected_coal_time_within_frag) > generations_since_fragmentation:
        weight_within = False
        expected_coal_time_within_frag = 0.0
        for gen in range(1, int(round(generations_since_fragmentation)) + 1):
            log_prob_coal = (math.log((n2 - 1.0) / n2) * (gen - 1)) + math.log(1.0 / n2)
            expected_coal_time_within_frag += math.exp(log_prob_coal + math.log(gen))

    # Get expected coalescence time of two gene copies condition on them NOT
    # coalescing within the fragment population
    expected_coal_time_between_frags = get_expected_coal_time_between_fragments(
            generations_since_fragmentation = generations_since_fragmentation,
            effective_pop_size_of_ancestor = effective_pop_size_of_ancestor)

    # Now average coal times weighted by the probability of coal versus no coal
    # within fragment
    prob_of_no_coal_within_frag = ((n2 - 1.0) / n2) ** generations_since_fragmentation
    prob_of_coal_within_frag = 1.0 - prob_of_no_coal_within_frag

    within_wt = 1.0
    if weight_within:
        within_wt = prob_of_coal_within_frag
    expected_coal = ((prob_of_no_coal_within_frag * expected_coal_time_between_frags) +
            (within_wt * expected_coal_time_within_frag))

    return expected_coal

def get_expected_divergence_within_fragments(
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000,
        effective_pop_size_of_fragment = 100,
        mutation_rate = 1e-8):
    expected_coal_time = get_expected_coal_time_within_fragments(
            generations_since_fragmentation = generations_since_fragmentation,
            effective_pop_size_of_ancestor = effective_pop_size_of_ancestor,
            effective_pop_size_of_fragment = effective_pop_size_of_fragment)
    return 2.0 * expected_coal_time * mutation_rate

def get_expected_coal_time(
        number_of_fragments = 5,
        number_of_genomes_per_fragment = 1,
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000,
        effective_pop_size_of_fragment = 100):
    if number_of_fragments < 2:
        e_coal_between = None
    else:
        e_coal_between = get_expected_coal_time_between_fragments(
                generations_since_fragmentation = generations_since_fragmentation,
                effective_pop_size_of_ancestor = effective_pop_size_of_ancestor)
    if number_of_genomes_per_fragment < 2:
        e_coal_within = None
    else:
        e_coal_within = get_expected_coal_time_within_fragments(
                generations_since_fragmentation = generations_since_fragmentation,
                effective_pop_size_of_ancestor = effective_pop_size_of_ancestor,
                effective_pop_size_of_fragment = effective_pop_size_of_fragment)

    if (e_coal_between is None) and (e_coal_within is None):
        return None, None, None
    if e_coal_between is None:
        return e_coal_within, e_coal_between, e_coal_within
    if e_coal_within is None:
        return e_coal_between, e_coal_between, e_coal_within

    n_comps_per_frag = (number_of_genomes_per_fragment * (number_of_genomes_per_fragment - 1.0)) / 2.0
    n_comps_within_frags = n_comps_per_frag * number_of_fragments
    ngenomes = number_of_fragments * number_of_genomes_per_fragment
    n_comps = (ngenomes * (ngenomes - 1)) / 2.0
    p_within = n_comps_within_frags / n_comps
    p_between = 1.0 - p_within

    e_coal = ((p_between * e_coal_between) +
            (p_within * e_coal_within))
    return e_coal, e_coal_between, e_coal_within


def get_expected_divergence(
        number_of_fragments = 5,
        number_of_genomes_per_fragment = 1,
        generations_since_fragmentation = 10.0,
        effective_pop_size_of_ancestor = 1000,
        effective_pop_size_of_fragment = 100,
        mutation_rate = 1e-7):

    m = 2.0 * mutation_rate
    e_coal, e_coal_between, e_coal_within = get_expected_coal_time(
            number_of_fragments = number_of_fragments,
            number_of_genomes_per_fragment = number_of_genomes_per_fragment,
            generations_since_fragmentation = generations_since_fragmentation,
            effective_pop_size_of_ancestor = effective_pop_size_of_ancestor,
            effective_pop_size_of_fragment = effective_pop_size_of_fragment)
    if (e_coal_between is None) and (e_coal_within is None):
        return None, None, None
    if e_coal_between is None:
        return m * e_coal_within, None, m * e_coal_within
    if e_coal_within is None:
        return m * e_coal_between, m * e_coal_between, None 
    return m * e_coal, m * e_coal_between, m * e_coal_within 


@contextmanager
def pool_context(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def run_mspi_simulations(frag_model, number_of_processes = None,
        number_of_replicates = 50):
    if number_of_processes is None:
        number_of_processes = multiprocessing.cpu_count()
    batch_size = int(math.ceil(number_of_replicates / float(number_of_processes)))
    args = ((frag_model.get_random_copy(), batch_size) for i in range(number_of_processes))
    with pool_context(processes = number_of_processes) as pool:
        results = pool.map(_mspi_sim_unpack, args)
    pi = []
    pi_between = []
    pi_within = []
    for p, p_a, p_w in results:
        pi.extend(p)
        pi_between.extend(p_a)
        pi_within.extend(p_w)
    return pi, pi_between, pi_within

def _mspi_sim(frag_model, number_of_replicates = 1):
    return frag_model.mspi_simulate(
            number_of_replicates = number_of_replicates)

def _mspi_sim_unpack(args):
    return _mspi_sim(*args)


class FragmentationModel(object):
    ms_segsite_pattern = re.compile(r"^segsites:\s+(?P<segsites>\d+)\s*$")
    ms_positions_pattern = re.compile(r"^positions:\s+(?P<positions>[0-9. ]+)$")
    ms_char_pattern = re.compile(r"^(?P<characters>[01]+)$")

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
        
        div, div_b, div_w = get_expected_divergence(
                number_of_fragments = self._number_of_fragments,
                number_of_genomes_per_fragment = self._sample_size,
                generations_since_fragmentation = self._generations_since_fragmentation,
                effective_pop_size_of_ancestor = self._ancestral_population_size,
                effective_pop_size_of_fragment = self._fragment_population_size,
                mutation_rate = self._mutation_rate)
        self._expected_div = div
        self._expected_div_between_fragments = div_b
        self._expected_div_within_fragments = div_w

        # self._population_configs = []
        # self._population_splits = []
        # self._fragmentation_size_change = None
        # self._migration_rate_matrix = []

        # for i in range(self._number_of_fragments):
        #     self._population_configs.append(
        #             msprime.PopulationConfiguration(
        #                 sample_size = self._sample_size,
        #                 initial_size = self._fragment_population_size,
        #                 growth_rate = 0.0))
        #     if i > 0:
        #         self._population_splits.append(
        #                 msprime.MassMigration(
        #                     time = self._generations_since_fragmentation,
        #                     source = i,
        #                     destination = 0,
        #                     proportion = 1.0))
        # self._fragmentation_size_change = msprime.PopulationParametersChange(
        #         time = self._generations_since_fragmentation,
        #         initial_size = self._ancestral_population_size,
        #         growth_rate = 0.0,
        #         population_id = 0)
        # self._migration_rate_matrix = [
        #         [0.0 for i in range(self._number_of_fragments)
        #                 ] for j in range(self._number_of_fragments)
        #         ]
        # for i in range(self._number_of_fragments):
        #     for j in range(self._number_of_fragments):
        #         if i != j:
        #             self._migration_rate_matrix[i][j] = self._migration_rate

        self._ms_path = check_external_tool("ms")
        self._mspi_path = check_external_tool("mspi")


    def _get_sim_seed(self):
        return self._rng.randint(1, 999999999)

    def get_random_copy(self):
        return FragmentationModel(
                seed = self._get_sim_seed(),
                number_of_fragments = self._number_of_fragments,
                number_of_genomes_per_fragment = self._sample_size,
                generations_since_fragmentation = self._generations_since_fragmentation,
                effective_pop_size_of_fragment = self._fragment_population_size,
                effective_pop_size_of_ancestor = self._ancestral_population_size,
                mutation_rate = self._mutation_rate,
                migration_rate = self._migration_rate)

    # def simulate(self, number_of_replicates):
    #     return msprime.simulate(
    #             sample_size = None, # will be determined from pop configs
    #             Ne = 1.0, # reference effective pop size,
    #             length = 1,
    #             recombination_rate = 0.0,
    #             mutation_rate = self._mutation_rate,
    #             population_configurations = self._population_configs,
    #             migration_matrix = self._migration_rate_matrix,
    #             demographic_events = self._population_splits + [self._fragmentation_size_change],
    #             random_seed = self._get_sim_seed(),
    #             num_replicates = number_of_replicates)

    # def sample_pi(self, number_of_replicates):
    #     return (x.get_pairwise_diversity() for x in self.simulate(
    #             number_of_replicates))

    # def sample_pi_within(self, number_of_replicates):
    #     divs = []
    #     for sim in self.simulate(number_of_replicates):
    #         for pop in sim.populations():
    #             div = sim.get_pairwise_diversity(samples = sim.samples(pop.id))
    #             divs.append(div)
    #     return divs


    def get_ms_command(self, locus_length = 1, number_of_replicates = 1):
        current_theta = 4.0 * self.fragment_population_size * self.mutation_rate * locus_length
        theta_change = self.ancestral_population_size / float(self.fragment_population_size)
        time = self.generations_since_fragmentation / (4.0 * self.fragment_population_size)
        migration = 4.0 * self.fragment_population_size * self.migration_rate
        ms_args = [
                self._ms_path,
                str(self.sample_size * self.number_of_fragments),
                str(number_of_replicates),
                "-t",
                str(current_theta),
                "-I",
                str(self.number_of_fragments),
                ]
        ms_args.extend(str(self.sample_size) for i in range(self.number_of_fragments))
        ms_args.append(str(migration * (self.number_of_fragments - 1)))
        for i in range(2, self.number_of_fragments + 1):
            ms_args.extend(["-ej", str(time), str(i), "1"])
        ms_args.extend(["-en", str(time), "1", str(theta_change)])
        ms_args.extend(["-eM", str(time), "0.0"])
        ms_args.extend(["-seeds", str(self._get_sim_seed())])
        return ms_args

    def mspi_simulate(self, number_of_replicates = 1):
        cmd = self.get_ms_command(
                locus_length = 1,
                number_of_replicates = number_of_replicates)
        pi, pi_b, pi_w = mspi.mspi_run_sims(cmd)
        assert len(pi) == number_of_replicates
        assert len(pi) == len(pi_b)
        assert len(pi) == len(pi_w)
        return pi, pi_b, pi_w

    def mspi_simulate_subprocess(self, locus_length = 1, number_of_replicates = 1):
        if not self._mspi_path:
            raise OSError(
                    "The mspi executable is not available for running simulations.\n"
                    "Please install mspi and try again.")
        sout = subprocess.PIPE
        serr = subprocess.PIPE
        cmd = self.get_ms_command(
                locus_length = locus_length,
                number_of_replicates = number_of_replicates)
        cmd[0] = self._mspi_path
        try:
            process = subprocess.Popen(
                    cmd,
                    stdout = sout,
                    stderr = serr,
                    shell = False)
        except OSError as e:
            _LOG.error("Problem running ms executable")
            raise e
        pi, pi_b, pi_w = self.parse_mspi_output(stream = process.stdout,
                locus_length = locus_length)
        exit_code = process.wait()
        assert len(pi) == number_of_replicates
        assert len(pi) == len(pi_b)
        assert len(pi) == len(pi_w)
        return pi, pi_b, pi_w

    def parse_mspi_output(self, stream, locus_length):
        results = json.loads("".join(stream.readlines()[2:]))
        pi = [x / float(locus_length) for x in results["pi"]]
        pi_b = [x / float(locus_length) for x in results["pi_b"]]
        pi_w = [x / float(locus_length) for x in results["pi_w"]]
        return pi, pi_b, pi_w


    def ms_simulate(self, locus_length = 1, number_of_replicates = 1):
        if not self._ms_path:
            raise OSError(
                    "The ms executable is not available for running simulations.\n"
                    "Please install ms and try again. You can find ms here:\n"
                    "http://home.uchicago.edu/rhudson1/source/mksamples.html")
        sout = subprocess.PIPE
        serr = subprocess.PIPE
        cmd = self.get_ms_command(
                locus_length = locus_length,
                number_of_replicates = number_of_replicates)
        try:
            process = subprocess.Popen(
                    cmd,
                    stdout = sout,
                    stderr = serr,
                    shell = False)
        except OSError as e:
            _LOG.error("Problem running ms executable")
            raise e
        pi, pi_b, pi_w = self.parse_ms_output(stream = process.stdout,
                locus_length = locus_length)
        # stdout, stderr = process.communicate()
        exit_code = process.wait()
        # print(stderr)
        # print(stdout)
        assert len(pi) == number_of_replicates
        assert len(pi) == len(pi_b)
        assert len(pi) == len(pi_w)
        return pi, pi_b, pi_w

    def parse_ms_output(self, stream, locus_length):
        pi, pi_b, pi_w = [], [], []
        for line in stream:
            if line.strip() == "//":
                p, b, w = self.process_ms_replicate_output(stream, locus_length)
                pi.append(p)
                pi_b.append(b)
                pi_w.append(w)
        return pi, pi_b, pi_w

    def process_ms_replicate_output(self, stream, locus_length):
        segsites_line = stream.next()
        m = self.ms_segsite_pattern.match(segsites_line)
        if not m:
            raise Exception("Unexpected first line of ms replicate {0!r}".format(
                    segsites_line.strip()))
        num_seg_sites = int(m.group("segsites"))
        if num_seg_sites < 1:
            return 0.0, 0.0, 0.0
        position_line = stream.next()
        char_matrix = []
        for pop_index in range(self.number_of_fragments):
            for sample_index in range(self.sample_size):
                char_line = stream.next()
                char_match = self.ms_char_pattern.match(char_line)
                if not char_match:
                    raise Exception("Problem parsing ms character line {0!r}".format(
                            char_line.strip()))
                char_matrix.append((pop_index, char_line.strip()))
        assert len(char_matrix) == self.sample_size * self.number_of_fragments
        pi, pi_b, pi_w = self.calculate_pi_from_ms_characters(char_matrix,
                locus_length)
        return pi, pi_b, pi_w

    def calculate_pi_from_ms_characters(cls, character_matrix, locus_length):
        num_comps = 0
        num_comps_within = 0
        num_comps_between = 0
        sum_diffs = 0
        sum_diffs_within = 0
        sum_diffs_between = 0
        for i, (c1, c2) in enumerate(itertools.combinations(character_matrix, 2)):
            pop_index1, chars1 = c1
            pop_index2, chars2 = c2
            assert len(chars1) == len(chars2)
            ndiffs = sum(1 for a, b in zip(chars1, chars2) if a != b)
            sum_diffs += ndiffs
            num_comps += 1
            if pop_index1 == pop_index2:
                sum_diffs_within += ndiffs
                num_comps_within += 1
            else:
                sum_diffs_between += ndiffs
                num_comps_between += 1
        assert num_comps == i + 1
        assert num_comps == (num_comps_within + num_comps_between)
        pi = (sum_diffs / float(num_comps)) / float(locus_length)
        if num_comps_within < 1:
            pi_w = 0.0
        else:
            pi_w = (sum_diffs_within / float(num_comps_within)) / float(locus_length)
        if num_comps_between < 1:
            pi_b = 0.0
        else:
            pi_b = (sum_diffs_between / float(num_comps_between)) / float(locus_length)
        return pi, pi_b, pi_w

    def get_prob_of_no_coal_within_frag(self):
        """
        Get the probability that two gene copies sampled from the same fragment
        fail to coalesce within the fragment:
            = (1 - 1/2N)^t = ((2N - 1)/2N)^t
        """
        n2 = self.fragment_population_size * 2.0
        return ((n2 - 1) / n2) ** self.generations_since_fragmentation

    def _get_expected_div(self):
        return self._expected_div
    expected_divergence = property(_get_expected_div)

    def _get_expected_div_between_fragments(self):
        return self._expected_div_between_fragments
    expected_divergence_between_fragments = property(_get_expected_div_between_fragments)

    def _get_expected_div_within_fragments(self):
        return self._expected_div_within_fragments
    expected_divergence_within_fragments = property(_get_expected_div_within_fragments)

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
            number_of_fragments = 5,
            number_of_genomes_per_fragment = 1,
            effective_pop_size_of_fragment = 100,
            effective_pop_size_of_ancestor = 1000,
            mutation_rate = 1e-8,
            migration_rate = 0.0,
            number_of_simulations_per_sample = 1000,
            number_of_processes = 2,
            locus_length = 1):
        self._seed = seed
        self._rng = random.Random()
        self._rng.seed(self._seed)
        self._years_to_sample = [float(y) for y in years_since_fragmentation_to_sample]
        self._generation_time = float(generation_time)
        self._number_of_simulations = number_of_simulations_per_sample
        self._number_of_processes = number_of_processes
        self._locus_length = locus_length
        assert len(self._years_to_sample) > 0

        self._number_of_fragments = number_of_fragments
        self._sample_size = number_of_genomes_per_fragment
        self._fragment_population_size = effective_pop_size_of_fragment
        self._ancestral_population_size = effective_pop_size_of_ancestor
        self._mutation_rate = mutation_rate
        self._migration_rate = migration_rate

        self.fragmentation_models = []
        self.pi_samples = []
        self.pi_between_samples = []
        for year in self._years_to_sample:
            frag_model = FragmentationModel(
                    seed = self._get_model_seed(),
                    number_of_fragments = self._number_of_fragments,
                    number_of_genomes_per_fragment = self._sample_size,
                    generations_since_fragmentation = year / self._generation_time,
                    effective_pop_size_of_fragment = self._fragment_population_size,
                    effective_pop_size_of_ancestor = self._ancestral_population_size,
                    mutation_rate = self._mutation_rate,
                    migration_rate = self._migration_rate)
            self.fragmentation_models.append(frag_model)

    def run_simulations(self):
        for frag_model in self.fragmentation_models:
            pi, pi_b, pi_w = run_mspi_simulations(frag_model,
                    number_of_processes = self._number_of_processes,
                    number_of_replicates = self._number_of_simulations)
            pi_summary = stats.SampleSummarizer(pi)
            pi_b_summary = stats.SampleSummarizer(pi_b)
            self.pi_samples.append(pi_summary)
            self.pi_between_samples.append(pi_b_summary)

    def _get_model_seed(self):
        return self._rng.randint(1, 999999999)

    def _get_seed(self):
        return self._seed
    seed = property(_get_seed)

    def _get_years_to_sample(self):
        return self._years_to_sample
    years = property(_get_years_to_sample)

    def _get_number_of_simulations(self):
        return self._number_of_simulations
    number_of_simulations = property(_get_number_of_simulations)

    def _get_generation_time(self):
        return self._generation_time
    generation_time = property(_get_generation_time)

    def _get_number_of_fragments(self):
        return self._number_of_fragments

    number_of_fragments = property(_get_number_of_fragments)

    def _get_sample_size(self):
        return self._sample_size

    sample_size = property(_get_sample_size)

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
