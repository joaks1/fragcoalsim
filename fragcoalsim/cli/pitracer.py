#! /usr/bin/env python

"""
CLI program for tracing pi following fragmentation.
"""

import os
import sys
import random
import argparse

import fragcoalsim

def main(argv = sys.argv):
    fragcoalsim.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--number-of-fragments',
            action = 'store',
            default = 2,
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Number of fragments. Default: 2'))
    parser.add_argument('-g', '--number-of-sampled-gene-copies',
            action = 'store',
            default = 1,
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Number of sampled gene copies per fragment. Default: 1'))
    parser.add_argument('-r', '--number-of-replicates',
            action = 'store',
            default = 10000,
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Number of coalescent simulations to run per taxon and per '
                    'sampled year. Default: 10000'))
    parser.add_argument('-y', '--years-to-sample',
            nargs = '?',
            default = [0.0, 25.0, 50.0, 75.0, 100.0],
            type = fragcoalsim.argparse_utils.arg_is_nonnegative_float,
            help = ('Years to sample diversity following fragmentation.'))
    parser.add_argument('--ancestral-pop-sizes',
            nargs = '?',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Effective population size before fragmentation for each '
                    'taxon.'))
    parser.add_argument('--fragment-pop-sizes',
            nargs = '?',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Effective fragement population size for each taxon.'))
    parser.add_argument('--mutation-rates',
            nargs = '?',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Mutation rate for each taxon.'))
    parser.add_argument('--migration-rates',
            nargs = '?',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Migration rate for each taxon.'))
    parser.add_argument('--generation-times',
            nargs = '?',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Generation time of each taxon.'))
    parser.add_argument('-l', '--labels',
            nargs = '?',
            default = [],
            type = str,
            help = ('Labels for each taxon.'))
    parser.add_argument('--seed',
            action = 'store',
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Seed for random number generator.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('--force',
            action = 'store_true',
            help = ('Overwrite any existing output files. By default, an error '
                    'is thrown if an output path exists.'))
    parser.add_argument('--x-label',
            action = 'store',
            type = str,
            default = "Time",
            help = ('Label for the X-axis. Default: \'Time\'.'))
    parser.add_argument('--y-label',
            action = 'store',
            type = str,
            default = "Diversity",
            help = ('Label for the Y-axis. Default: \'Diversity\'.'))
    parser.add_argument('--no-plot',
            action = 'store_true',
            help = ('Skip plotting; only report summary table.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    rng = random.Random()
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    rng.seed(args.seed)

    prefix = args.prefix
    if len(prefix.split(os.path.sep)) < 2:
        prefix = os.path.join(os.curdir, prefix)

    svg_path = prefix + "pytrace-plot.svg"
    output_dir = os.path.dirname(svg_path)
    if not output_dir:
        output_dir = os.curdir

    if not args.force:
        for p in [svg_path]:
            if os.path.exists(p):
                raise Exception(
                        "\nERROR: File {0!r} already exists.\n"
                        "Use \'-p/--prefix\' option to specify a different prefix,\n"
                        "or the \'--force\' option to overwrite existing "
                        "files.".format(p))

    number_of_taxa = max([
            len(args.labels),
            len(args.mutation_rates),
            len(args.generation_times),
            len(args.migration_rates),
            len(args.ancestral_pop_sizes),
            len(args.fragment_pop_sizes),
            ])
    if number_of_taxa < 1:
        number_of_taxa = 1

    if not args.labels:
        args.labels = [str(i + 1) for i in range(number_of_taxa)]
    else:
        assert len(args.labels) == number_of_taxa, (
                "The length of taxon-specific arguments must match")

    if not args.mutation_rates:
        args.mutation_rates = [1e-7 for i in range(number_of_taxa)]
    else:
        assert len(args.mutation_rates) == number_of_taxa, (
                "The length of taxon-specific arguments must match")

    if not args.migration_rates:
        args.migration_rates = [0.0 for i in range(number_of_taxa)]
    else:
        assert len(args.migration_rates) == number_of_taxa, (
                "The length of taxon-specific arguments must match")

    if not args.generation_times:
        args.generation_times = [1.0 for i in range(number_of_taxa)]
    else:
        assert len(args.generation_times) == number_of_taxa, (
                "The length of taxon-specific arguments must match")

    if not args.ancestral_pop_sizes:
        args.ancestral_pop_sizes = [10000.0 for i in range(number_of_taxa)]
    else:
        assert len(args.ancestral_pop_sizes) == number_of_taxa, (
                "The length of taxon-specific arguments must match")

    if not args.fragment_pop_sizes:
        args.fragment_pop_sizes = [10000.0 for i in range(number_of_taxa)]
    else:
        assert len(args.fragment_pop_sizes) == number_of_taxa, (
                "The length of taxon-specific arguments must match")

    sys.stdout.write("label\tyears\tmean_pi\te_frag_div\tmutation_rate\tgeneration_time\tancestral_pop_size\tfragment_pop_size\tmigration_rate\n")
    tracers = []
    for i in range(number_of_taxa):
        pi_tracer = fragcoalsim.frag.FragmentationDiversityTracker(
                seed = rng.randint(1, 999999999),
                years_since_fragmentation_to_sample = args.years_to_sample,
                generation_time = args.generation_times[i],
                number_of_simulations_per_sample = args.number_of_replicates,
                number_of_fragments = args.number_of_fragments,
                number_of_genomes_per_fragment = args.number_of_sampled_gene_copies,
                effective_pop_size_of_fragment = args.fragment_pop_sizes[i],
                effective_pop_size_of_ancestor = args.ancestral_pop_sizes[i],
                mutation_rate = args.mutation_rates[i],
                migration_rate = args.migration_rate[i])
        tracers.append(pi_tracer)
        for j, y in enumerate(pi_tracer.years):
            sys.stdout.write("{label}\t{years}\t{mean_pi}\t{e_frag_div}\t{mutation_rate}\t{generation_time}\t{ancestral_pop_size}\t{fragment_pop_size}\t{migration_rate}\n".format(
                    label = args.labels[i],
                    years = y,
                    mean_pi = pi_tracer.pi_samples[j],
                    e_frag_div = pi_tracer.fragmentation_models[j].expected_divergence_between_fragments,
                    mutation_rate = pi_tracer.mutation_rate,
                    generation_time = pi_tracer.generation_time,
                    ancestral_pop_size = pi_tracer.ancestral_population_size,
                    fragment_pop_size = pi_tracer.fragment_population_size,
                    migration_rate = pi_tracer.migration_rate))

    if args.no_plot:
        sys.exit(0)

if __name__ == "__main__":
    main()
