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

    parser.add_argument('-y', '--years-to-sample',
            nargs = '+',
            type = fragcoalsim.argparse_utils.arg_is_nonnegative_float,
            help = ('Years to sample diversity following fragmentation.'))
    parser.add_argument('-s', '--seed',
            action = 'store',
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Seed for random number generator.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('-f', '--force',
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
                        "or the \'-f/--force\' option to overwrite existing "
                        "files.".format(p))


    pi_tracer = fragcoalsim.frag.FragmentationDiversityTracker(
            seed = rng.randint(1, 999999999),
            years_since_fragmentation_to_sample = args.years_to_sample,
            generation_time = 1.0,
            number_of_simulations_per_sample = 10000,
            number_of_fragments = 5,
            number_of_genomes_per_fragment = 1,
            sequence_length = 1000,
            effective_pop_size_of_fragment = 100,
            effective_pop_size_of_ancestor = 1000,
            mutation_rate = 1e-8,
            migration_rate = 0.0)

    print(", ".join(str(x.mean) for x in pi_tracer.pi_samples))

    if args.no_plot:
        sys.exit(0)

if __name__ == "__main__":
    main()
