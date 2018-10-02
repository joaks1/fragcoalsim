#! /usr/bin/env python

"""
CLI program for tracing pi following fragmentation.
"""

import os
import sys
import random
import argparse
import multiprocessing
import math

try:
    import matplotlib as mpl
    mpl_available = True
except:
    mpl_available = False

if mpl_available:
    # Use TrueType (42) fonts rather than Type 3 fonts
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    tex_font_settings = {
            "text.usetex": True,
            "font.family": "sans-serif",
            # "font.serif": [
            #         "Computer Modern Roman",
            #         "Times",
            #         ],
            # "font.sans-serif": [
            #         "Computer Modern Sans serif",
            #         "Helvetica",
            #         ],
            # "font.cursive": [
            #         "Zapf Chancery",
            #         ],
            # "font.monospace": [
            #         "Computer Modern Typewriter",
            #         "Courier",
            #         ],
            "text.latex.preamble" : [
                    "\\usepackage[T1]{fontenc}",
                    "\\usepackage[cm]{sfmath}",
                    ]
    }

    mpl.rcParams.update(tex_font_settings)

    import matplotlib.pyplot as plt
    from matplotlib import gridspec

import fragcoalsim
import fragcoalsim.frag

palette = [
(57  / 255.0, 115 / 255.0, 124 / 255.0), # teal
(184 / 255.0, 90  / 255.0, 13  / 255.0), # auburn
(50  / 255.0, 162 / 255.0, 81  / 255.0), # green
(60  / 255.0, 183 / 255.0, 204 / 255.0), # blue
(255 / 255.0, 127 / 255.0, 15  / 255.0), # orange
(255 / 255.0, 217 / 255.0, 74  / 255.0), # yellow
]

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
    parser.add_argument('-y', '--years-to-sample',
            nargs = '+',
            default = [0.0, 25.0, 50.0, 75.0, 100.0],
            type = fragcoalsim.argparse_utils.arg_is_nonnegative_float,
            help = ('Years to sample diversity following fragmentation.'))
    parser.add_argument('--ancestral-pop-sizes',
            nargs = '+',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Effective population size before fragmentation for each '
                    'taxon.'))
    parser.add_argument('--fragment-pop-sizes',
            nargs = '+',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Effective fragement population size for each taxon.'))
    parser.add_argument('--mutation-rates',
            nargs = '+',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Mutation rate for each taxon.'))
    parser.add_argument('--migration-rates',
            nargs = '+',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_nonnegative_float,
            help = ('Migration rate for each taxon.'))
    parser.add_argument('--generation-times',
            nargs = '+',
            default = [],
            type = fragcoalsim.argparse_utils.arg_is_positive_float,
            help = ('Generation time of each taxon.'))
    parser.add_argument('-l', '--labels',
            nargs = '+',
            default = [],
            type = str,
            help = ('Labels for each taxon.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('--force',
            action = 'store_true',
            help = ('Overwrite any existing output files. By default, an error '
                    'is thrown if an output path exists.'))
    parser.add_argument('--no-plot',
            action = 'store_true',
            help = ('Skip plotting; only report summary table.'))
    parser.add_argument('-r', '--number-of-replicates',
            action = 'store',
            default = 0,
            type = fragcoalsim.argparse_utils.arg_is_nonnegative_int,
            help = ('Number of coalescent simulations to run per taxon and per '
                    'sampled year. Default: 0 (no simulations; only report '
                    'expected diversity'))
    parser.add_argument('--np',
            action = 'store',
            default = multiprocessing.cpu_count(),
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Number of of processes to use during simulations. '
                    'Default: One per available CPU. Only relevant if '
                    'simulations are being run.'))
    parser.add_argument('--seed',
            action = 'store',
            type = fragcoalsim.argparse_utils.arg_is_positive_int,
            help = ('Seed for random number generator. Only relevant if '
                    'simulations are being run.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if not mpl_available:
        if not args.no_plot:
            sys.stderr.write("NOTE: Could not import matplotlib, "
                    "so turning plotting options off.\n")
        args.no_plot = True

    rng = random.Random()
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    rng.seed(args.seed)

    running_simulations = False
    if args.number_of_replicates > 0:
        running_simulations = True

    prefix = args.prefix
    if len(prefix.split(os.path.sep)) < 2:
        prefix = os.path.join(os.curdir, prefix)

    expected_div_path = prefix + "pitracer-plot-expected-div.pdf"
    sim_div_path = prefix + "pitracer-plot-sim-div.pdf"
    overlay_path = prefix + "pitracer-plot-overlay.pdf"
    output_dir = os.path.dirname(expected_div_path)
    if not output_dir:
        output_dir = os.curdir

    if not args.force:
        for p in [expected_div_path, sim_div_path, overlay_path]:
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

    if running_simulations:
        sys.stdout.write("label\tyears\te_pi\tmean_pi\te_pi_between\tmean_pi_between\tmutation_rate\tgeneration_time\tancestral_pop_size\tfragment_pop_size\tmigration_rate\tnfrags\tnsamples\n")
    else:
        sys.stdout.write("label\tyears\te_pi\te_pi_between\tmutation_rate\tgeneration_time\tancestral_pop_size\tfragment_pop_size\tmigration_rate\tnfrags\tnsamples\n")
    tracers = []
    for i in range(number_of_taxa):
        pi_tracer = fragcoalsim.frag.FragmentationDiversityTracker(
                seed = rng.randint(1, 999999999),
                years_since_fragmentation_to_sample = args.years_to_sample,
                generation_time = args.generation_times[i],
                number_of_fragments = args.number_of_fragments,
                number_of_genomes_per_fragment = args.number_of_sampled_gene_copies,
                effective_pop_size_of_fragment = args.fragment_pop_sizes[i],
                effective_pop_size_of_ancestor = args.ancestral_pop_sizes[i],
                mutation_rate = args.mutation_rates[i],
                migration_rate = args.migration_rates[i],
                number_of_simulations_per_sample = args.number_of_replicates,
                number_of_processes = args.np,
                locus_length = 1)
        if running_simulations:
            pi_tracer.run_simulations()
        tracers.append(pi_tracer)
        for j, y in enumerate(pi_tracer.years):
            if running_simulations:
                sys.stdout.write("{label}\t{years}\t{e_pi}\t{mean_pi}\t{e_pi_between}\t{mean_pi_between}\t{mutation_rate}\t{generation_time}\t{ancestral_pop_size}\t{fragment_pop_size}\t{migration_rate}\t{nfrags}\t{nsamples}\n".format(
                        label = args.labels[i],
                        years = y,
                        e_pi = pi_tracer.fragmentation_models[j].expected_divergence,
                        mean_pi = pi_tracer.pi_samples[j].mean,
                        e_pi_between = pi_tracer.fragmentation_models[j].expected_divergence_between_fragments,
                        mean_pi_between = pi_tracer.pi_between_samples[j].mean,
                        mutation_rate = pi_tracer.mutation_rate,
                        generation_time = pi_tracer.generation_time,
                        ancestral_pop_size = pi_tracer.ancestral_population_size,
                        fragment_pop_size = pi_tracer.fragment_population_size,
                        migration_rate = pi_tracer.migration_rate,
                        nfrags = pi_tracer.number_of_fragments,
                        nsamples = pi_tracer.sample_size,
                        ))
            else:
                sys.stdout.write("{label}\t{years}\t{e_pi}\t{e_pi_between}\t{mutation_rate}\t{generation_time}\t{ancestral_pop_size}\t{fragment_pop_size}\t{migration_rate}\t{nfrags}\t{nsamples}\n".format(
                        label = args.labels[i],
                        years = y,
                        e_pi = pi_tracer.fragmentation_models[j].expected_divergence,
                        e_pi_between = pi_tracer.fragmentation_models[j].expected_divergence_between_fragments,
                        mutation_rate = pi_tracer.mutation_rate,
                        generation_time = pi_tracer.generation_time,
                        ancestral_pop_size = pi_tracer.ancestral_population_size,
                        fragment_pop_size = pi_tracer.fragment_population_size,
                        migration_rate = pi_tracer.migration_rate,
                        nfrags = pi_tracer.number_of_fragments,
                        nsamples = pi_tracer.sample_size,
                        ))

    if args.no_plot:
        sys.exit(0)

    # Expected divergence plot
    plt.close('all')
    fig = plt.figure(figsize = (4, 3))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    for i, pi_tracer in enumerate(tracers):
        if i < len(palette):
            line_color = palette[i]
        else:
            line_color = palette[-1]
        line, = ax.plot(
                pi_tracer.years,
                [m.expected_divergence for m in pi_tracer.fragmentation_models],
                label = "{0} expected $\\pi$".format(args.labels[i]),
                )
        plt.setp(line,
                marker = 'o',
                markerfacecolor = line_color,
                markeredgecolor = line_color,
                markersize = 3.5,
                linestyle = '--',
                linewidth = 2,
                color = line_color,
                zorder = 200,
                rasterized = False)
        line, = ax.plot(
                pi_tracer.years,
                [m.expected_divergence_between_fragments for m in pi_tracer.fragmentation_models],
                label = "{0} expected $\\pi_F$".format(args.labels[i]),
                )
        plt.setp(line,
                linestyle = '-',
                linewidth = 2,
                color = line_color,
                alpha = 0.5,
                zorder = 100,
                rasterized = False)
    xlabel_text = ax.set_xlabel("Years since fragmentation")
    ylabel_text = ax.set_ylabel("Diversity")
    line_types = []
    line_colors = []
    for ls in ['--', '-']:
        line_types.append(ax.plot([],[], color = "black", ls = ls)[0])
    for i in range(number_of_taxa):
        line_colors.append(ax.plot([],[], color = palette[i], ls = '-')[0])
    lines = ax.get_lines()
    legend1 = plt.legend(
            line_colors,
            args.labels,
            loc = "center right",
            ncol = 1,
            edgecolor = "none",
            )
    plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 1.02), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
    legend2 = plt.legend(
            line_types, 
            [
                    "Overall ($\\pi$)",
                    "Among fragments ($\\pi_F$)"
            ],
            loc = "upper left",
            edgecolor = "none",
            )
    ax.add_artist(legend1)
    gs.update(left = 0.14, right = 0.995, bottom = 0.14, top = 0.995)
    plt.savefig(expected_div_path)

    if running_simulations:
        # Simulation plot
        plt.close('all')
        fig = plt.figure(figsize = (4, 3))
        gs = gridspec.GridSpec(1, 1,
                wspace = 0.0,
                hspace = 0.0)
        ax = plt.subplot(gs[0, 0])
        for i, pi_tracer in enumerate(tracers):
            if i < len(palette):
                line_color = palette[i]
            else:
                line_color = palette[-1]
            line, = ax.plot(
                    pi_tracer.years,
                    [s.mean for s in pi_tracer.pi_samples],
                    label = "{0} mean $\\pi$".format(args.labels[i]),
                    )
            plt.setp(line,
                    marker = 'o',
                    markerfacecolor = line_color,
                    markeredgecolor = line_color,
                    markersize = 3.5,
                    linestyle = '--',
                    linewidth = 2,
                    color = line_color,
                    zorder = 200,
                    rasterized = False)
            line, = ax.plot(
                    pi_tracer.years,
                    [s.mean for s in pi_tracer.pi_between_samples],
                    label = "{0} mean $\\pi_F$".format(args.labels[i]),
                    )
            plt.setp(line,
                    linestyle = '-',
                    linewidth = 2,
                    color = line_color,
                    alpha = 0.5,
                    zorder = 100,
                    rasterized = False)
        xlabel_text = ax.set_xlabel("Years since fragmentation")
        ylabel_text = ax.set_ylabel("Diversity")
        line_types = []
        line_colors = []
        for ls in ['--', '-']:
            line_types.append(ax.plot([],[], color = "black", ls = ls)[0])
        for i in range(number_of_taxa):
            line_colors.append(ax.plot([],[], color = palette[i], ls = '-')[0])
        lines = ax.get_lines()
        legend1 = plt.legend(
                line_colors,
                args.labels,
                loc = "center right",
                ncol = 1,
                edgecolor = "none",
                )
        plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 1.02), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
        legend2 = plt.legend(
                line_types, 
                [
                        "Overall ($\\pi$)",
                        "Among fragments ($\\pi_F$)"
                ],
                loc = "upper left",
                edgecolor = "none",
                )
        ax.add_artist(legend1)
        gs.update(left = 0.14, right = 0.995, bottom = 0.14, top = 0.995)
        plt.savefig(sim_div_path)


if __name__ == "__main__":
    main()
