#!/usr/bin/env bash

# generation time for mouse lemurs is 3-4.5 years, lice its 28 days, rna
# viruses can produce 100,000 copies in 10 hours

# Viral Ne based on:
#       Bedford, T., S. Cobey, and M. Pascual. 2011. Strength and tempo of
#       selection reveaeled in viral gene genealogies. BMC Evolutionary Biology
#       11:220. https://doi.org/10.1186/1471-2148-11-220.

# Viral mu based on:
#       Bedford, T., S. Cobey, and M. Pascual. 2011. Strength and tempo of
#       selection reveaeled in viral gene genealogies. BMC Evolutionary Biology
#       11:220. https://doi.org/10.1186/1471-2148-11-220.
# and
#       https://www.nature.com/nrmicro/journal/v11/n5/box/nrmicro3003_BX1.html

# Viral generation time based on:
#       Bedford, T., S. Cobey, and M. Pascual. 2011. Strength and tempo of
#       selection reveaeled in viral gene genealogies. BMC Evolutionary Biology
#       11:220. https://doi.org/10.1186/1471-2148-11-220.

# Louse Ne loosely based on:
#       Toups, M. A., A. Kitchen, J. E. Light, and D. L. Reed. 2011. Origin of
#       clothing lice indicates early clothing use by anatomically modern
#       humans in Africa. Molecular Biology and Evolution 28(1):29--32.
#       https://doi.org/10.1093/molbev/msq234.
 
pitracer --number-of-fragments 5 \
    --number-of-sampled-gene-copies 10 \
    --number-of-replicates 100000 \
    --years-to-sample 0.0 20.0 40.0 60.0 80.0 100.0 \
    --labels Lemurs Lice Viruses \
    --ancestral-pop-sizes 50000 200000 4000 \
    --fragment-pop-sizes 1000 100000 4000 \
    --mutation-rates 1e-8 2e-8 1e-5 \
    --generation-times 4.0 0.08 0.011 \
    --migration-rates 0.0 0.0 0.0 \
    --seed 29478952 \
    --force
