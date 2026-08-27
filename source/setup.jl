# The standard first lines for each analysis script. Include this file like
# this, from a script in analysis/<subdirectory>/:
#
#     include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
#
# This file gives you three variables:
#   REPO      The absolute path of the repository root. Use it for output paths.
#   ss        The cached baseline steady state under current law.
#   dep_fit   The cached path of the ratio of retirees to workers, with one
#             value per projection year. It is fitted to the SSA Trustees Report
#             baseline. Refer to docs/PARAMETERS.md, keyword dep_path.
#
# A new calculation of `ss` takes about 1 minute. It is in the cache so that
# every analysis starts from the same baseline. Refer to docs/EXTENDING.md to
# find when you must build the cache again.

const REPO = normpath(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "functions_pwbm_w_spouse.jl"))

using JLD2
@load joinpath(REPO, "cached_objects.jld2") ss dep_fit

# The analysis scripts write spreadsheets to output/. XLSX does not make
# directories, so make them here.
mkpath(joinpath(REPO, "output", "new_policies"))
