using Test
using Tidejinks
using Dates
import SPICE

t = DateTime(1993, 1, 1, 1) + Second(times[1])

kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

@test isfile(kernel_meta_file)

