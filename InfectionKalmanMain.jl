module InfectionKalmanMain

using CSVFiles
using DataFrames

include("InfectionKalman.jl")

import .InfectionKalman
export res, n, r, s, x, pk, pkk

if length(ARGS) == 0
    fdt = "2020-10-12"
    loc = "36"
else
    fdt = ARGS[1]
    loc = ARGS[2]
end

pdata = DataFrame(load(string("initial-pars--", fdt, "--", loc, ".csv")))
cdata = DataFrame(load(string("data--", fdt, "--", loc, ".csv")))

res, n, r, s, x, pk, pkk = InfectionKalman.fit(cdata, pdata; show_trace = true, detailed_results = true, time_limit = 600)
pdata.minimizer = res.minimizer

save(string("minimizer--", fdt, "--", loc, ".csv"), pdata)
#save(string("hessian--", fdt, "--", loc, ".csv"), DataFrame(h))

end