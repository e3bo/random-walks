module InfectionKalmanTest



if length(ARGS) == 0
    fdt = "2020-10-12"
    loc = "36"
else
    fdt = ARGS[1]
    loc = ARGS[2]
end

pdata = DataFrame(load(string("initial-pars--", fdt, "--", loc, ".csv")))
cdata = DataFrame(load(string("data--", fdt, "--", loc, ".csv")))

save(string("minimizer--", fdt, "--", loc, ".csv"), pdata)

end