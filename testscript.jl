Pkg.activate(".")
Pkg.add(url="https://github.com/HedwigNordlinder/CodonMolecularEvolution.jl")
using MolecularEvolution, CodonMolecularEvolution, Distributions

function flu_Ne(t; 
    baseline_Ne = 1000,     # baseline infected population
    seasonal_amplitude = 0.8, # how much seasonal variation (0-1)
    peak_time = 0.25,       # winter peak (fraction of year)
    period = 1.0            # yearly cycle
)
    # Seasonal pattern: peaks in winter, troughs in summer
    seasonal_factor = 1 + seasonal_amplitude * cos(2π * (t/period - peak_time))
    return baseline_Ne * seasonal_factor
end
function hospital_sampling_rate(t;
    base_rate = 0.1,
    seasonal_strength = 0.5,  # how much seasonal variation
    peak_time = 0.25,
    period = 1.0
)
    # exp(cos()) is always positive and gives nice seasonal variation
    seasonal_factor = exp(seasonal_strength * cos(2π * (t/period - peak_time)))
    return base_rate * seasonal_factor
end

α_distribution = Gamma(10,0.1)
β_distribution = Exponential(1)
ntaxa = 100

nsites = 100
diversifying_sites = 5

# This is the main simulation function

result = simulate_k_diversifying_sites(ntaxa,flu_Ne, hospital_sampling_rate, α_distribution, β_distribution, 
                                        nsites, diversifying_sites, CodonMolecularEvolution.demo_nucmat, 
                                        CodonMolecularEvolution.demo_f3x4)

save_simulation_data(result)                                    
