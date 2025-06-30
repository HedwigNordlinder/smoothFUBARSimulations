using Pkg
Pkg.activate(".")
Pkg.rm("CodonMolecularEvolution")
Pkg.add(url="https://github.com/HedwigNordlinder/CodonMolecularEvolution.jl")
using MolecularEvolution, CodonMolecularEvolution, Distributions, Plots, Phylo, DataFrames, CSV
include("plot_simulation_data.jl")

α_distribution = Gamma(10,0.1)
β_distribution = Exponential(1)

function create_simulation_parameter_csv(output_file::String = "simulation_parameters.csv")
    # Parameter combinations to explore
    #ntaxa_values = [10, 20, 50, 100, 200]
    ntaxa_values = [10,20]
    # nsites_values = [100, 500, 1000]
    nsites_values = [50,100]
    # diversifying_sites_values = [1, 10, 50]
    diversifying_sites_values = [0]
    scenarios = ["LogisticScenario()", "SeasonalScenario(;sin_divisor=1.0)"]
    
    # Static parameters
    static_params = (
        nucleotide_model = "default",
        f3x4_model = "default",
        target_normalisation = 60.0
    )
    
    # Create all combinations
    combinations = collect(Iterators.product(ntaxa_values, nsites_values, diversifying_sites_values, scenarios))
    
    # Create DataFrame
    df = DataFrame(
        scenario_name = String[],
        ntaxa = Int[],
        nsites = Int[],
        coalescence_scenario = String[],
        rate_sampler = String[],
        nucleotide_model = String[],
        f3x4_model = String[],
        target_normalisation = Float64[]
    )
    
    # Generate rows
    for (ntaxa, nsites, div_sites, scenario) in combinations
        scenario_name = "$(lowercase(replace(scenario, "()" => "")))_t$(ntaxa)_s$(nsites)_d$(div_sites)"
        
        # Create the rate sampler specification in Julia code
        rate_sampler = "DiversifyingSitesSampler(UnivariateRateSampler(Gamma(10,0.1), Exponential(1)), $div_sites, $nsites)"
        
        push!(df, (
            scenario_name = scenario_name,
            ntaxa = ntaxa,
            nsites = nsites,
            coalescence_scenario = scenario,
            rate_sampler = rate_sampler,
            nucleotide_model = static_params.nucleotide_model,
            f3x4_model = static_params.f3x4_model,
            target_normalisation = static_params.target_normalisation
        ))
    end
    
    # Write to CSV
    CSV.write(output_file, df)
    
    println("Created $(nrow(df)) simulation scenarios in $output_file")
    println("Combinations: $(length(ntaxa_values)) ntaxa × $(length(nsites_values)) nsites × $(length(diversifying_sites_values)) diversifying_sites × $(length(scenarios)) scenarios")
    
    return df
end

# Run it
df = create_simulation_parameter_csv("simulation_parameters.csv")

# Run the batch simulation
run_simulation_batch("simulation_parameters.csv", "simulations")
run_fubar_benchmark("simulations/",[DirichletFUBAR(), SKBDIFUBAR(), FIFEFUBAR()])
generate_roc_curves("simulations/")
collect_global_values("simulations/")