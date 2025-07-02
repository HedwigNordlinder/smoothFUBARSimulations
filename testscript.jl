using Pkg
Pkg.activate(".")
#Pkg.rm("CodonMolecularEvolution")
#Pkg.develop(path="/home/hedwig/Documents/Julia/karolinska/summer/CodonMolecularEvolution.jl")
Pkg.add(url="https://github.com/HedwigNordlinder/CodonMolecularEvolution.jl")
using MolecularEvolution, CodonMolecularEvolution, Distributions, Plots, Phylo, DataFrames, CSV
include("plot_simulation_data.jl")

α_distribution = Gamma(10, 0.1)
β_distribution = Exponential(1)

function create_simulation_parameter_csv(output_file::String="simulation_parameters.csv")
    # Parameter combinations to explore
    #ntaxa_values = [20,100,500]
    ntaxa_values = [20, 100]
    #nsites_values = [100,500,2500]
    nsites_values = [100, 200]

    #diversifying_sites_values = [0.01,0.05,0.1,0.2]
    diversifying_sites_values = [0.01]
    scenarios = [
        "LogisticScenario()",
        "StandardLadderScenario()"
    ]
    #normalisations = [1.0,5.0,10.0,20.0] 
    normalisations = [1.0]
    #λs = [0.5,1.0,2.0]
    λs = [1.0,2.0]
    # Static parameters
    static_params = (
        nucleotide_model="default",
        f3x4_model="default",
        #target_normalisation = 60.0,
        n_replicates=1
    )

    # Create all combinations
    combinations = collect(Iterators.product(ntaxa_values, nsites_values, diversifying_sites_values, scenarios, normalisations, λs))

    # Create DataFrame
    df = DataFrame(
        scenario_name=String[],
        ntaxa=Int[],
        nsites=Int[],
        coalescence_scenario=String[],
        rate_sampler=String[],
        nucleotide_model=String[],
        f3x4_model=String[],
        target_normalisation=Float64[],
        n_replicates=Int64[]
    )

    # Generate rows
    for (ntaxa, nsites, div_sites, scenario, normalisation, λ) in combinations
        scenario_name = "$(lowercase(replace(scenario, "()" => "")))_t$(ntaxa)_s$(nsites)_d$(div_sites)_n$(normalisation)_l$(λ)"


        # Create the rate sampler specification in Julia code
        rate_sampler = "DiversifyingSitesSampler(UnivariateRateSampler(Gamma(10,0.1), Exponential($λ)), $div_sites, $nsites)"
        push!(df, (
            scenario_name=scenario_name,
            ntaxa=ntaxa,
            nsites=nsites,
            coalescence_scenario=scenario,
            rate_sampler=rate_sampler,
            nucleotide_model=static_params.nucleotide_model,
            f3x4_model=static_params.f3x4_model,
            target_normalisation=normalisation,
            n_replicates=static_params.n_replicates
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
#underlying_sampler = UnivariateRateSampler(Gamma(10,0.1),Exponential(1))
#sampler = DiversifyingSitesSampler(underlying_sampler, 1, 100)
#CodonMolecularEvolution.serialize_sampler_to_dict(sampler)