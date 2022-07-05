using Agents
using LinearAlgebra
using Random
using Distributions
using StatsBase
using SharedArrays
using CSV
using DataFrames

##

mutable struct Cell <: AbstractAgent
    id::Int
    pos::Dims{2}
    isStem::Bool
    isQuiescent::Bool
    p_proliferate::Float64
    p_migrate::Float64
    p_symmetry::Float64
    p_dedifferentiation::Float64
    p_capacity::Int
    resistance::Float64
end

Base.@kwdef mutable struct CellProperties
    p_mitotic_death::Float64 = 0.01         #probability of dying when replicate
    p_max::Int = 15                         #maximum proliferation capacity
    stem_cell_resistance::Float64 = 0.1376  #resistance of CSC. Smaller is more resist
    cancer_cell_resistance::Float64 = 1.0   #resistance of non-CSC.
    irradiation::Bool = false               #irradiation status
    is_targeted::Bool = false               #irradiation status
    dose_per_fraction::Float64 = 2.0        #radiation dosage per fraction
    therapy_time = 1:24:24*30               #fractionation schedule
    step_counter::Int = 0                   #1 step = 1 hour
    p_mutation::Float64 = 0.0001            #mutation rate
end


"""
Create a list of the 8 neighboring positions
"""
function build_neighbor_pos_list(pos, boundary, n=1)
    x, y = pos
    min_x = max(x - n, 1)
    max_x = min(x + n, boundary[1])
    min_y = max(y - n, 1)
    max_y = min(y + n, boundary[2])
    neighbor = Iterators.product(min_x:max_x, min_y:max_y)

    return neighbor
end

"""
Return an empty position from 8 neighboring position if there is any
"""
function locate_empty_neighbor_position(agent::Cell, model, n=1)
    pos = agent.pos
    boundary = size(model.space)
    neighbor_pos_list = build_neighbor_pos_list(pos, boundary, n)
    empty_positions = [pos for pos in neighbor_pos_list if isempty(pos, model)]
    if length(empty_positions) > 0
        idx_empty_pos = rand(1:length(empty_positions))
        empty_position = empty_positions[idx_empty_pos]
    else
        empty_position = Dims{2}((-1, -1))
    end
    return empty_position
end

"""
Calculate the survival fraction
"""
function radiation_sf(agent, dosage, is_targeted;
    quiescent_resistance=0.5,
    non_quiescent_resistance=1.0,
    targeted_resistance=0.5,
    α=0.3859, β=0.01148)
    σ = agent.isQuiescent ? quiescent_resistance : non_quiescent_resistance # radioresistance  quiescent vs non-quiescent cell
    if is_targeted && agent.isStem
        resistance = σ * (agent.resistance + targeted_resistance)
    else
        resistance = σ * agent.resistance
    end


    exponent = -resistance * (α * dosage + β * dosage^2)
    return exp(exponent)
end

"""
Calculate the survival fraction
"""
function radiation_sf(agent, model::ABM;
    quiescent_resistance=0.5,
    non_quiescent_resistance=1.0,
    targeted_resistance=0.5,
    α=0.3859, β=0.01148)
    dosage = model.dose_per_fraction
    is_targeted = model.is_targeted
    return radiation_sf(agent, dosage, is_targeted;
        quiescent_resistance=quiescent_resistance,
        non_quiescent_resistance=non_quiescent_resistance,
        targeted_resistance=targeted_resistance,
        α=0.3859, β=0.01148)
end



"""
Simulate the cell replication. If the replication succeeded, it would added an agent
"""
function replication!(agent::Cell, new_pos, model)

    if rand(model.rng) < agent.p_dedifferentiation      # cell plasticity change cell stemness phenotype
        agent.isStem = ~(agent.isStem)
        ∇resistance = model.cancer_cell_resistance - model.stem_cell_resistance
        if agent.isStem
            agent.resistance = agent.resistance - ∇resistance  # if change phenotype to CSC, decrease the lambda
        else
            agent.resistance = agent.resistance + ∇resistance  # if change phenotype to CC, increase the lambda
        end
    end

    is_stem = agent.isStem
    agent.isQuiescent = false   # set isQuiescent to false because the cell replicate
    is_quiescent = false
    dist = Bernoulli(model.p_mutation)
    proliferation_sd = 0.01
    migration_sd = 0.05
    symmetry_sd = 0.005
    dedifferentiation_sd = 0.005
    resistance_sd = 0.05

    # mutation
    proliferation = rand(model.rng, dist) ?
                    rand(model.rng, Normal(agent.p_proliferate, proliferation_sd)) :
                    agent.p_proliferate
    migration = rand(model.rng, dist) ?
                rand(model.rng, Normal(agent.p_migrate, migration_sd)) :
                agent.p_migrate
    symmetry = rand(model.rng, dist) ?
               rand(model.rng, Normal(agent.p_symmetry, symmetry_sd)) :
               agent.p_symmetry
    dedifferentiation = rand(model.rng, dist) ?
                        rand(model.rng, Normal(agent.p_dedifferentiation, dedifferentiation_sd)) :
                        agent.p_dedifferentiation


    if is_stem
        daughter_is_stem = rand(model.rng) < agent.p_symmetry # symmetric division
        capacity = agent.p_capacity
        if daughter_is_stem
            resistance = rand(model.rng, dist) ?
                         rand(model.rng, Normal(agent.resistance, resistance_sd)) :
                         agent.resistance
        else
            resistance = rand(model.rng, dist) ?
                         rand(model.rng, Normal(model.cancer_cell_resistance, resistance_sd)) :
                         model.cancer_cell_resistance
        end
        add_agent!(new_pos,
            model,
            daughter_is_stem,
            is_quiescent,
            proliferation,
            migration,
            symmetry,
            dedifferentiation,
            capacity,
            resistance,
        )
    else
        resistance = rand(model.rng, dist) ?
                     rand(model.rng, Normal(agent.resistance, resistance_sd)) :
                     agent.resistance
        agent.p_capacity = agent.p_capacity - 1
        capacity = agent.p_capacity
        add_agent!(new_pos,
            model,
            agent.isStem,
            is_quiescent,
            proliferation,
            migration,
            symmetry,
            dedifferentiation,
            capacity,
            resistance,
        )
    end
end

"""
Check whether the agent has negative parameter value
"""
function is_deleterious_mutation(agent)
    return any([agent.p_proliferate <= 0,
        agent.p_migrate <= 0,
        agent.p_symmetry <= 0,
        agent.p_dedifferentiation <= 0,
        agent.p_capacity <= 0,
        agent.resistance <= 0,])
end

function agent_step!(agent::Cell, model)
    if is_deleterious_mutation(agent)                   # remove unviable cell with negative parameters
        kill_agent!(agent, model)
        return Nothing
    end
    new_pos = locate_empty_neighbor_position(agent, model)
    new_position_exist = new_pos[1] == -1 ? false : true
    proliferate = rand(model.rng) < agent.p_proliferate
    if new_position_exist && proliferate
        if agent.p_capacity < 1 || rand(model.rng) < model.p_mitotic_death
            kill_agent!(agent, model)
            return Nothing
        else
            replication!(agent, new_pos, model)
        end
    elseif new_position_exist && rand(model.rng) < agent.p_migrate
        move_agent!(agent, new_pos, model)
    elseif !new_position_exist
        agent.isQuiescent = true
    end
    if model.irradiation && rand(model.rng) < (1 - radiation_sf(agent, model)) # probability of death in irradiation
        kill_agent!(agent, model)
        return Nothing
    end
    return Nothing
end


"""
Initiate the model. Populate the tumor with the predefined number of the cell or time step.
"""
function model_initiation(;
    N=1000,                     #square domain dimension
    p_proliferate=1 / 24,          # cell proliferate once a day (0.0417)
    p_migrate=15 / 24,             # cell move 150 μm per day (0.625)
    p_symmetry=0.01,
    p_dedifferentiation=0.01,
    p_mitotic_death=0.01,         #probability of dying when replicate
    p_max=10,                         #maximum proliferation capacity
    stem_cell_resistance=0.1376,  #resistance of CSC. Smaller is more resist
    cancer_cell_resistance=1.0,   #resistance of non-CSC.
    p_mutation=0.0001,            #mutation rate
    seed=1
)
    properties = CellProperties(
        p_mitotic_death=p_mitotic_death,
        p_max=p_max,
        stem_cell_resistance=stem_cell_resistance,
        cancer_cell_resistance=cancer_cell_resistance,
        p_mutation=p_mutation,
    )
    x, y = div(N, 2), div(N, 2)
    initial_pos = [
        (x, y)
    ]
    space = GridSpace((N, N), periodic=false)
    rng = Random.MersenneTwister(seed)
    model = ABM(
        Cell,
        space,
        properties=properties,
        rng=rng
    )
    p_capacity = model.p_max
    resistance = model.stem_cell_resistance
    for pos in initial_pos
        add_agent!(pos,
            model,
            true,
            false,
            p_proliferate,
            p_migrate,
            p_symmetry,
            p_dedifferentiation,
            p_capacity,
            resistance,
        )
    end
    return model
end

"""
Create the n replicates of the initial tumors
"""
function models_initiation(;
    N=1000,                     #square domain dimension
    n_step=20_000,                  #number of step after first agent added,
    p_proliferate=1 / 24,          # cell proliferate once a day (0.0417)
    p_migrate=15 / 24,             # cell move 150 μm per day (0.625)
    p_symmetry=0.01,
    p_dedifferentiation=0.01,
    p_mitotic_death=0.01,         #probability of dying when replicate
    p_max=15,                         #maximum proliferation capacity
    stem_cell_resistance=0.1376,  #resistance of CSC. Smaller is more resist
    cancer_cell_resistance=1.0,   #resistance of non-CSC.
    p_mutation=0.1,            #mutation rate
    n_cell=50000,
    n_models=10
)
    models = [model_initiation(
        N=N,                     #square domain dimension
        p_proliferate=p_proliferate,          # cell proliferate once a day (0.0417)
        p_migrate=p_migrate,             # cell move 150 μm per day (0.625)
        p_symmetry=p_symmetry,
        p_dedifferentiation=p_dedifferentiation,
        p_mitotic_death=p_mitotic_death,         #probability of dying when replicate
        p_max=p_max,                         #maximum proliferation capacity
        stem_cell_resistance=stem_cell_resistance,  #resistance of CSC. Smaller is more resist
        cancer_cell_resistance=cancer_cell_resistance,   #resistance of non-CSC.
        p_mutation=p_mutation,            #mutation rate
        seed=i
    ) for i in 1:n_models]
    stop_function(m, s) = nagents(m) > n_cell || s > n_step
    ensemblerun!(models, agent_step!, dummystep, stop_function)
    return models
end



"""
Step rule for model in control group
"""
function control_step!(m)
    m.step_counter += 1
end


"""
Step rule for model in treated group
"""
function therapy_step!(m)
    m.step_counter += 1
    m.irradiation = m.step_counter in m.therapy_time
end

function targeted_therapy_step!(m)
    m.step_counter += 1
    m.irradiation = m.step_counter in m.therapy_time
    m.is_targeted = true
end

"""
Cell color for ABM visualization
"""
function cell_color(agent, pmax::Int)
    if agent.isStem
        return colormap("Reds", pmax + 1)[agent.p_capacity+1]
    else
        return colormap("Blues", pmax + 1)[agent.p_capacity+1]
    end
end

as(agent) = agent.isStem ? 2 : 1 # cell size for ABM visualization
ac(agent) = agent.isStem ? :red : :blue # cell color for ABM visualization
is_stem(a) = a.isStem
adata = [:p_proliferate, :p_migrate, :p_symmetry, :p_dedifferentiation, :resistance, :p_capacity] # data collected when run the model
stem_cell_proportion(m) = count(is_stem, allagents(m)) / nagents(m) # count the proportion of the stem cell
stem_cell_count(m) = count(is_stem, allagents(m)) # count the number of the stem cell
