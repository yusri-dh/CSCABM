using Distributed
addprocs(8, exeflags="--project")
using Colors

##
@everywhere include("csc.jl")
params_mut_05 = Dict(
    :N => 1000,
    :p_mutation => 0.5,
    :p_max => 15,
    :n_cell => 50_000,
)

params_mut_01 = Dict(
    :N => 1000,
    :p_mutation => 0.1,
    :p_max => 15,
    :n_cell => 50_000,
)

params_mut_03 = Dict(
    :N => 1000,
    :p_mutation => 0.3,
    :p_max => 15,
    :n_cell => 50_000,
)
n_models = 20

##
Random.seed!(185)
# create 20 initial tumor with p_mut= 0.5
models_mut_05 = models_initiation(;
    n_models=n_models,
    params_mut_05...)

# create 20 initial tumor with p_mut= 0.1
models_mut_01 = models_initiation(;
    n_models=n_models,
    params_mut_01...)

# create 20 initial tumor with p_mut= 0.3
models_mut_03 = models_initiation(;
    n_models=n_models,
    params_mut_03...)
const n_step = 1200

dir_name = "/media/yusri/data/project/data_agent_4"


##
# convenient function to run model in control group
function control_condition(models;
    path="control",
    n_step=200,
)
    isdir(path) || mkdir(path)
    adata = [:p_proliferate,
        :p_migrate,
        :p_symmetry,
        :p_dedifferentiation,
        :resistance,
        :p_capacity,
        :isStem,
        :isQuiescent]
    @distributed for i in eachindex(models)
        model = deepcopy(models[i])
        adf, _ = run!(model,
            agent_step!,
            control_step!,
            n_step,
            when=0:12:n_step,
            adata=adata)
        CSV.write(joinpath(path, "$i.csv"), adf)
    end
    return Nothing
end

# control p_mut=0.5
control_condition(models_mut_05; n_step=n_step,
    path=joinpath(dir_name, "control_mut_05"))

# control p_mut=0.1
control_condition(models_mut_01; n_step=n_step,
    path=joinpath(dir_name, "control_mut_01"))

# control p_mut=0.3
control_condition(models_mut_03; n_step=n_step,
    path=joinpath(dir_name, "control_mut_03"))
##
# convenient function to run model in treated group
function therapy_condition(models, agent_step!, model_step!;
    path="therapy",
    dose_per_fraction=2.0,
    therapy_time=1:24:24*30,
    p_mutation=0.5,
    n_step=200, seed=18,
    kwargs...)
    Random.seed!(seed)
    isdir(path) || mkdir(path)
    adata = [:p_proliferate,
        :p_migrate,
        :p_symmetry,
        :p_dedifferentiation,
        :resistance,
        :p_capacity,
        :isStem,
        :isQuiescent]
    @distributed for i in eachindex(models)
        model = deepcopy(models[i])
        model.dose_per_fraction = dose_per_fraction
        model.p_mutation = p_mutation
        model.therapy_time = therapy_time
        adf, _ = run!(model,
            agent_step!,
            model_step!,
            n_step,
            when=0:12:n_step,
            adata=adata,
        )
        CSV.write(joinpath(path, "$i.csv"), adf)
    end
    return Nothing
end
##

# conventional radiotherapy p_mut=0.5
therapy_condition(models_mut_05, agent_step!, therapy_step!;
    path=joinpath(dir_name, "therapy_2Gy_24h_30days_mut_05"),
    dose_per_fraction=2.0,
    therapy_time=1:24:24*30,
    p_mutation=0.5,
    n_step=n_step, seed=1)

# hyperfractionated radiotherapy p_mut=0.5
therapy_condition(models_mut_05, agent_step!, therapy_step!;
    path=joinpath(dir_name, "therapy_1Gy_12h_30days_mut_05"),
    dose_per_fraction=1.0,
    therapy_time=1:12:24*30,
    p_mutation=0.5,
    n_step=n_step, seed=2)

# hypofractionated radiotherapy p_mut=0.5
therapy_condition(models_mut_05, agent_step!, therapy_step!;
    path=joinpath(dir_name, "therapy_4Gy_48h_30days_mut_05"),
    dose_per_fraction=4.0,
    therapy_time=1:48:24*30,
    p_mutation=0.5,
    n_step=n_step, seed=3)


# targeted radiotherapy p_mut=0.5
therapy_condition(models_mut_05, agent_step!, targeted_therapy_step!;
    path=joinpath(dir_name, "targeted_therapy_2Gy_24h_30days_mut_05"),
    dose_per_fraction=2.0,
    therapy_time=1:24:24*30,
    p_mutation=0.5,
    n_step=n_step, seed=4)



##


##

# conventional radiotherapy p_mut=0.1
therapy_condition(models_mut_01, agent_step!, therapy_step!;
    path=joinpath(dir_name, "therapy_2Gy_24h_30days_mut_01"),
    dose_per_fraction=2.0,
    therapy_time=1:24:24*30,
    p_mutation=0.1,
    n_step=n_step, seed=5)

# conventional radiotherapy p_mut=0.3
therapy_condition(models_mut_03, agent_step!, therapy_step!;
    path=joinpath(dir_name, "therapy_2Gy_24h_30days_mut_03"),
    dose_per_fraction=2.0,
    therapy_time=1:24:24*30,
    p_mutation=0.3,
    n_step=n_step, seed=6)
##
