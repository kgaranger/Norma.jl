# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using YAML

include("simulation_def.jl")
include("model.jl")
include("time_integrator.jl")
include("solver.jl")
include("schwarz.jl")

function create_simulation(params::Dict{String,Any}, name::String)
    params["name"] = name
    sim_type = params["type"]
    if sim_type == "single"
        sim = SingleDomainSimulation(params)
        create_bcs(sim)
        return sim
    elseif sim_type == "multi"
        sim = MultiDomainSimulation(params)
        create_bcs(sim)
        return sim
    else
        error("Unknown type of simulation: ", sim_type)
    end
end

function create_simulation(input_file::String)
    println("Reading simulation file: ", input_file)
    params = YAML.load_file(input_file; dicttype=Dict{String,Any})
    return create_simulation(params, input_file)
end

function create_bcs(sim::SingleDomainSimulation)
    boundary_conditions = create_bcs(sim.params)
    for bc in boundary_conditions
        if isa(bc, SMDirichletInclined) || isa(bc, SMContactSchwarzBC)
            sim.model.inclined_support = true
            break
        end
    end
    sim.model.boundary_conditions = boundary_conditions
    
end

function create_bcs(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        create_bcs(subsim)
    end
    pair_schwarz_bcs(sim)
end

function SingleDomainSimulation(params::Dict{String,Any})
    name = params["name"]
    input_mesh_file = params["input mesh file"]
    output_mesh_file = params["output mesh file"]
    rm(output_mesh_file, force=true)
    input_mesh = Exodus.ExodusDatabase(input_mesh_file, "r")
    Exodus.copy(input_mesh, output_mesh_file)
    output_mesh = Exodus.ExodusDatabase(output_mesh_file, "rw")
    params["output_mesh"] = output_mesh
    params["input_mesh"] = input_mesh
    model = create_model(params)
    integrator = create_time_integrator(params,model)
    solver = create_solver(params,model)
    failed = false
    return SingleDomainSimulation(name, params, integrator, solver, model, failed)
end

function MultiDomainSimulation(params::Dict{String,Any})
    name = params["name"]
    domain_names = params["domains"]
    subsims = Vector{SingleDomainSimulation}()
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    same_step = get(params, "same time step for domains", false)
    exodus_interval = get(params, "Exodus output interval", 1)
    csv_interval = get(params, "CSV output interval", 0)
    sim_type = "none"
    subsim_name_index_map = Dict{String,Int64}()
    subsim_index = 1
    for domain_name ∈ domain_names
        println("Reading subsimulation file: ", domain_name)
        subparams = YAML.load_file(domain_name; dicttype=Dict{String,Any})
        subparams["name"] = domain_name
        subparams["time integrator"]["initial time"] = initial_time
        subparams["time integrator"]["final time"] = final_time
        if same_step == true
            subparams["time integrator"]["time step"] = time_step
        end
        subparams["Exodus output interval"] = exodus_interval
        subparams["CSV output interval"] = csv_interval
        subsim = SingleDomainSimulation(subparams)
        params[domain_name] = subsim.params
        subsim_type =
            get_analysis_type(subsim.integrator) * " " * subparams["model"]["type"]
        if sim_type == "none"
            sim_type = subsim_type
        elseif subsim_type ≠ sim_type && subsim_type != "dynamic linear opinf rom"
            error("Multidomain subdomains must all have the same physics")
        end
        push!(subsims, subsim)
        subsim_name_index_map[domain_name] = subsim_index
        subsim_index += 1
    end
    params["subdomains type"] = sim_type
    schwarz_controller = create_schwarz_controller(params)
    failed = false
    sim = MultiDomainSimulation(
        name,
        params,
        schwarz_controller,
        subsims,
        subsim_name_index_map,
        failed
    )
    for subsim ∈ sim.subsims
        subsim.params["global_simulation"] = sim
    end
    return sim
end

function get_block_connectivity(mesh::ExodusDatabase, blk_id::Integer)
    _, num_elems, num_nodes, _, _, _ = Exodus.read_block_parameters(mesh, Int32(blk_id))
    conn = Exodus.read_block_connectivity(mesh, Int32(blk_id), num_elems * num_nodes)
    return reshape(conn, (num_elems, num_nodes))
end
