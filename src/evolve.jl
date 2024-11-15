function evolve(sim::SingleDomainSimulation)
    watch_keep_time(sim)
    apply_ics(sim)
    apply_bcs(sim)
    initialize(sim)
    initialize_writing(sim)
    write_step(sim)
    while true
        advance_time(sim)
        if stop_evolve(sim) == true
            break
        end
        watch_keep_time(sim)
        apply_bcs(sim)
        advance(sim)
        write_step(sim)
    end
    finalize_writing(sim)
end

function evolve(sim::MultiDomainSimulation)
    watch_keep_time(sim)
    apply_ics(sim)
    initialize(sim)
    initialize_writing(sim)
    write_step(sim)
    while true
        advance_time(sim)
        if stop_evolve(sim) == true
            break
        end
        watch_keep_time(sim)
        advance(sim)
        write_step(sim)
    end
    finalize_writing(sim)
end

function stop_evolve(sim::SingleDomainSimulation)
    return sim.integrator.time > sim.integrator.final_time
end

function stop_evolve(sim::MultiDomainSimulation)
    return sim.schwarz_controller.time > sim.schwarz_controller.final_time
end

function advance(sim::SingleDomainSimulation)
    solve(sim)
end

function solve_contact(sim::MultiDomainSimulation)
    if sim.schwarz_controller.active_contact == true
        schwarz(sim)
    else
        advance_independent(sim)
    end
end

function advance(sim::MultiDomainSimulation)
    if sim.schwarz_controller.schwarz_contact == false
        schwarz(sim)
        return
    end
    save_stop_solutions(sim)
    solve_contact(sim)
    if sim.failed == true
        return
    end
    was_in_contact = sim.schwarz_controller.active_contact
    detect_contact(sim)
    if sim.schwarz_controller.active_contact ≠ was_in_contact
        if was_in_contact == true
            println("Contact release detected, redoing control step")
        else
            println("Contact initiation detected, redoing control step")
        end
        restore_stop_solutions(sim)
        solve_contact(sim)
    end
end

function apply_ics(sim::SingleDomainSimulation)
    apply_ics(sim.params, sim.model)
end

function apply_ics(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        apply_ics(subsim)
    end
end

function apply_bcs(sim::SingleDomainSimulation)
    apply_bcs(sim.model)
end

function apply_bcs(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        apply_bcs(subsim)
    end
end

function initialize(sim::SingleDomainSimulation)
    initialize(sim.integrator, sim.solver, sim.model)
end

function initialize(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        apply_bcs(subsim)
        initialize(subsim)
    end
    detect_contact(sim)
end

function solve(sim::SingleDomainSimulation)
    solve(sim.integrator, sim.solver, sim.model)
    sim.failed = sim.model.failed
end

function solve(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        solve(subsim)
        if subsim.failed == true
            sim.failed = true
            return
        end
    end
end

function initialize_writing(sim::SingleDomainSimulation)
    initialize_writing(sim.params, sim.integrator, sim.model)
end

function initialize_writing(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        initialize_writing(subsim)
    end
end

function write_step(sim::SingleDomainSimulation)
    write_step(sim.params, sim.integrator, sim.model)
end

function write_step(sim::MultiDomainSimulation)
    time = sim.schwarz_controller.time
    stop = sim.schwarz_controller.stop
    for subsim ∈ sim.subsims
        subsim.integrator.time = subsim.model.time = time
        subsim.integrator.stop = stop
        write_step(subsim)
    end
end

function finalize_writing(sim::SingleDomainSimulation)
    finalize_writing(sim.params)
end

function finalize_writing(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        finalize_writing(subsim)
    end
end

using Printf

function watch_keep_time(sim::SingleDomainSimulation)
    synchronize(sim)
    stop = sim.integrator.stop
    initial_time = sim.integrator.time - sim.integrator.time_step
    final_time = sim.integrator.time
    if stop == 0
        @printf("Initializing run at stop 0 with time = %6.2e\n", final_time)
    else
        @printf(
            "Advancing from stop %d with time = %6.2e to stop %d with time = %6.2e\n",
            stop - 1,
            initial_time,
            stop,
            final_time
        )
    end
end

function watch_keep_time(sim::MultiDomainSimulation)
    synchronize(sim)
    stop = sim.schwarz_controller.stop
    initial_time = sim.schwarz_controller.time - sim.schwarz_controller.time_step
    final_time = sim.schwarz_controller.time
    if stop == 0
        @printf("Initializing run at stop 0 with time = %6.2e\n", final_time)
    else
        @printf(
            "Advancing from stop %d with time = %6.2e to stop %d with time = %6.2e\n",
            stop - 1,
            initial_time,
            stop,
            final_time
        )
    end
end

function synchronize(sim::SingleDomainSimulation)
    sim.model.time = sim.integrator.time
end

function synchronize(sim::MultiDomainSimulation)
    time = sim.schwarz_controller.prev_time
    stop = 0
    for subsim ∈ sim.subsims
        subsim.integrator.time = subsim.model.time = time
        subsim.integrator.stop = stop
    end
end

function advance_time(sim::SingleDomainSimulation)
    next_time = round(sim.integrator.time + sim.integrator.time_step; digits = 12)
    sim.integrator.time = sim.model.time = next_time
    sim.integrator.stop += 1
end

function advance_time(sim::MultiDomainSimulation)
    sim.schwarz_controller.prev_time = sim.schwarz_controller.time
    next_time =
        round(sim.schwarz_controller.time + sim.schwarz_controller.time_step, digits = 12)
    sim.schwarz_controller.time = next_time
    sim.schwarz_controller.stop += 1
end
