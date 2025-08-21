# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using LinearAlgebra
using Random
using Exodus

Random.seed!(0)

tetra4_a = 1
tetra4_t = sqrt(3) / 2 * tetra4_a
tetra4_h = sqrt(2.0 / 3.0) * tetra4_a
reg_tetra4_coords = [
    -tetra4_t/3 -tetra4_t/3 2*tetra4_t/3 0
    tetra4_a/2 -tetra4_a/2 0.0 0
    0.0 0.0 0.0 tetra4_h
]
tetra4_conn = [1 2 3 4]'
tetra4_base = [1, 2, 3]
tetra4_back_edge = [1, 2]
tetra4_first_node = [1]

hex8_a = 1
reg_hex8_coords = [
    -hex8_a/2 hex8_a/2 hex8_a/2 -hex8_a/2 -hex8_a/2 hex8_a/2 hex8_a/2 -hex8_a/2
    -hex8_a/2 -hex8_a/2 hex8_a/2 hex8_a/2 -hex8_a/2 -hex8_a/2 hex8_a/2 hex8_a/2
    -hex8_a/2 -hex8_a/2 -hex8_a/2 -hex8_a/2 hex8_a/2 hex8_a/2 hex8_a/2 hex8_a/2
]
hex8_conn = [1 2 3 4 5 6 7 8]'
hex8_bottom_face = [1, 2, 3, 4]
hex8_bottom_front_edge = [2, 3]
hex8_first_node = [1]
hex8_top_face = [5, 6, 7, 8]
hex8_top_front_edge = [6, 7]
hex8_top_back_edge = [5, 8]

@testset "tetra4_smooth_reference" begin
    for (sx, sy, sz) in [(1, 1, 1),
                         (1, 1/2, 2),
                         (1/2, 4, 1/2)]
        scaled_tetra4_coords = diagm([sx, sy, sz]) * reg_tetra4_coords
        ref_tetra4_eqvol = Norma.create_smooth_reference("equal volume", Norma.TETRA4, scaled_tetra4_coords);
        @test size(ref_tetra4_eqvol) == (3, 4)
        a1 = norm(ref_tetra4_eqvol[:, 1] - ref_tetra4_eqvol[:, 2])
        a2 = norm(ref_tetra4_eqvol[:, 1] - ref_tetra4_eqvol[:, 3])
        a3 = norm(ref_tetra4_eqvol[:, 1] - ref_tetra4_eqvol[:, 4])
        a4 = norm(ref_tetra4_eqvol[:, 2] - ref_tetra4_eqvol[:, 3])
        a5 = norm(ref_tetra4_eqvol[:, 2] - ref_tetra4_eqvol[:, 4])
        a6 = norm(ref_tetra4_eqvol[:, 3] - ref_tetra4_eqvol[:, 4])
        @test a1 ≈ tetra4_a atol = 1.0e-12
        @test a1 ≈ a2 atol = 1.0e-12
        @test a1 ≈ a3 atol = 1.0e-12
        @test a1 ≈ a4 atol = 1.0e-12
        @test a1 ≈ a5 atol = 1.0e-12
        @test a1 ≈ a6 atol = 1.0e-12
    end

    ref_tetra4_eqedl = Norma.create_smooth_reference("average edge length", Norma.TETRA4, reg_tetra4_coords);
    @test size(ref_tetra4_eqedl) == (3, 4)
    a1 = norm(ref_tetra4_eqedl[:, 1] - ref_tetra4_eqedl[:, 2])
    a2 = norm(ref_tetra4_eqedl[:, 1] - ref_tetra4_eqedl[:, 3])
    a3 = norm(ref_tetra4_eqedl[:, 1] - ref_tetra4_eqedl[:, 4])
    a4 = norm(ref_tetra4_eqedl[:, 2] - ref_tetra4_eqedl[:, 3])
    a5 = norm(ref_tetra4_eqedl[:, 2] - ref_tetra4_eqedl[:, 4])
    a6 = norm(ref_tetra4_eqedl[:, 3] - ref_tetra4_eqedl[:, 4])
    @test a1 ≈ tetra4_a atol = 1.0e-12
    @test a1 ≈ a2 atol = 1.0e-12
    @test a1 ≈ a3 atol = 1.0e-12
    @test a1 ≈ a4 atol = 1.0e-12
    @test a1 ≈ a5 atol = 1.0e-12
    @test a1 ≈ a6 atol = 1.0e-12

    rand_tests = 10
    for _ = 1:rand_tests
        random_tetra4 = reg_tetra4_coords + Random.randn(3, 4) * 0.1 * tetra4_a
        random_tetra4_vol = det(random_tetra4[:, 2:4] .- random_tetra4[:, 1]) / 6.0
        random_tetra4_edl =
            norm(random_tetra4[:, 1] - random_tetra4[:, 2]) +
            norm(random_tetra4[:, 1] - random_tetra4[:, 3]) +
            norm(random_tetra4[:, 1] - random_tetra4[:, 4]) +
            norm(random_tetra4[:, 2] - random_tetra4[:, 3]) +
            norm(random_tetra4[:, 3] - random_tetra4[:, 4]) +
            norm(random_tetra4[:, 4] - random_tetra4[:, 2])

        ref_tetra4_eqvol = Norma.create_smooth_reference("equal volume", Norma.TETRA4, random_tetra4);
        a_eqvol = norm(ref_tetra4_eqvol[:, 1] - ref_tetra4_eqvol[:, 2])
        @test random_tetra4_vol ≈ a_eqvol^3 / 6.0 / sqrt(2.0) atol = 1.0e-12

        ref_tetra4_eqedl =
            Norma.create_smooth_reference("average edge length", Norma.TETRA4, random_tetra4);
        a_eqedl = norm(ref_tetra4_eqedl[:, 1] - ref_tetra4_eqedl[:, 2])
        @test random_tetra4_edl ≈ 6 * a_eqedl atol = 1.0e-12

        @test a_eqedl >= a_eqvol
    end
end

@testset "hex8_smooth_reference" begin
    for (sx, sy, sz) in [(1, 1, 1),
                         (1, 1/2, 2),
                         (1/2, 4, 1/2)]
        scaled_hex8_coords = diagm([sx, sy, sz]) * reg_hex8_coords
        ref_hex8_eqvol = Norma.create_smooth_reference("equal volume", Norma.HEX8, scaled_hex8_coords);
        @test size(ref_hex8_eqvol) == (3, 8)
        a1 = norm(ref_hex8_eqvol[:, 1] - ref_hex8_eqvol[:, 2])
        a3 = norm(ref_hex8_eqvol[:, 1] - ref_hex8_eqvol[:, 4])
        a2 = norm(ref_hex8_eqvol[:, 1] - ref_hex8_eqvol[:, 5])
        a4 = norm(ref_hex8_eqvol[:, 2] - ref_hex8_eqvol[:, 3])
        a5 = norm(ref_hex8_eqvol[:, 2] - ref_hex8_eqvol[:, 6])
        a6 = norm(ref_hex8_eqvol[:, 3] - ref_hex8_eqvol[:, 4])
        a7 = norm(ref_hex8_eqvol[:, 3] - ref_hex8_eqvol[:, 7])
        a8 = norm(ref_hex8_eqvol[:, 4] - ref_hex8_eqvol[:, 8])
        a9 = norm(ref_hex8_eqvol[:, 5] - ref_hex8_eqvol[:, 6])
        a10 = norm(ref_hex8_eqvol[:, 5] - ref_hex8_eqvol[:, 8])
        a11 = norm(ref_hex8_eqvol[:, 6] - ref_hex8_eqvol[:, 7])
        a12 = norm(ref_hex8_eqvol[:, 7] - ref_hex8_eqvol[:, 8])
        @test a1 ≈ hex8_a atol = 1.0e-12
        @test a1 ≈ a2 atol = 1.0e-12
        @test a1 ≈ a3 atol = 1.0e-12
        @test a1 ≈ a4 atol = 1.0e-12
        @test a1 ≈ a5 atol = 1.0e-12
        @test a1 ≈ a6 atol = 1.0e-12
        @test a1 ≈ a7 atol = 1.0e-12
        @test a1 ≈ a8 atol = 1.0e-12
        @test a1 ≈ a9 atol = 1.0e-12
        @test a1 ≈ a10 atol = 1.0e-12
        @test a1 ≈ a11 atol = 1.0e-12
        @test a1 ≈ a12 atol = 1.0e-12
    end

    ref_hex8_eqedl = Norma.create_smooth_reference("average edge length", Norma.HEX8, reg_hex8_coords);
    @test size(ref_hex8_eqedl) == (3, 8)
    a1 = norm(ref_hex8_eqedl[:, 1] - ref_hex8_eqedl[:, 2])
    a3 = norm(ref_hex8_eqedl[:, 1] - ref_hex8_eqedl[:, 4])
    a2 = norm(ref_hex8_eqedl[:, 1] - ref_hex8_eqedl[:, 5])
    a4 = norm(ref_hex8_eqedl[:, 2] - ref_hex8_eqedl[:, 3])
    a5 = norm(ref_hex8_eqedl[:, 2] - ref_hex8_eqedl[:, 6])
    a6 = norm(ref_hex8_eqedl[:, 3] - ref_hex8_eqedl[:, 4])
    a7 = norm(ref_hex8_eqedl[:, 3] - ref_hex8_eqedl[:, 7])
    a8 = norm(ref_hex8_eqedl[:, 4] - ref_hex8_eqedl[:, 8])
    a9 = norm(ref_hex8_eqedl[:, 5] - ref_hex8_eqedl[:, 6])
    a10 = norm(ref_hex8_eqedl[:, 5] - ref_hex8_eqedl[:, 8])
    a11 = norm(ref_hex8_eqedl[:, 6] - ref_hex8_eqedl[:, 7])
    a12 = norm(ref_hex8_eqedl[:, 7] - ref_hex8_eqedl[:, 8])
    @test a1 ≈ hex8_a atol = 1.0e-12
    @test a1 ≈ a2 atol = 1.0e-12
    @test a1 ≈ a3 atol = 1.0e-12
    @test a1 ≈ a4 atol = 1.0e-12
    @test a1 ≈ a5 atol = 1.0e-12
    @test a1 ≈ a6 atol = 1.0e-12
    @test a1 ≈ a7 atol = 1.0e-12
    @test a1 ≈ a8 atol = 1.0e-12
    @test a1 ≈ a9 atol = 1.0e-12
    @test a1 ≈ a10 atol = 1.0e-12
    @test a1 ≈ a11 atol = 1.0e-12
    @test a1 ≈ a12 atol = 1.0e-12

    rand_tests = 10
    for _ = 1:rand_tests
        # We can't just perturb all the vertices of the initial hexahedron, as it
        # would make some faces non-planar. Instead, all but two opposed vertices
        # (1 and 7) are perturbed. The remaining vertices are then projected onto
        # the planes defined by the other six vertices using Cramer's rule.
        random_hex8 = copy(reg_hex8_coords)
        random_hex8[:, [2, 3, 4, 5, 6, 8]] .+= Random.randn(3, 6) * 0.1 * hex8_a
        u265 = cross(random_hex8[:, 5] - random_hex8[:, 2],
                      random_hex8[:, 6] - random_hex8[:, 2])
        d265 = dot(u265, random_hex8[:, 2])
        u243 = cross(random_hex8[:, 3] - random_hex8[:, 2],
                      random_hex8[:, 4] - random_hex8[:, 2])
        d243 = dot(u243, random_hex8[:, 2])
        u458 = cross(random_hex8[:, 8] - random_hex8[:, 4],
                      random_hex8[:, 5] - random_hex8[:, 4])
        d458 = dot(u458, random_hex8[:, 4])
        random_hex8[:, 1] = (d265*cross(u243, u458) +
                             d243*cross(u458, u265) +
                             d458*cross(u265, u243)) /
                            dot(u265, cross(u243, u458))
        u348 = cross(random_hex8[:, 8] - random_hex8[:, 3],
                    random_hex8[:, 4] - random_hex8[:, 3])
        d348 = dot(u348, random_hex8[:, 3])
        u586 = cross(random_hex8[:, 6] - random_hex8[:, 5],
                      random_hex8[:, 8] - random_hex8[:, 5])
        d586 = dot(u586, random_hex8[:, 5])
        u236 = cross(random_hex8[:, 6] - random_hex8[:, 2],
                      random_hex8[:, 3] - random_hex8[:, 2])
        d236 = dot(u236, random_hex8[:, 2])
        random_hex8[:, 7] = (d348*cross(u586, u236) +
                             d586*cross(u236, u348) +
                             d236*cross(u348, u586)) /
                            dot(u348, cross(u586, u236))

        # The implementation of the equal volume smoothing method relies on dividing
        # a hexahedron into five tetrahedra, so in the test we compute the volume
        # with a different method by dividing the hexahedron into three pyramids
        # First pyramid: base 2, 3, 7, 6 apex 1
        # Second pyramid: base 3, 4, 8, 7 apex 1
        # Third pyramid: base 5, 6, 7, 8 apex 1
        # pyr_1_base_area = det(random_hex8[:, (2, 7)] .- random_hex8[:, (3, 6)]) / 2.0
        pyr_1_volume = dot(cross(random_hex8[:, 7] - random_hex8[:, 2],
                                 random_hex8[:, 3] - random_hex8[:, 6]),
                           random_hex8[:, 1] - random_hex8[:, 2]) / 6.0
        pyr_2_volume = dot(cross(random_hex8[:, 8] - random_hex8[:, 3],
                                 random_hex8[:, 4] - random_hex8[:, 7]),
                            random_hex8[:, 1] - random_hex8[:, 3]) / 6.0
        pyr_3_volume = dot(cross(random_hex8[:, 7] - random_hex8[:, 5],
                                 random_hex8[:, 6] - random_hex8[:, 8]),
                            random_hex8[:, 1] - random_hex8[:, 6]) / 6.0
        random_hex8_vol = pyr_1_volume + pyr_2_volume + pyr_3_volume

        random_hex8_edl =
            norm(random_hex8[:, 1] - random_hex8[:, 2]) +
            norm(random_hex8[:, 1] - random_hex8[:, 4]) +
            norm(random_hex8[:, 1] - random_hex8[:, 5]) +
            norm(random_hex8[:, 2] - random_hex8[:, 3]) +
            norm(random_hex8[:, 2] - random_hex8[:, 6]) +
            norm(random_hex8[:, 3] - random_hex8[:, 4]) +
            norm(random_hex8[:, 3] - random_hex8[:, 7]) +
            norm(random_hex8[:, 4] - random_hex8[:, 8]) +
            norm(random_hex8[:, 5] - random_hex8[:, 6]) +
            norm(random_hex8[:, 5] - random_hex8[:, 8]) +
            norm(random_hex8[:, 6] - random_hex8[:, 7]) +
            norm(random_hex8[:, 7] - random_hex8[:, 8])

        ref_hex8_eqvol = Norma.create_smooth_reference("equal volume", Norma.HEX8, random_hex8);
        a_eqvol = norm(ref_hex8_eqvol[:, 1] - ref_hex8_eqvol[:, 2])
        @test random_hex8_vol ≈ a_eqvol^3 atol = 1.0e-12

        ref_hex8_eqedl =
            Norma.create_smooth_reference("average edge length", Norma.HEX8, random_hex8);
        a_eqedl = norm(ref_hex8_eqedl[:, 1] - ref_hex8_eqedl[:, 2])
        @test random_hex8_edl ≈ 12 * a_eqedl atol = 1.0e-12

        @test a_eqedl >= a_eqvol
    end


end

base_params = Dict{String,Any}(
    "type" => "single",
    "model" => Dict{String,Any}(
        "material" => Dict{String,Any}(
            "elastic" => Dict{String,Any}(
                "model" => "seth-hill",
                "m" => 2,
                "n" => 2,
                "density" => 1.0,
            ),
            "blocks" => Dict{String,Any}("block" => "elastic"),
        ),
        "type" => "mesh smoothing",
    ),
    "Exodus output interval" => 0,
    "CSV output interval" => 0,
    "time integrator" => Dict{String,Any}(
        "type" => "quasi static",
        "initial time" => 0.0,
        "final time" => 1.0,
        "time step" => 1.0e-1,
    ),
    "solver" => Dict{String,Any}(
        "step" => "steepest descent",
        "type" => "steepest descent",
        "minimum iterations" => 1,
        "maximum iterations" => 64,
        "absolute tolerance" => 1.0e-8,
        "relative tolerance" => 1.0e-12,
        "step length" => 1.0e-3,
        "use line search" => true,
        "line search backtrack factor" => 0.5,
        "line search decrease factor" => 1.0e-04,
        "line search maximum iterations" => 16,
    ),
)


input_mesh_file = "tetra4_smoothing.g"
output_mesh_file = "tetra4_smoothing.e"

@testset "tetra4_smoothing" begin
    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by perturbating the xy coordinates of the top vertex,
    # corresponding to a pure shear deformation.
    # The resulting tetrahedron's volume is unchanged, so the smoothing
    # procedure should find the original coordinates of the top vertex when the
    # base is fixed, using a deviatoric energy term only.
    top_xy_disp = Random.rand(2) * tetra4_a * 0.1
    tetra4_coords =
        reg_tetra4_coords + [
            0.0 0.0 0.0 top_xy_disp[1]
            0.0 0.0 0.0 top_xy_disp[2]
            0.0 0.0 0.0 0.0
        ]

    node_sets = Dict{String,Vector}("base" => tetra4_base)

    num_dim, num_nodes = size(tetra4_coords)
    num_elems = size(tetra4_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        tetra4_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end

        tetra4_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tetra4_init)
        write_coordinates(tetra4_exo, tetra4_coords)
        write_block(tetra4_exo, 1, "TETRA4", Matrix{Int32}(tetra4_conn))
        write_name(tetra4_exo, Block, 1, "block")
        for (i, (node_set_name, node_set_nodes)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(node_set_nodes))
            write_set(tetra4_exo, node_set)
            write_name(tetra4_exo, node_set, node_set_name)
        end

        close(tetra4_exo)

        tetra4_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "shear_tetra4_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        tetra4_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tetra4_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        tetra4_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(tetra4_test_params)

        @test sim.integrator.displacement ≈ vec(reg_tetra4_coords - tetra4_coords) atol = tetra4_a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end

    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by perturbating the z coordinates of the top vertex,
    # corresponding to a uniaxial z deformation. An additional uniaxial x displacement
    # is applied to the top and front vertices to preserve the volume of the tetrahedron.
    top_z_disp = Random.rand() * tetra4_a * 0.1
    front_x_disp = tetra4_t * tetra4_h / (tetra4_h + top_z_disp) - tetra4_t
    top_x_disp = front_x_disp / 3
    tetra4_coords =
        reg_tetra4_coords + [
            0.0 0.0 front_x_disp top_x_disp
            0.0 0.0 0.0 0
            0.0 0.0 0.0 top_z_disp
        ]

    node_sets = Dict{String,Vector}("base" => tetra4_base, "back edge" => tetra4_back_edge)

    num_dim, num_nodes = size(tetra4_coords)
    num_elems = size(tetra4_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        tetra4_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end

        tetra4_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tetra4_init)
        write_coordinates(tetra4_exo, tetra4_coords)
        write_block(tetra4_exo, 1, "TETRA4", Matrix{Int32}(tetra4_conn))
        write_name(tetra4_exo, Block, 1, "block")
        for (i, (node_set_name, node_set_nodes)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(node_set_nodes))
            write_set(tetra4_exo, node_set)
            write_name(tetra4_exo, node_set, node_set_name)
        end

        close(tetra4_exo)

        tetra4_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "bulk_tetra4_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}(
                            "node set" => "back edge",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "back edge",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        tetra4_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tetra4_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        tetra4_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(tetra4_test_params)

        @test sim.integrator.displacement ≈ vec(reg_tetra4_coords - tetra4_coords) atol = tetra4_a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end
    
    # This test creates a tetrahedron from a regular tetrahedron with a base
    # triangle in the xy-plane by applying equal uniaxial xy deformation gradients
    # to all vertices and then adjusting the z coordinate of the top vertex to
    # preserve the total edge length of the tetrahedron.
    # Using the average edge length smoothing method, the smoothing procedure
    # should retrieve the original tetrahedron coordinates.

    xy_strain_init = Random.rand() * 0.1
    tetra4_coords =
        diagm([1+xy_strain_init, 1+xy_strain_init, 1]) * (reg_tetra4_coords .- reg_tetra4_coords[:, 1]) .+ 
        reg_tetra4_coords[:, 1]
    total_edge_l = norm(tetra4_coords[:, 1] - tetra4_coords[:, 2]) +
                   norm(tetra4_coords[:, 1] - tetra4_coords[:, 3]) +
                   norm(tetra4_coords[:, 1] - tetra4_coords[:, 4]) +
                   norm(tetra4_coords[:, 2] - tetra4_coords[:, 3]) +
                   norm(tetra4_coords[:, 2] - tetra4_coords[:, 4]) +
                   norm(tetra4_coords[:, 3] - tetra4_coords[:, 4])
    xyz_scale = 6*tetra4_a / total_edge_l
    tetra4_coords = xyz_scale * (tetra4_coords .- tetra4_coords[:, 1]) .+ tetra4_coords[:, 1]

    node_sets = Dict{String,Vector}("base" => tetra4_base, "back edge" => tetra4_back_edge, "first node" => tetra4_first_node)

    num_dim, num_nodes = size(tetra4_coords)
    num_elems = size(tetra4_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        tetra4_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end

        tetra4_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", tetra4_init)
        write_coordinates(tetra4_exo, tetra4_coords)
        write_block(tetra4_exo, 1, "TETRA4", Matrix{Int32}(tetra4_conn))
        write_name(tetra4_exo, Block, 1, "block")
        for (i, (node_set_name, node_set_nodes)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(node_set_nodes))
            write_set(tetra4_exo, node_set)
            write_name(tetra4_exo, node_set, node_set_name)
        end

        close(tetra4_exo)

        tetra4_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "bulk_tetra4_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}(
                            "node set" => "back edge",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "first node",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "base",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        tetra4_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        tetra4_test_params["model"]["material"]["elastic"]["bulk modulus"] = 1.0
        tetra4_test_params["model"]["smooth reference"] = "average edge length"

        sim = Norma.run(tetra4_test_params)

        @test sim.integrator.displacement ≈ vec(reg_tetra4_coords - tetra4_coords) atol = tetra4_a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end
end

@testset "hex8_smoothing" begin
    # This test creates a hexahedron from a cube by shifting its top face along
    # the xy plane, preserving its total volume.
    top_xy_disp = Random.rand(2) * hex8_a * 0.1
    hex8_coords = copy(reg_hex8_coords)
    hex8_coords[1:2, hex8_top_face] .+= top_xy_disp

    node_sets = Dict{String,Vector}("bottom_face" => hex8_bottom_face)

    num_dim, num_nodes = size(hex8_coords)
    num_elems = size(hex8_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        hex8_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end

        hex8_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", hex8_init)
        write_coordinates(hex8_exo, hex8_coords)
        write_block(hex8_exo, 1, "HEX8", Matrix{Int32}(hex8_conn))
        write_name(hex8_exo, Block, 1, "block")
        for (i, (node_set_name, node_set_nodes)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(node_set_nodes))
            write_set(hex8_exo, node_set)
            write_name(hex8_exo, node_set, node_set_name)
        end

        close(hex8_exo)

        hex8_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "shear_hex8_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}(
                            "node set" => "bottom_face",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "bottom_face",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "bottom_face",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        hex8_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        hex8_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        hex8_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(hex8_test_params)

        @test sim.integrator.displacement ≈ vec(reg_hex8_coords - hex8_coords) atol = hex8_a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end

    # This test creates a hexadron by shifting the top front edge and back front
    # edge along z in opposite directions with same magnitude, preserving its
    # total volume.
    top_front_z_disp = Random.rand() * hex8_a * 0.1
    hex8_coords = copy(reg_hex8_coords)
    hex8_coords[3, hex8_top_front_edge] .+= top_front_z_disp
    hex8_coords[3, hex8_top_back_edge] .-= top_front_z_disp

    node_sets = Dict{String,Vector}("bottom_face" => hex8_bottom_face)

    num_dim, num_nodes = size(hex8_coords)
    num_elems = size(hex8_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        hex8_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end

        hex8_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", hex8_init)
        write_coordinates(hex8_exo, hex8_coords)
        write_block(hex8_exo, 1, "HEX8", Matrix{Int32}(hex8_conn))
        write_name(hex8_exo, Block, 1, "block")
        for (i, (node_set_name, node_set_nodes)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(node_set_nodes))
            write_set(hex8_exo, node_set)
            write_name(hex8_exo, node_set, node_set_name)
        end

        close(hex8_exo)


        hex8_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "shear_hex8_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}(
                            "node set" => "bottom_face",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "bottom_face",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "bottom_face",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        hex8_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        hex8_test_params["model"]["material"]["elastic"]["bulk modulus"] = 0.0
        hex8_test_params["model"]["smooth reference"] = "equal volume"

        sim = Norma.run(hex8_test_params)

        @test sim.integrator.displacement ≈ vec(reg_hex8_coords - hex8_coords) atol = hex8_a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end
    
    # This test creates a hexahedron by applying uniaxial strain deformation to
    # the top face along the x and y directions. The top face height is then
    # adjusted to preserve the total edge length of the hexahedron.
    # Using the average edge length smoothing method, the smoothing procedure
    # should retrieve the original hexahdron coordinates.

    top_xy_strains = Random.rand(2) * 0.1
    hex8_coords = copy(reg_hex8_coords)
    hex8_coords[1:2, hex8_top_face] .*= (top_xy_strains .+ 1.0)
    dz = sqrt((1-top_xy_strains[1]/2-top_xy_strains[2]/2)^2 - (top_xy_strains[1]^2 + top_xy_strains[2]^2)/4.0) * hex8_a
    hex8_coords[3, hex8_top_face] .= dz - hex8_a / 2.0

    node_sets = Dict{String,Vector}("bottom face" => hex8_bottom_face, "bottom front edge" => hex8_bottom_front_edge, "first node" => hex8_first_node)

    num_dim, num_nodes = size(hex8_coords)
    num_elems = size(hex8_conn, 2)
    num_elem_blks = 1
    num_node_sets = length(node_sets)
    num_side_sets = 0

    try
        hex8_init = Initialization{Int32}(
            Int32(num_dim),
            Int32(num_nodes),
            Int32(num_elems),
            Int32(num_elem_blks),
            Int32(num_node_sets),
            Int32(num_side_sets),
        )

        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end

        hex8_exo = ExodusDatabase{Int32,Int32,Int32,Float64}(input_mesh_file, "w", hex8_init)
        write_coordinates(hex8_exo, hex8_coords)
        write_block(hex8_exo, 1, "HEX8", Matrix{Int32}(hex8_conn))
        write_name(hex8_exo, Block, 1, "block")
        for (i, (node_set_name, node_set_nodes)) in enumerate(node_sets)
            node_set = NodeSet(i, Vector{Int32}(node_set_nodes))
            write_set(hex8_exo, node_set)
            write_name(hex8_exo, node_set, node_set_name)
        end

        close(hex8_exo)

        hex8_test_params = merge(
            base_params,
            Dict{String,Any}(
                "name" => "bulk_hex8_smoothing",
                "input mesh file" => input_mesh_file,
                "output mesh file" => output_mesh_file,
                "boundary conditions" => Dict{String,Any}(
                    "Dirichlet" => [
                        Dict{String,Any}(
                            "node set" => "bottom front edge",
                            "component" => "x",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "first node",
                            "component" => "y",
                            "function" => "0.0",
                        ),
                        Dict{String,Any}(
                            "node set" => "bottom face",
                            "component" => "z",
                            "function" => "0.0",
                        ),
                    ],
                ),
            ),
        )
        hex8_test_params["model"]["material"]["elastic"]["shear modulus"] = 1.0
        hex8_test_params["model"]["material"]["elastic"]["bulk modulus"] = 1.0
        hex8_test_params["model"]["smooth reference"] = "average edge length"

        sim = Norma.run(hex8_test_params)

        @test sim.integrator.displacement ≈ vec(reg_hex8_coords - hex8_coords) atol = hex8_a*1.0e-6
    finally
        if isfile(input_mesh_file)
            rm(input_mesh_file)
        end
        if isfile(output_mesh_file)
            rm(output_mesh_file)
        end
    end
end
