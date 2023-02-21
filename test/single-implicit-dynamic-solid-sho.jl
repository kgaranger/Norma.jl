cp("../examples/single/implicit-dynamic-solid/sho/sho.yaml", "sho.yaml", force=true)
cp("../examples/single/implicit-dynamic-solid/sho/sho.g", "sho.g", force=true)
integrator, solver, model = Norma.run("sho.yaml")
rm("sho.yaml")
rm("sho.g")
rm("sho.e")
displacement = integrator.displacement[end]
velocity = integrator.velocity[end]
acceleration = integrator.acceleration[end]
@testset "single-implicit-dynamic-solid-sho" begin
    @test displacement ≈ 0.0 atol = 2.0e-03
    @test velocity ≈ -1.0 rtol = 1.0e-05
    @test acceleration ≈ 0.0 atol = 2.0e-03
end