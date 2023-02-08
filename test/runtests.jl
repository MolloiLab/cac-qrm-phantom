using DrWatson, Test
@quickactivate "cac-qrm-phantom"

# Here you include files using `srcdir`
include(srcdir("masks.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "masks tests" begin
    @test 1 + 1 == 2
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
