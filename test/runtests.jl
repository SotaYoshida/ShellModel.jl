using ShellModel
using Test

@testset "ShellModel.jl" begin
    main("../snts/x_mass.snt","Be8",[[0,1]];mdimmode=true)
end
