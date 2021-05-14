include("./src/ShellModel.jl")
using .ShellModel

# similar to samplerun()
function run()
    @time main("./snts/x_mass.snt","Be8",1,[];is_show=true);println("\n")
    @time main("./snts/usdb.snt","Si28",10,[];is_show=true);println("\n")
    #@time main("./snts/usdb.snt","Mg24",2,[];is_show=true);println("\n")
    @time main("./snts/gxpf1a.snt","Cr48",10,[];is_show=true);println("\n")
end
run()
