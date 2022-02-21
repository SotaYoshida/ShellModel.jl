include("./src/ShellModel.jl")
using .ShellModel

# similar to samplerun()
function run()
    ### for compilation
    println("For JIT compilation...")
    @time main("./snts/x_mass.snt","Be8",[[0,5]];is_show=true)
    @time main("./snts/x_mass.snt","Be8",[[0,5]];is_block=true,q=3)
    println("Done: JIT compilation")
    ### input
    tJ = 0 
    num_ev_target = 10
    sntf = "./snts/usdb.snt"; target_nuc = "Si28"
    
    println("\n### sample (a). 10-lowest states for $target_nuc")    
    JPNs = [10]
    @time main(sntf,target_nuc,[10];is_show=true)

    println("\n### sample (b). 10-lowest states with J=0")    
    @time main(sntf,target_nuc,[[0,10]])

    println("\n### sample (c). calc. EC estimates for 10 J=$tJ states ")    
    Hs = ["./snts/usdb.snt"]
    tJNs = [ [tJ,num_ev_target] ]
    solveEC(Hs,target_nuc,tJNs;
            verbose=true, exact_logf = "runlog_sigma1_valid.dat")    

    println("\n### sample (d).  (b). with the preprocessing")
    println("Block Lanczos method with n_block=4")
    @time main(sntf,target_nuc,[[0,4]];is_block=true,q=4)
    println("\n+ preprocessing")
    @time main(sntf,target_nuc,[[0,4]];is_block=true,q=4,
               in_wf="./appwavs/"*target_nuc*"_usdb_j"*string(tJ)*".appwav")

end
run()
