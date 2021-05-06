include("./src/shellmodel_main.jl")
include("./src/lanczos_methods.jl")
include("./src/transit.jl")
include("./src/input_int_snt.jl")
include("./src/eigenvector_continuation.jl")

"""
Arguments for main: digonalize the model-space Hamiltonian
sntf:       input interaction file (.snt)
target_nuc: target nucleus
num_ev:     number of eigenstates to be evaluated
target_J    J*2 for the target states 0=>0, 1=>2, 3/2=>3, ...
(Optional)
save_wav=false  # to specify the option to save wavefunctions
is_show = true  # to show elapsed time & allocations
lm = 100        # number of Lanczos vectors
ls = 15         # number of vectors to be used for Thick-Restart
tol= 1.e-6      # tolerance for convergence check in the Lanczos methods
"""
function samplerun(num_ev=10)
    ### for compilation
    println("For JIT compilation...")
    main("./snts/x_mass.snt","Be8",5,[0])
    main("./snts/x_mass.snt","Be8",5,[0];is_block=true,q=3)
    println("Done: JIT compilation")

    ### input
    tJ = 0 # target J
    tJs =[ tJ ]   
    num_ev = 4; n_block = 4
    num_ev_target = 10
    sntf = "./snts/usdb.snt"; target_nuc = "Si28"
    ###
    
    println("\n### sample (a). 10-lowest states for $target_nuc")    
    @time main(sntf,target_nuc,10,[])

    println("\n### sample (b). 10-lowest states with J=$tJs")    
    @time main(sntf,target_nuc,10,tJs)

    println("\n### sample (c). calc. EC estimates for $num_ev_target J=$tJs states ")    
    Hs = ["./snts/usdb.snt"]
    tJNs = [ [tJ,num_ev_target] ]
    solveEC(Hs,target_nuc,tJNs;
            verbose=true, exact_logf = "runlog_sigma1_valid.dat")    

    println("\n### sample (d).  (b). with the preprocessing")
    println("Block Lanczos method with n_block=$n_block")
    @time main(sntf,target_nuc,num_ev,tJs;is_block=true,q=n_block)
    println("\n+ preprocessing")
    @time main(sntf,target_nuc,num_ev,tJs;is_block=true,q=n_block,
               in_wf="./appwavs/"*target_nuc*"_usdb_j"*string(tJ)*".appwav")

end
samplerun()



