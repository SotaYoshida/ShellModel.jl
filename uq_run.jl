include("./src/ShellModel.jl")
using .ShellModel
#include("./src/lanczos_methods.jl")
#include("./src/transit.jl")
#include("./src/input_int_snt.jl")
#include("./src/eigenvector_continuation.jl")

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
function run(num_ev=10)    
    Hs = ["./snts/usdb.snt"]
    intf = "./snts/random_input/usdb.int"

    target_nuc = "Si28"
    tJNs = [ [0,2],[2,2],[4,1],[6,4],[8,3]]

    Erefs = [ [-135.70000,-130.72100], [-127.37200,-126.20400],[-133.92100],
              [-129.42400,-127.90100,-127.11200,-125.49100],
              [-131.08200,-128.81200,-126.53600] ]

    Eexact = [ [-135.861,-131.024], [ -127.669,-126.278],[-133.929],
               [-129.531,-127.852,-126.826,-126.443],
               [-131.254,-128.856,-126.471] ]

    errors = [[0.06764981972995088,0.11807824551232216],[0.14625584538718783,0.28744120626794256],
              [0.05034078995433333],
              [0.10487511604722499,0.23692392369497384, 0.21450202957453257, 0.31933638080241167],
              [0.06909094120663895,0.14593063441913046, 0.26155402077326073]]
    #bnum = 10000 ; itnum = 100000 + bnum
    bnum = 2000; itnum = 10000 + bnum
    #bnum = 10; itnum = 100 + bnum
    #bnum = 50000; itnum = 500000
    solveEC_UQ(Hs,target_nuc,tJNs,Erefs,errors;
               intf=intf,itnum_MCMC=itnum,burnin=bnum,var_M=6.e-2,
               num_replica=1)
    #intf=intf,itnum_MCMC=10000,burnin=1000,var_M=3.e-2)
    #plot_from_history(tJNs,Erefs)
end
run()



