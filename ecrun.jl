include("./src/ShellModel.jl")
using .ShellModel

##########################
### solve EC for approx. eigenpairs
##########################
"""
calc_moment=true   # calc. mu&Q-moment 
write_appwav0false # to write .appwav you need to make sample w.f.
"""
function solve(write_appwav=false)
    sntpath = "snts/random_input/validsnts_sigma1/sdshl/"
    Hs = [  sntpath*"tmp_$i"*".snt" for i = 0:99 ]
    #tpath = "/Users/Sota/Desktop/develop_ShellModel.jl/"
    tpath = "/Users/sotauu/Desktop/develop_ShellModel.jl/"
    elogf = "runlog_sigma1_valid.dat"
        
    num_ev_target = 1
    verbose = false # to show all Energy(EC)
    tJNs = [ [tJ,num_ev_target] for tJ in 0:2:8]
    for target_nuc in ["Mg24","Al26","Si28"]
        solveEC(Hs,target_nuc,tJNs;wpath=tpath,verbose=verbose,exact_logf=elogf)
    end

    tJNs = [ [tJ,num_ev_target] for tJ in 1:2:9]
    for target_nuc in ["Mg25"]
        solveEC(Hs,target_nuc,tJNs;wpath=tpath,verbose=verbose,exact_logf=elogf)
    end    
end

solve()
