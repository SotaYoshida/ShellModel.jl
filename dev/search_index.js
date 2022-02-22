var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ShellModel","category":"page"},{"location":"#ShellModel","page":"Home","title":"ShellModel","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ShellModel.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ShellModel]","category":"page"},{"location":"#ShellModel.HbitT1-NTuple{6, Any}","page":"Home","title":"ShellModel.HbitT1","text":"function HbitT1(p_sps::Array{Array{Int64,1}},n_sps::Array{Array{Int64,1}},\n            mstates_p::Array{Array{Int64,1},1},mstates_n::Array{Array{Int64,1},1},\n            labels::Array{Array{Array{Int64,1},1},1},TBMEs::Array{Array{Float64,1}})\n\nmake bit representation of T=1 (proton-proton&neutron-neutron) interactions for each {m_z}\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.Hbitpn","page":"Home","title":"ShellModel.Hbitpn","text":"Hbitpn(p_sps::Array{Array{Int64,1}},n_sps::Array{Array{Int64,1}},\n       mstates_p::Array{Array{Int64,1},1},mstates_n::Array{Array{Int64,1},1},\n       labels::Array{Array{Int64,1}},TBMEs::Array{Float64,1},zeroME=false)\n\nmake bit representation of T=0 (proton-neutron) interactions for each {m_z}\n\n\n\n\n\n","category":"function"},{"location":"#ShellModel.all_perm!-Tuple{Int64, Int64, Array{Vector{Bool}, N} where N}","page":"Home","title":"ShellModel.all_perm!","text":"all_perm!(ln::Int64,num_valence::Int64,occs::Array{Array{Bool,1}})\n\nmake all possible permutation of 'bits'\n\nExample: If 2 protons and 1 neutron are in the 0p-shell space, valence orbits(0p1/2,0p3/2) => -1/2, 1/2, -3/2, -1/2, 1/2, 3/2\n\nconfigurations are represented like:\n\nproton: 000011, 000101, ..., 110000\n\nneutron: 000001, 000010, ..., 100000\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.def_mstates-Tuple{Any, Any}","page":"Home","title":"ShellModel.def_mstates","text":"defmstates(psps,n_sps)\n\nto define the single particle states specified by m_z\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.embedding_crsp-NTuple{6, Any}","page":"Home","title":"ShellModel.embedding_crsp","text":"To store correspondance between index for vectors w/ and w/o truncation\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.main-Tuple{Any, Any, Any}","page":"Home","title":"ShellModel.main","text":"main(sntf,targetnuc,targetJNs;      savewav=false,q=1,isblock=false,isshow=false,numhistory=3,lm=100,ls=20,tol=1.e-6,      inwf=\"\",mdimmode=false,calcmoment = false,gfactors = [1.0,0.0,5.586,-3.826],effcharge=[1.5,0.5])\n\nDigonalize the model-space Hamiltonian \n\nArguments:\n\nsntf:       path to input interaction file (.snt fmt)\ntarget_nuc: target nucleus\nnum_ev:     number of eigenstates to be evaluated\ntarget_J:   target total J (specified by e.g. [0]). Set to [] if you want lowest states with any J.   Note that J should be doubled (J=0=>[0], J=1/2=>[1], J=1=>[2],...) \n\nOptional arguments:\n\nq=1              block size for Block-Lanczos methods \nis_block=false   whether or not to use Block algorithm \nsave_wav=false   whether or not to save wavefunction file \nis_show = true   to show elapsed time & allocations \nlm = 100         number of Lanczos vectors to store \nls = 15          number of vectors to be used for Thick-Restart \ntol= 1.e-6       tolerance for convergence check in the Lanczos method \nin_wf=\"\"      path to initial w.f. (for preprocessing) \nmdimmode=false   true => calculate only the M-scheme dimension\ncalc_moment=false  true => calculate mu&Q moments \ngfactors=[1.0,0.0,5.586,-3.826] angular momentum and spin g-factors \neffcgarge=[1.5,0.5] effective charges \n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.occ-Tuple{Any, Any, Any, Bool}","page":"Home","title":"ShellModel.occ","text":"occ function to embedding truncated vector to full space one\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.occ-Tuple{Any, Any, Any}","page":"Home","title":"ShellModel.occ","text":"occ(modelspace,Mtot,truncation_params)\n\nprepare bit representations of proton/neutron Slater determinants => pbits/nbits\n\nMps/Mns: total M for proton/neutron \"blocks\"\n\nFor 6Li in the p shell and M=0, Mps = [-3,-1,1,3] & Mns = [3,1,-1,-3]   blocks => [ (Mp,Mn)=(-3,3),(Mp,Mn)=(-1,1),...]  \n\ntdims: array of cumulative number of M-scheme dimensions for \"blocks\"  \n\ntdims =[ # of possible configurations of (-3,3),            # of possible configurations of (-1,1),...]  \n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.readsnt-Tuple{Any, Any, Any}","page":"Home","title":"ShellModel.readsnt","text":"readsnt(sntf,Anum)\n\nTo read interaction file in \".snt\" format.\n\nsntf: path to the interaction file\nAnum: mass number (used for \"scaling\" of TBMEs)\n\nnote: Note\nThe current version supports only \"properly ordered\" interaction file in .snt format.A .snt file can be ordered to be a<=b,c<=d,a<=c for V(abcd;J) by the Python script \"ShellModel.jl/src/makeorderedsnt.py\"(, which will be replaced by Julia implementation...).\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.samplerun-Tuple{}","page":"Home","title":"ShellModel.samplerun","text":"samplerun()\n\nSample scripts to calculate\n\n(a) 10 lowest states of 28Si in sd shell\n(b) 10 lowest states of 28Si with J=0\n(c) EC estimates of 10 lowest J=0 states of 28Si\n(d) (b) with the preprocessing\n\nnote: Note\nTo run samplerun(), you need a copy of ShellModel.jl in your environment. Try e.g.,git clone https://github.com/SotaYoshida/ShellModel.jland then, execute samplerun() in the repository:>julia using ShellModel\n>julia samplerun()or julia -t 12 sample_run.jl\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.solveEC-Tuple{Any, Any, Any}","page":"Home","title":"ShellModel.solveEC","text":"solveEC(Hs,target_nuc,tJNs;\n        write_appwav=false,verbose=false,calc_moment=true,wpath=\"./\",is_show=false,\n        gfactors = [1.0,0.0,5.586,-3.826],effcharge=[1.5,0.5],exact_logf=\"\")\n\nTo solve EC (generalized eigenvalue problem) to approximate the eigenpairs for a given interaction.\n\nH vecv = lambda N vecv\n\nTransition densities and overlap matrix for H and N are read from \"tdmat/\" directory (To be changed to more flexible)\n\nArguments:\n\nHs:array of paths to interaction files (.snt fmt)\ntarget_nuc: target nucleus\ntJNs:array of target total J (doubled) and number of eigenstates to be evaluated   e.g., [ [0,5], [2,5] ], when you want five lowest J=0 and J=1 states.\n\nOptional arguments:\n\nwrite_appwav=false:write out the approximate wavefunction\nverbose=false:to print (stdout) approx. energies for each interaction\ncalc_moment=true: to calculate mu&Q moments\nwpath=\"./\": path to sample eigenvectors to construct approx. w.f.\nis_show=false: to show TimerOutput\ngfactors=[1.0,0.0,5.586,-3.826]: g-factors to evaluate moments\neffcharge=[1.5,0.5]:effective charges to evaluate moments\n\nOptional arguments for author's own purpose\n\nexact_logf=\"\":path to logfile for E(EC) vs E(Exact) plot\n\nnote: Note\nAll the effective interactions must be in \"the same order\" and must be consistent with interaction file from which the transition density matrices were made.\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.tf_truncation_occ!-Tuple{Any, Any, Any}","page":"Home","title":"ShellModel.tf_truncation_occ!","text":"truncation scheme by specifying occupation numbers for a given [n,l,j,tz]\n\n\n\n\n\n","category":"method"},{"location":"#ShellModel.transit_main-NTuple{5, Any}","page":"Home","title":"ShellModel.transit_main","text":"transit_main(sntf,target_nuc,jl2,jr2,in_wfs;\n             num_ev_l=100,num_ev_r=100,q=1,is_block=false,is_show=true,\n             calc_EM=true,gfactors=[1.0,0.0,5.586,-3.826],eff_charge=[1.5,0.5])\n\nCalculate the M1&E2 transitions for two wavefunctions\n\nArguments\n\nsntf:path to the interaction file\ntarget_nuc:target nucleus in string e.g., \"Si28\"\njl2:J*2 for the left w.f.\njr2:J*2 for the right w.f.\nin_wfs:[\"path to left wf\",\"path to right wf\"]\n\nOptional arguments\n\nnum_ev_l=100:upper limit of the number of eigenvectors for the left w.f.\nnum_ev_r=100:upper limit of the number of eigenvectors for the right w.f.\nis_show=true:to display the TimerOutput\ngfactors=[1.0,0.0,5.586,-3.826]:g factors [glp,gln,gsp,gsn]\neff_charge=[1.5,0.5]:effective charges [ep,en]\n\nOptional arguments (not used, but for future developments)\n\nq=1:block size\nis_block=false:to use Block algorithm\n\n\n\n\n\n","category":"method"}]
}
