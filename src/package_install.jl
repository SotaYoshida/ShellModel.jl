using Pkg

ts = ["Arpack",
      "Combinatorics",
      "Distributions",
      "FLoops",
      "LaTeXStrings",
      "LinearAlgebra",
      "Printf",
      "PyCall",
      "Random",
      "SIMD",
      "StatsBase",
      "TimerOutputs",
      "ThreadPools",
      "WignerSymbols"]
for tmp in ts
    Pkg.add(tmp)
end
