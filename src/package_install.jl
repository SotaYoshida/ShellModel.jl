using Pkg

ts = ["Arpack",
      "Combinatorics",
      "Distributions",
      "LaTeXStrings",
      "LinearAlgebra",
      "Printf",
      "PyCall",
      "Random",
      "SIMD",
      "StatsBase",
      "TimerOutputs",
      "WignerSymbols"]
for tmp in ts
    Pkg.add(tmp)
end
