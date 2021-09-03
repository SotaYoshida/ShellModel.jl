using Pkg

ts = ["Combinatorics",
      "Distributions",
      "KrylovKit",
      "LaTeXStrings",
      "LinearAlgebra",
      "Printf",
      "Random",
      "StatsBase",
      "TimerOutputs",
      "WignerSymbols"]
for tmp in ts
    Pkg.add(tmp)
end
