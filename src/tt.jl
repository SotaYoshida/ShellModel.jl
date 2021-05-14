struct ifph
    i::Int64
    f::Int64
    phase::Bool
end

struct jump1b
    bi::Int64
    bf::Int64
    tifph::Array{ifph,1}
end

ar = jump1b[]
bi=2;bf=10;tjs = [ifph(1,2,false),ifph(3,2,true)]
push!(ar,jump1b(bi,bf,tjs))
bi=5;bf=6;tjs = [ifph(11,22,false),ifph(413,42,true)]
push!(ar,jump1b(bi,bf,tjs))
println("ar $ar")

println("ar[2].tifph ",ar[2].tifph)
push!(ar[2].tifph,ifph(11,0,false))
println("ar[2].tifph ",ar[2].tifph)

