# Copyright 2021 Hironori Uchikawa
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Galois Field over 2^8
"""
module GF256

export GF, bvec

const prim_poly = 285 # 1+x^2+x^3+x^4+x^8
const gf_size = 256
const gf_size_m1 = gf_size-1

struct GF
    x::Int64
end

"""
transform a GF symbol to corresponding binary vector
>> bvec(GF(3)) == [1,1,0,0,0,0,0,0]
"""
bvec(a::GF) = digits(Int64, a.x, base=2, pad=8)

Base.:(==)(a::GF, b::GF)::Bool = a.x == b.x

Base.:(+)(a::GF) = a
Base.:(+)(a::GF, b::GF)::GF = GF(a.x ⊻ b.x)

Base.:-(a::GF) = a
Base.:-(a::GF, b::GF)::GF = GF(a.x ⊻ b.x)

Base.:⊻(a::GF) = a
Base.:⊻(a::GF, b::GF)::GF = GF(a.x ⊻ b.x)

"""
makepowertable()

Retrun a list of α^p over GF(256).
"""
function makepowertable()::Array{Int64}
    t = ones(Int64, gf_size)
    v = 1
    for i in 2:gf_size
        v <<= 1
        if v > (gf_size_m1)
            v ⊻= prim_poly
        end
        t[i] = v
    end
    return t
end

"""
powertable[p] is α^p over GF(256)
"""
const powertable = Dict{Int64, Int64}(zip(0:(gf_size_m1), makepowertable()))

"""
logarithm table for GF(256).
logtable[p] means log(p)
"""
const logtable = Dict{Int64, Int64}(zip(makepowertable(), 0:(gf_size-2)))

function Base.log(a::GF)::Int64
    if a.x == 0
        return Inf
    end
    try
        return logtable[a.x]
    catch e
        println("GF(256) symbol must be less than 256!")
        bt = backtrace()
        msg = sprint(showerror, e, bt)
        println(msg)
    end
end

function makemultable()::Array{Int64, 2}
    t = zeros(Int64, gf_size_m1, gf_size_m1)
    for i in 1:gf_size_m1
        for j in 1:gf_size_m1
            t[i,j] = powertable[(logtable[i] + logtable[j]) % (gf_size_m1)]
        end
    end
    return t
end

const multable = makemultable()

function mul(a::GF, b::GF)::GF
    if a.x == 0 | b.x == 0
        return GF(0)
    end
    try 
        return GF(multable[a.x,b.x])
    catch e
        println("GF(256) symbol must be less than 256!")
        bt = backtrace()
        msg = sprint(showerror, e, bt)
        println(msg)
    end
end

Base.:*(a::GF, b::GF)::GF = mul(a, b)

function makeinvtable()::Array{Int64}
    t = ones(Int64, gf_size)
    t[1] = 0
    for i in 3:gf_size
        t[i] = powertable[gf_size_m1-logtable[i-1]]
    end
    return t
end

const invtable = Dict{Int64, Int64}(zip(0:(gf_size_m1), makeinvtable()))

function Base.inv(a::GF)::GF
    try
        return GF(invtable[a.x])
    catch e
        println("GF(256) symbol must be less than 256!")
        bt = backtrace()
        msg = sprint(showerror, e, bt)
        println(msg)
    end
end

function makepowtable()::Array{Int64, 2}
    t = ones(Int64, gf_size_m1, gf_size_m1)
    for i in 2:gf_size_m1
        for j in 1:gf_size_m1
            t[i,j] = powertable[(logtable[i] * j) % (gf_size_m1)]
        end
    end
    return t
end

const powtable = makepowtable()

function pow(a::GF, b::Int64)::GF
    if a.x == 0
        return GF(0)
    end
    @assert a.x < gf_size "GF(256) symbol must be less than 256!"
    b %= (gf_size_m1)
    if b == 0
        return GF(1)
    end
    # don't need to think about -b because α^-b is evaluated as inv(α)^b based on literal_pow
    return GF(powtable[a.x,b])
end

Base.:^(a::GF, b::Int64)::GF = pow(a, b)

const α = GF(2)

end # module

