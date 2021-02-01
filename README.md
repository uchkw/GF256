# GF256
Manipulation module for Galois Field over $2^8$

```julia
julia> a = GF(2)
GF(2)

julia> b = GF(2)
GF(2)

julia> c = GF(3)
GF(3)

julia> a == b
true

julia> a == c
false

julia> bvec(c)
8-element Array{Int64,1}:
 1
 1
 0
 0
 0
 0
 0
 0

julia>  a + c
GF(1)

julia>  c + c
GF(0)

julia>  log(c)
25

julia> a * a
GF(4)

julia> c ^ 5
GF(51)

julia> inv(c)
GF(244)

julia> c * GF(244)
GF(1)
```
