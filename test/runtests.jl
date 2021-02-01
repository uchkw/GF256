using Revise
using GF256
using Test

function check_inv()::Bool
    for i = 1:255
        if GF(1) != (GF(i) * inv(GF(i)))
            return false
        end
    end
    return true
end

a = GF(2)
b = GF(2)
c = GF(3)

@test a == GF(2)
@test (a == b) == true
@test (a == c) == false
@test bvec(c) == [1,1,0,0,0,0,0,0]
@test (a + c) == GF(1)
@test (c + c) == GF(0)
@test log(c) == 25
@test (a * a) == GF(4)
@test (GF(7) * GF(22)) == GF(98)
@test check_inv()
@test inv(GF(0)) == GF(0)
@test (c ^ 5) == GF(51)
@test (a ^ 300) == GF(193)
@test (a ^ -257) == GF(71)