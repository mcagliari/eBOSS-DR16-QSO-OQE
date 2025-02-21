#As the code was unstable at very low σk, I attached the Taylor expansions to I2, I4, I6 and I8.
#I could improve with:
#    1. Higher order expansions
#    2. Smooth transition
function compute_I_a(a::Int64, σ_fog, ks)
    if !(a in (0,2,4,6,8))
        throw(DomainError(a, "a must be 0, 2, 4, 6, or 8"))
    end

    kσ = ks .* σ_fog
    kσp = (kσ).^2 .+ 2
    arctankσ = √2 .* atan.(√2 ./ 2 .* kσ)

    if a == 0 #I0
        I0 = 2 ./ kσp .+ arctankσ ./ kσ
        return I0
    elseif a == 2 #I2
        I2 = 2 .* (arctankσ ./ kσ.^3 .- 2 ./ kσp ./ kσ.^2)
        sel = kσ .< 1e-5
        I2[sel] .= (2 ./ 3 .- 2 .* σ_fog^2 .* ks[sel].^2 ./ 5 .+ 3 .* σ_fog^4 .* ks[sel].^4 ./ 14 
                    .- σ_fog^6 .* ks[sel].^6 ./ 9 .+ 5 .* σ_fog^8 .* ks[sel].^8 ./ 88)
        return I2
    elseif a == 4 #I4
        I4 = 4 .* (2 .* (1 ./ kσp .+ 1) ./ kσ.^4 - 3 .* arctankσ ./ kσ.^5)
        sel = kσ .< 5.6e-3
        I4[sel] .= (2 ./ 5 .- 2 .* σ_fog^2 .* ks[sel].^2 ./ 7 .+ σ_fog^4 .* ks[sel].^4 ./ 6
                    .- σ_fog^6 .* ks[sel].^6 ./ 11 .+ 5 .* σ_fog^8 .* ks[sel].^8 ./ 104)
        return I4
    elseif a == 6 #I6
        I6 = 8 ./ 3 .* ((kσ.^2 .- 6 ./ kσp .- 12) ./ kσ.^6 .+ 15 .* arctankσ ./ kσ.^7)
        sel = kσ .< 2.4e-2
        I6[sel] .= (2 ./ 7 .- 2 .* σ_fog^2 .* ks[sel].^2 ./ 9 .+ 3 .* σ_fog^4 .* ks[sel].^4 ./ 22
                    .- σ_fog^6 .* ks[sel].^6 ./ 13 .+ σ_fog^8 .* ks[sel].^8 ./ 24)
        return I6
    elseif a == 8 #I8
        I8 = 8 ./ 15 .* (3 ./ kσ.^4 .- 20 ./ kσ.^6 .+ 60 .* (1 ./ kσp .+ 3) ./ kσ.^8 .- 210 .* arctankσ ./ kσ.^9) 
        sel = kσ .< 9e-2
        I8[sel] .= (2 ./ 9 .- 2 .* σ_fog^2 .* ks[sel].^2 ./ 11 .+ 3 .* σ_fog^4 .* ks[sel].^4 ./ 26 
                    .- σ_fog^6 .* ks[sel].^6 ./ 26 .+ 5 .* σ_fog^8 .* ks[sel].^8 ./ 136)
        return I8
    end
end