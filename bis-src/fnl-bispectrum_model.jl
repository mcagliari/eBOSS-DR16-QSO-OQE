using Distributions
using Turing
using Tullio


include("I_abc.jl")
#include("../fnl-model.jl")

#F2 s^2 G2
#in @model get_Is 

function F₂(k1, k2, k3)
    μ = μij(k1, k2, k3)
    return 5 / 7 .+ 0.5 .* (k1 ./ k2 + k2 ./ k1) .* μ .+ 2 ./ 7 .* μ.^2
end

function G₂(k1, k2, k3)
    μ = μij(k1, k2, k3)
    return 3 / 7 .+ 0.5 .* (k1 ./ k2 + k2 ./ k1) .* μ .+ 4 ./ 7 .* μ.^2
end

function s₂(k1, k2 , k3)
    return μij(k1, k2, k3).^2 .- 1 / 3
end

function Bℓ(k, Pₘ, Minv, Iabc, ℓ)
    #attualmente scritto considerando che k, Pₘ, e Minv sono delle matrici bidimensinali, ogni colonna corrisponde a uno tra k1, k2 , k3 (inpratica niente mixing matrix per ora)
    I(α, β, γ) = read_Is_matrix(α, β, γ, Iabc, ℓ)

    k1, k2, k3 = k[:,1], k[:,2], k[:,3]

    μ12 = μij(k1, k2, k3)
    μ23 = μij(k2, k3, k1)
    μ31 = μij(k3, k1, k2)

    F₂12 = F₂(k1, k2, k3)
    F₂23 = F₂(k2, k3, k1)
    F₂31 = F₂(k3, k1, k2)

    G₂12 = G₂(k1, k2, k3)
    G₂23 = G₂(k2, k3, k1)
    G₂31 = G₂(k3, k1, k2)

    s₂12 = s₂(k1, k2, k3)
    s₂23 = s₂(k2, k3, k1)
    s₂31 = s₂(k3, k1, k2)
    
    Pₘ1, Pₘ2, Pₘ3 = Pₘ[:,1], Pₘ[:,2], Pₘ[:,3]
    Pₘ12 = Pₘ1 .* Pₘ2
    Pₘ23 = Pₘ2 .* Pₘ3
    Pₘ31 = Pₘ3 .* Pₘ1
    
    #NOTE: in the Pk code α_k = 2 * δ_c * Minv
    M_1, M_2, M_3 = Minv[:,1], Minv[:,2], Minv[:,3]

    #Terms without fNL and without FoG (see Eqs. 75, 76, 79, 80, 81, and 82 in the notes)
    Bb₁3 = 2 .* I(0,0,0) .* (F₂12 .* Pₘ12 .+ F₂23 .* Pₘ23 .+ F₂31 .* Pₘ31)

    Bb₁2b₂ = I(0,0,0) .* (Pₘ12 .+ Pₘ23 .+ Pₘ31)

    Bb₁2bₛ = 2 .* I(0,0,0) .* (s₂12 .* Pₘ12 .+ s₂23 .* Pₘ23 .+ s₂31 .* Pₘ31)

    Bb₁2f = 2 .* (((I(2,0,0) .+ I(0,2,0)) .* F₂12 .+ I(0,0,2) .* G₂12) .* Pₘ12 .+
                  ((I(0,2,0) .+ I(0,0,2)) .* F₂23 .+ I(2,0,0) .* G₂23) .* Pₘ23 .+
                  ((I(0,0,2) .+ I(2,0,0)) .* F₂31 .+ I(0,2,0) .* G₂31) .* Pₘ31)

    Bb₁3f = -1 .* ((k3 ./ k1 .* I(1,0,1) .+ k3 ./ k2 .* I(0,1,1)) .* Pₘ12 .+
                   (k1 ./ k2 .* I(1,1,0) .+ k1 ./ k3 .* I(1,0,1)) .* Pₘ23 .+
                   (k2 ./ k3 .* I(0,1,1) .+ k2 ./ k1 .* I(1,1,0)) .* Pₘ31)

    Bb₁2f2 = -1 .* ((k3 ./ k1 .* (I(3,0,1) .+ 2 .* I(1,2,1)) .+ k3 ./ k2 .* (I(0,3,1) .+ 2 .* I(2,1,1))) .* Pₘ12 .+
                    (k1 ./ k2 .* (I(1,3,0) .+ 2 .* I(1,1,2)) .+ k1 ./ k3 .* (I(1,0,3) .+ 2 .* I(1,2,1))) .* Pₘ23 .+
                    (k2 ./ k3 .* (I(0,1,3) .+ 2 .* I(2,1,1)) .+ k2 ./ k1 .* (I(3,1,0) .+ 2 .* I(1,1,2))) .* Pₘ31)

    Bb₁b₂f = (I(2,0,0) .+ I(0,2,0)) .* Pₘ12 .+
             (I(0,2,0) .+ I(0,0,2)) .* Pₘ23 .+
             (I(0,0,2) .+ I(2,0,0)) .* Pₘ31

    Bb₁bₛf = 2 .* ((I(2,0,0) .+ I(0,2,0)) .* s₂12 .* Pₘ12 .+ 
                   (I(0,2,0) .+ I(0,0,2)) .* s₂23 .* Pₘ23 .+
                   (I(0,0,2) .+ I(2,0,0)) .* s₂31 .* Pₘ31)

    Bb₁f2 = 2 .* ((I(2,2,0) .* F₂12 .+ (I(2,0,2) .+ I(0,2,2)) .* G₂12) .* Pₘ12 .+ 
                 (I(0,2,2) .* F₂23 .+ (I(2,2,0) .+ I(2,0,2)) .* G₂23) .* Pₘ23 .+
                 (I(2,0,2) .* F₂31 .+ (I(0,2,2) .+ I(2,2,0)) .* G₂31) .* Pₘ31)

    Bb₁f3 = -1 .* ((k3 ./ k1 .* (2 .* I(3,2,1) .+ I(1,4,1)) .+ k3 ./ k2 .* (2 .* I(2,3,1) .+ I(4,1,1))) .* Pₘ12 .+
                   (k1 ./ k2 .* (2 .* I(1,3,2) .+ I(1,1,4)) .+ k1 ./ k3 .* (2 .* I(1,2,3) .+ I(1,4,1))) .* Pₘ23 .+
                   (k2 ./ k3 .* (2 .* I(2,1,3) .+ I(4,1,1)) .+ k2 ./ k1 .* (2 .* I(3,1,2) .+ I(1,1,4))) .* Pₘ31)

    Bb₂f2 = I(2,2,0) .* Pₘ12 .+
            I(0,2,2) .* Pₘ23 .+
            I(2,0,2) .* Pₘ31

    Bbₛf2 = 2 .* (I(2,2,0) .* s₂12 .* Pₘ12 .+
                  I(0,2,2) .* s₂23 .* Pₘ23 .+
                  I(2,0,2) .* s₂31 .* Pₘ31)

    Bf3 = 2 .* I(2,2,2) .* (G₂12 .* Pₘ12 .+ G₂23 .* Pₘ23 .+ G₂31 .* Pₘ31)

    Bf4 = -1 .* ((k3 ./ k1 .* I(3,4,1) .+ k3 ./ k2 .* I(4,3,1)) .* Pₘ12 .+
                 (k1 ./ k2 .* I(1,3,4) .+ k1 ./ k3 .* I(1,4,3)) .* Pₘ23 .+
                 (k2 ./ k3 .* I(4,1,3) .+ k2 ./ k1 .* I(3,1,4)) .* Pₘ31)
    
    #Terms with fNL and without FoG (see Eqs. 74, 75, 76, 77, 78, 79, 80, 81, and 82 in the notes)
    #Terms proportional to bϕ (see Eqs. 75, 76, 77, 79, 80, 81, and 82 in the notes)
    Bb₁2bϕfNL = 2 .* I(0,0,0) .* (((M_1 .+ M_2) .* F₂12 .+ 0.5 .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12) .* Pₘ12 .+
                                  ((M_2 .+ M_3) .* F₂23 .+ 0.5 .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23) .* Pₘ23 .+
                                  ((M_3 .+ M_1) .* F₂31 .+ 0.5 .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31) .* Pₘ31)

    Bb₁bϕffNL = 2 .* ((I(2,0,0) .* (M_2 .* F₂12 .+ 0.5 .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12) .+
                       I(0,2,0) .* (M_1 .* F₂12 .+ 0.5 .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12) .+
                       I(0,0,2) .* (M_1 .+ M_2) .* G₂12) .* Pₘ12 .+
                      (I(0,2,0) .* (M_3 .* F₂23 .+ 0.5 .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23) .+
                       I(0,0,2) .* (M_2 .* F₂23 .+ 0.5 .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23) .+
                       I(2,0,0) .* (M_2 .+ M_3) .* G₂23) .* Pₘ23 .+
                      (I(0,0,2) .* (M_1 .* F₂31 .+ 0.5 .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31) .+
                       I(2,0,0) .* (M_3 .* F₂31 .+ 0.5 .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31) .+
                       I(0,2,0) .* (M_3 .+ M_1) .* G₂31) .* Pₘ31)

    Bb₁2bϕffNL = -1 .* ((I(1,0,1) .* k3 ./ k1 .* (M_1 .+ 2 .* M_2) .+ I(0,1,1) .* k3 ./ k2 .* (2 .* M_1 .+ M_2)) .* Pₘ12 .+
                        (I(1,1,0) .* k1 ./ k2 .* (M_2 .+ 2 .* M_3) .+ I(1,0,1) .* k1 ./ k3 .* (2 .* M_2 .+ M_3)) .* Pₘ23 .+
                        (I(0,1,1) .* k2 ./ k3 .* (M_3 .+ 2 .* M_1) .+ I(1,1,0) .* k2 ./ k1 .* (2 .* M_3 .+ M_1)) .* Pₘ31)

    Bb₁bϕf2fNL = -2 .* ((k3 ./ k1 .* (I(3,0,1) .* M_2 .+ I(1,2,1) .* (M_1 .+ M_2)) .+
                         k3 ./ k2 .* (I(0,3,1) .* M_1 .+ I(2,1,1) .* (M_1 .+ M_2))) .* Pₘ12 .+
                        (k1 ./ k2 .* (I(1,3,0) .* M_3 .+ I(1,1,2) .* (M_2 .+ M_3)) .+
                         k1 ./ k3 .* (I(1,0,3) .* M_2 .+ I(1,2,1) .* (M_2 .+ M_3))) .* Pₘ23 .+
                        (k2 ./ k3 .* (I(0,1,3) .* M_1 .+ I(2,1,1) .* (M_3 .+ M_1)) .+
                         k2 ./ k1 .* (I(3,1,0) .* M_3 .+ I(1,1,2) .* (M_3 .+ M_1))) .* Pₘ31)

    Bb₁b₂bϕfNL = I(0,0,0) .* ((M_1 .+ M_2) .* Pₘ12 .+
                              (M_2 .+ M_3) .* Pₘ23 .+
                              (M_3 .+ M_1) .* Pₘ31)

    Bb₁bₛbϕfNL = 2 .* I(0,0,0) .* ((M_1 .+ M_2) .* s₂12 .* Pₘ12 .+
                                  (M_2 .+ M_3) .* s₂23 .* Pₘ23 .+
                                  (M_3 .+ M_1) .* s₂31 .* Pₘ31)

    Bb₂bϕffNL = (I(2,0,0) .* M_2 .+ I(0,2,0) .* M_1) .* Pₘ12 .+
                (I(0,2,0) .* M_3 .+ I(0,0,2) .* M_2) .* Pₘ23 .+
                (I(0,0,2) .* M_1 .+ I(2,0,0) .* M_3) .* Pₘ31

    BbₛbϕffNL = 2 .* ((I(2,0,0) .* M_2 .+ I(0,2,0) .* M_1) .* s₂12 .* Pₘ12 .+
                      (I(0,2,0) .* M_3 .+ I(0,0,2) .* M_2) .* s₂23 .* Pₘ23 .+
                      (I(0,0,2) .* M_1 .+ I(2,0,0) .* M_3) .* s₂31 .* Pₘ31)

    Bbϕf2fNL = 2 .* ((I(2,2,0) .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12 .* 0.5 .+ (I(2,0,2) .* M_2 .+ I(0,2,2) .* M_1) .* G₂12) .* Pₘ12 .+
                     (I(0,2,2) .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23 .* 0.5 .+ (I(2,2,0) .* M_3 .+ I(2,0,2) .* M_2) .* G₂23) .* Pₘ23 .+
                     (I(2,0,2) .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31 .* 0.5 .+ (I(0,2,2) .* M_1 .+ I(2,2,0) .* M_3) .* G₂31) .* Pₘ31)

    Bbϕf3fNL = -1 .* ((k3 ./ k1 .* (2 .* I(3,2,1) .* M_2 .+ I(1,4,1) .* M_1) .+ k3 ./ k2 .* (2 .* I(2,3,1) .* M_1 .+ I(4,1,1) .* M_2)) .* Pₘ12 .+
                      (k1 ./ k2 .* (2 .* I(1,3,2) .* M_3 .+ I(1,1,4) .* M_2) .+ k1 ./ k3 .* (2 .* I(1,2,3) .* M_2 .+ I(1,4,1) .* M_3)) .* Pₘ23 .+
                      (k2 ./ k3 .* (2 .* I(2,1,3) .* M_1 .+ I(4,1,1) .* M_3) .+ k2 ./ k1 .* (2 .* I(3,1,2) .* M_3 .+ I(1,1,4) .* M_1)) .* Pₘ31)
    
    #Terms proportional to bϕδ (see Eq. 78 in the notes)
    Bb₁2bϕδfNL = I(0,0,0) .* ((M_1 .+ M_2) .* Pₘ12 .+ 
                              (M_2 .+ M_3) .* Pₘ23 .+  
                              (M_3 .+ M_1) .* Pₘ31) #I missed a broadcasting here!

    Bb₁bϕδffNL = ((I(2,0,0) .+ I(0,2,0)) .* (M_1 .+ M_2) .* Pₘ12 .+
                  (I(0,2,0) .+ I(0,0,2)) .* (M_2 .+ M_3) .* Pₘ23 .+
                  (I(0,0,2) .+ I(2,0,0)) .* (M_3 .+ M_1) .* Pₘ31)

    Bbϕδf2fNL = (I(2,2,0) .* (M_1 .+ M_2) .* Pₘ12 .+
                 I(0,2,2) .* (M_2 .+ M_3) .* Pₘ23 .+
                 I(2,0,2) .* (M_3 .+ M_1) .* Pₘ31)

    #Terms coming from Z_1(k1)*Z_1(k2)*Z_1(k3) (see Eq. 74 in the notes)
    B123 = Pₘ12 .* M_1 .* M_2 ./ M_3 .+ Pₘ23 ./ M_1 .* M_2 .* M_3 .+ Pₘ31 .* M_1 ./ M_2 .* M_3

    Bb₁3fNL = 2 .* I(0,0,0) .* B123

    Bb₁2ffNL = 2 .* (I(2,0,0) .+ I(0,2,0) .+ I(0,0,2)) .* B123

    Bb₁f2fNL = 2 .* (I(2,2,0) .+ I(2,0,2) .+ I(0,2,2)) .* B123

    Bf3fNL = 2 .* I(2,2,2) .* B123

    #Stochastic terms without FoG and assuming α₂=0 and α₃=-1 (see Eq. 83 in the notes)
    #NOTE: the terms needs to be multiplied by 1/n_bar
    Bb₁2α₁ = I(0,0,0) .* (Pₘ1 .+ Pₘ2 .+ Pₘ3)

    Bb₁bϕα₁fNL = 2 .* I(0,0,0) .* (M_1 .* Pₘ1 .+ M_2 .* Pₘ2 .+ M_3 .* Pₘ3)

    Bb₁fα₁ = I(2,0,0) .* Pₘ1 .+ I(0,2,0) .* Pₘ2 .+ I(0,0,2) .* Pₘ3

    Bbϕfα₁fNL = I(2,0,0) .* M_1 .* Pₘ1 .+ I(0,2,0) .* M_2 .* Pₘ2 .+ I(0,0,2) .* M_3 .* Pₘ3

    #Terms without fNL with FoG (see Eqs. 53, 54, 55, 58, and 59 in the notes) these terms come from 2 Z₁Z₁Z₂
    #Terms proportional to c₁ (see Eqs. 53, 54, 55, and 58 in the notes)

    kNL = 0.3
    κ1 = k1 ./ kNL
    κ2 = k2 ./ kNL
    κ3 = k3 ./ kNL

    Bb₁b₂c₁ = -1 .* ((κ1.^2 .* I(2,0,0) .+  κ2.^2 .* I(0,2,0)) .* Pₘ12 .+
                     (κ2.^2 .* I(0,2,0) .+  κ3.^2 .* I(0,0,2)) .* Pₘ23 .+
                     (κ3.^2 .* I(0,0,2) .+  κ1.^2 .* I(2,0,0)) .* Pₘ31)

    Bb₁2c₁ = -2 .* ((κ1.^2 .* I(2,0,0) .+  κ2.^2 .* I(0,2,0)) .* F₂12 .* Pₘ12 .+
                    (κ2.^2 .* I(0,2,0) .+  κ3.^2 .* I(0,0,2)) .* F₂23 .* Pₘ23 .+
                    (κ3.^2 .* I(0,0,2) .+  κ1.^2 .* I(2,0,0)) .* F₂31 .* Pₘ31)

    Bb₁bₛc₁ = -2 .* ((κ1.^2 .* I(2,0,0) .+  κ2.^2 .* I(0,2,0)) .* s₂12 .* Pₘ12 .+
                     (κ2.^2 .* I(0,2,0) .+  κ3.^2 .* I(0,0,2)) .* s₂23 .* Pₘ23 .+
                     (κ3.^2 .* I(0,0,2) .+  κ1.^2 .* I(2,0,0)) .* s₂31 .* Pₘ31)

    Bb₁fc₁ = -2 .* (((κ1.^2 .* I(2,0,2) .+ κ2.^2 .* I(0,2,2)) .* G₂12 .+ I(2,2,0) .* (κ1.^2 .+ κ2.^2) .* F₂12) .* Pₘ12 .+ 
                    ((κ2.^2 .* I(2,2,0) .+ κ3.^2 .* I(2,0,2)) .* G₂23 .+ I(0,2,2) .* (κ2.^2 .+ κ3.^2) .* F₂23) .* Pₘ23 .+
                    ((κ3.^2 .* I(0,2,2) .+ κ1.^2 .* I(2,2,0)) .* G₂31 .+ I(2,0,2) .* (κ3.^2 .+ κ1.^2) .* F₂31) .* Pₘ31)

    Bb₁2fc₁ = (k3 ./ k1 .* (κ1.^2 .* I(3,0,1) .+ 2 .* κ2.^2 .* I(1,2,1)) .+ k3 ./ k2 .* (2 .* κ1.^2 .* I(2,1,1) .+ κ2.^2 .* I(0,3,1))) .* Pₘ12 .+
              (k1 ./ k2 .* (κ2.^2 .* I(1,3,0) .+ 2 .* κ3.^2 .* I(1,1,2)) .+ k1 ./ k3 .* (2 .* κ2.^2 .* I(1,2,1) .+ κ3.^2 .* I(1,0,3))) .* Pₘ23 .+
              (k2 ./ k3 .* (κ3.^2 .* I(0,1,3) .+ 2 .* κ1.^2 .* I(2,1,1)) .+ k2 ./ k1 .* (2 .* κ3.^2 .* I(1,1,2) .+ κ1.^2 .* I(3,1,0))) .* Pₘ31 # (53), (54). (55)

    Bb₁f2c₁ = 2 .* ((k3 ./ k1 .* ((κ1.^2 .+ κ2.^2) .* I(3,2,1) .+ κ2.^2 .* I(1,4,1)) .+ k3 ./ k2 .* (κ1.^2 .* I(4,1,1) .+ (κ1.^2 .+ κ2.^2) .* I(2,3,1))) .* Pₘ12 .+
                    (k1 ./ k2 .* ((κ2.^2 .+ κ3.^2) .* I(1,3,2) .+ κ3.^2 .* I(1,1,4)) .+ k1 ./ k3 .* (κ2.^2 .* I(1,4,1) .+ (κ2.^2 .+ κ3.^2) .* I(1,2,3))) .* Pₘ23 .+
                    (k2 ./ k3 .* ((κ3.^2 .+ κ1.^2) .* I(2,1,3) .+ κ1.^2 .* I(4,1,1)) .+ k2 ./ k1 .* (κ3.^2 .* I(1,1,4) .+ (κ3.^2 .+ κ1.^2) .* I(3,1,2))) .* Pₘ31) #(53), (54), (55), (58)

    Bb₂fc₁ = -1 .* ((κ1.^2 .+ κ2.^2) .* I(2,2,0) .* Pₘ12 .+
                    (κ2.^2 .+ κ3.^2) .* I(0,2,2) .* Pₘ23 .+
                    (κ3.^2 .+ κ1.^2) .* I(2,0,2) .* Pₘ31)

    Bbₛfc₁ = -2 .* ((κ1.^2 .+ κ2.^2) .* I(2,2,0) .* s₂12 .* Pₘ12 .+
                    (κ2.^2 .+ κ3.^2) .* I(0,2,2) .* s₂23 .* Pₘ23 .+
                    (κ3.^2 .+ κ1.^2) .* I(2,0,2) .* s₂31 .* Pₘ31)

    Bf2c₁ = -2 .* I(2,2,2) .* ((κ1.^2 .+ κ2.^2) .* G₂12 .* Pₘ12 .+
                               (κ2.^2 .+ κ3.^2) .* G₂23 .* Pₘ23 .+
                               (κ3.^2 .+ κ1.^2) .* G₂31 .* Pₘ31)

    Bf3c₁ = (k3 ./ k1 .* (κ1.^2 .+ 2 .* κ2.^2) .* I(3,4,1) .+ k3 ./ k2 .* (2 .* κ1.^2 .+ κ2.^2) .* I(4,3,1)) .* Pₘ12 .+
            (k1 ./ k2 .* (κ2.^2 .+ 2 .* κ3.^2) .* I(1,3,4) .+ k1 ./ k3 .* (2 .* κ2.^2 .+ κ3.^2) .* I(1,4,3)) .* Pₘ23 .+
            (k2 ./ k3 .* (κ3.^2 .+ 2 .* κ1.^2) .* I(4,1,3) .+ k2 ./ k1 .* (2 .* κ3.^2 .+ κ1.^2) .* I(3,1,4)) .* Pₘ31 #(53), (54), (58)

    #Terms proportional to c₁ⁿ (see Eqs. 55, 58, and 59)

    Bb₂c₁2 = κ1.^2 .* κ2.^2 .* I(2,2,0) .* Pₘ12 .+
             κ2.^2 .* κ3.^2 .* I(0,2,2) .* Pₘ23 .+
             κ3.^2 .* κ1.^2 .* I(2,0,2) .* Pₘ31

    Bb₁c₁2 = 2 .* (κ1.^2 .* κ2.^2 .* I(2,2,0) .* F₂12 .* Pₘ12 .+
                   κ2.^2 .* κ3.^2 .* I(0,2,2) .* F₂23 .* Pₘ23 .+
                   κ3.^2 .* κ1.^2 .* I(2,0,2) .* F₂31 .* Pₘ31)

    Bbₛc₁2 = 2 .* (κ1.^2 .* κ2.^2 .* I(2,2,0) .* s₂12 .* Pₘ12 .+
                   κ2.^2 .* κ3.^2 .* I(0,2,2) .* s₂23 .* Pₘ23 .+
                   κ3.^2 .* κ1.^2 .* I(2,0,2) .* s₂31 .* Pₘ31)

    Bfc₁2 = 2 .* I(2,2,2) .* (κ1.^2 .* κ2.^2 .* G₂12 .* Pₘ12 .+
                              κ2.^2 .* κ3.^2 .* G₂23 .* Pₘ23 .+
                              κ3.^2 .* κ1.^2 .* G₂31 .* Pₘ31)

    Bb₁fc₁2 = -1 .* ((k3 ./ k1 .* (2 .* κ1.^2 .* κ2.^2 .* I(3,2,1) .+ κ2.^4 .* I(1,4,1)) .+ k3 ./ k2 .* (κ1.^4 .* I(4,1,1) .+ 2 .* κ1.^2 .* κ2.^2 .* I(2,3,1))) .* Pₘ12 .+
                     (k1 ./ k2 .* (2 .* κ2.^2 .* κ3.^2 .* I(1,3,2) .+ κ3.^4 .* I(1,1,4)) .+ k1 ./ k3 .* (κ2.^4 .* I(1,4,1) .+ 2 .* κ2.^2 .* κ3.^2 .* I(1,2,3))) .* Pₘ23 .+
                     (k2 ./ k3 .* (2 .* κ3.^2 .* κ1.^2 .* I(2,1,3) .+ κ1.^4 .* I(4,1,1)) .+ k2 ./ k1 .* (κ3.^4 .* I(1,1,4) .+ 2 .* κ3.^2 .* κ1.^2 .* I(3,1,2))) .* Pₘ31) #(55), (59)

    Bf2c₁2 = -1 .* ((k3 ./ k1 .* (2 .* κ1.^2 .* κ2.^2 .+ κ2.^4) .* I(3,4,1) .+ k3 ./ k2 .* (κ1.^4 .+ 2 .* κ1.^2 .* κ2.^2) .* I(4,3,1)) .* Pₘ12 .+
                    (k1 ./ k2 .* (2 .* κ2.^2 .* κ3.^2 .+ κ3.^4) .* I(1,3,4) .+ k1 ./ k3 .* (κ2.^4 .+ 2 .* κ2.^2 .* κ3.^2) .* I(1,4,3)) .* Pₘ23 .+
                    (k2 ./ k3 .* (2 .* κ3.^2 .* κ1.^2 .+ κ1.^4) .* I(4,1,3) .+ k2 ./ k1 .* (κ3.^4 .+ 2 .* κ3.^2 .* κ1.^2) .* I(3,1,4)) .* Pₘ31) #(58), (59)

    Bfc₁3 = (k3 ./ k1 .* κ1.^2 .* κ2.^4 .* I(3,4,1) .+ k3 ./ k2 .* κ1.^4 .* κ2.^2 .* I(4,3,1)) .* Pₘ12 .+
            (k1 ./ k2 .* κ2.^2 .* κ3.^4 .* I(1,3,4) .+ k1 ./ k3 .* κ2.^4 .* κ3.^2 .* I(1,4,3)) .* Pₘ23 .+
            (k2 ./ k3 .* κ3.^2 .* κ1.^4 .* I(4,1,3) .+ k2 ./ k1 .* κ3.^4 .* κ1.^2 .* I(3,1,4)) .* Pₘ31        

    #Terms with fNL and FoG (see Eqs. 53, 54, 55, 57, 58, 59, 60, 62, and 63 in the notes)
    #Terms proportional to c₁ and from 2 Z₁Z₁Z₂ (see Eqs. 53, 54, 55, 57, and 58 in the notes)
    
    Bb₁bϕδfNLc₁ = -1 .* ((κ1.^2 .* I(2,0,0) .+ κ2.^2 .* I(0,2,0)) .* (M_1 + M_2) .* Pₘ12 .+
                         (κ2.^2 .* I(0,2,0) .+ κ3.^2 .* I(0,0,2)) .* (M_2 + M_3) .* Pₘ23 .+
                         (κ3.^2 .* I(0,0,2) .+ κ1.^2 .* I(2,0,0)) .* (M_3 + M_1) .* Pₘ31)

    Bb₁bϕfNLc₁ = -2 .* (((κ1.^2 .* I(2,0,0) .+ κ2.^2 .* I(0,2,0)) .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12 .* 0.5 .+ 
                          (κ1.^2 .* I(2,0,0) .* M_2 .+ κ2.^2 .* I(0,2,0) .* M_1) .* F₂12) .* Pₘ12 .+
                        ((κ2.^2 .* I(0,2,0) .+ κ3.^2 .* I(0,0,2)) .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23 .* 0.5 .+ 
                          (κ2.^2 .* I(0,2,0) .* M_3 .+ κ3.^2 .* I(0,0,2) .* M_2) .* F₂23) .* Pₘ23 .+
                        ((κ3.^2 .* I(0,0,2) .+ κ1.^2 .* I(2,0,0)) .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31 .* 0.5 .+ 
                          (κ3.^2 .* I(0,0,2) .* M_1 .+ κ1.^2 .* I(2,0,0) .* M_3) .* F₂31) .* Pₘ31)

    Bb₁bϕffNLc₁ = 2 .* ((k3 ./ k1 .* (κ2.^2 .* I(1,2,1) .* (M_1 .+ M_2) .+ κ1.^2 .* I(3,0,1) .* M_2) .+ 
                         k3 ./ k2 .* (κ1.^2 .* I(2,1,1) .* (M_1 .+ M_2) .+ κ2.^2 .* I(0,3,1) .* M_1)) .* Pₘ12 .+
                        (k1 ./ k2 .* (κ3.^2 .* I(1,1,2) .* (M_2 .+ M_3) .+ κ2.^2 .* I(1,3,0) .* M_3) .+ 
                         k1 ./ k3 .* (κ2.^2 .* I(1,2,1) .* (M_2 .+ M_3) .+ κ3.^2 .* I(1,0,3) .* M_2)) .* Pₘ23 .+
                        (k2 ./ k3 .* (κ1.^2 .* I(2,1,1) .* (M_3 .+ M_1) .+ κ3.^2 .* I(0,1,3) .* M_1) .+ 
                         k2 ./ k1 .* (κ3.^2 .* I(1,1,2) .* (M_3 .+ M_1) .+ κ1.^2 .* I(3,1,0) .* M_3)) .* Pₘ31)

    BbϕδffNLc₁ = -1 .* ((κ1.^2 .+ κ2.^2) .* I(2,2,0) .* (M_1 .+ M_2) .* Pₘ12 .+
                        (κ2.^2 .+ κ3.^2) .* I(0,2,2) .* (M_2 .+ M_3) .* Pₘ23 .+
                        (κ3.^2 .+ κ1.^2) .* I(2,0,2) .* (M_3 .+ M_1) .* Pₘ31)

    BbϕffNLc₁ = -2 .* (((κ1.^2 .* I(2,0,2) .* M_2 .+ κ2.^2 .* I(0,2,2) .* M_1) .* G₂12 .+ 
                        (κ1.^2 .+ κ2.^2) .* I(2,2,0) .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12 .* 0.5) .* Pₘ12 .+
                       ((κ2.^2 .* I(2,2,0) .* M_3 .+ κ3.^2 .* I(2,0,2) .* M_2) .* G₂23 .+ 
                        (κ2.^2 .+ κ3.^2) .* I(0,2,2) .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23 .* 0.5) .* Pₘ23 .+
                       ((κ3.^2 .* I(0,2,2) .* M_1 .+ κ1.^2 .* I(2,2,0) .* M_3) .* G₂31 .+ 
                        (κ3.^2 .+ κ1.^2) .* I(2,0,2) .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31 .* 0.5) .* Pₘ31) #(57), (58)

    Bbϕf2fNLc₁ = 2 .* ((k3 ./ k1 .* ((κ1.^2 .+ κ2.^2) .* I(3,2,1) .* M_2 .+ κ2.^2 .* I(1,4,1) .* M_1) .+ 
                        k3 ./ k2 .* ((κ1.^2 .+ κ2.^2) .* I(2,3,1) .* M_1 .+ κ1.^2 .* I(4,1,1) .* M_2)) .* Pₘ12 .+
                       (k1 ./ k2 .* ((κ2.^2 .+ κ3.^2) .* I(1,3,2) .* M_3 .+ κ3.^2 .* I(1,1,4) .* M_2) .+ 
                        k1 ./ k3 .* ((κ2.^2 .+ κ3.^2) .* I(1,2,3) .* M_2 .+ κ2.^2 .* I(1,4,1) .* M_3)) .* Pₘ23 .+
                       (k2 ./ k3 .* ((κ3.^2 .+ κ1.^2) .* I(2,1,3) .* M_1 .+ κ1.^2 .* I(4,1,1) .* M_3) .+ 
                        k2 ./ k1 .* ((κ3.^2 .+ κ1.^2) .* I(3,1,2) .* M_3 .+ κ3.^2 .* I(1,1,4) .* M_1)) .* Pₘ31) #(53), (54), (57), (58)

    Bb₂bϕfNLc₁ = -1 .* ((κ1.^2 .* I(2,0,0) .* M_2 .+ κ2.^2 .* I(0,2,0) .* M_1) .* Pₘ12 .+
                        (κ2.^2 .* I(0,2,0) .* M_3 .+ κ3.^2 .* I(0,0,2) .* M_2) .* Pₘ23 .+
                        (κ3.^2 .* I(0,0,2) .* M_1 .+ κ1.^2 .* I(2,0,0) .* M_3) .* Pₘ31)

    BbₛbϕfNLc₁ = -2 .* ((κ1.^2 .* I(2,0,0) .* M_2 .+ κ2.^2 .* I(0,2,0) .* M_1) .* s₂12 .* Pₘ12 .+
                        (κ2.^2 .* I(0,2,0) .* M_3 .+ κ3.^2 .* I(0,0,2) .* M_2) .* s₂23 .* Pₘ23 .+
                        (κ3.^2 .* I(0,0,2) .* M_1 .+ κ1.^2 .* I(2,0,0) .* M_3) .* s₂31 .* Pₘ31)

    #Terms proportional to c₁ⁿ and from 2 Z₁Z₁Z₂ (see Eqs. 57 and 59 in the notes)

    BbϕδfNLc₁2 = (κ1.^2 .* κ2.^2 .* I(2,2,0) .* (M_1 .+ M_2) .* Pₘ12 .+
                  κ2.^2 .* κ3.^2 .* I(0,2,2) .* (M_2 .+ M_3) .* Pₘ23 .+
                  κ3.^2 .* κ1.^2 .* I(2,0,2) .* (M_3 .+ M_1) .* Pₘ31)

    BbϕfNLc₁2 = (κ1.^2 .* κ2.^2 .* I(2,2,0) .* (k1 ./ k2 .* M_1 .+ k2 ./ k1 .* M_2) .* μ12 .* Pₘ12 .+
                 κ2.^2 .* κ3.^2 .* I(0,2,2) .* (k2 ./ k3 .* M_2 .+ k3 ./ k2 .* M_3) .* μ23 .* Pₘ23 .+
                 κ3.^2 .* κ1.^2 .* I(2,0,2) .* (k3 ./ k1 .* M_3 .+ k1 ./ k3 .* M_1) .* μ31 .* Pₘ31)

    BbϕffNLc₁2 = -1 .* ((k3 ./ k1 .* (2 .* κ1.^2 .* κ2.^2 .* I(3,2,1) .* M_2 .+ κ2.^4 .* I(1,4,1) .* M_1) .+
                         k3 ./ k2 .* (2 .* κ1.^2 .* κ2.^2 .* I(2,3,1) .* M_1 .+ κ1.^4 .* I(4,1,1) .* M_2)) .* Pₘ12 .+
                        (k1 ./ k2 .* (2 .* κ2.^2 .* κ3.^2 .* I(1,3,2) .* M_3 .+ κ3.^4 .* I(1,1,4) .* M_2) .+
                         k1 ./ k3 .* (2 .* κ2.^2 .* κ3.^2 .* I(1,2,3) .* M_2 .+ κ2.^4 .* I(1,4,1) .* M_3)) .* Pₘ23 .+
                        (k2 ./ k3 .* (2 .* κ3.^2 .* κ1.^2 .* I(2,1,3) .* M_1 .+ κ1.^4 .* I(4,1,1) .* M_3) .+
                         k2 ./ k1 .* (2 .* κ3.^2 .* κ1.^2 .* I(3,1,2) .* M_3 .+ κ3.^4 .* I(1,1,4) .* M_1)) .* Pₘ31) #(57), (59)
    
    #Terms from Z₁Z₁Z₁ (see Eqs. 60, 62, and 63 in the notes)

    Bb₁2fNLc₁ = -2 .* (κ1.^2 .* I(2,0,0) .+ κ2.^2 .* I(0,2,0) .+ κ3.^2 .* I(0,0,2)) .* B123

    Bb₁ffNLc₁ = -2 .* ((κ1.^2 .+ κ2.^2) .* I(2,2,0) .+ 
                       (κ2.^2 .+ κ3.^2) .* I(0,2,2) .+
                       (κ3.^2 .+ κ1.^2) .* I(2,0,2)) .* B123

    Bf2fNLc₁ = -2 .* (κ1.^2 .+ κ2.^2 .+ κ3.^2) .* I(2,2,2) .* B123

    Bb₁fNLc₁2 = 2 .* (κ1.^2 .* κ2.^2 .* I(2,2,0) .+
                      κ2.^2 .* κ3.^2 .* I(0,2,2) .+
                      κ3.^2 .* κ1.^2 .* I(2,0,2)) .* B123

    BffNLc₁2 = 2 .* (κ1.^2 .* κ2.^2 .+ κ2.^2 .* κ3.^2 .+ κ3.^2 .* κ1.^2) .* I(2,2,2) .* B123

    BfNLc₁3 = -2 .* κ1.^2 .* κ2.^2 .* κ3.^2 .* I(2,2,2) .* B123

    #Stochastic terms with FoG and assuming α₂=0 and α₃=-1 (see Eq. 17 in the notes)
    #NOTE: the terms needs to be multiplied by 1/n_bar

    Bb₁α₁c₁ = -1 .* (κ1.^2 .* I(2,0,0) .* Pₘ1 .+ κ2.^2 .* I(0,2,0) .* Pₘ2 .+ κ3.^2 .* I(0,0,2) .* Pₘ3)

    Bbϕα₁fNLc₁ = -1 .* (κ1.^2 .* I(2,0,0) .* M_1 .* Pₘ1 .+ κ2.^2 .* I(0,2,0) .* M_2 .* Pₘ2 .+ κ3.^2 .* I(0,0,2) .* M_3 .* Pₘ3)

    return [Bb₁3, Bb₁2b₂, Bb₁2bₛ, Bb₁2f, Bb₁3f, Bb₁2f2, Bb₁b₂f, Bb₁bₛf, Bb₁f2, Bb₁f3, Bb₂f2, Bbₛf2, Bf3, Bf4, #Terms without fNL and without FoG
            Bb₁2bϕfNL, Bb₁bϕffNL, Bb₁2bϕffNL, Bb₁bϕf2fNL, Bb₁b₂bϕfNL, Bb₁bₛbϕfNL, Bb₂bϕffNL, BbₛbϕffNL, Bbϕf2fNL, Bbϕf3fNL, #Terms  (proportional to bϕ)
            Bb₁2bϕδfNL, Bb₁bϕδffNL, Bbϕδf2fNL, #Terms with fNL and without FoG (proportional to bϕδ)
            Bb₁3fNL, Bb₁2ffNL, Bb₁f2fNL, Bf3fNL, #Terms coming from Z_1(k1)*Z_1(k2)*Z_1(k3) without FoG
            Bb₁2α₁, Bb₁bϕα₁fNL, Bb₁fα₁, Bbϕfα₁fNL, #Stochastic terms without FoG
            Bb₁b₂c₁, Bb₁2c₁, Bb₁bₛc₁, Bb₁fc₁, Bb₁2fc₁, Bb₁f2c₁, Bb₂fc₁, Bbₛfc₁, Bf2c₁, Bf3c₁, #Terms without fNL with FoG from 2 Z₁Z₁Z₂ (proportional to c₁)
            Bb₂c₁2, Bb₁c₁2, Bbₛc₁2, Bfc₁2, Bb₁fc₁2, Bf2c₁2, Bfc₁3, #Terms without fNL with FoG from 2 Z₁Z₁Z₂ (proportional to c₁ⁿ)
            Bb₁bϕδfNLc₁, Bb₁bϕfNLc₁, Bb₁bϕffNLc₁, BbϕδffNLc₁, BbϕffNLc₁, Bbϕf2fNLc₁, Bb₂bϕfNLc₁, BbₛbϕfNLc₁, #Terms with fNL and FoG from 2 Z₁Z₁Z₂ (proportional to c₁)
            BbϕδfNLc₁2, BbϕfNLc₁2, BbϕffNLc₁2, #Terms with fNL and FoG from 2 Z₁Z₁Z₂ (proportional to c₁ⁿ)
            Bb₁2fNLc₁, Bb₁ffNLc₁, Bf2fNLc₁, Bb₁fNLc₁2, BffNLc₁2, BfNLc₁3, #Terms from Z₁Z₁Z₁
            Bb₁α₁c₁, Bbϕα₁fNLc₁] #Stochastic terms with FoG

    
end

function compute_B0(k, Pₘ, Minv, Iabc, bias, f, fNL, Psn)
    b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂ = bias

    bias_vector = [b₁*b₁*b₁, b₁*b₁*b₂, b₁*b₁*bₛ, b₁*b₁*f, b₁*b₁*b₁*f, b₁*b₁*f*f, b₁*b₂*f, b₁*bₛ*f, b₁*f*f, b₁*f*f*f, b₂*f*f, bₛ*f*f, f*f*f, f*f*f*f, #Terms without fNL and without FoG
    b₁*b₁*bϕ*fNL, b₁*bϕ*f*fNL, b₁*b₁*bϕ*f*fNL, b₁*bϕ*f*f*fNL, b₁*b₂*bϕ*fNL, b₁*bₛ*bϕ*fNL, b₂*bϕ*f*fNL, bₛ*bϕ*f*fNL, bϕ*f*f*fNL, bϕ*f*f*f*fNL, #Terms  (proportional to bϕ)
    b₁*b₁*bϕδ*fNL, b₁*bϕδ*f*fNL, bϕδ*f*f*fNL, #Terms with fNL and without FoG (proportional to bϕδ)
    b₁*b₁*b₁*fNL, b₁*b₁*f*fNL, b₁*f*f*fNL, f*f*f*fNL, #Terms coming from Z_1(k1)*Z_1(k2)*Z_1(k3) without FoG
    b₁*b₁*(1+α₁)*Psn, b₁*bϕ*(1+α₁)*fNL*Psn, b₁*f*(1+α₁)*Psn, bϕ*f*(1+α₁)*fNL*Psn, #Stochastic terms without FoG ####NOTA SERVE TERMINE APLHA2!
    b₁*b₂*c₁, b₁*b₁*c₁, b₁*bₛ*c₁, b₁*f*c₁, b₁*b₁*f*c₁, b₁*f*f*c₁, b₂*f*c₁, bₛ*f*c₁, f*f*c₁, f*f*f*c₁, #Terms without fNL with FoG from 2 Z₁Z₁Z₂ (proportional to c₁)
    b₂*c₁*c₁, b₁*c₁*c₁, bₛ*c₁*c₁, f*c₁*c₁, b₁*f*c₁*c₁, f*f*c₁*c₁, f*c₁*c₁*c₁, #Terms without fNL with FoG from 2 Z₁Z₁Z₂ (proportional to c₁ⁿ)
    b₁*bϕδ*fNL*c₁, b₁*bϕ*fNL*c₁, b₁*bϕ*f*fNL*c₁, bϕδ*f*fNL*c₁, bϕ*f*fNL*c₁, bϕ*f*f*fNL*c₁, b₂*bϕ*fNL*c₁, bₛ*bϕ*fNL*c₁, #Terms with fNL and FoG from 2 Z₁Z₁Z₂ (proportional to c₁)
    bϕδ*fNL*c₁*c₁, bϕ*fNL*c₁*c₁, bϕ*f*fNL*c₁*c₁, #Terms with fNL and FoG from 2 Z₁Z₁Z₂ (proportional to c₁ⁿ)
    b₁*b₁*fNL*c₁, b₁*f*fNL*c₁, f*f*fNL*c₁, b₁*fNL*c₁*c₁, f*fNL*c₁*c₁, fNL*c₁*c₁*c₁, #Terms from Z₁Z₁Z₁
    b₁*(1+α₁)*c₁*Psn, bϕ*(1+α₁)*fNL*c₁*Psn, (1+α₂)*Psn*Psn] #Stochastic terms with FoG + Pns²

    terms = Bℓ(k, Pₘ, Minv, Iabc, 0)
    nk = length(k[:,1])
    Bα₂ = ones(nk)
    terms = vcat(terms, [Bα₂])
    B₀terms = reduce(vcat, terms')

    @tullio B0[j] := bias_vector[i] * B₀terms[i,j]

    return B0

end

function compute_B0(B0_terms, bias, f, fNL, Psn)
    b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂ = bias

    bias_vector = [b₁*b₁*b₁, b₁*b₁*b₂, b₁*b₁*bₛ, b₁*b₁*f, b₁*b₁*b₁*f, b₁*b₁*f*f, b₁*b₂*f, b₁*bₛ*f, b₁*f*f, b₁*f*f*f, b₂*f*f, bₛ*f*f, f*f*f, f*f*f*f, #Terms without fNL and without FoG
    b₁*b₁*bϕ*fNL, b₁*bϕ*f*fNL, b₁*b₁*bϕ*f*fNL, b₁*bϕ*f*f*fNL, b₁*b₂*bϕ*fNL, b₁*bₛ*bϕ*fNL, b₂*bϕ*f*fNL, bₛ*bϕ*f*fNL, bϕ*f*f*fNL, bϕ*f*f*f*fNL, #Terms  (proportional to bϕ)
    b₁*b₁*bϕδ*fNL, b₁*bϕδ*f*fNL, bϕδ*f*f*fNL, #Terms with fNL and without FoG (proportional to bϕδ)
    b₁*b₁*b₁*fNL, b₁*b₁*f*fNL, b₁*f*f*fNL, f*f*f*fNL, #Terms coming from Z_1(k1)*Z_1(k2)*Z_1(k3) without FoG
    b₁*b₁*(1+α₁)*Psn, b₁*bϕ*(1+α₁)*fNL*Psn, b₁*f*(1+α₁)*Psn, bϕ*f*(1+α₁)*fNL*Psn, #Stochastic terms without FoG ####NOTA SERVE TERMINE APLHA2!
    b₁*b₂*c₁, b₁*b₁*c₁, b₁*bₛ*c₁, b₁*f*c₁, b₁*b₁*f*c₁, b₁*f*f*c₁, b₂*f*c₁, bₛ*f*c₁, f*f*c₁, f*f*f*c₁, #Terms without fNL with FoG from 2 Z₁Z₁Z₂ (proportional to c₁)
    b₂*c₁*c₁, b₁*c₁*c₁, bₛ*c₁*c₁, f*c₁*c₁, b₁*f*c₁*c₁, f*f*c₁*c₁, f*c₁*c₁*c₁, #Terms without fNL with FoG from 2 Z₁Z₁Z₂ (proportional to c₁ⁿ)
    b₁*bϕδ*fNL*c₁, b₁*bϕ*fNL*c₁, b₁*bϕ*f*fNL*c₁, bϕδ*f*fNL*c₁, bϕ*f*fNL*c₁, bϕ*f*f*fNL*c₁, b₂*bϕ*fNL*c₁, bₛ*bϕ*fNL*c₁, #Terms with fNL and FoG from 2 Z₁Z₁Z₂ (proportional to c₁)
    bϕδ*fNL*c₁*c₁, bϕ*fNL*c₁*c₁, bϕ*f*fNL*c₁*c₁, #Terms with fNL and FoG from 2 Z₁Z₁Z₂ (proportional to c₁ⁿ)
    b₁*b₁*fNL*c₁, b₁*f*fNL*c₁, f*f*fNL*c₁, b₁*fNL*c₁*c₁, f*fNL*c₁*c₁, fNL*c₁*c₁*c₁, #Terms from Z₁Z₁Z₁
    b₁*(1+α₁)*c₁*Psn, bϕ*(1+α₁)*fNL*c₁*Psn, (1+α₂)*Psn*Psn] #Stochastic terms with FoG + Pns²

    #@tullio B0[j] := bias_vector[i] * B0_terms[i,j]
    B0 = vecmat_tullio(bias_vector, B0_terms)

    return B0
end

#Turing model

@model function B_qso(data, k, P0_model, Minv, Iabc, f, Psn, iΓ, p)
    
    δ_c = 1.686

    #per ora ho scritto dei prior solo per far funzionare il codice
    fNL ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    b₂ ~ Uniform(-4, 4)
    bₛ ~ Uniform(-4, 4)
    bϕ = 2 * δ_c * (b₁ - p) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁ - 1)) - b₁ + 1) #usiamo Universal relation
    c₁ ~ Uniform(-1, 1)
    α₁ ~ Uniform(-1.5, 1.5)
    α₂ ~ Uniform(-1.5, 1.5) #forse messo a 0

    #Bispectrum monopole
    bias = [b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂]
    prediction = compute_B0(k, P0_model, Minv, Iabc, bias, f, fNL, Psn)
    prediction_recats = iΓ * prediction

    data ~ MvNormal(prediction_recats, I)

    return nothing

end

@model function B_qso(data, k, P0_model, Minv, Iabc, f, Psn, iΓ, p, WBric)
    
    δ_c = 1.686

    #per ora ho scritto dei prior solo per far funzionare il codice
    fNL ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    b₂ ~ Uniform(-4, 4)
    bₛ ~ Uniform(-4, 4)
    bϕ = 2 * δ_c * (b₁ - p) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁ - 1)) - b₁ + 1) #usiamo Universal relation
    c₁ ~ Uniform(-1, 1)
    α₁ ~ Uniform(-1, 1)
    α₂ ~ Uniform(-1, 1) #forse messo a 0

    #Bispectrum monopole

    bias = [b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂]
    B0 = compute_B0(k, P0_model, Minv, Iabc, bias, f, fNL, Psn)

    #RIC contribution (GIC contribution in P0_model)
    prediction = B0 .- B0 .* WBric
    prediction_recats = iΓ * prediction

    data ~ MvNormal(prediction_recats, I)

    return nothing

end

@model function B_qso(data, B0_terms, f, Psn, iΓ, p)
    
    δ_c = 1.686

    #per ora ho scritto dei prior solo per far funzionare il codice
    fNL ~ Uniform(-500, 500)
    b₁ ~ Normal(2.3, 0.5)#Uniform(0.1, 6)
    b₂ ~ Normal(-3.7,5)#Uniform(-4, 4)
    bₛ ~ Normal(-0.7,0.6)#Uniform(-4, 4)
    bϕ = 2 * δ_c * (b₁ - p) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁ - 1)) - b₁ + 1) #usiamo Universal relation 2 * (δ_c * (b₂ - 25 / 3 * bₛ + 13 / 21 * (b₁ - p)) - b₁ + p)
    c₁ = 0 #~ Uniform(-1, 1)
    α₁ ~ Uniform(-1.5, 1.5)
    α₂ ~ Uniform(-1.5, 1.5) #forse messo a 0

    #Bispectrum monopole
    bias = [b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂]
    prediction = compute_B0(B0_terms, bias, f, fNL, Psn)
    prediction_recats = iΓ * prediction

    data ~ MvNormal(prediction_recats, I)

    return nothing

end

@model function B_qso(data, B0_terms, BGIC_terms, f, Psn, iΓ, p, WBric)
    
    δ_c = 1.686

    #per ora ho scritto dei prior solo per far funzionare il codice
    fNL ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6) #Normal(2.28, 2.)#
    b₂ ~ Uniform(-4, 4)
    bₛ ~ Uniform(-4, 4)
    bϕ = 2 * δ_c * (b₁ - p) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁ - 1)) - b₁ + 1) #usiamo Universal relation
    c₁ ~ Uniform(-1, 1)
    α₁ ~ Uniform(-1.5, 1.5)
    α₂ ~ Uniform(-1.5, 1.5) #forse messo a 0

    #Bispectrum monopole

    bias = [b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂]
    B0 = compute_B0(B0_terms, bias, f, fNL, Psn)
    BGIC = compute_B0(BGIC_terms, bias, f, fNL, Psn)

    #RIC contribution (GIC contribution in P0_model)
    prediction = BGIC .- B0 .* WBric
    prediction_recats = iΓ * prediction

    data ~ MvNormal(prediction_recats, I)

    return nothing

end

##fnl=0
#@model function B_qso(data, B0_terms, f, Psn, iΓ, p; fNL=0)
#    
#    δ_c = 1.686
#
#    #per ora ho scritto dei prior solo per far funzionare il codice
#    fNL = 0
#    b₁ ~ Uniform(0.1, 6)
#    b₂ = 0.412-2.143*b₁+0.929*b₁^2+0.008*b₁^3 #~ Uniform(-1, 1)
#    bₛ ~ Uniform(-4, 4)
#    bϕ = 0 #2 * δ_c * (b₁ - p) #usiamo Universal relation
#    bϕδ = 0 #2 * (δ_c * (b₂ + 13 / 21 * (b₁ - p)) - b₁ + p) #usiamo Universal relation
#    c₁ ~ Uniform(-1, 1)
#    α₁ ~ Uniform(-1.5, 1.5)
#    α₂ ~ Uniform(-1.5, 1.5) #forse messo a 0
#
#    #Bispectrum monopole
#    bias = [b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂]
#    prediction = compute_B0(B0_terms, bias, f, fNL, Psn)
#    prediction_recats = iΓ * prediction
#
#    data ~ MvNormal(prediction_recats, I)
#
#    return nothing
#
#end
#
#@model function B_qso(data, B0_terms, BGIC_terms, f, Psn, iΓ, p, WBric; fNL=false)
#    
#    δ_c = 1.686
#
#    #per ora ho scritto dei prior solo per far funzionare il codice
#    fNL = 0 #~ Uniform(-500, 500)
#    b₁ ~ Uniform(0.1, 6)
#    b₂ ~ Uniform(-4, 4)
#    bₛ ~ Uniform(-4, 4)
#    bϕ = 0 #2 * δ_c * (b₁ - p) #usiamo Universal relation
#    bϕδ = 0 #2 * (δ_c * (b₂ + 13 / 21 * (b₁ - p)) - b₁ + p) #usiamo Universal relation
#    c₁ ~ Uniform(-1, 1)
#    α₁ ~ Uniform(-1.5, 1.5)
#    α₂ ~ Uniform(-1.5, 1.5) #forse messo a 0
#
#    #Bispectrum monopole
#
#    bias = [b₁, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂]
#    B0 = compute_B0(B0_terms, bias, f, fNL, Psn)
#    BGIC = compute_B0(BGIC_terms, bias, f, fNL, Psn)
#
#    #RIC contribution (GIC contribution in P0_model)
#    prediction = BGIC .- B0 .* WBric
#    prediction_recats = iΓ * prediction
#
#    data ~ MvNormal(prediction_recats, I)
#
#    return nothing
#
#end

#join anaylses
@model function B_qso(dataN, dataS, p, B0_termsN, B0_termsS, fN, fS, PsnN, PsnS, iΓN, iΓS)
    
    δ_c = 1.686

    #prior:
    fNL ~ Uniform(-500, 500)

    b₁N ~ Uniform(0.1, 6)
    b₂N ~ Uniform(-4, 4)
    bₛN ~ Uniform(-4, 4)
    c₁N ~ Uniform(-1, 1)
    α₁N ~ Uniform(-1.5, 1.5)
    α₂N ~ Uniform(-1.5, 1.5) #forse messo a 0

    b₁S ~ Uniform(0.1, 6)
    b₂S ~ Uniform(-4, 4)
    bₛS ~ Uniform(-4, 4)
    c₁S ~ Uniform(-1, 1)
    α₁S ~ Uniform(-1.5, 1.5)
    α₂S ~ Uniform(-1.5, 1.5) #forse messo a 0

    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    bϕN = 2 * δ_c * (b₁N - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδN = bϕN + 2 * (δ_c * (b₂N - 8 / 21 * (b₁N - 1)) - b₁N + 1) #a caso #usiamo Universal relation
    bϕS = 2 * δ_c * (b₁S - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδS = bϕS + 2 * (δ_c * (b₂S - 8 / 21 * (b₁S - 1)) - b₁S + 1) #a caso #usiamo Universal relation

    #parameter vectors
    biasN = [b₁N, b₂N, bₛN, bϕN, bϕδN, c₁N, α₁N, α₂N] #for B
    biasS = [b₁S, b₂S, bₛS, bϕS, bϕδS, c₁S, α₁S, α₂S] #for B

    #B likelihood
    predictionBN = compute_B0(B0_termsN, biasN, fN, fNL, PsnN)
    predictionBS = compute_B0(B0_termsS, biasS, fS, fNL, PsnS)

    predictionBN_recast = iΓN * predictionBN
    predictionBS_recast = iΓS * predictionBS

    dataN ~ MvNormal(predictionBN_recast, I)
    dataS ~ MvNormal(predictionBS_recast, I)

    return nothing

end

@model function B_qso(dataN, dataS, p, B0_termsN, B0_termsS, BGIC_termsN, BGIC_termsS, fN, fS, PsnN, PsnS, iΓN, iΓS, WBricN, WBricS)
    
    δ_c = 1.686

    #prior:
    fNL ~ Uniform(-500, 500)
    
    b₁N ~ Uniform(0.1, 6)
    b₂N ~ Uniform(-4, 4)
    bₛN ~ Uniform(-4, 4)
    c₁N ~ Uniform(-1, 1)
    α₁N ~ Uniform(-1.5, 1.5)
    α₂N ~ Uniform(-1.5, 1.5) #forse messo a 0

    b₁S ~ Uniform(0.1, 6)
    b₂S ~ Uniform(-4, 4)
    bₛS ~ Uniform(-4, 4)
    c₁S ~ Uniform(-1, 1)
    α₁S ~ Uniform(-1.5, 1.5)
    α₂S ~ Uniform(-1.5, 1.5) #forse messo a 0

    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    bϕN = 2 * δ_c * (b₁N - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδN = bϕN + 2 * (δ_c * (b₂N - 8 / 21 * (b₁N - 1)) - b₁N + 1) #a caso #usiamo Universal relation
    bϕS = 2 * δ_c * (b₁S - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδS = bϕS + 2 * (δ_c * (b₂S - 8 / 21 * (b₁S - 1)) - b₁S + 1) #a caso #usiamo Universal relation

    #parameter vectors
    biasN = [b₁N, b₂N, bₛN, bϕN, bϕδN, c₁N, α₁N, α₂N] #for B
    biasS = [b₁S, b₂S, bₛS, bϕS, bϕδS, c₁S, α₁S, α₂S] #for B

    #B likelihood
    B0N = compute_B0(B0_termsN, biasN, fN, fNL, PsnN)
    B0S = compute_B0(B0_termsS, biasS, fS, fNL, PsnS)
    BGICN = compute_B0(BGIC_termsN, biasN, fN, fNL, PsnN)
    BGICS = compute_B0(BGIC_termsS, biasS, fS, fNL, PsnS)

    #RIC contribution (GIC contribution in P0_model)
    predictionBN = BGICN .- B0N .* WBricN
    predictionBS = BGICS .- B0S .* WBricS

    predictionBN_recast = iΓN * predictionBN
    predictionBS_recast = iΓS * predictionBS

    dataN ~ MvNormal(predictionBN_recast, I)
    dataS ~ MvNormal(predictionBS_recast, I)

    return nothing

end