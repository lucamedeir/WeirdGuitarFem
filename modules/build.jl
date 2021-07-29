module Builders

    using FastGaussQuadrature
    include("./interpol.jl")
    using .Interpol

    export buildAe,buildCe,buildGlobalA,buildGlobalC,buildGlobalMLeft,buildGlobalMRight,buildLM

    function buildAe(n) #=Local matrix Aᵢⱼᵉ=∫Φ'ᵢ(x)Φ'ⱼ(x)dx in Ωᵉ=#
        x,w = gausslegendre(n)

        a1 = sum(w .* pa.(x) .* pa.(x) )
        a2 = sum(w .* pa.(x) .* pb.(x) )
        a3 = sum(w .* pb.(x) .* pa.(x) )
        a4 = sum(w .* pb.(x) .* pb.(x) )

        return [a1 a2; a3 a4]
    end

    function buildCe(n) #=Local matrix Cᵢⱼᵉ=∫Φᵢ(x)Φⱼ(x)dx in Ωᵉ=#

        x,w = gausslegendre(n)

        c1 = sum(w .* dpa.(x) .* dpa.(x) )
        c2 = sum(w .* dpa.(x) .* dpb.(x) )
        c3 = sum(w .* dpb.(x) .* dpa.(x) )
        c4 = sum(w .* dpb.(x) .* dpb.(x) )  

        return [c1 c2; c3 c4]
    end

    function buildLM(Nn)
        E = 1:Nn
        LM = [0 0]
        for e in E
            LM = [LM;[e e+1]]
        end
        return LM[2:end-1,:] #=Vector indentifying the Nn-1 finite elements=#
    end

    function buildGlobalA(X,LM,Nn) #=Contructing global matrix Aᵢⱼ=∫Φ'ᵢ(x)Φ'ⱼ(x)dx=#
        he = X[2]-X[1] #=Spatial step Δx=#
        A = zeros((Nn,Nn))
        Ae = buildAe(2)
        for e in 1:(Nn-1)
            globalP = LM[e,:]
            localP = globalP.-(e-1)
            g_i = globalP[1]
            g_j = globalP[2]
            l_i = localP[1]
            l_j = localP[2]
            A[g_i,g_i] += 2*Ae[l_i,l_i]/he 
            A[g_i,g_j] += 2*Ae[l_i,l_j]/he
            A[g_j,g_i] += 2*Ae[l_j,l_i]/he
            A[g_j,g_j] += 2*Ae[l_j,l_j]/he
        end
        return A
    end

    function buildGlobalC(X,LM,Nn) #=Contructing global matrix Cᵢⱼ=∫Φᵢ(x)Φⱼ(x)dx=#
        he = X[2]-X[1] #=Spatial step Δx=#
        C = zeros((Nn,Nn))
        Ce = buildCe(2)
        for e in 1:(Nn-1)
            globalP = LM[e,:]
            localP = globalP.-(e-1)
            g_i = globalP[1]
            g_j = globalP[2]
            l_i = localP[1]
            l_j = localP[2]
            C[g_i,g_i] += he*Ce[l_i,l_i]/2 
            C[g_i,g_j] += he*Ce[l_i,l_j]/2
            C[g_j,g_i] += he*Ce[l_j,l_i]/2
            C[g_j,g_j] += he*Ce[l_j,l_j]/2
        end
        return C
    end

    function buildGlobalMLeft(X,LM,Nn) #=Constructing matrix on the left-hand side of (11) (see LaTeX)=#
        Cl = buildGlobalC(X,LM,Nn) #=Index l for "left"=#
        Cl[1,1] = 1 #=Imposing boundary conditions=#
        Cl[1,2] = 0
        Cl[2,1] = 0
        Cl[end,end] = 1
        Cl[end,end-1] = 0
        Cl[end-1,end] = 0
        return [Cl zeros((Nn,Nn)); zeros((Nn,Nn)) Cl] #Ml
        
    end

    function buildGlobalMRight(X,LM,Nn) #=Constructing matrix on the right-hand side of (11) (see LaTeX)=#
        c = 1 #=speed of sound=#
        Cr = buildGlobalC(X,LM,Nn) #=Index r for "right"=#
        Cr[1,1] = 0
        Cr[1,2] = 0
        Cr[2,1] = 0
        Cr[end,end] = 0
        Cr[end,end-1] = 0
        Cr[end-1,end] = 0
        Ar = buildGlobalA(X,LM,Nn) #=Index r for "right"=#
        Ar[1,1] = 0
        Ar[1,2] = 0
        Ar[2,1] = 0
        Ar[end,end] = 0
        Ar[end,end-1] = 0
        Ar[end-1,end] = 0
        return [zeros((Nn,Nn)) Cr; -(c^2).*Ar zeros((Nn,Nn))] #Mr
    end

end