module Interpol

    using SymPy
    x,a,b = Sym("x,a,b")

    export pa,pb,dpa,dpb

    #=First-order Lagrange polynomials for a general finite element Ω (on the interval [-1,1])=#
    pa(x) = subs(((x-b)/(a-b)),(a,-1),(b,1)) #=Negative angular coefficient=#
    pb(x) = subs(((x-a)/(b-a)),(a,-1),(b,1)) #=Positive angular coefficient=#
    dpa = diff(pa(x),x) #=Derivatives=#
    dpb = diff(pb(x),x)

    
    function buildAe() #=Local matrix Aᵢⱼᵉ=∫Φ'ᵢ(x)Φ'ⱼ(x)dx in Ωᵉ=#
        
        a1 = integrate(dpa*dpa,(x,-1,1))
        a2 = integrate(dpa*dpb,(x,-1,1))
        a3 = integrate(dpb*dpa,(x,-1,1))
        a4 = integrate(dpb*dpb,(x,-1,1))
        return [a1 a2; a3 a4]
    end


    function buildCe() #=Local matrix Cᵢⱼᵉ=∫Φᵢ(x)Φⱼ(x)dx in Ωᵉ=#
        c1 = integrate(pa(x)*pa(x),(x,-1,1))
        c2 = integrate(pa(x)*pb(x),(x,-1,1))
        c3 = integrate(pb(x)*pa(x),(x,-1,1))
        c4 = integrate(pb(x)*pb(x),(x,-1,1))
        return [c1 c2; c3 c4]
    end

end