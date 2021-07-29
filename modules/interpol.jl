module Interpol

    export pa,pb,dpa,dpb

    #=First-order Lagrange polynomials for a general finite element Ω (on the interval [-1,1])=#
    pa(x) = 1/2-x/2 #=Negative angular coefficient=#
    pb(x) = 1/2+x/2 #=Positive angular coefficient=#
    #=Derivatives=#
    dpa(x) = -1/2 
    dpb(x) = 1/2

end