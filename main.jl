#=Variables used in the construction of the local matrices=#
using Plots
using LinearAlgebra

include("./modules/build.jl")
using .Builders

Nn = 20 #=Number of nodes=#
L = 1 #=Length of relaxed string=#
X = LinRange(0,L,Nn)
E = 1:Nn
LM = buildLM(Nn)

Ml = buildGlobalMLeft(X,LM,Nn) #=Constructing matrix on the left-hand side of (11) (see LaTeX)=#
Mr = buildGlobalMRight(X,LM,Nn) #=Constructing matrix on the right-hand side of (11) (see LaTeX)=#


M = inv(Ml)*Mr #=Constructing matrix M (see equation (19) in LaTeX)=#

begin #=Constructing the Crank-Nicolson matrix CK used for time evolution=#
    T = LinRange(0,2,100) #=Time discretization=#
    Δt = 2/(100-1) #=Time step=#
    CK = inv(I-M.*(Δt/2))*(I+M.*(Δt/2)) 
end

begin #=Initial condition=#
    function gaussian(x,μ,σ)#=Gaussian packet with mean μ and variance σ²=#
        exp(-(x-μ)^2/(2σ)^2)
    end
    u0 = gaussian.(X,L/2,L/20).-gaussian(0,L/2,L/20) #=Initial position=#
    p0 = zeros(Nn) #=Initial velocity=#
end

begin #=Time evolution=#
    sol = zeros(length(T),2*Nn)
    sol[1,:] += [u0;p0]
    for t in 1:length(T)-1
        sol[t+1,:] += CK*sol[t,:]
    end
end

i = 5

p = plot(X,sol[i,1:Nn],ylims=(-1,1),label="Guitar String")
xlabel!("x (m)")
ylabel!("y (m)")

savefig(p,"output.png")


function energy(sol,T,C,A,c) #=Calculates total energy and plots relative error=#
    U = sol[:,1:Nn] #=Position=#
    P = sol[:,Nn+1:end] #=Velocity=#
    E = zeros(length(T))
    E0 = 1/2*1/c^2*(transpose(P[1,:])*C*P[1,:])+1/2*(transpose(U[1,:])*A*U[1,:]) #=Initial energy=#
    E[1] += E0
    for t in 2:length(T)
        E[t] += 1/2*1/c^2*(transpose(P[t,:])*C*P[t,:])+1/2*(transpose(U[t,:])*A*U[t,:]) #=Energy=#
    end
    Erel = abs.((E.-E0)./E0) #=Relative error=#
    plot(T,Erel,title="Relative Error for Conserved Energy",label="")
    xlabel!("Time (s)")
    ylabel!("\$\\frac{|E(t)-E(0)|}{|E(0)|}\$",guidefontsize=10)
end
