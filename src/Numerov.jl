using LinearAlgebra
using Unitful
using PhysicalConstants

function get_wvn_prefactor(x_unit, mass)

    h = PhysicalConstants.CODATA2018.PlanckConstant
    c = PhysicalConstants.CODATA2018.SpeedOfLightInVacuum

    prefact = uconvert(x_unit^2 * u"cm^-1", h / (8* π^2 * mass * c) )

    return ustrip(prefact)
end

"""
    numerov(Vf, N, xmin, xmax, prefactor)

Implementation of the Matrix Numerov methods as described by Pillai, 
Goglio, and Walker (http://dx.doi.org/10.1119/1.4748813)

"""
function numerov(Vf, N, xmin, xmax, prefactor)

    xvals = range(start=xmin, stop=xmax, length=N)

    d = xvals[2] - xvals[1]

    Im = zeros(N,N)
    Ip = zeros(N,N)

    Im[diagind(Im, -1)] .= 1.0
    Im[diagind(Ip, 1)] .= 1.0

    # Build A, B and V matrices
    A = (Im -2 .* I + Ip) / d^2
    B = (Im + 10 .* I + Ip)/12
    V = diagm(Vf.(xvals))

    # Build Hamiltonian
    H = -prefactor .* inv(B)*A + V

    # Get eigenvalues
    e, C = eigen(H)
end 
