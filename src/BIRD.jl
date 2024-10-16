module BIRD
using PhysicalConstants.CODATA2018
using Unitful

export @u_str, ν2ω, ω2ν

const e = CODATA2018.ElementaryCharge # Elementary charge in C, use to convert from Coulumb to e
const ε0 = ustrip(u"Å^-1 * J^-1", CODATA2018.VacuumElectricPermittivity/e^2) # ε0 in e^2 / J / Å
const ħ = ustrip(u"J*s", CODATA2018.PlanckConstant/(2π)) # ħ in J*s
const k = ustrip(u"J*K^-1", CODATA2018.BoltzmannConstant) #k_b in J/K
const c = ustrip(u"Å/s", CODATA2018.SpeedOfLightInVacuum) # c in Å/s
const c_cm = ustrip(u"cm/s", CODATA2018.SpeedOfLightInVacuum) # c in Å/s
const eV2cm = 8065.54393734921

"""
    ν2ω(ν)

Convert the wavenumber ν (cm⁻¹) to angular frequency ω (s⁻¹).
"""
function ν2ω(ν)
    return 2π*c_cm* ν
end

"""
    ω2ν(ν)

Convert the angular frequency ω (s⁻¹) to wavenumbers ν (cm⁻¹).
"""
function ω2ν(ω)
    return ω / (2π*c_cm)
end

include("Numerov.jl")
include("Morse.jl")
include("DensityofStates.jl")
include("Rates.jl")

end # module BIRD
