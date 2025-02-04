using LaTeXStrings

function plot_dipole()

    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"Bond length ($\AA$)", ylabel=L"Dipole ($e\cdot \AA$)")

    # Fit Morse
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    scatter!(ax, rvals, dips)
    rs = rvals[1]:0.01:rvals[end]
    lines!(ax, rs, Dfunc.(rs))

    fig

end

function overtones_osc_str()

    # Fit Morse
    data = [parse(Float64, x) for x in readlines(joinpath(@__DIR__, "Fedorov2014.txt"))]
    rvals = data[1:18]

    debye2au = 1/2.5417464519
    au2eA = 0.529177
    dips = au2eA .* debye2au .* data[37:54]
    Dfunc = BSplineInterpolation(dips, rvals, 3, :ArcLen, :Average, extrapolate=true)

    # Anharmonic constant was obtained from the dissociation energy
    nali = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=257.4u"cm^-1", νχ=2.3273u"cm^-1", 
        re=2.889u"Å", De=0.882388u"eV")

    E = transition_energy_matrix(nali)
    V = transition_dipole_matrix(Dfunc, nali, rvals=0.0:0.01:50, overtones=true)

    DOSp = get_DOS_matrix(E, L=26.5*1e4, g=100, νM=get_transition_wvn(nali, 26, 54)+10^-6, regime=2)
    DOSc = get_DOS_matrix(E, L=26.5*1e4, regime=1)

    J1 = get_transport_matrix(E, V, DOSp, 700)
    J2 = get_transport_matrix(E, V, DOSp, 800)

    diff = abs.(J1 .- J2)
    println(findmax(diff))

    #Vovert = [abs2(V[i,55]) for i in 1:55]

    #BF(e, T) = (1 / (exp(BIRD.ħ*e/(BIRD.k*T)) - 1))

    #fig = Figure()
    #x = Axis(fig[1,1])
    #for T in [400, 600, 800]
    #    lines!(ax, 0:54, [BF(e,T) for e in E[end,:]])
    #end

    fig, ax, hm = heatmap(J1 ./ J2)
    Colorbar(fig[:,end+1], hm)
    fig
end

function madeupmolecule(;precision=Float64)

    # Use HF dipole function
    dipolef(r) = 0.4541416928679838 * r * exp(-0.081616*r^4) 

    # Anharmonic constant was obtained from the dissociation energy
    mol = Morse(mA=22.98u"u", mB=6.941u"u", 
        ν=450.0u"cm^-1", νχ=39.245u"cm^-1", 
        re=2.889u"Å", De=0.16u"eV")

    E = transition_energy_matrix(mol)
    V = transition_dipole_matrix(dipolef, mol, rvals=0.0:0.01:50, overtones=true)

    DOS0 = get_DOS_matrix(E, regime=0)
    J0 = get_transport_matrix(E, V, DOS0, 800)

    println(eigmin(J0))
    
    fig = Figure(fontsize=20)
    ax = Axis(fig[1,1], xlabel=L"$i \rightarrow 5\;$ overtone resonant to $\omega_M$", ylabel=L"k_p/k_0", xticks=0:5, title="Maximum BIRD enhancement")

    for T in [400, 600, 800]
        J0 = precision.(get_transport_matrix(E, V, DOS0, T))
        #k0 = partialschur(J0, nev=1, which=:SR)[1].eigenvalues[1]
        k0 = eigmin(J0)
        println(k0)

        x = zeros(mol.nmax+1)
        for i = 0:mol.nmax
            w = get_transition_wvn(mol, i, mol.nmax)
            println("$i -> 5 -- $w")
            DOSp = get_DOS_matrix(E, L=26.5*1e4, g=100, νM=w + (10^(-6)), regime=2)
            Jp = precision.(get_transport_matrix(E, V, DOSp, T))
            kp = eigmin(Jp)
            println(kp)
            #kp = partialschur(Jp, nev=1, which=:SR)[1].eigenvalues[1]
            x[i+1] = kp/k0
        end
        lines!(ax, 0:mol.nmax, x, label="$T")
        scatter!(ax, 0:mol.nmax, x, label="$T")
    end
    
    axislegend(ax, L"T\;(K)", merge=true)
    fig


end