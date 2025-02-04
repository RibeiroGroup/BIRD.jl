function get_pes(name)

    if name == "LiH"
        path1 = joinpath(@__DIR__, "../../QMcalcs/LiH/dip/")
        path2 = joinpath(@__DIR__, "../../QMcalcs/LiH/pes/")
        return BIRD.read_molecule(joinpath(path1, "rvals.txt"), joinpath(path2, "energies.txt"), joinpath(path1, "dip.txt"), 50.0, 0.8740712756579776;N=2000)

    elseif name == "NaH"
        path1 = joinpath(@__DIR__, "../../QMcalcs/NaH/dip")
        path2 = joinpath(@__DIR__, "../../QMcalcs/NaH/pes")
        return BIRD.read_molecule(joinpath(path1, "rvals.txt"), joinpath(path2, "energies.txt"), joinpath(path1, "dips.txt"), 50.0, 0.9656600641501301;N=2000)

    elseif name == "NaLi"
        path1 = joinpath(@__DIR__, "../../QMcalcs/NaLi/dip")
        path2 = joinpath(@__DIR__, "../../QMcalcs/NaLi/pes")
        return BIRD.read_molecule(joinpath(path1, "rvals.txt"), joinpath(path2, "energies.txt"), joinpath(path1, "dips.txt"), 50.0, 5.331369422182236;N=2000)

    elseif name == "MgS"
        path1 = joinpath(@__DIR__, "../../QMcalcs/MgS/dip")
        path2 = joinpath(@__DIR__, "../../QMcalcs/MgS/pes")
        return BIRD.read_molecule(joinpath(path2, "rvals.txt"), joinpath(path2, "energies.txt"), joinpath(path1, "dips.txt"), 50.0, 13.825435958843356;N=2000)

    elseif name == "LiNa"

        rvals = [
        1.7500
        2.0000
        2.2500
        2.5000
        2.8959
        3.5000
        3.7500
        4.0000
        4.2500
        4.5000
        5.0000
        5.5000
        6.0000
        7.0000
        8.0000
        10.0000
        12.0000
        15.0000]

        energies = [
        -169.65053237456
        -169.67999419305
        -169.69857809774
        -169.70886968455
        -169.71372076627
        -169.70793116788
        -169.70393087648
        -169.69987664653
        -169.69609523782
        -169.69277686308
        -169.68778773590
        -169.68478123653
        -169.68313793070
        -169.68183960364
        -169.68148263066
        -169.68132971248
        -169.68130192951
        -169.68129223629
        ]
        energies = 219474.63136320*(energies .- minimum(energies))

        debye2au = 1/2.5417464519
        au2eA = 0.529177
        d = au2eA .* debye2au .* [
        -0.631	
        -0.512	
        -0.469	
        -0.477	
        -0.535	
        -0.626	
        -0.638	
        -0.624	
        -0.583	
        -0.519	
        -0.355	
        -0.208	
        -0.111	
        -0.029	
        -0.009	
        -0.003	
        -0.003	
        -0.003	
        ] 

        return BIRD.PES(rvals, energies, d, 50.0, 5.331369422182236, N=2000)

    else
        error("No $name molecule found.")
    end
end