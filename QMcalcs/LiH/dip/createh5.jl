using HDF5

rvals = 0.5:0.1:6.6
dips = zeros(length(rvals))

# Open the file for reading
open("dip.txt", "r") do file
    # Read and process each line
    i = 1
    for line in eachline(file)
        dips[i] = parse(Float64, line)
        i += 1
    end
end

h5write("dips.h5", "rvals", collect(rvals))
h5write("dips.h5", "dips", dips)
