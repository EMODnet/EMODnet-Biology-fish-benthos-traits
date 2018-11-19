using DIVAnd
using PyPlot
using NCDatasets
using Missings
using Interpolations
using Plots

if VERSION >= v"0.7"
    using Random
    using DelimitedFiles
    using Statistics
    using Printf
    using FileIO
else
    using Compat: @info, @warn, range, cat
end

include("../src/emodnet_bio_grid.jl");

"""
```julia-repl
obslon, obslat, obsyear, g1, g2, g3, g4 = read_fish(filename)
```
Read the coordinates, the year and the abundance of each of the 4 types of fish
"""
function read_fish(filename::String)
    data,header = readdlm(filename,',',header = true)
    header = header[:]
    if "year" in header
        @info "Working on a temporal data file"
        stationname = Vector{String}(data[:,findfirst(header .== "samp")]);
        obsyear = Vector{Int32}(data[:,findfirst(header .== "year")]);
    else
        @info "Working on a spatial data file"
        obsyear = undef;
    end;

    obslon = Vector{Float64}(data[:,findfirst(header .== "x")]);
    obslat = Vector{Float64}(data[:,findfirst(header .== "y")]);

    g1 = Vector{Float64}(data[:,findfirst(header .== "g1")]);
    g2 = Vector{Float64}(data[:,findfirst(header .== "g2")]);
    g3 = Vector{Float64}(data[:,findfirst(header .== "g3")]);
    g4 = Vector{Float64}(data[:,findfirst(header .== "g4")]);

    @info "Number of data points: $(length(g1))"

    return obslon, obslat, obsyear, g1, g2, g3, g4
end

"""
```julia-repl
obslon, obslat, obsyear, g1, g2, g3, g4 = read_fish_specific(filename)
```
Read the coordinates, the year and the abundance of each of the 4 types of fish
for the specific cases
"""
function read_fish_specific(filename::String)
    """
    Header:
    "samp","x","y","Entelurus aequoreus","Molva molva","Sprattus sprattus","Squalus acanthias"
    """
    data,header = readdlm(filename,',',header = true)
    header = header[:]

    obslon = Vector{Float64}(data[:,findfirst(header .== "x")]);
    obslat = Vector{Float64}(data[:,findfirst(header .== "y")]);

    stationname = Vector{String}(data[:,findfirst(header .== "samp")]);

    g1 = Vector{Float64}(data[:,findfirst(header .== "Entelurus aequoreus")]);
    g2 = Vector{Float64}(data[:,findfirst(header .== "Molva molva")]);
    g3 = Vector{Float64}(data[:,findfirst(header .== "Sprattus sprattus")]);
    g4 = Vector{Float64}(data[:,findfirst(header .== "Squalus acanthias")]);

    return obslon, obslat, g1, g2, g3, g4
end

"""
```julia-repl
f1, f2, f3 = make_fish_analysis(obslon, obslat, g1, g2, g3, g4)
```
Perform
1. data transformation,
2. DIVAnd interpolation and
3. inverse transformation
"""
function make_fish_analysis(obslon, obslat, g1, g2, g3, g4)

    # Transformed fields
    g1log = log.(g1 .+ 1);
    g2log = log.(g2 .+ 1);
    g3log = log.(g3 .+ 1);
    g4log = log.(g4 .+ 1);

    # Perform analysis using the selected reference field
    @time fi1,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g1log .- mean(g1log),len,epsilon2,alphabc=2);

    @time fi2,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g2log .- mean(g2log),len,epsilon2,alphabc=2);

    @time fi3,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g3log .- mean(g3log),len,epsilon2,alphabc=2);

    @time fi4,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g4log .- mean(g4log),len,epsilon2,alphabc=2);

    # Tranform back and relative fields
    f1ori = exp.(fi1 .+ mean(g1log)) .- 1;
    f2ori = exp.(fi2 .+ mean(g2log)) .- 1;
    f3ori = exp.(fi3 .+ mean(g3log)) .- 1;
    f4ori = exp.(fi4 .+ mean(g4log)) .- 1;
    f1ori[f1ori.<0] .= 0.
    f2ori[f2ori.<0] .= 0.
    f3ori[f3ori.<0] .= 0.
    f4ori[f4ori.<0] .= 0.
    totalfield = f1ori + f2ori + f3ori + f4ori;

    f1rel, f2rel, f3rel, f4rel = f1ori./totalfield, f2ori./totalfield, f3ori./totalfield, f4ori./totalfield;

    return f1rel, f2rel, f3rel, f4rel, totalfield;
end;

"""
```julia-repl
err1, err2, err3 = compute_fish_error(obslon, obslat, g1, g2, g3, g4)
```
Compute the error fields using the Clever Poor Man's Estimate.
Data are first transformed using a log function.
"""
function compute_fish_error(obslon, obslat, g1, g2, g3, g4)
    g1log = log.(g1 .+ 1);
    g2log = log.(g2 .+ 1);
    g3log = log.(g3 .+ 1);
    g4log = log.(g3 .+ 1);
    g1_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g1log .- mean(g1log),len,epsilon2,alphabc=2);
    g2_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g2log .- mean(g2log),len,epsilon2,alphabc=2);
    g3_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g3log .- mean(g3log),len,epsilon2,alphabc=2);
    g4_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g4log .- mean(g4log),len,epsilon2,alphabc=2);

    return g1_err, g2_err, g3_err, g4_err
end

"""
```julia-repl
add_mask(bx, by, b)
```
Add the land-sea mask to the plot with a grey (0.5) color
"""
function add_mask(bx, by, b)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0], colors = [[.5,.5,.5]])
end

"""
```julia-repl
plot_fish_data(obslonS, obslatS, g1, g2, g3, g4)
```
Make a scatter of the data values, one subplot per fish category
"""
function plot_fish_data(obslonS, obslatS, g1, g2, g3, g4)

    function makesubplot(g, ax)
        ax[:tick_params]("both",labelsize=6)
        scat1 = PyPlot.scatter(obslonS, obslatS, s=.01, c=g)
        add_mask(bx, by, b)
        gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
        colorbar(scat1)[:ax][:tick_params](labelsize=8)
        add_mask(bx, by, b)
    end

    figure("fish_data")
    ax1 = subplot(2,2,1)
    makesubplot(g1, ax1)
    ax2 = subplot(2,2,2)
    makesubplot(g2, ax2)
    ax3 = subplot(2,2,3)
    makesubplot(g3, ax3)
    ax4 = subplot(2,2,4)
    makesubplot(g4, ax4)
end

"""
```julia-repl
make_scatter_grid_specific_fish(fields, titlelist)
```
Make the plot of the data stored in `fields` and use the `titlelist`
for the title of the subplots.
"""
function make_scatter_grid_specific_fish(fields, titlelist)

    nfield = length(fields)
    if nfield == 3
        nx = 1
        ny = 3
        figure("fish_analysis_spec", figsize=(12,8))
    elseif nfield == 4
        nx = 2
        ny = 2
        figure("fish_analysis_spec", figsize=(8,8))
    end

    for i = 1:nfield
        ax = subplot(nx,ny,i)
        title(titlelist[i], fontsize=8)
        ax[:tick_params]("both",labelsize=6)
        PyPlot.plot(obslon[g1 .== 1], obslat[g1 .== 1], "ro", markersize=1, label="Presence")
        PyPlot.plot(obslon[g1 .== 0], obslat[g1 .== 0], "go", markersize=.1, label="Absence")
        gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
        add_mask(bx, by, b)
        if i == 1
            legend()
        end
    end
end

"""
```julia-repl
plot_fish_results(gridlon, gridlat, f1, f2, f3, f4, bx, by, b)
```
Make the plot of the 4 interpolated fields f1, f2, f3, f4
and overlay the mask based on the bathymetry defined by bx, by and b
"""
function plot_fish_results(gridlon, gridlat, f1, f2, f3, f4, bx, by, b)

    function makesubplot(fi, ax)
        ax[:tick_params]("both",labelsize=6)
        pcm1 = PyPlot.pcolormesh(gridlon, gridlat, permutedims(fi, [2,1]), vmin=0, vmax=1.)
        add_mask(bx, by, b)
        gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
        colorbar(pcm1)[:ax][:tick_params](labelsize=8)
    end

    figure()
    ax1 = subplot(2,2,1)
    makesubplot(f1, ax1)
    ax2 = subplot(2,2,2)
    makesubplot(f2, ax2)
    ax3 = subplot(2,2,3)
    makesubplot(f3, ax3)
    ax4 = subplot(2,2,4)
    makesubplot(f4, ax4)

end

"""
```julia-repl
make_scatter_grid_specific_fish(fields, titlelist)
```
Make the plot of the 4 interpolated fields stored in `fields` and use `titlelist`
for the titles of the subplots.
"""
function make_scatter_grid_specific_fish(fields, titlelist; vmin=0., vmax=1.0)
    nfield = length(fields)
    if nfield == 3
        nx = 1
        ny = 3
        figure("fish_analysis_spec", figsize=(12,8))
    elseif nfield == 4
        nx = 2
        ny = 2
        figure("fish_analysis_spec", figsize=(8,8))
    end

    for i = 1:nfield
        ax = subplot(nx,ny,i)
        title(titlelist[i], fontsize=8)
        ax[:tick_params]("both",labelsize=6)
        pcm = PyPlot.pcolormesh(gridlonFish, gridlatFish, permutedims(fields[i], [2,1]), vmin=vmin, vmax=vmax)
        gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
        add_mask(bx, by, b)
        colorbar(pcm,shrink=.7)[:ax][:tick_params](labelsize=8)
    end
end

"""
```julia-repl
write_fish_nc(filename, gridlon, gridlat, field1, field2, field3, field4,
             err1, err2, err3, err4)
```
Write the interpolated and error fields for the 4 fish categories in the
netCDF file `filename`
"""
function write_fish_nc(filename::String, gridlon, gridlat,
        field1::Array, field2::Array, field3::Array, field4::Array,
        err1::Array, err2::Array, err3::Array, err4::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);

        # Define a global attribute
        ds.attrib["title"] = "Interpolated Fish"

        # Define the variables and coordinates
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))

        # Attributes
        lat.attrib["long_name"] = "Latitude";
        lat.attrib["standard_name"] = "latitude";
        lat.attrib["units"] = "degrees_north";

        lon.attrib["long_name"] = "Longitude";
        lon.attrib["standard_name"] = "longitude";
        lon.attrib["units"] = "degrees_east";

        # Interpolated fields
        g1 = defVar(ds,"g1",Float64,("lon","lat"))
        g2 = defVar(ds,"g2",Float64,("lon","lat"))
        g3 = defVar(ds,"g3",Float64,("lon","lat"))
        g4 = defVar(ds,"g4",Float64,("lon","lat"))

        # Error fields
        g1_err = defVar(ds,"g1_err",Float64,("lon","lat"))
        g2_err = defVar(ds,"g2_err",Float64,("lon","lat"))
        g3_err = defVar(ds,"g3_err",Float64,("lon","lat"))
        g4_err = defVar(ds,"g4_err",Float64,("lon","lat"))

        # Fill the coord vectors and the fields
        lon[:] = gridlon;
        lat[:] = gridlat;

        g1[:,:] = field1;
        g2[:,:] = field2;
        g3[:,:] = field3;
        g4[:,:] = field4;

        g1_err[:,:] = err1;
        g2_err[:,:] = err2;
        g3_err[:,:] = err3;
        g4_err[:,:] = err4;

    end

end;


"""
```julia-repl
write_fish_time_nc(filename, gridlon, gridlat, gridtime,
                    field1, field2, field3, field4,
                    err1, err2, err3, err4)
```
Write a netCDF file containing the gridded and error fields for the
temporal fish product
"""
function write_fish_time_nc(filename::String, gridlon, gridlat, gridtime,
        field1::Array, field2::Array, field3::Array, field4::Array,
        err1::Array, err2::Array, err3::Array, err4::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);
        ntimes = length(gridtime);

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);
        defDim(ds,"time",ntimes);

        # Define a global attribute
        ds.attrib["title"] = "Interpolated Fish - temporal product"

        # Define the variables and coordinates
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))
        time = defVar(ds,"time",Float32,("time",))

        # Attributes
        lat.attrib["long_name"] = "Latitude";
        lat.attrib["standard_name"] = "latitude";
        lat.attrib["units"] = "degrees_north";

        lon.attrib["long_name"] = "Longitude";
        lon.attrib["standard_name"] = "longitude";
        lon.attrib["units"] = "degrees_east";

        time.attrib["long_name"] = "Time";
        time.attrib["standard_name"] = "time";
        time.attrib["units"] = "years";

        # Interpolated fields
        g1 = defVar(ds,"g1",Float64,("lon","lat","time"))
        g2 = defVar(ds,"g2",Float64,("lon","lat","time"))
        g3 = defVar(ds,"g3",Float64,("lon","lat","time"))
        g4 = defVar(ds,"g4",Float64,("lon","lat","time"))

        # Error fields
        g1_err = defVar(ds,"g1_err",Float64,("lon","lat","time"))
        g2_err = defVar(ds,"g2_err",Float64,("lon","lat","time"))
        g3_err = defVar(ds,"g3_err",Float64,("lon","lat","time"))
        g4_err = defVar(ds,"g4_err",Float64,("lon","lat","time"))

        # Fill the coord vectors and the fields
        lon[:] = gridlon;
        lat[:] = gridlat;
        time[:] = gridtime;

        g1[:,:,:] = field1;
        g2[:,:,:] = field2;
        g3[:,:,:] = field3;
        g4[:,:,:] = field4;

        g1_err[:,:,:] = err1;
        g2_err[:,:,:] = err2;
        g3_err[:,:,:] = err3;
        g4_err[:,:,:] = err4;
    end

end;


"""
```julia-repl
write_fish_specific_nc(filename, gridlon, gridlat, gridtime,
                    field1, field2, field3, field4,
                    err1, err2, err3, err4)
```
Write a netCDF file containing the gridded and error fields for the
specific fish product.
"""
function write_fish_specific_nc(filename::String, gridlon, gridlat,
        field1::Array, field2::Array, field3::Array, field4::Array,
        err1::Array, err2::Array, err3::Array, err4::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);

        # Define a global attribute
        ds.attrib["title"] = "Interpolated Fish"

        # Define the variables and coordinates
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))

        # Attributes
        lat.attrib["long_name"] = "Latitude";
        lat.attrib["standard_name"] = "latitude";
        lat.attrib["units"] = "degrees_north";

        lon.attrib["long_name"] = "Longitude";
        lon.attrib["standard_name"] = "longitude";
        lon.attrib["units"] = "degrees_east";

        # Interpolated fields
        g1 = defVar(ds,"g1",Float64,("lon","lat"))
        g2 = defVar(ds,"g2",Float64,("lon","lat"))
        g3 = defVar(ds,"g3",Float64,("lon","lat"))
        g4 = defVar(ds,"g4",Float64,("lon","lat"))

        g1.attrib["long_name"] = "Entelurus aequoreus";
        g2.attrib["long_name"] = "Molva molva";
        g3.attrib["long_name"] = "Sprattus sprattus";
        g4.attrib["long_name"] = "Squalus acanthias";

        # Error fields
        g1_err = defVar(ds,"g1_err",Float64,("lon","lat"))
        g2_err = defVar(ds,"g2_err",Float64,("lon","lat"))
        g3_err = defVar(ds,"g3_err",Float64,("lon","lat"))
        g4_err = defVar(ds,"g4_err",Float64,("lon","lat"))

        # Fill the coord vectors and the fields
        lon[:] = gridlon;
        lat[:] = gridlat;

        g1[:,:] = field1;
        g2[:,:] = field2;
        g3[:,:] = field3;
        g4[:,:] = field4;

        g1_err[:,:] = err1;
        g2_err[:,:] = err2;
        g3_err[:,:] = err3;
        g4_err[:,:] = err4;

    end

end;


"""
```julia-repl
plot_fish_results_time(years,gridlon, gridlat, f1, bx, by, b)
```
Create the general plot for the fish products.
One plot per fish type, all the years on the plot.
"""
function plot_fish_results_time(years,gridlon, gridlat, f1, bx, by, b)

    yearinit = collect(years) .-1;
    yearinit[1] = yearinit[2];
    yearend = collect(years) .+ 1;
    yearend[end] = yearend[end-1]

    ffig = figure("fish_time",figsize=(12,12))
    for ii = 1:size(f1)[3]
        ax = subplot(5,4,ii)
        title("$(yearinit[ii])-$(yearend[ii])", fontsize=8)
        ax[:tick_params]("both",labelsize=6)
        pcm = PyPlot.pcolormesh(gridlon, gridlat, permutedims(f1[:,:,ii], [2,1]), vmin=0, vmax=1.)
        add_mask(bx, by, b)
        gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
        colorbar(pcm)[:ax][:tick_params](labelsize=8)
    end
end
