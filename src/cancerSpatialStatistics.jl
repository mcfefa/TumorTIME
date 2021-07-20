#=

	TODO:
	- User needs to provide:
		CSV file to be read
		String of marker names?
		Cutoff (minimum number of pairwise distances)
		Output desired

	- Read CSV file
	- Compute cell center locations
	- Normalize to unit box
	- Compute pair-wise distances


=#

using CSV,DataFrames,Parameters

# Used to split a tuple by components
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function readFile(filename::String)

	# Check if the file exists
	if !isfile(filename)
		error("$filename not found.")
	end

	# return the data as a dataframe
	df = CSV.read(filename,DataFrame)

	# convert all headers to lowercase to avoid parsing errors
	rename!(df,lowercase.(names(df)))

	return df

end

function normalizeCoordinates(v::AbstractArray{T}) where T <: Real
	"""
	Normalize so coordinates live on the unit box

		v:	vector of numeric values

	"""

	min_v,max_v = minimum(v),maximum(v)

	return (v .- min_v)/(max_v - min_v)
end


function computeCellCenterLocations(df::DataFrame,normalize::Bool = true)

	# We calculate the centroid of x,y coordinates
	x_center,y_center = (df.xmin .+ df.xmax)/2.0, (df.ymin .+ df.ymax)/2.0

	# If normalization is desired (default is true)
	if normalize
		x_center = normalizeCoordinates(x_center)
		y_center = normalizeCoordinates(y_center)
	end

	return x_center, y_center
end

function extractCellTypes(df::DataFrame)

	# Dictionary which will contain info on indices of cell types.
	cellTypes = Dict{String,Vector{Int}}()

	#=
		Loop over the column names and find which ones end with the word
		positive (FIXME??). Then we fill the dictionary with the indices which
		indicate that the marker is present at that cell.
	=#
	for name in names(df)
		if occursin(r"positive$",name)	# regexp checks end of string
			firstWord = split(name)[1]		# FIXME??
			get!(cellTypes, firstWord) do 	# function to fill the dictionary
				findall(x->x==1,df[:,name])
			end
		end
	end

	# We need the complement of the cellTypes set to find which remaining cells
	# are stroma and tumor. Note that multiple markers can show up for the same
	# cell and so we need to remove unique ones.
	allTypes = unique(sort(vcat(values(cellTypes)...)))

	# We assume (or hope) that tumor and stroma are in the data frame.
	# This can potentially be the users responsibility.
	for name in ["Tumor","Stroma"] ### USER-DEFINED?
		get!(cellTypes, lowercase(name)) do
			setdiff(findall(x->x==name,df[:,end]),allTypes)
		end
	end

	x_center, y_center = computeCellCenterLocations(df)

	cellLocations = Dict{String,Any}()

	for (key,val) in cellTypes
		get!(cellLocations, key) do
			(x_center[val],y_center[val])
		end
	end

	return cellTypes,cellLocations

end

function computeDistances(locations::Dict)

	# We will save the intra and inter-pairwise distances to individual 
	# dictionaries
	intradist = Dict{String,Vector{Float64}}()
	interdist = Dict{String,Vector{Float64}}()

	# Save enumerator
	enumKeys = enumerate(keys(cellLocations))

	# Calculate inter- and intra-type pair-wise distances
	for (i,k1) in enumKeys
		intradist[k1] = intraPairwiseDistances(cellLocations[k1]...)
		for (j,k2) in enumKeys
			if i > j
				interdist[k1*"/"*k2] = interPairwiseDistances(
										cellLocations[k1]...,
										cellLocations[k2]...)
			end
		end
	end

	return interdist,intradist

end

function intraPairwiseDistances(x::T,y::T) where {T<:AbstractArray}

	@assert length(x) == length(y)
	N = length(x)
	[ sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) for i in 1:N for j in 1:N if i<j]

end

function interPairwiseDistances(x1::T,y1::T,x2::T,y2::T) where {T<:AbstractArray}

	@assert length(x1) == length(y1) && length(x2) == length(y2)
	N1,N2 = length(x1),length(x2)
	[ sqrt((x1[i]-x2[j])^2 + (y1[i]-y2[j])^2) for i in 1:N1 for j in 1:N2 ]

end

function nullpdf(r::Real)
	if r<=1
		val = 2r*(pi - 4*r + r^2)
	else
		val = 2*r*(2*(atan(1/sqrt(r^2 - 1)) - atan(sqrt(r^2 - 1)) - 1) + 4*sqrt(r^2 - 1) - r^2)
	end
	return val
end

function makePlots(dist::Dict;len::Int=1800,wid::Int=2400,
	biomarkerColors::Symbol=:rainbow,
	classifierColors::Vector{String}=["lightgray","pink"],
	)

	# Generate figure handle
	fig = Figure(resolution = (len, wid))

	# Create empty lists for storing axes and figure info
	ax,f = [],[]
	
	c = cgrad(:rainbow,15,categorical=true);
	c = [color.(classifierColors),c...];
	row,col = 1,1
	for n = 1 : length(dist)
		axtemp = fig[row,col] = Axis(fig,title=titles[n],
		xlabel = "Distance", ylabel = "Density")
		ftemp=hist!(axtemp,dist[n],normalization=:pdf,color=c[n],
		strokewidth=1,strokecolor=:black)
		lines!(axtemp,x,map(nullpdf,x))
		push!(ax,axtemp);push!(f,ftemp);
		if col > 1
			hideydecorations!(axtemp,grid=false)
		end
		if col < 3
			col += 1
		else
			col = 1
			row += 1
		end
	end
	
	linkaxes!(ax...)
	
	return fig
	
	end

# function makePlots(plotinfo::Dict,cutoff::Int=0;
# 	xrotation::Int=0,α::Real=1.0)

# 	# Histogram will plot only if there are at least 1 element in the set
# 	title,dist = unzip([(key,val) for (key,val) in plotinfo
# 	if length(val)>cutoff])

# 	numPlots = length(title)

# 	if numPlots < 5
# 		layout = (numPlots,1)
# 	else
# 		sqrootplt=sqrt(numPlots)
# 		mygrid = convert.(Int,(ceil(sqrootplt),floor(sqrootplt)))
# 	end
	
# 	myplot=plot(dist,
# 	normalize=:pdf,
# 	bins=:sturges,
# 	layout=mygrid,
# 	title=permutedims(title),
# 	label=nothing,
# 	seriestype=:barhist,
# 	xrotation=xrotation,
# 	tickfontsize=α,
# 	titlefontsize=α,
# 	yticks=0:0.5:2
# 	)

# 	x = 0:0.01:sqrt(2)

# 	plot!([(x,map(nullpdf,x)) for i in 1:numPlots],
# 	layout=mygrid,
# 	label=nothing
# 	)

# 	return myplot

# end


function main(user_par=nothing)

	# Grab the parameters from the struct
	if user_par !== nothing
		@unpack (inputfile,
		classifierLabels,
		numPairCutoff,
		savePlots,
		saveData
		) = user_par
	else
		inputfile = "SP09-1997 A8_[52341,12388].tif_11881_job18745.object_results.csv"
	end

	#= For now, we keep the df, but we could just call it directly to
	save memory if needed =#
	df = readFile(inputfile)

	# Get inter and intra cellular distances
	interdist,intradist = computeDistances(df)

	return df,interdist,intradist


end



#

#plotting method?
# histogram([data for data in dist if length(data)>50],normalize=:pdf,bins=:scott)

# plot([data for data in intradist if length(data)>50],normalize=:pdf,bins=:sturges,layout=(3,1),label=nothing,seriestype=:barhist,xrotation=45,xtickfontsize=6,ytickfontsize=6,seriescolor=[:green :red :blue])