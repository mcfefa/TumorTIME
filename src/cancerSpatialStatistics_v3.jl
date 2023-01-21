#=

	TODO:
	- User needs to provide:
		CSV file to be read
		String of marker names?
		Cutoff (minimum number of pairwise distances)
		Output desired

	- Read CSV file
	- Compute cell center location
	- Normalize to unit box
	- Compute pair-wise distances


=#

using CSV, DataFrames, Parameters, Statistics, StatsBase, HypothesisTests, Cairo, CairoMakie, Colors, Makie, Distributions, QuadGK, EmpiricalDistributions, Dates, Query, Missings, Distances;
import Distributions: cdf
# Used to split a tuple by components
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function readFile(filename::String)

	# Check if the file exists
	if !isfile(filename)
		error("$filename not found.")
	end

	# return the data as a dataframe
	df = CSV.read(filename,DataFrame; header=true)

	# convert all headers to lowercase to avoid parsing errors
	for i in 1:size(names(df))[1]
		validity = isvalid(names(df)[i])
		if validity == true
			# rename!(df, Dict(:i => "A", :x => "X"))
			rename!(df, Dict(names(df)[i] => lowercase(names(df)[i])))
		else
			# replace(str,r => "" )
			fixedName = replace(names(df)[i], "\xb5m" => "micro")
			fixedName2 = replace(fixedName, "\xb2" => "sq")
			rename!(df, Dict(names(df)[i] => lowercase(fixedName2)))
		end

	end

	# rename!(df,lowercase.(filter(isvalid.(names(df),names(df)))))
	# replaced this with for loop above, mu in header was throwing error with lowercase()

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

function computeDistances(df::DataFrame)

	# We will save the intra and inter-pairwise distances to individual 
	# dictionaries
	intradist = Dict{String,Vector{Float64}}()
	interdist = Dict{String,Vector{Float64}}()

	cellTypes,cellLocations = extractCellTypes(df)

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



function my_cdf(r::Real)
	I, _ = quadgk(nullpdf, 0, r)
	return I
end

struct customDist <: ContinuousUnivariateDistribution

end

function cdf(::customDist, x::Float64)
	my_cdf(x)
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
	x = (0:0.01:sqrt(2))
	for (n,(k,v)) in enumerate(dist)
		axtemp = fig[row,col] = Axis(fig,title=titles[n],
		xlabel = "Distance", ylabel = "Density")
		ftemp=hist!(axtemp,v,normalization=:pdf,color=c[n],
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

#function makePlots(plotinfo::Dict,cutoff::Int=0;
 #	xrotation::Int=0,α::Real=1.0)
 #
 #	# Histogram will plot only if there are at least 1 element in the set
 #	title,dist = unzip([(key,val) for (key,val) in plotinfo
 #	if length(val)>cutoff])
 #
 #	numPlots = length(title)
 #
 #	if numPlots < 5
 #		layout = (numPlots,1)
 #	else
 #		sqrootplt=sqrt(numPlots)
 #		mygrid = convert.(Int,(ceil(sqrootplt),floor(sqrootplt)))
 #	end
 #
 #	myplot=plot(dist,
 #	normalize=:pdf,
 #	bins=:sturges,
 #	layout=mygrid,
 #	title=permutedims(title),
 #	label=nothing,
 #	seriestype=:barhist,
 #	xrotation=xrotation,
 #	tickfontsize=α,
 #	titlefontsize=α,
 #	yticks=0:0.5:2
 #	)

 #	x = 0:0.01:sqrt(2)

 #	plot!([(x,map(nullpdf,x)) for i in 1:numPlots],
 #	layout=mygrid,
 #	label=nothing
 #	)

 #	return myplot

 #end


function main(user_par=nothing)

	#dir = "C:\\Users\\camara.casson\\OneDrive - University of Florida\\Desktop\\TumorTIME\\src"
	dir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor"
	# Grab the parameters from the struct
	if user_par !== nothing
		show(user_par)

		@unpack (inputfile,
		classifierLabels,
		numPairCutoff,
		savePlots,
		saveData
		) = user_par

		#inputfile = user_par
	else
		inputfile = "SP09-1997 A8_[52341,12388].tif_11881_job18745.object_results.csv"
	end

	#inputfile = user_par
	# show(inputfile)
	#= For now, we keep the df, but we could just call it directly to save memory if needed =#
	df = readFile(string(dir,"\\",inputfile))	
	#show(df)

	# Get inter and intra cellular distances
	interdist,intradist = computeDistances(df)

	

	return df,interdist,intradist,inputfile


end


function TumorTIMEPipeline(directory1, file, marker, panelName, panelLoc)

	# # directory1 is a directory to pull files from 

    # # file wants a string name of the file that has the single-cell measurements in it
	#file = "P1a_SP04-4722_[37094,12061].tif_94266_job58256.object_results.csv"

	#marker wants a vector of strings that are the markers in the panel
	#marker = ["tumor","stroma","CD68","CD163","CD206","PD-L1"]

	#panelName wants a string that identifies which panel subfolder ("Panel-1" or "Panel-2") to use
	#panelLoc wants a string that identifies which location folder (tumor, stroma, interface, normal) to use

	# identifies the patient name based on the SP number
	println(file)
	indName = string(split(file,"_")[2])
	

	# ## identifies the patient name with the location on the slide
	splitname = string(split(file,".")[1])


	df, interdist, intradist, inputfile = main((inputfile=file,classifierLabels=marker,numPairCutoff=5,savePlots=false,saveData=false))

	stats = [:mean,:median,:std,:var,:kurtosis,:skewness]
	interdist_stats = Dict(k => NamedTuple(stat => eval(stat)(v) for stat in stats) for (k,v) in interdist if !isempty(v))
	intradist_stats = Dict(k => NamedTuple(stat => eval(stat)(v) for stat in stats) for (k,v) in intradist if !isempty(v))

	interdist_stats_df = DataFrame(name=[], mean=[], median=[], std=[], var=[], kurtosis=[], skewness=[])
	for (k,v) in interdist_stats
		append!(interdist_stats_df, hcat(DataFrame(name = k), DataFrame([v])))
	end

	function MakeDistributions(data)
		result = Dict()
		for (k,v) in data 
			if !isempty(v)
				Hist = fit(Histogram, v)
				dist = UvBinnedDist(Hist)
				result[k] = dist
			end
		end
		return result
	end

    ## creates distribution based on interdistances calculated from slide information
    interdist_distr = MakeDistributions(interdist)
    println(interdist_distr)

    ##Creates distribution based on theoretical pdf equation
    ## add in theoretical PDF calculation here with range related to min and max of interdist in the sample
    
	#NamesofInterdistP1=
	NamesOfInterdistP2=["cd68/cd20", "nucleus/pd-l1", "cd20/stroma", "cd68/stroma", "nucleus/tumor", "cd163/nucleus", "cd20/nucleus", "cd206/stroma", "cd206/tumor", "cd206/cd163", "pd-l1/tumor", "cd206/cd20", "pd-l1/stroma", "cd163/pd-l1", 
	"cd163/stroma", "cd68/nucleus", "cd20/tumor", "cd206/pd-l1", "nucleus/stroma", "cd20/pd-l1", "cd163/tumor", "cd206/cd68", "cd163/cd20", "tumor/stroma", "cd68/pd-l1", "cd206/nucleus", "cd68/cd163", "cd68/tumor"]
   

    for i in eachindex(NamesOfInterdistP2)
		if !isempty(interdist[NamesOfInterdistP2[i]])
    		min1= minimum(interdist[NamesOfInterdistP2[i]])
    		println(min1)
    		max1= maximum(interdist[NamesOfInterdistP2[i]])
    		println(max1)
    		size1=length(interdist[NamesOfInterdistP2[i]])
   			println(size1)
    		range1= LinRange(min1, max1, size1)
   			println(range1)		
				for j in range1
				TheoreticalPDF = nullpdf(range1[j])
				end
		end
	end
end
TheoreticalPDF = nullpdf(min1:incr:max1)

	# function RunCramerVonMisesTest(data,stats)
	# 	CVMTest_col = DataFrame(CVM_result = [], CVM_p=[])
	# 	for (k,v) in data
	# 		if !isempty(v)
	# 			res1 = 

# 	function RunKSTest(dataEmp, dataNull, stats)
#         ## dataEmp <DEFINE> 
#         ## dataNull <DEFINE>
#         ## stats is the dataframe with distribution statistcs that gets appended 

# 		kstests_col = DataFrame(ks_result = [], ks_p = [])
# 		for (k,v) in data
# 			if !isempty(v)
# 				res = ExactOneSampleKSTest(v, customDist())
# 				append!(kstests_col, DataFrame(ks_result=res.δ, ks_p = pvalue(res)))
# 			end
# 		end
# 		stats=hcat(stats,kstests_col)
# 		return stats
# 	end

# 	KSResults = RunKSTest(interdist,interdist_stats_df)

# 	# function RunADTest(data, stats)
# 	# 	adtests_col = DataFrame(ad_result = [], ad_p = [])
# 	# 	for (k,v) in data
# 	# 		if !isempty(v)
# 	# 			res1 = OneSampleADTest(v, customDist())
# 	# 			append!(adtests_col, DataFrame(ad_result = res1.A², ad_p=(pvalue(res1))))
# 	# 		end
# 	# 	end
# 	# 	stats=hcat(stats,adtests_col)
# 	# 	return stats
# 	# end

# 	# ADResults=RunADTest(interdist,KSResults)

# 	# add a column to statistics DF that has patient name 
# 	KSResults[:,:patient] .= indName 

# 	# write interdistance statistcs to a CSV file 
# 	@time CSV.write(string(splitname,".interdistances_stats",Dates.today(),".csv"), KSResults)

# 	coltypesS = Any[Float64 for i=1:35]
# 	coltypesS[1]=String
# 	coltypesS[2]=String

# 	statsdir="C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results and Analysis\\"
# 	@time KSResults1 = CSV.read(string(statsdir,"All-Clinical-Data-for-Interdistances-at-Panel1-Tumor2022-10-27.csv"),DataFrame,types=coltypesS)
# 	KSResults = KSResults1[:,1:12]

# 	println(KSResults[1:5,:])

# 	# loading clinical information
# 	println("INFO: Time to load clinical data ")
# 	rawdatadirbase = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data"
# 	rawdatadir = string(rawdatadirbase,"\\", panelName, "\\", panelLoc)
# 	clindir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data" 
# 	coltypes1 = Any[String for i=1:30]
# 	coltypes1[4]=Union{Missing, String}   ## EGFR norm reads
# 	coltypes1[5]=Union{Missing, String}   ## EGFR alpha reads
# 	coltypes1[6]=Union{Missing, String}   ## EGFR beta reads
# 	coltypes1[7]=Union{Missing, String}   ## EGFR gamma reads
# 	coltypes1[10]=Float64    ## age at diagnosis
# 	coltypes1[11]=Float64    ## age at surgery
# 	coltypes1[14]=Union{Missing, Int}        ## grade
# 	coltypes1[15]=Float64    ## size
# 	coltypes1[22]=Union{Missing, String}   ## RFS
# 	coltypes1[24]=Union{Missing, String}   ## CauseOfDeath
# 	coltypes1[25]=Int        ## event of death
# 	coltypes1[26]=Int 		 ## overall survival
# 	coltypes1[27]=Union{Missing, Float64}	   ## OS-IT
# 	coltypes1[28]=Union{Missing, Float64}	   ## OS-TT
# 	coltypes1[29]=Union{Missing, Float64}	   ## IMDC
# 	coltypes1[30]=Union{Missing, String}   ## sutent neoadvant
# 	@time clinical_raw = CSV.read(string(clindir,"\\ccRCC-TIME-clinical-2022-09-27.csv"), DataFrame, header=1, types=coltypes1);
	
# 	# show(ADResults)
# 	# show(clinical_raw)

# 	joinkey(i) = (i.patient)

# 	#println(joinkey(KSResults))
# 	#println(joinkey(clinical_raw))

# 	## !!!!!!!!!!!!!!!!!!!!!
# 	## double check that all clinical variables are present 
# 	@time DFFinal = @from i in KSResults begin
# 		@join j in clinical_raw on joinkey(i) equals joinkey(j)
# 		@select {i.patient, 
# 				i.name, 
# 				i.mean, 
# 				i.median, 
# 				i.std, i.var, 
# 				i.kurtosis, i.skewness, i.ks_result, i.ks_p, 
# 				i.ad_result, i.ad_p, 
# 				EGFRnormReads=j.EGFRnormReads, 
# 				EGFRalphaReads=j.EGFRalphaReads,
# 				EGFRbetaReads=j.EGFRbetaReads,
# 				EGFRgammaReads=j.EGFRgammaReads,
# 				Gender=j.Gender, Race=j.Race, 
# 				AgeDiagnosis=j.AgeDiagnosis, AgeSurgery=j.AgeSurgery, 
# 				Histology=j.Histology, Laterality=j.Laterality, 
# 				Grade=j.Grade, Size=j.Size, 
# 				SarcomatoidStatus=j.SarcomatoidStatus, 
# 				RhabdoidStatus=j.RhabdoidStatus, 
# 				pT=j.pT, pN=j.pN, pM=j.pM, 
# 				CytoreductiveSurgStatus=j.CytoreductiveSurg, 
# 				RFS=j.RFS, OS=j.OverallSurvival, CauseOfDeath=j.CauseOfDeath, 
# 				DeathStatus=j.EventDeath
# 		}
# 		@collect DataFrame;
# 	end

# 	# write merged statistics and clinical data for a single patient to a file 
# 	@time CSV.write(string(directory1, "\\", splitname,".interdistances_stats+clin",Dates.today(),".csv"), DFFinal)
# 	println(string(splitname, " interdistance statistics and clinical data has been saved to CSV"))

# end

NamesOfInterdistP2=["cd68/cd20", "nucleus/pd-l1", "cd20/stroma", "cd68/stroma", "nucleus/tumor", "cd163/nucleus", "cd20/nucleus", "cd206/stroma", "cd206/tumor", "cd206/cd163", "pd-l1/tumor", "cd206/cd20", "pd-l1/stroma", "cd163/pd-l1", 
"cd163/stroma", "cd68/nucleus", "cd20/tumor", "cd206/pd-l1", "nucleus/stroma", "cd20/pd-l1", "cd163/tumor", "cd206/cd68", "cd163/cd20", "tumor/stroma", "cd68/pd-l1", "cd206/nucleus", "cd68/cd163", "cd68/tumor"]

# file1 = "P1a_SP04-4722_[37094,12061].tif_94266_job58256.object_results.csv"

## Panel 2 Markers
markerPanel = ["tumor","stroma","CD68","CD163","CD206","PD-L1"]

directory1 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor"
## Reads all files in the directory
##    note, this is looking at the interface location
ReadingFiles = readdir("C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor")

for i in 1:2 #size(ReadingFiles)[1]-4
	filenametemp = ReadingFiles[i]
	TumorTIMEPipeline(directory1, filenametemp, markerPanel, "Panel-2", "tumor")
end

# function ConcatFiles(directory)

# 	#Read files in results from pipeline folder  
# 	ReadingFiles2 = readdir(directory)
# 	println(ReadingFiles2)

# 	coltypesL = Any[String for i=1:34]
# 	coltypesL[1]=String
# 	coltypesL[2]=String
# 	coltypesL[13]=Union{Missing, String}
# 	coltypesL[14]=Union{Missing, String}
# 	coltypesL[15]=Union{Missing, String}
# 	coltypesL[16]=Union{Missing, String}
# 	coltypesL[17]=String
# 	coltypesL[18]=String
# 	coltypesL[19]=String
# 	coltypesL[20]=String
# 	coltypesL[21]=String
# 	coltypesL[22]=String
# 	coltypesL[23]= Union{Missing, Float64}
# 	coltypesL[24]=Union{Missing, Float64}
# 	coltypesL[25]=String
# 	coltypesL[26]=String
# 	coltypesL[27]=String
# 	coltypesL[28]=Union{Missing, String}
# 	coltypesL[29]=String
# 	coltypesL[30]=String
# 	coltypesL[31]=Union{Missing, String}
# 	longDF = CSV.read(string(directory,"\\",ReadingFiles2[1]), DataFrame, header=1, types=coltypesL)#Dict(:col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int, :col31=>Union{Missing,Int,String}))
# 	# add a column to individual patient and loc DF that has the filename to be able to pull out location later if needed
# 	longDF[:,:filename] .= ReadingFiles2[1] 

# 	#println(names(longDF)[31])
# 	println(longDF)

# 	#Create for loop that runs through files of stats and clincial data 
# 	for i in 2:size(ReadingFiles2)[1]-1
# 		filenametemp2 = ReadingFiles2[i]
# 		println(filenametemp2)

# 		#tmpdf = CSV.read(string(directory,"\\",ReadingFiles2[i]), DataFrame, header=1, types=Dict(:col31=>Union{Missing,Int,String}, :col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int))
# 		coltypes = Any[String for i=1:34]
# 		coltypes[1]=String
# 		coltypes[2]=String
# 		coltypes[13]=Union{Missing, String}
# 		coltypes[14]=Union{Missing, String}
# 		coltypes[15]=Union{Missing, String}
# 		coltypes[16]=Union{Missing, String}
# 		coltypes[17]=String
# 		coltypes[18]=String
# 		coltypesL[19]=String
# 		coltypesL[20]=String
# 		coltypes[21]=String
# 		coltypes[22]=String
# 		coltypes[23]=Union{Missing, Float64}
# 		coltypes[24]=Union{Missing, Float64}
# 		coltypes[25]=String
# 		coltypes[26]=String
# 		coltypes[27]=String
# 		coltypes[28]=Union{Missing, String}
# 		coltypes[29]=String
# 		coltypes[30]=String
# 		coltypes[31]=Union{Missing, String}
# 		#tmpdf = CSV.read(string(directory,"\\",filenametemp2), DataFrame, header=1, types=Dict(:col31=>Union{Missing,String}, :col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int))
# 		tmpdf = CSV.read(string(directory,"\\",filenametemp2), DataFrame, header=1, types=coltypes)
# 		tmpdf[:,:filename] .= ReadingFiles2[i] 

# 		append!(longDF,tmpdf)
# 	end

# 	# write out CSV of full dataframe --- longDF
# 	@time CSV.write(string(directory1,"\\All-Clinical-Data-for-Interdistances-at-Panel2-Tumor",Dates.today(),".csv"), longDF)

# 	# Query to identify individual pairs and then save those as files 
# 	#Loop over each of the different names
# 	listNames = unique(longDF[:,:name])

# 	for n in 1:length(listNames)
# 		Query1= @from i in longDF begin
# 			@where i.name == listNames[n]
# 			@select i #{pair=i.name, i.filename}
# 			@collect DataFrame
# 		end
# 		replace(listNames[n], "/" => "+") 
# 		@time CSV.write(string(directory1,"\\","Query for ", replace(listNames[n], "/" => "+"), Dates.today(),".csv"), Query1)
# 	end




# end

# individdir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-1\\tumor\\Results from Pipeline\\_CutOff3\\"
# #individdir2="C:\\Users\\camara.casson\\OneDrive - University of Florida\\Desktop\\TumorTIME\\src"
# ConcatFiles(individdir)




# #plotting method?
# # histogram([data for data in dist if length(data)>50],normalize=:pdf,bins=:scott)

# # plot([data for data in intradist if length(data)>50],normalize=:pdf,bins=:sturges,layout=(3,1),label=nothing,seriestype=:barhist,xrotation=45,xtickfontsize=6,ytickfontsize=6,seriescolor=[:green :red :blue])