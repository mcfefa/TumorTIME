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

using CSV, DataFrames, Parameters, Statistics, StatsBase, HypothesisTests, Cairo, CairoMakie, Colors, Makie, Distributions, QuadGK, EmpiricalDistributions, Dates, Query, Missings, Distances, StatsPlots;
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

# function makePlots(plotinfo::Dict,cutoff::Int=0;
#  	xrotation::Int=0,α::Real=1.0)
 
#  	# Histogram will plot only if there are at least 1 element in the set
#  	title,dist = unzip([(key,val) for (key,val) in plotinfo
#  	if length(val)>cutoff])
 
#  	numPlots = length(title)
 
#  	if numPlots < 5
#  		layout = (numPlots,1)
#  	else
#  		sqrootplt=sqrt(numPlots)
#  		mygrid = convert.(Int,(ceil(sqrootplt),floor(sqrootplt)))
#  	end
 
#  	myplot=plot(dist,
#  	normalize=:pdf,
#  	bins=:sturges,
#  	layout=mygrid,
#  	title=permutedims(title),
#  	label=nothing,
#  	seriestype=:barhist,
#  	xrotation=xrotation,
#  	tickfontsize=α,
#  	titlefontsize=α,
#  	yticks=0:0.5:2
#  	)

#  	x = 0:0.01:sqrt(2)

#  	plot!([(x,map(nullpdf,x)) for i in 1:numPlots],
#  	layout=mygrid,
#  	label=nothing
#  	)

#  	return myplot

end


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


function TumorTIMEPipeline_Inter(directory1, file, marker, panelName, panelLoc)

	# # directory1 is a directory to pull files from 
	## directory1 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor"

    # # file wants a string name of the file that has the single-cell measurements in it
	#file = "P2a_SP04-4722_[41516,9786].tif_95219_job59171.object_results.csv"

	#marker wants a vector of strings that are the markers in the panel
	#marker = ["tumor","stroma","CD68","CD163","CD206","PD-L1"]

	#panelName wants a string that identifies which panel subfolder ("Panel-1" or "Panel-2") to use
		## panelName = "Panel-2"
	#panelLoc wants a string that identifies which location folder (tumor, stroma, interface, normal) to use
		## panleLoc = "tumor"

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

	#Plotting data with StatsPlots to visualize distributions
	# histogram(interdist["cd68/stroma"])
	# dist = A["cd68/stroma"]._bin_pdf
	# StatsPlots.scatter(dist, leg=false)
	# StatsPlots.bar!(dist, func=cdf, alpha=0.3)

	function MakeDistributions(data)
		#Creates distributions fitting with a maximum likelihood estimator to a DiscreteNonParametric distribution
		#Creates a dictionary with distributions for each marker pair (non-empty)
		result = Dict()
		for (k,v) in data 
			if !isempty(v)
				Dist = fit_mle(DiscreteNonParametric,v); 
				result[k] = Dist
			end
		end
		return result
	end
    ## creates distribution based on interdistances calculated from slide information
    interdist_distr = MakeDistributions(interdist)
    println(interdist_distr)
	#Plotting the PDF for empirical and theoretical distribution
	#StatsPlots.scatter(hcat(TheoreticalPDFs["cd68/tumor"]...)[1,:],hcat(TheoreticalPDFs["cd68/tumor"]...)[2,:])
	#StatsPlots.bar!(z,alpha=0.3)

	#OLD, Don't need this
	# interdist_distr2 = Dict{String,Vector{Float64}}()
    # for i in keys(interdist_distr)
	# 	interdist_distr[i]._bin_probmass
	# 	interdist_distr2[i] = interdist_distr[i]._bin_probmass
	# end

    ##Creates distribution based on theoretical pdf equation
    ## add in theoretical PDF calculation here with range related to min and max of interdist in the sample
    
   #This prints out all the keys individually 
    # for i in keys(interdist)
    #     println(i)
    # end
   
    #creates a blank dictionary for the ranges of the interdistances at random distances
	InterdistRanges = Dict{String,Vector{Float64}}()
    #creates a blank dictionary for the distribution of the ranges of the interdistances
    TheoreticalPDFs = Dict{String,AbstractVector}()
    for i in keys(interdist_distr)
        #if !isempty(interdist[i])
            TheoreticalPDF = [];
			R = range(minimum(interdist[i]),maximum(interdist[i]),length=length(interdist[i]))
			for s in R
                # println("entered for loop")
			    push!(TheoreticalPDF,[s,nullpdf(s)])
			end
            TheoreticalPDFs[i]=TheoreticalPDF
	end
	#Plotting the nullpdf
	#StatsPlots.scatter(hcat(TheoreticalPDFs["cd68/stroma"]...)[1,:],hcat(TheoreticalPDFs["cd68/stroma"]...)[2,:])
    directory2 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results-from-Pipeline"

    # #write merged theoretical ranges to a file 
 	# @time CSV.write(string(directory2, "\\", splitname,".interdistances_ranges",Dates.today(),".csv"), InterdistRanges)
 	# println(string(splitname, " interdistance theoretical ranges have been saved to CSV"))

    #write merged theoretical PDFs to a file 
    @time CSV.write(string(directory2, "\\", splitname,".interdistances_PDFs",Dates.today(),".csv"), TheoreticalPDFs, bufsize=128*1024*1024)
    println(string(splitname, " interdistance theoretical distributions have been saved to CSV"))

	#calculate CDF from nullpdf using integration
	TheoreticalCDFs = Dict{String,AbstractVector}()
    for i in keys(interdist_distr)
        #if !isempty(interdist[i])
            TheoreticalCDF = [];
			R = range(minimum(interdist[i]),maximum(interdist[i]),length=length(interdist[i]))
			for s in R
                # println("entered for loop")
			    push!(TheoreticalCDF,[s,my_cdf(s)])
			end
            TheoreticalCDFs[i]=TheoreticalCDF
	end
#StatsPlots.scatter(hcat(TheoreticalCDFs["cd68/stroma"]...)[1,:],hcat(TheoreticalCDFs["cd68/stroma"]...)[2,:])

	#Calculate CDF of empirical distribution
	function MakeCDFs(data)
		#Creates a CDF distribution for empirical data
		#Creates a dictionary with distributions for each marker pair (non-empty)
		result = Dict()
		for (k,v) in data 
			if !isempty(v)
				Dist = ecdf(v) 
				result[k] = Dist
			end
		end
		return result
	end
    ## creates distribution based on interdistances calculated from slide information
    interdist_CDF = MakeCDFs(interdist)
    
	
	EmpiricalCDFs = Dict{String,AbstractVector}()
    for i in keys(interdist_distr)
            EmpiricalCDF = [];
			Interdist_Sort = sort(interdist[i])
			#println("Show Interdist_Sort")
			for s in 1:length(interdist_CDF[i].sorted_values)
                # println("entered for loop")
			    push!(EmpiricalCDF,[Interdist_Sort[s],interdist_CDF[i](Interdist_Sort[s])])
			end
            EmpiricalCDFs[i]=EmpiricalCDF
	end
	#StatsPlots.scatter!(hcat(EmpiricalCDFs["cd68/stroma"]...)[1,:],hcat(EmpiricalCDFs["cd68/stroma"]...)[2,:])
	
	function MakeDNPTheo(dataNull)
		DNPTheos = Dict{String,DiscreteNonParametric}()
		for i in keys(dataNull)
			dist1 = fit_mle(DiscreteNonParametric,hcat(dataNull[i]...))
			DNPTheos[i] = dist1
		end
		return DNPTheos
	end
	C = MakeDNPTheo(TheoreticalCDFs)


	function RunKSTest(dataEmp, dataNull, stats)
        ## dataEmp is the empirical data with the made distributions found in EmpiricalCDFs
        ## dataNull the theoretical distribution found inside TheoreticalCDFs
        ## stats is the dataframe with distribution statistcs that gets appended 

		kstests_col = DataFrame(ks_result = [], ks_p = [])
		for i in keys(dataNull) 
				A = hcat(dataEmp[i]...)[2,:]
				# x1 = vcat([s[1] for i in dataEmp])
				# x2 = vcat([s[1] for i in dataNull])
				# y1 = vcat([s[2] for i in dataEmp])
				# y2 = vcat([s[2] for i in dataNull])
				res = ExactOneSampleKSTest(A,C[i])#(dataEmp[i],dataNull[i],Any,Any) #(dataEmp[i][s][1], dataEmp[i][s][2],dataNull[i][s][1], dataNull[i][s][2]) #[x1,y1],[x2,y2]
				append!(kstests_col, DataFrame(ks_result=res.δ, ks_p = pvalue(res)))
		end

		stats=hcat(stats,kstests_col)
		return stats
	end

	KSResults = RunKSTest(EmpiricalCDFs,TheoreticalCDFs,interdist_stats_df)
	display(KSResults)

	 ADResults=RunADTest(EmpiricalCDFs,TheoreticalCDFs,KSResults)

	# function RunCramerVonMisesTest(dataEmp, dataNull, stats)
	# 	CVMTest_col = DataFrame(CVM_result = [], CVM_p=[])
	# 	for i in keys(dataNull) 	
	# 		dist3 = fit(DiscreteNonParametric,dataNull[i])
	# 		for j in 
	# 	CVM_Results =  

	#  ADResults=RunADTest(EmpiricalCDFs,TheoreticalCDFs,KSResults)

	# function RunCramerVonMisesTest(dataEmp, dataNull, stats)
	# 	CVMTest_col = DataFrame(CVM_result = [], CVM_p=[])
	#	for i in keys(dataNull) 	
	#		dist3 = fit(DiscreteNonParametric,dataNull[i])
	# 		for j in 
	#	CVM_Results =  

	println("Starting KL Divergence")
	function KullbackLeibler(dataEmp,dataNull,stats)
		kld_col = DataFrame(kld_result= [])
		for i in keys(dataNull)
			A = hcat(dataEmp[i]...)[2,:]
			B = hcat(dataNull[i]...)[2,:]
			kld_result = kl_divergence(A,B)
			append!(kld_col, DataFrame(kld_result=kld_result))
		end
 
		stats=hcat(stats,kld_col)
		return stats
	end
	KLDResults = KullbackLeibler(EmpiricalCDFs,TheoreticalCDFs,KSResults)
	println("Finished KL Divergence")
	
	println("Starting Chebyshev Distance")
	function Chebyshev(dataEmp,dataNull,stats)
		cbs_col = DataFrame(cbs_result= [])
		for i in keys(dataNull)
			A = hcat(dataEmp[i]...)[1,:]
			B = hcat(dataNull[i]...)[1,:]
			cbs_result = chebyshev(A,B)
			append!(cbs_col, DataFrame(cbs_result=cbs_result))
		end
 
		stats=hcat(stats,cbs_col)
		return stats
	end
	CBSResults = Chebyshev(EmpiricalCDFs,TheoreticalCDFs,KLDResults)
	println("Finished Chebyshev Distance")

	println("Starting Jensen-Shannon Divergence")
	function JensenShannon(dataEmp,dataNull,stats)
		js_col = DataFrame(js_result= [])
		for i in keys(dataNull)
			A = hcat(dataEmp[i]...)[2,:]
			B = hcat(dataNull[i]...)[2,:]
			js_result = js_divergence(A,B)
			append!(js_col, DataFrame(js_result=js_result))
		end
 
		stats=hcat(stats,js_col)
		return stats
	end
	JSResults = JensenShannon(EmpiricalCDFs,TheoreticalCDFs,CBSResults)
	println("Finished Jensen-Shannon Divergence")

 	# add a column to statistics DF that has patient name 
 	JSResults[:,:patient] .= indName 

 	# write interdistance statistcs to a CSV file 
	 statsdir="C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\stats\\"
 	@time CSV.write(string(statsdir, "\\", splitname,".interdistances_stats",Dates.today(),".csv"), JSResults)

	coltypesS = Any[Float64 for i=1:35]
 	coltypesS[1]=String
 	coltypesS[2]=String

	statsdir="C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\"
	@time JSResults1 = CSV.read(string(statsdir,"All-Clinical-Data-for-Interdistances-at-Panel1-Tumor2022-10-27.csv"),DataFrame,types=coltypesS)
	JSResults = JSResults1[:,1:12]

	println(JSResults[1:5,:])

	# loading clinical information
	println("INFO: Time to load clinical data ")
	rawdatadirbase = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data"
	rawdatadir = string(rawdatadirbase,"\\", panelName, "\\", panelLoc)
	clindir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results-from-Pipeline\\_CutOff5\\Stats" 
	coltypes1 = Any[String for i=1:30]
	coltypes1[4]=Union{Missing, String}   ## EGFR norm reads
	coltypes1[5]=Union{Missing, String}   ## EGFR alpha reads
	coltypes1[6]=Union{Missing, String}   ## EGFR beta reads
	coltypes1[7]=Union{Missing, String}   ## EGFR gamma reads
	coltypes1[10]=Float64    ## age at diagnosis
	coltypes1[11]=Float64    ## age at surgery
	coltypes1[14]=Union{Missing, Int}        ## grade
	coltypes1[15]=Float64    ## size
	coltypes1[22]=Union{Missing, String}   ## RFS
	coltypes1[24]=Union{Missing, String}   ## CauseOfDeath
	coltypes1[25]=Int        ## event of death
	coltypes1[26]=Int 		 ## overall survival
	coltypes1[27]=Union{Missing, Float64}	   ## OS-IT
	coltypes1[28]=Union{Missing, Float64}	   ## OS-TT
	coltypes1[29]=Union{Missing, Float64}	   ## IMDC
	coltypes1[30]=Union{Missing, String}   ## sutent neoadvant
	@time clinical_raw = CSV.read(string(clindir,"\\ccRCC-TIME-clinical-2022-09-27.csv"), DataFrame, header=1, types=coltypes1);
	
	# show(ADResults)
	# show(clinical_raw)

	joinkey(i) = (i.patient)

	#println(joinkey(KSResults))
	#println(joinkey(clinical_raw))

	## !!!!!!!!!!!!!!!!!!!!!
	## double check that all clinical variables are present 
	@time DFFinal = @from i in JSResults begin
		@join j in clinical_raw on joinkey(i) equals joinkey(j)
		@select {i.patient, 
				i.name, 
				i.mean, 
				i.median, 
				i.std, i.var, 
				i.kurtosis, i.skewness, i.ks_result, i.ks_p, 
				i.ad_result, i.ad_p, 
				EGFRnormReads=j.EGFRnormReads, 
				EGFRalphaReads=j.EGFRalphaReads,
				EGFRbetaReads=j.EGFRbetaReads,
				EGFRgammaReads=j.EGFRgammaReads,
				Gender=j.Gender, Race=j.Race, 
				AgeDiagnosis=j.AgeDiagnosis, AgeSurgery=j.AgeSurgery, 
				Histology=j.Histology, Laterality=j.Laterality, 
				Grade=j.Grade, Size=j.Size, 
				SarcomatoidStatus=j.SarcomatoidStatus, 
				RhabdoidStatus=j.RhabdoidStatus, 
				pT=j.pT, pN=j.pN, pM=j.pM, 
				CytoreductiveSurgStatus=j.CytoreductiveSurg, 
				RFS=j.RFS, OS=j.OverallSurvival, CauseOfDeath=j.CauseOfDeath, 
				DeathStatus=j.EventDeath
		}
		@collect DataFrame;
	end

	directory3 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results-from-Pipeline\\_CutOff5\\interdist"
	# write merged statistics and clinical data for a single patient to a file 
	@time CSV.write(string(directory3, "\\", splitname,".interdistances_stats+clin",Dates.today(),".csv"), DFFinal)
	println(string(splitname, " interdistance statistics and clinical data has been saved to CSV"))

end

NamesOfInterdistP2=["cd68/cd20", "nucleus/pd-l1", "cd20/stroma", "cd68/stroma", "nucleus/tumor", "cd163/nucleus", "cd20/nucleus", "cd206/stroma", "cd206/tumor", "cd206/cd163", "pd-l1/tumor", "cd206/cd20", "pd-l1/stroma", "cd163/pd-l1", 
"cd163/stroma", "cd68/nucleus", "cd20/tumor", "cd206/pd-l1", "nucleus/stroma", "cd20/pd-l1", "cd163/tumor", "cd206/cd68", "cd163/cd20", "tumor/stroma", "cd68/pd-l1", "cd206/nucleus", "cd68/cd163", "cd68/tumor"]

# file1 = "P1a_SP04-4722_[37094,12061].tif_94266_job58256.object_results.csv"

## Panel 2 Markers
markerPanel = ["tumor","stroma","CD68","CD163","CD206","PD-L1"]

directory1 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor"
## Reads all files in the directory
##    note, this is looking at the interface location
ReadingFiles = readdir(directory1)

for i in 1:2#size(ReadingFiles)[1]-4
	filenametemp = ReadingFiles[i]
	TumorTIMEPipeline(directory1, filenametemp, markerPanel, "Panel-2", "tumor")
end

function ConcatFiles_Inter(directory)

	#Read files in results from pipeline folder  
	ReadingFiles2 = readdir(directory)
	println(ReadingFiles2)

 	coltypesL = Any[String for i=1:35]
 	coltypesL[1]=String
 	coltypesL[2]=String
 	coltypesL[14]=Union{Missing, String}
 	coltypesL[15]=Union{Missing, String}
 	coltypesL[16]=Union{Missing, String}
 	coltypesL[17]=Union{Missing, String}
 	coltypesL[18]=String
 	coltypesL[19]=String
 	coltypesL[20]=String
 	coltypesL[21]=String
 	coltypesL[22]=String
 	coltypesL[23]=String
 	coltypesL[24]= Union{Missing, Float64}
 	coltypesL[25]=Union{Missing, Float64}
 	coltypesL[26]=String
 	coltypesL[27]=String
 	coltypesL[28]=String
 	coltypesL[29]=Union{Missing, String}
 	coltypesL[30]=String
 	coltypesL[31]=String
 	coltypesL[32]=Union{Missing, String}
 	coltypesL[34]=Union{Missing, String}
 	longDF = CSV.read(string(directory,"\\",ReadingFiles2[1]), DataFrame, header=1, types=coltypesL)#Dict(:col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int, :col31=>Union{Missing,Int,String}))
 	# add a column to individual patient and loc DF that has the filename to be able to pull out location later if needed
 	longDF[:,:filename] .= ReadingFiles2[1] 

	#println(names(longDF)[31])
	println(longDF)

 	#Create for loop that runs through files of stats and clincial data 
 	for i in 1:size(ReadingFiles2)[1]-1
 		filenametemp2 = ReadingFiles2[i]
 		println(filenametemp2)

 		#tmpdf = CSV.read(string(directory,"\\",ReadingFiles2[i]), DataFrame, header=1, types=Dict(:col31=>Union{Missing,Int,String}, :col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int))
 		coltypes = Any[String for i=1:35]
 		coltypes[1]=String
 		coltypes[2]=String
 		coltypes[14]=Union{Missing, String}
 		coltypes[15]=Union{Missing, String}
 		coltypes[16]=Union{Missing, String}
 		coltypes[17]=Union{Missing, String}
 		coltypes[18]=String
 		coltypes[19]=String
 		coltypesL[20]=String
 		coltypesL[21]=String
 		coltypes[22]=String
 		coltypes[23]=String
 		coltypes[24]=Union{Missing, Float64}
 		coltypes[25]=Union{Missing, Float64}
 		coltypes[26]=String
 		coltypes[27]=String
 		coltypes[28]=String
 		coltypes[29]=Union{Missing, String}
 		coltypes[30]=String
 		coltypes[31]=String
 		coltypes[32]=Union{Missing, String}
 		coltypes[34]=Union{Missing, String}
 		#tmpdf = CSV.read(string(directory,"\\",filenametemp2), DataFrame, header=1, types=Dict(:col31=>Union{Missing,String}, :col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int))
 		tmpdf = CSV.read(string(directory,"\\",filenametemp2), DataFrame, header=1, types=coltypes)
 		tmpdf[:,:filename] .= ReadingFiles2[i] 

		append!(longDF,tmpdf)
	end
	directory4 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\"
	# write out CSV of full dataframe --- longDF
	@time CSV.write(string(directory4,"\\All-Clinical-Data-for-Interdistances-at-Panel2-Tumor",Dates.today(),".csv"), longDF)

	# Query to identify individual pairs and then save those as files 
	#Loop over each of the different names
	listNames = unique(longDF[:,:name])

	directory5 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\interdist\\"
 	for n in 1:length(listNames)
 		Query1= @from i in longDF begin
 			@where i.name == listNames[n]
 			@select i #{pair=i.name, i.filename}
 			@collect DataFrame
 		end
 		replace(listNames[n], "/" => "+") 
 		@time CSV.write(string(directory5,"\\","Query for (interdist)", replace(listNames[n], "/" => "+"), Dates.today(),".csv"), Query1)
 	end
  end

  #individdir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results from Pipeline\\_CutOff5\\interdist\\"
 # #individdir2="C:\\Users\\camara.casson\\OneDrive - University of Florida\\Desktop\\TumorTIME\\src"
  ConcatFiles1(directory3)


 function TumorTIMEPipeline_Intra(directory1, file, marker, panelName, panelLoc)

	# # directory1 is a directory to pull files from 
	## directory1 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor"

    # # file wants a string name of the file that has the single-cell measurements in it
	#file = "P2a_SP04-4722_[41516,9786].tif_95219_job59171.object_results.csv"

	#marker wants a vector of strings that are the markers in the panel
	#marker = ["tumor","stroma","CD68","CD163","CD206","PD-L1"]

	#panelName wants a string that identifies which panel subfolder ("Panel-1" or "Panel-2") to use
		## panelName = "Panel-2"
	#panelLoc wants a string that identifies which location folder (tumor, stroma, interface, normal) to use
		## panleLoc = "tumor"

	# identifies the patient name based on the SP number
	println(file)
	indName = string(split(file,"_")[2])
	

	# ## identifies the patient name with the location on the slide
	splitname = string(split(file,".")[1])


	df, interdist, intradist, inputfile = main((inputfile=file,classifierLabels=marker,numPairCutoff=5,savePlots=false,saveData=false))

	stats = [:mean,:median,:std,:var,:kurtosis,:skewness]
	interdist_stats = Dict(k => NamedTuple(stat => eval(stat)(v) for stat in stats) for (k,v) in interdist if !isempty(v))
	intradist_stats = Dict(k => NamedTuple(stat => eval(stat)(v) for stat in stats) for (k,v) in intradist if !isempty(v))

	intradist_stats_df = DataFrame(name=[], mean=[], median=[], std=[], var=[], kurtosis=[], skewness=[])
	for (k,v) in intradist_stats
		append!(intradist_stats_df, hcat(DataFrame(name = k), DataFrame([v])))
	end

	#Plotting data with StatsPlots to visualize distributions
	# histogram(interdist["cd68/stroma"])
	# dist = A["cd68/stroma"]._bin_pdf
	# StatsPlots.scatter(dist, leg=false)
	# StatsPlots.bar!(dist, func=cdf, alpha=0.3)

	function MakeDistributions(data)
		#Creates distributions fitting with a maximum likelihood estimator to a DiscreteNonParametric distribution
		#Creates a dictionary with distributions for each marker pair (non-empty)
		result = Dict()
		for (k,v) in data 
			if !isempty(v)
				Dist = fit_mle(DiscreteNonParametric,v); 
				result[k] = Dist
			end
		end
		return result
	end
    ## creates distribution based on interdistances calculated from slide information
    intradist_distr = MakeDistributions(intradist)
    println(intradist_distr)
	#Plotting the PDF for empirical and theoretical distribution
	#StatsPlots.scatter(hcat(TheoreticalPDFs["cd68/tumor"]...)[1,:],hcat(TheoreticalPDFs["cd68/tumor"]...)[2,:])
	#StatsPlots.bar!(z,alpha=0.3)

	#OLD, Don't need this
	# interdist_distr2 = Dict{String,Vector{Float64}}()
    # for i in keys(interdist_distr)
	# 	interdist_distr[i]._bin_probmass
	# 	interdist_distr2[i] = interdist_distr[i]._bin_probmass
	# end

    ##Creates distribution based on theoretical pdf equation
    ## add in theoretical PDF calculation here with range related to min and max of interdist in the sample
    
   #This prints out all the keys individually 
    # for i in keys(interdist)
    #     println(i)
    # end
   
    #creates a blank dictionary for the ranges of the interdistances at random distances
	IntradistRanges = Dict{String,Vector{Float64}}()
    #creates a blank dictionary for the distribution of the ranges of the interdistances
    TheoreticalPDFs = Dict{String,AbstractVector}()
    for i in keys(intradist_distr)
        #if !isempty(interdist[i])
            TheoreticalPDF = [];
			R = range(minimum(intradist[i]),maximum(intradist[i]),length=length(intradist[i]))
			for s in R
                # println("entered for loop")
			    push!(TheoreticalPDF,[s,nullpdf(s)])
			end
            TheoreticalPDFs[i]=TheoreticalPDF
	end
	#Plotting the nullpdf
	#StatsPlots.scatter(hcat(TheoreticalPDFs["cd68/stroma"]...)[1,:],hcat(TheoreticalPDFs["cd68/stroma"]...)[2,:])
    directory2 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results-from-Pipeline"

    # #write merged theoretical ranges to a file 
 	# @time CSV.write(string(directory2, "\\", splitname,".interdistances_ranges",Dates.today(),".csv"), InterdistRanges)
 	# println(string(splitname, " interdistance theoretical ranges have been saved to CSV"))

    #write merged theoretical PDFs to a file 
    @time CSV.write(string(directory2, "\\", splitname,".intradistances_PDFs",Dates.today(),".csv"), TheoreticalPDFs, bufsize=128*1024*1024)
    println(string(splitname, " intradistance theoretical distributions have been saved to CSV"))

	#calculate CDF from nullpdf using integration
	TheoreticalCDFs = Dict{String,AbstractVector}()
    for i in keys(intradist_distr)
        #if !isempty(interdist[i])
            TheoreticalCDF = [];
			R = range(minimum(intradist[i]),maximum(intradist[i]),length=length(intradist[i]))
			for s in R
                # println("entered for loop")
			    push!(TheoreticalCDF,[s,my_cdf(s)])
			end
            TheoreticalCDFs[i]=TheoreticalCDF
	end
#StatsPlots.scatter(hcat(TheoreticalCDFs["cd68/stroma"]...)[1,:],hcat(TheoreticalCDFs["cd68/stroma"]...)[2,:])

	#Calculate CDF of empirical distribution
	function MakeCDFs(data)
		#Creates a CDF distribution for empirical data
		#Creates a dictionary with distributions for each marker pair (non-empty)
		result = Dict()
		for (k,v) in data 
			if !isempty(v)
				Dist = ecdf(v) 
				result[k] = Dist
			end
		end
		return result
	end
    ## creates distribution based on interdistances calculated from slide information
    intradist_CDF = MakeCDFs(intradist)
    
	
	EmpiricalCDFs = Dict{String,AbstractVector}()
    for i in keys(intradist_distr)
            EmpiricalCDF = [];
			Intradist_Sort = sort(intradist[i])
			#println("Show Interdist_Sort")
			for s in 1:length(intradist_CDF[i].sorted_values)
                # println("entered for loop")
			    push!(EmpiricalCDF,[Intradist_Sort[s],intradist_CDF[i](Intradist_Sort[s])])
			end
            EmpiricalCDFs[i]=EmpiricalCDF
	end
	#StatsPlots.scatter!(hcat(EmpiricalCDFs["cd68/stroma"]...)[1,:],hcat(EmpiricalCDFs["cd68/stroma"]...)[2,:])
	
	function MakeDNPTheo(dataNull)
		DNPTheos = Dict{String,DiscreteNonParametric}()
		for i in keys(dataNull)
			dist1 = fit_mle(DiscreteNonParametric,hcat(dataNull[i]...))
			DNPTheos[i] = dist1
		end
		return DNPTheos
	end
	C = MakeDNPTheo(TheoreticalCDFs)


	function RunKSTest(dataEmp, dataNull, stats)
        ## dataEmp is the empirical data with the made distributions found in EmpiricalCDFs
        ## dataNull the theoretical distribution found inside TheoreticalCDFs
        ## stats is the dataframe with distribution statistcs that gets appended 

		kstests_col = DataFrame(ks_result = [], ks_p = [])
		for i in keys(dataNull) 
				A = hcat(dataEmp[i]...)[2,:]
				# x1 = vcat([s[1] for i in dataEmp])
				# x2 = vcat([s[1] for i in dataNull])
				# y1 = vcat([s[2] for i in dataEmp])
				# y2 = vcat([s[2] for i in dataNull])
				res = ExactOneSampleKSTest(A,C[i])#(dataEmp[i],dataNull[i],Any,Any) #(dataEmp[i][s][1], dataEmp[i][s][2],dataNull[i][s][1], dataNull[i][s][2]) #[x1,y1],[x2,y2]
				append!(kstests_col, DataFrame(ks_result=res.δ, ks_p = pvalue(res)))
		end

		stats=hcat(stats,kstests_col)
		return stats
	end

	KSResults = RunKSTest(EmpiricalCDFs,TheoreticalCDFs,intradist_stats_df)
	display(KSResults)

	# function RunADTest(dataEmp, dataNull, stats)
	# 	adtests_col = DataFrame(ad_result = [], ad_p = [])
	# 	for i in keys(dataNull)
	# 		A = hcat(dataEmp[i]...)[2,:]
	# 		res1 = OneSampleADTest(A,C[i])
	# 		append!(adtests_col, DataFrame(ad_result = res1.A², ad_p=(pvalue(res1))))
	# 	end
		
	# 	stats=hcat(stats,adtests_col)
	# 	return stats
	# end

	#  ADResults=RunADTest(EmpiricalCDFs,TheoreticalCDFs,KSResults)

	# function RunCramerVonMisesTest(dataEmp, dataNull, stats)
	# 	CVMTest_col = DataFrame(CVM_result = [], CVM_p=[])
	#	for i in keys(dataNull) 	
	#		dist3 = fit(DiscreteNonParametric,dataNull[i])
	# 		for j in 
	#	CVM_Results =  

	println("Starting KL Divergence")
	function KullbackLeibler(dataEmp,dataNull,stats)
		kld_col = DataFrame(kld_result= [])
		for i in keys(dataNull)
			A = hcat(dataEmp[i]...)[2,:]
			B = hcat(dataNull[i]...)[2,:]
			kld_result = kl_divergence(A,B)
			append!(kld_col, DataFrame(kld_result=kld_result))
		end
 
		stats=hcat(stats,kld_col)
		return stats
	end
	KLDResults = KullbackLeibler(EmpiricalCDFs,TheoreticalCDFs,KSResults)
	println("Finished KL Divergence")
	
	println("Starting Chebyshev Distance")
	function Chebyshev(dataEmp,dataNull,stats)
		cbs_col = DataFrame(cbs_result= [])
		for i in keys(dataNull)
			A = hcat(dataEmp[i]...)[1,:]
			B = hcat(dataNull[i]...)[1,:]
			cbs_result = chebyshev(A,B)
			append!(cbs_col, DataFrame(cbs_result=cbs_result))
		end
 
		stats=hcat(stats,cbs_col)
		return stats
	end
	CBSResults = Chebyshev(EmpiricalCDFs,TheoreticalCDFs,KLDResults)
	println("Finished Chebyshev Distance")


	println("Starting Jensen-Shannon Divergence")
	function JensenShannon(dataEmp,dataNull,stats)
		js_col = DataFrame(js_result= [])
		for i in keys(dataNull)
			A = hcat(dataEmp[i]...)[2,:]
			B = hcat(dataNull[i]...)[2,:]
			js_result = js_divergence(A,B)
			append!(js_col, DataFrame(js_result=js_result))
		end
 
		stats=hcat(stats,js_col)
		return stats
	end
	JSResults = JensenShannon(EmpiricalCDFs,TheoreticalCDFs,CBSResults)
	println("Finished Jensen-Shannon Divergence")


 	# add a column to statistics DF that has patient name 
 	JSResults[:,:patient] .= indName 

  	# write interdistance statistcs to a CSV file 
	  statsdir="C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\stats\\"
	  @time CSV.write(string(statsdir, "\\", splitname,".intradistances_stats",Dates.today(),".csv"), JSResults)

	coltypesS = Any[Float64 for i=1:35]
 	coltypesS[1]=String
 	coltypesS[2]=String

	
	#@time JSResults1 = CSV.read(string(statsdir,"All-Clinical-Data-for-Intradistances-at-Panel1-Tumor2022-10-27.csv"),DataFrame,types=coltypesS)
	JSResults1 = JSResults[:,1:12]

	println(JSResults1[1:5,:])

	# loading clinical information
	println("INFO: Time to load clinical data ")
	rawdatadirbase = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data"
	rawdatadir = string(rawdatadirbase,"\\", panelName, "\\", panelLoc)
	clindir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results-from-Pipeline\\_CutOff5\\Stats" 
	coltypes1 = Any[String for i=1:30]
	coltypes1[4]=Union{Missing, String}   ## EGFR norm reads
	coltypes1[5]=Union{Missing, String}   ## EGFR alpha reads
	coltypes1[6]=Union{Missing, String}   ## EGFR beta reads
	coltypes1[7]=Union{Missing, String}   ## EGFR gamma reads
	coltypes1[10]=Float64    ## age at diagnosis
	coltypes1[11]=Float64    ## age at surgery
	coltypes1[14]=Union{Missing, Int}        ## grade
	coltypes1[15]=Float64    ## size
	coltypes1[22]=Union{Missing, String}   ## RFS
	coltypes1[24]=Union{Missing, String}   ## CauseOfDeath
	coltypes1[25]=Int        ## event of death
	coltypes1[26]=Int 		 ## overall survival
	coltypes1[27]=Union{Missing, Float64}	   ## OS-IT
	coltypes1[28]=Union{Missing, Float64}	   ## OS-TT
	coltypes1[29]=Union{Missing, Float64}	   ## IMDC
	coltypes1[30]=Union{Missing, String}   ## sutent neoadvant
	@time clinical_raw = CSV.read(string(rawdatadirbase,"\\ccRCC-TIME-clinical-2022-09-27.csv"), DataFrame, header=1, types=coltypes1);
	
	# show(JSResults)
	# show(clinical_raw)

	joinkey(i) = (i.patient)

	#println(joinkey(KSResults))
	#println(joinkey(clinical_raw))

	## !!!!!!!!!!!!!!!!!!!!!
	## double check that all clinical variables are present 
	@time DFFinal = @from i in JSResults begin
		@join j in clinical_raw on joinkey(i) equals joinkey(j)
		@select {i.patient, 
				i.name, 
				i.mean, 
				i.median, 
				i.std, i.var, 
				i.kurtosis, i.skewness, i.ks_result, i.ks_p, 
				# i.ad_result, i.ad_p, 
				i.kld_result, i.cbs_result, i.js_result,
				EGFRnormReads=j.EGFRnormReads, 
				EGFRalphaReads=j.EGFRalphaReads,
				EGFRbetaReads=j.EGFRbetaReads,
				EGFRgammaReads=j.EGFRgammaReads,
				Gender=j.Gender, Race=j.Race, 
				AgeDiagnosis=j.AgeDiagnosis, AgeSurgery=j.AgeSurgery, 
				Histology=j.Histology, Laterality=j.Laterality, 
				Grade=j.Grade, Size=j.Size, 
				SarcomatoidStatus=j.SarcomatoidStatus, 
				RhabdoidStatus=j.RhabdoidStatus, 
				pT=j.pT, pN=j.pN, pM=j.pM, 
				CytoreductiveSurgStatus=j.CytoreductiveSurg, 
				RFS=j.RFS, OS=j.OverallSurvival, CauseOfDeath=j.CauseOfDeath, 
				DeathStatus=j.EventDeath
		}
		@collect DataFrame;
	end

	directory6 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results-from-Pipeline\\_Cutoff5\\intradist"
	# write merged statistics and clinical data for a single patient to a file 
	@time CSV.write(string(directory6, "\\", splitname,".intradistances_stats+clin",Dates.today(),".csv"), DFFinal)
	println(string(splitname, " intradistance statistics and clinical data has been saved to CSV"))

end

# NamesOfInterdistP2=["cd68/cd20", "nucleus/pd-l1", "cd20/stroma", "cd68/stroma", "nucleus/tumor", "cd163/nucleus", "cd20/nucleus", "cd206/stroma", "cd206/tumor", "cd206/cd163", "pd-l1/tumor", "cd206/cd20", "pd-l1/stroma", "cd163/pd-l1", 
# "cd163/stroma", "cd68/nucleus", "cd20/tumor", "cd206/pd-l1", "nucleus/stroma", "cd20/pd-l1", "cd163/tumor", "cd206/cd68", "cd163/cd20", "tumor/stroma", "cd68/pd-l1", "cd206/nucleus", "cd68/cd163", "cd68/tumor"]

# # file1 = "P1a_SP04-4722_[37094,12061].tif_94266_job58256.object_results.csv"

## Panel 2 Markers
markerPanel = ["tumor","stroma","CD68","CD163","CD206","PD-L1"]

directory1 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor"
## Reads all files in the directory
##    note, this is looking at the interface location
ReadingFiles = readdir(directory1)

for i in 1:size(ReadingFiles)[1]-4
	filenametemp = ReadingFiles[i]
	TumorTIMEPipeline_Intra(directory1, filenametemp, markerPanel, "Panel-2", "tumor")
end

function ConcatFiles_Intra(directory)

	#Read files in results from pipeline folder  
	ReadingFiles2 = readdir(directory)
	println(ReadingFiles2)

	coltypesL = Any[String for i=1:35]
	coltypesL[1]=String
	coltypesL[2]=String
	coltypesL[14]=Union{Missing, String}
	coltypesL[15]=Union{Missing, String}
	coltypesL[16]=Union{Missing, String}
	coltypesL[17]=Union{Missing, String}
	coltypesL[18]=String
	coltypesL[19]=String
	coltypesL[20]=String
	coltypesL[21]=String
	coltypesL[22]=String
	coltypesL[23]=String
	coltypesL[24]= Union{Missing, Float64}
	coltypesL[25]=Union{Missing, Float64}
	coltypesL[26]=String
	coltypesL[27]=String
	coltypesL[28]=String
	coltypesL[29]=Union{Missing, String}
	coltypesL[30]=String
	coltypesL[31]=String
	coltypesL[32]=Union{Missing, String}
	coltypesL[34]=Union{Missing, String}
	longDF = CSV.read(string(directory,"\\",ReadingFiles2[1]), DataFrame, header=1, types=coltypesL)#Dict(:col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int, :col31=>Union{Missing,Int,String}))
	# add a column to individual patient and loc DF that has the filename to be able to pull out location later if needed
	longDF[:,:filename] .= ReadingFiles2[1] 

	#println(names(longDF)[31])
	println(longDF)

	#Create for loop that runs through files of stats and clincial data 
	for i in 1:size(ReadingFiles2)[1]-1
		filenametemp2 = ReadingFiles2[i]
		println(filenametemp2)

		#tmpdf = CSV.read(string(directory,"\\",ReadingFiles2[i]), DataFrame, header=1, types=Dict(:col31=>Union{Missing,Int,String}, :col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int))
		coltypes = Any[String for i=1:35]
		coltypes[1]=String
		coltypes[2]=String
		coltypes[14]=Union{Missing, String}
		coltypes[15]=Union{Missing, String}
		coltypes[16]=Union{Missing, String}
		coltypes[17]=Union{Missing, String}
		coltypes[18]=String
		coltypes[19]=String
		coltypesL[20]=String
		coltypesL[21]=String
		coltypes[22]=String
		coltypes[23]=String
		coltypes[24]=Union{Missing, Float64}
		coltypes[25]=Union{Missing, Float64}
		coltypes[26]=String
		coltypes[27]=String
		coltypes[28]=String
		coltypes[29]=Union{Missing, String}
		coltypes[30]=String
		coltypes[31]=String
		coltypes[32]=Union{Missing, String}
		coltypes[34]=Union{Missing, String}
		#tmpdf = CSV.read(string(directory,"\\",filenametemp2), DataFrame, header=1, types=Dict(:col31=>Union{Missing,String}, :col13=>Int, :col14=>Int, :col15=>Int, :col16=>Int, :col23=>Int))
		tmpdf = CSV.read(string(directory,"\\",filenametemp2), DataFrame, header=1, types=coltypes)
		tmpdf[:,:filename] .= ReadingFiles2[i] 

		append!(longDF,tmpdf)
	end
	directory7 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\intradist\\"
	# write out CSV of full dataframe --- longDF
	@time CSV.write(string(directory4,"\\All-Clinical-Data-for-Intradistances-at-Panel2-Tumor",Dates.today(),".csv"), longDF)

	# Query to identify individual pairs and then save those as files 
	#Loop over each of the different names
	listNames = unique(longDF[:,:name])

	for n in 1:length(listNames)
		Query1= @from i in longDF begin
			@where i.name == listNames[n]
			@select i #{pair=i.name, i.filename}
			@collect DataFrame
		end
		replace(listNames[n], "/" => "+") 
		@time CSV.write(string(directory7,"\\","Query for (intradist)", replace(listNames[n], "/" => "+"), Dates.today(),".csv"), Query1)
	end
 end

 #individdir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\data\\Panel-2\\tumor\\Results from Pipeline\\_Cutoff5\\intradist\\"
# #individdir2="C:\\Users\\camara.casson\\OneDrive - University of Florida\\Desktop\\TumorTIME\\src"
 ConcatFiles2(directory6)

function RunStatistics_Inter(directory)
	QueryDirectory = directory
	JuliaStatsDir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\interdist\\Julia-Stats"
	ReadingFiles3 = readdir(directory)
    println(ReadingFiles3)

	#1 - Unequal Variance Test and Mann Whitney of KS and Gender
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            # Convert the "Gender" column to lowercase
            query_file.Gender = lowercase.(query_file.Gender)
            
            selected_columns = select(query_file, :ks_result, :Gender)
            grouped_data = groupby(selected_columns, :Gender)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(ks_result = pvalue_test1, Gender = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and Gender on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(ks_result = pvalue_test2, Gender = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and Gender on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#2 - Unequal Variance Test and Mann Whitney of chebyshev and Gender
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            # Convert the "Gender" column to lowercase
            query_file.Gender = lowercase.(query_file.Gender)
            
            selected_columns = select(query_file, :cbs_result, :Gender)
            grouped_data = groupby(selected_columns, :Gender)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(cbs_result = pvalue_test1, Gender = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of Chebyshev Test and Gender on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(cbs_result = pvalue_test2, Gender = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and Gender on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

		#3 - Unequal Variance Test and Mann Whitney of KS and Race
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
				
				# Convert the "Gender" column to lowercase
				query_file.Race = lowercase.(query_file.Race)
				
				selected_columns = select(query_file, :ks_result, :Race)
				grouped_data = groupby(selected_columns, :Race)
				grouped_grouped_data = vcat(grouped_data[2], grouped_data[3], grouped_data[4])
				grouped_data_DataFrame = DataFrame(grouped_data)
				grouped_grouped_data_DataFrame = DataFrame(grouped_grouped_data)
				#println(grouped_data_DataFrame)
				#println(grouped_grouped_data)
	
				result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_grouped_data[!,1])
				pvalue_test1 = pvalue(result1)
				#println(pvalue_test1)
				n1=length(grouped_data[1][!,1])
				n2 = length(grouped_grouped_data[!,1])
				DoF = n1+n2-2
				#println(DoF)
				test_info = DataFrame(ks_result = pvalue_test1, Race = DoF)
				Test_Info = vcat(grouped_data_DataFrame, test_info)
				println(Test_Info)
				@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and Race on ", file, Dates.today(),".csv"), Test_Info)
				
				result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_grouped_data[!,1])
				pvalue_test2 = pvalue(result2)
				#println(pvalue_test2)
				testinfo = DataFrame(ks_result = pvalue_test2, Race = ~)
				TestInfo = vcat(grouped_data_DataFrame, testinfo)
				println(TestInfo)
				@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and Race on ", file, Dates.today(),".csv"), TestInfo)
			end
	end

	#4 - Unequal Variance Test and Mann Whitney of chebyshev and Race
	for file in ReadingFiles3
		if occursin(".csv", file)
			filepath = joinpath(directory, file)
			query_file = CSV.File(filepath) |> DataFrame
			column_names = strip.(names(query_file))
			println(column_names)
			
			# Convert the "Gender" column to lowercase
			query_file.Race = lowercase.(query_file.Race)
			
			selected_columns = select(query_file, :cbs_result, :Race)
			grouped_data = groupby(selected_columns, :Race)
			grouped_grouped_data = vcat(grouped_data[2], grouped_data[3], grouped_data[4])
			grouped_data_DataFrame = DataFrame(grouped_data)
			grouped_grouped_data_DataFrame = DataFrame(grouped_grouped_data)
			#println(grouped_data_DataFrame)
			#println(grouped_grouped_data)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_grouped_data[!,1])
			pvalue_test1 = pvalue(result1)
			#println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_grouped_data[!,1])
			DoF = n1+n2-2
			#println(DoF)
			test_info = DataFrame(cbs_result = pvalue_test1, Race = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			println(Test_Info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of Chebyshev Test and Race on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_grouped_data[!,1])
			pvalue_test2 = pvalue(result2)
			#println(pvalue_test2)
			testinfo = DataFrame(cbs_result = pvalue_test2, Race = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and Race on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#5 - Pearson and Spearman correlations of KS and age at diagnosis
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            selected_columns = select(query_file, :ks_result, :AgeDiagnosis)
			println(selected_columns)
            # grouped_data = groupby(selected_columns, :AgeDiagnosis)
			# grouped_data_DataFrame = DataFrame(grouped_data)
			# println(grouped_data_DataFrame)

			result1 = Statistics.cor(selected_columns[!,1],selected_columns[!,2])
			correlation1 = result1
			println(correlation1)
				if length(selected_columns[!,1]) == length(selected_columns[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(cor(selected_columns[!,1],selected_columns[!,2]))) * sqrt(length(selected_columns[!,1]) - 3))
					println(pvalue)
					test_info = DataFrame(ks_result = correlation1, AgeDiagnosis = pvalue)
					Test_Info = vcat(selected_columns, test_info)
					@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of KS Test and Age at Diagnosis on ", file, Dates.today(),".csv"), Test_Info)
				else
					error("x and y have different lengths")
				end
			
			result2 = StatsBase.corspearman(selected_columns[!,1],selected_columns[!,2])
			correlation2 = result2
			println(correlation2)
				if length(selected_columns[!,1]) == length(selected_columns[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(selected_columns[!,1],selected_columns[!,2]))) * sqrt(length(selected_columns[!,1]) - 3))
					println(pvalue)
					testinfo = DataFrame(ks_result = correlation2, AgeDiagnosis = pvalue)
					TestInfo = vcat(selected_columns, testinfo)
					@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of KS Test and Age at Diagnosis on ", file, Dates.today(),".csv"), TestInfo)
				else
					error("x and y have different lengths")
				end
		end
	end

	#6 - Pearson and Spearman correlations of chebyshev and age at diagnosis
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            selected_columns = select(query_file, :cbs_result, :AgeDiagnosis)
			println(selected_columns)
            # grouped_data = groupby(selected_columns, :AgeDiagnosis)
			# grouped_data_DataFrame = DataFrame(grouped_data)
			# println(grouped_data_DataFrame)

			result1 = Statistics.cor(selected_columns[!,1],selected_columns[!,2])
			correlation1 = result1
			println(correlation1)
				if length(selected_columns[!,1]) == length(selected_columns[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(cor(selected_columns[!,1],selected_columns[!,2]))) * sqrt(length(selected_columns[!,1]) - 3))
					println(pvalue)
					test_info = DataFrame(cbs_result = correlation1, AgeDiagnosis = pvalue)
					Test_Info = vcat(selected_columns, test_info)
					@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of Chebyshev Test and Age at Diagnosis on ", file, Dates.today(),".csv"), Test_Info)

				else
					error("x and y have different lengths")
				end
			
			result2 = StatsBase.corspearman(selected_columns[!,1],selected_columns[!,2])
			correlation2 = result2
			println(correlation2)
				if length(selected_columns[!,1]) == length(selected_columns[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(selected_columns[!,1],selected_columns[!,2]))) * sqrt(length(selected_columns[!,1]) - 3))
					println(pvalue)
					testinfo = DataFrame(cbs_result = correlation2, AgeDiagnosis = pvalue)
					TestInfo = vcat(selected_columns, testinfo)
					@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of Chebyshev Test and Age at Diagnosis on ", file, Dates.today(),".csv"), TestInfo)
				else
					error("x and y have different lengths")
				end
		end
	end

	#7 -Pearson and Spearman correlations of KS and size
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            selected_columns = select(query_file, :ks_result, :Size)
			println(selected_columns)
			cleaned_data = dropmissing(selected_columns)
            # grouped_data = groupby(selected_columns, :AgeDiagnosis)
			# grouped_data_DataFrame = DataFrame(grouped_data)
			# println(grouped_data_DataFrame)

			result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
			correlation1 = result1
			println(correlation1)
				if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
					println(pvalue)
					test_info = DataFrame(ks_result = correlation1, Size = pvalue)
					Test_Info = vcat(cleaned_data, test_info)
					@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of KS Test and Size on ", file, Dates.today(),".csv"), Test_Info)
				else
					error("x and y have different lengths")
				end
			
			result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
			correlation2 = result2
			println(correlation2)
				if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
					println(pvalue)
					testinfo = DataFrame(ks_result = correlation2, Size = pvalue)
					TestInfo = vcat(cleaned_data, testinfo)
					@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of KS Test and Size on ", file, Dates.today(),".csv"), TestInfo)
				else
					error("x and y have different lengths")
				end
		end
	end

	#8 - Pearson and Spearman correlations of chebyshev and size
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            selected_columns = select(query_file, :cbs_result, :Size)
			println(selected_columns)
			cleaned_data = dropmissing(selected_columns)
            # grouped_data = groupby(selected_columns, :AgeDiagnosis)
			# grouped_data_DataFrame = DataFrame(grouped_data)
			# println(grouped_data_DataFrame)

			result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
			correlation1 = result1
			println(correlation1)
				if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
					println(pvalue)
					test_info = DataFrame(cbs_result = correlation1, Size = pvalue)
					Test_Info = vcat(cleaned_data, test_info)
					@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of Chebyshev Test and Size on ", file, Dates.today(),".csv"), Test_Info)
				else
					error("x and y have different lengths")
				end
			
			
			result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
			correlation2 = result2
			println(correlation2)
				if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
					pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
					println(pvalue)
					testinfo = DataFrame(cbs_result = correlation2, Size = pvalue)
					TestInfo = vcat(cleaned_data, testinfo)
					@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of Chebyshev Test and Size on ", file, Dates.today(),".csv"), TestInfo)
				else
					error("x and y have different lengths")
				end
		end
	end

	# #9 - Unequal Variance Test and Mann Whitney of KS and Grade (pairwise)
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            selected_columns = select(query_file, :ks_result, :Grade)
			#println(selected_columns)
            grouped_data = groupby(selected_columns, :Grade)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)
			#println(grouped_data[1][!,1])
			
			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF1 = n1+n2-2
			println(DoF1)
			result2 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[3][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[3][!,1])
			DoF2 = n1+n2-2
			println(DoF2)
			result3 = UnequalVarianceTTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test3 = pvalue(result3)
			println(pvalue_test3)
			n1=length(grouped_data[2][!,1])
			n2 = length(grouped_data[3][!,1])
			DoF3 = n1+n2-2
			println(DoF2)
			test_info1 = DataFrame(ks_result = pvalue_test1, Grade = DoF1)
			Test_Info = vcat(grouped_data_DataFrame, test_info1)
			test_info2 = DataFrame(ks_result = pvalue_test2, Grade = DoF2)
			Test_Info2 = vcat(Test_Info, test_info2)
			test_info3 = DataFrame(ks_result = pvalue_test3, Grade = DoF3)
			Test_Info3 = vcat(Test_Info2,test_info3)
			println(Test_Info3)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and Grade on ", file, Dates.today(),".csv"), Test_Info3)
			
			result4 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test4 = pvalue(result4)
			println(pvalue_test4)
			result5 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[3][!,1])
			pvalue_test5 = pvalue(result5)
			println(pvalue_test5)
			result6 = MannWhitneyUTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test6 = pvalue(result6)
			println(pvalue_test6)
			testinfo4 = DataFrame(ks_result = pvalue_test4, Grade = ~)
			TestInfo4 = vcat(grouped_data_DataFrame, testinfo4)
			testinfo5 = DataFrame(ks_result = pvalue_test5, Grade = ~)
			TestInfo5 = vcat(TestInfo4, testinfo5)
			testinfo6 = DataFrame(ks_result = pvalue_test6, Grade = ~)
			TestInfo6 = vcat(TestInfo5, testinfo6)
			println(TestInfo6)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and Grade on ", file, Dates.today(),".csv"), TestInfo6)

	# 		 result1 = KruskalWallisTest(grouped_data[1][!,1],grouped_data[2][!,1],grouped_data[3][!,1])
	# 		 x=result1
	# 		 println(x)
	# 		 # Get the chi-square statistic
	# 		#chi_sq_stat = result1.statistc

	# 		# Get the degrees of freedom
	# 		#df = result1.df

	# 		#println(pvalue(result1))
	# 		y = pvalue(result1)
	# 		println(y)
	# 		PVALUE= DataFrame(ks_result = y, Grade = ~)
	# 		Test_Info = vcat(grouped_data_DataFrame, PVALUE)
	# 		#@time CSV.write(string(JuliaStatsDir,"\\","Kruskal Wallis of KS Test and Grade on ", file, Dates.today(),".csv"), Test_Info)

	# 		# Compute the p-value using the chi-square distribution
	# 		#p_value = ccdf(Chisq(df), chi_sq_stat)
	# 		 #pvalue = parse(Float64, match(r"one-sided p-value:\s+([\d.]+)", result1).captures[1])
	# 		 #pvalue = result1.pvalue
	# 		 #pvalue = ccdf(Chisq(length(grouped_data-1)), x)

	# 		 result2 = OneWayANOVATest(grouped_data[1][!,1],grouped_data[2][!,1],grouped_data[3][!,1])
	# 		 z = pvalue(result2)
	# 		 println(z)
	# 		 PVALUE2 = DataFrame(ks_result = z, Grade = ~)
	# 		 TestInfo = vcat(grouped_data_DataFrame, PVALUE2)
	# 		 #@time CSV.write(string(JuliaStatsDir,"\\","One Way Anova of KS Test and Grade on ", file, Dates.today(),".csv"), TestInfo)
	 	end
	 end

	# #10 - Unequal Variance Test and Mann Whitney of Chebyshev and Grade (pairwise)
	for file in ReadingFiles3
		if occursin(".csv", file)
			filepath = joinpath(directory, file)
			query_file = CSV.File(filepath) |> DataFrame
			column_names = strip.(names(query_file))
			println(column_names)
				
			selected_columns = select(query_file, :cbs_result, :Grade)
			#println(selected_columns)
			grouped_data = groupby(selected_columns, :Grade)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)
				
			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF1 = n1+n2-2
			println(DoF1)
			result2 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[3][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[3][!,1])
			DoF2 = n1+n2-2
			println(DoF2)
			result3 = UnequalVarianceTTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test3 = pvalue(result3)
			println(pvalue_test3)
			n1=length(grouped_data[2][!,1])
			n2 = length(grouped_data[3][!,1])
			DoF3 = n1+n2-2
			println(DoF2)
			test_info1 = DataFrame(cbs_result = pvalue_test1, Grade = DoF1)
			Test_Info = vcat(grouped_data_DataFrame, test_info1)
			test_info2 = DataFrame(cbs_result = pvalue_test2, Grade = DoF2)
			Test_Info2 = vcat(Test_Info, test_info2)
			test_info3 = DataFrame(cbs_result = pvalue_test3, Grade = DoF3)
			Test_Info3 = vcat(Test_Info2,test_info3)
			println(Test_Info3)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of Chebyshev Test and Grade on ", file, Dates.today(),".csv"), Test_Info3)
			
			result4 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test4 = pvalue(result4)
			println(pvalue_test4)
			result5 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[3][!,1])
			pvalue_test5 = pvalue(result5)
			println(pvalue_test5)
			result6 = MannWhitneyUTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test6 = pvalue(result6)
			println(pvalue_test6)
			testinfo4 = DataFrame(cbs_result = pvalue_test4, Grade = ~)
			TestInfo4 = vcat(grouped_data_DataFrame, testinfo4)
			testinfo5 = DataFrame(cbs_result = pvalue_test5, Grade = ~)
			TestInfo5 = vcat(TestInfo4, testinfo5)
			testinfo6 = DataFrame(cbs_result = pvalue_test6, Grade = ~)
			TestInfo6 = vcat(TestInfo5, testinfo6)
			println(TestInfo6)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and Grade on ", file, Dates.today(),".csv"), TestInfo6)
			# result1 = KruskalWallisTest(grouped_data[1][!,1],grouped_data[2][!,1],grouped_data[3][!,1])
			# x=result1
			# println(x)
			# # Get the chi-square statistic
			# #chi_sq_stat = result1.statistc
	
			# # Get the degrees of freedom
			# #df = result1.df
	
			# #println(pvalue(result1))
			# y = pvalue(result1)
			# println(y)
			# PVALUE= DataFrame(cbs_result = y, Grade = ~)
			# Test_Info = vcat(grouped_data_DataFrame, PVALUE)
			# @time CSV.write(string(JuliaStatsDir,"\\","Kruskal Wallis of Chebyshev Test and Grade on ", file, Dates.today(),".csv"), Test_Info)
	
			# # Compute the p-value using the chi-square distribution
			# #p_value = ccdf(Chisq(df), chi_sq_stat)
			# #pvalue = parse(Float64, match(r"one-sided p-value:\s+([\d.]+)", result1).captures[1])
			# #pvalue = result1.pvalue
			# #pvalue = ccdf(Chisq(length(grouped_data-1)), x)
	
			# result2 = OneWayANOVATest(grouped_data[1][!,1],grouped_data[2][!,1],grouped_data[3][!,1])
			# z = pvalue(result2)
			# println(z)
			# PVALUE2 = DataFrame(cbs_result = z, Grade = ~)
			# TestInfo = vcat(grouped_data_DataFrame, PVALUE2)
			# @time CSV.write(string(JuliaStatsDir,"\\","One Way Anova of Chebyshev Test and Grade on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#11 - Unequal Variance Test and Mann Whitney of KS and SarcomatoidStatus
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            
            selected_columns = select(query_file, :ks_result, :SarcomatoidStatus)
            grouped_data = groupby(selected_columns, :SarcomatoidStatus)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(ks_result = pvalue_test1, SarcomatoidStatus = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and SarcomatoidStatus on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(ks_result = pvalue_test2, SarcomatoidStatus = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and SarcomatoidStatus on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#12 - Unequal Variance Test and Mann Whitney of chebyshev and SarcomatoidStatus
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            
            selected_columns = select(query_file, :cbs_result, :SarcomatoidStatus)
            grouped_data = groupby(selected_columns, :SarcomatoidStatus)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(cbs_result = pvalue_test1, Gender = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of Chebyshev Test and SarcomatoidStatus on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(cbs_result = pvalue_test2, Gender = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and SarcomatoidStatus on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#13 - Unequal Variance Test and Mann Whitney of KS and RhabdoidStatus
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            
            selected_columns = select(query_file, :ks_result, :RhabdoidStatus)
            grouped_data = groupby(selected_columns, :RhabdoidStatus)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(ks_result = pvalue_test1, RhabdoidStatus = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and RhabdoidStatus on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(ks_result = pvalue_test2, RhabdoidStatus = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and RhabdoidStatus on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#14 - Unequal Variance Test and Mann Whitney of chebyshev and RhabdoidStatus
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            
            selected_columns = select(query_file, :cbs_result, :RhabdoidStatus)
            grouped_data = groupby(selected_columns, :RhabdoidStatus)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(cbs_result = pvalue_test1, RhabdoidStatus = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of Chebyshev Test and SarcomatoidStatus on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(cbs_result = pvalue_test2, RhabdoidStatus = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and RhabdoidStatus on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#15 - Unequal Variance Test and Mann Whitney of KS and Laterality
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            
            selected_columns = select(query_file, :ks_result, :Laterality)
            grouped_data = groupby(selected_columns, :Laterality)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(ks_result = pvalue_test1, Laterality = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and Laterality on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(ks_result = pvalue_test2, Laterality = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and Laterality on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

	#16 - Unequal Variance Test and Mann Whitney of chebyshev and Laterality
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            
            
            selected_columns = select(query_file, :cbs_result, :Laterality)
            grouped_data = groupby(selected_columns, :Laterality)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)

			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = pvalue(result1)
			println(pvalue_test1)
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF = n1+n2-2
			println(DoF)
			test_info = DataFrame(cbs_result = pvalue_test1, Laterality = DoF)
			Test_Info = vcat(grouped_data_DataFrame, test_info)
			@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of Chebyshev Test and Laterality on ", file, Dates.today(),".csv"), Test_Info)
			
			result2 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			testinfo = DataFrame(cbs_result = pvalue_test2, Laterality = ~)
			TestInfo = vcat(grouped_data_DataFrame, testinfo)
			println(TestInfo)
			@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and Laterality on ", file, Dates.today(),".csv"), TestInfo)
		end
	end

		#17 -Pearson and Spearman correlations of KS and EGFRnormReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
				
				selected_columns = select(query_file, :ks_result, :EGFRnormReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
				# grouped_data = groupby(selected_columns, :AgeDiagnosis)
				# grouped_data_DataFrame = DataFrame(grouped_data)
				# println(grouped_data_DataFrame)
	
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(ks_result = correlation1, EGFRnormReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of KS Test and EGFRnormReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
				
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(ks_result = correlation2, EGFRnormReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of KS Test and EGFRnormReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
			end
		end
	
		#18 - Pearson and Spearman correlations of chebyshev and EGFRnormReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
				
				selected_columns = select(query_file, :cbs_result, :EGFRnormReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
				# grouped_data = groupby(selected_columns, :EGFRnormReads)
				# grouped_data_DataFrame = DataFrame(grouped_data)
				# println(grouped_data_DataFrame)
	
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(cbs_result = correlation1, EGFRnormReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of Chebyshev Test and EGFRnormReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
				
				
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(cbs_result = correlation2, EGFRnormReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of Chebyshev Test and EGFRnormReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
			end
		end

		#19 -Pearson and Spearman correlations of KS and EGFRalphaReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
						
				selected_columns = select(query_file, :ks_result, :EGFRalphaReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
						# grouped_data = groupby(selected_columns, :EGFRalphaReads)
						# grouped_data_DataFrame = DataFrame(grouped_data)
						# println(grouped_data_DataFrame)
			
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(ks_result = correlation1, EGFRalphaReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of KS Test and EGFRalphaReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
						
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(ks_result = correlation2, EGFRalphaReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of KS Test and EGFRalphaReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
				end
		end
			
		#20 - Pearson and Spearman correlations of chebyshev and EGFRnormReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
						
				selected_columns = select(query_file, :cbs_result, :EGFRalphaReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
						# grouped_data = groupby(selected_columns, :EGFRalphaReads)
						# grouped_data_DataFrame = DataFrame(grouped_data)
						# println(grouped_data_DataFrame)
			
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(cbs_result = correlation1, EGFRalphaReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of Chebyshev Test and EGFRalphaReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
						
						
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(cbs_result = correlation2, EGFRalphaReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of Chebyshev Test and EGFRalphaReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
			end
		end

		#21 -Pearson and Spearman correlations of KS and EGFRbetaReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
						
				selected_columns = select(query_file, :ks_result, :EGFRbetaReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
						# grouped_data = groupby(selected_columns, :EGFRbetaReads)
						# grouped_data_DataFrame = DataFrame(grouped_data)
						# println(grouped_data_DataFrame)
			
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(ks_result = correlation1, EGFRbetaReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of KS Test and EGFRbetaReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
						
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(ks_result = correlation2, EGFRbetaReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of KS Test and EGFRbetaReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
				end
		end
			
		#22 - Pearson and Spearman correlations of chebyshev and EGFRnormReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
						
				selected_columns = select(query_file, :cbs_result, :EGFRbetaReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
						# grouped_data = groupby(selected_columns, :EGFRbetaReads)
						# grouped_data_DataFrame = DataFrame(grouped_data)
						# println(grouped_data_DataFrame)
			
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(cbs_result = correlation1, EGFRbetaReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of Chebyshev Test and EGFRbetaReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
						
						
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(cbs_result = correlation2, EGFRbetaReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of Chebyshev Test and EGFRbetaReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
			end
		end

		#23 -Pearson and Spearman correlations of KS and EGFRgammaReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
						
				selected_columns = select(query_file, :ks_result, :EGFRgammaReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
						# grouped_data = groupby(selected_columns, :EGFRgammaReads)
						# grouped_data_DataFrame = DataFrame(grouped_data)
						# println(grouped_data_DataFrame)
			
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(ks_result = correlation1, EGFRgammaReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of KS Test and EGFRgammaReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
						
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(ks_result = correlation2, EGFRgammaReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of KS Test and EGFRgammaReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
				end
		end
			
		#24 - Pearson and Spearman correlations of chebyshev and EGFRgammaReads
		for file in ReadingFiles3
			if occursin(".csv", file)
				filepath = joinpath(directory, file)
				query_file = CSV.File(filepath) |> DataFrame
				column_names = strip.(names(query_file))
				println(column_names)
						
				selected_columns = select(query_file, :cbs_result, :EGFRgammaReads)
				println(selected_columns)
				cleaned_data = dropmissing(selected_columns)
						# grouped_data = groupby(selected_columns, :EGFRgammaReads)
						# grouped_data_DataFrame = DataFrame(grouped_data)
						# println(grouped_data_DataFrame)
			
				result1 = Statistics.cor(cleaned_data[!,1],cleaned_data[!,2])
				correlation1 = result1
				println(correlation1)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(cor(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						test_info = DataFrame(cbs_result = correlation1, EGFRgammaReads = pvalue)
						Test_Info = vcat(cleaned_data, test_info)
						@time CSV.write(string(JuliaStatsDir,"\\","Pearson Correlation of Chebyshev Test and EGFRgammaReads on ", file, Dates.today(),".csv"), Test_Info)
					else
						error("x and y have different lengths")
					end
						
						
				result2 = StatsBase.corspearman(cleaned_data[!,1],cleaned_data[!,2])
				correlation2 = result2
				println(correlation2)
					if length(cleaned_data[!,1]) == length(cleaned_data[!,2])
						pvalue = 2 * ccdf(Normal(), atanh(abs(corspearman(cleaned_data[!,1],cleaned_data[!,2]))) * sqrt(length(cleaned_data[!,1]) - 3))
						println(pvalue)
						testinfo = DataFrame(cbs_result = correlation2, EGFRgammaReads = pvalue)
						TestInfo = vcat(cleaned_data, testinfo)
						@time CSV.write(string(JuliaStatsDir,"\\","Spearmann Correlation of Chebyshev Test and EGFRgammaReads on ", file, Dates.today(),".csv"), TestInfo)
					else
						error("x and y have different lengths")
					end
			end
		end
function RunStatistics_Inter(directory)
		QueryDirectory = directory
		JuliaStatsDir = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\interdist\\Julia-Stats"
		ReadingFiles3 = readdir(directory)
		println(ReadingFiles3)

	# 25 - Unequal Variance Test and Mann Whitney of KS and pT (pairwise)
	for file in ReadingFiles3
		if occursin(".csv", file)
            filepath = joinpath(directory, file)
            query_file = CSV.File(filepath) |> DataFrame
            column_names = strip.(names(query_file))
            println(column_names)
            query_file.pT = replace.(query_file.pT, r"\s+$" => "")

            selected_columns = select(query_file, :ks_result, :pT)
			#println(selected_columns)
            grouped_data = groupby(selected_columns, :pT)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data)
			#println(grouped_data[1][!,1])
			
			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = NaN 
			try
				pvalue_test1 = pvalue(result1)
				println(pvalue_test1)
			catch err
				pvalue_test1 = NaN
				println("Error occurred; pvalue = NaN")
			end
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF1 = n1+n2-2
			println(DoF1)
			result2 = UnequalVarianceTTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test2 = NaN
			try
				pvalue_test2 = pvalue(result2)
				println(pvalue_test2)
			catch err
				pvalue_test2 = NaN
				println("Error occurred; pvalue = NaN")
			end
			n1=length(grouped_data[2][!,1])
			n2 = length(grouped_data[3][!,1])
			DoF2 = n1+n2-2
			println(DoF2)
			result3 = UnequalVarianceTTest(grouped_data[3][!,1],grouped_data[4][!,1])
			pvalue_test3 = NaN
			try
				pvalue_test3 = pvalue(result3)
				println(pvalue_test3)
			catch err 
				pvalue_test3 = NaN
				println("Error occurred; PValue = NaN")
			end
			n1=length(grouped_data[3][!,1])
			n2 = length(grouped_data[4][!,1])
			DoF3 = n1+n2-2
			println(DoF3)
			result4 = UnequalVarianceTTest(grouped_data[4][!,1],grouped_data[5][!,1])
			pvalue_test4 = NaN
			try
				pvalue_test4 = pvalue(result4)
				println(pvalue_test4)
			catch err
				pvalue_test4 = NaN
				println("Error occurred; PValue = NaN")
			end
			n1=length(grouped_data[4][!,1])
			n2 = length(grouped_data[5][!,1])
			DoF4 = n1+n2-2
			println(DoF4)
			result5 = UnequalVarianceTTest(grouped_data[5][!,1],grouped_data[6][!,1])
			pvalue_test5 = NaN
			try
				pvalue_test5 = pvalue(result5)
				println(pvalue_test5)
			catch err
				pvalue_test5 = NaN
				println("Error occurred; PValue = NaN")
			end
			n1=length(grouped_data[5][!,1])
			n2 = length(grouped_data[6][!,1])
			DoF5 = n1+n2-2
			println(DoF5)
			result6 = UnequalVarianceTTest(grouped_data[6][!,1],grouped_data[7][!,1])
			pvalue_test6 = NaN
			try
				pvalue_test6 = pvalue(result6)
				println(pvalue_test6)
			catch err
				pvalue_test6 = NaN
				println("Error occured; pvalue = NaN")
			end
			n1=length(grouped_data[6][!,1])
			n2 = length(grouped_data[7][!,1])
			DoF6 = n1+n2-2
			println(DoF6)
			pvalue_test7 = NaN
			DoF7 = NaN
			try
				result7 = UnequalVarianceTTest(grouped_data[7][!,1],grouped_data[8][!,1])
				pvalue_test7 = pvalue(result7)
				println(pvalue_test7)
				n1=length(grouped_data[7][!,1])
				n2 = length(grouped_data[8][!,1])
				DoF7 = n1+n2-2
				println(DoF7)
			catch err
				pvalue_test7 = NaN
				println("Error occurred; PValue = NaN")
			end
			test_info1 = DataFrame(ks_result = pvalue_test1, pT = DoF1)
			Test_Info = vcat(grouped_data_DataFrame, test_info1)
			test_info2 = DataFrame(ks_result = pvalue_test2, pT = DoF2)
			Test_Info2 = vcat(Test_Info, test_info2)
			test_info3 = DataFrame(ks_result = pvalue_test3, pT = DoF3)
			Test_Info3 = vcat(Test_Info2,test_info3)
			test_info4 = DataFrame(ks_result = pvalue_test4, pT = DoF4)
			Test_Info4 = vcat(Test_Info3,test_info4)
			test_info5 = DataFrame(ks_result = pvalue_test5, pT = DoF5)
			Test_Info5 = vcat(Test_Info4,test_info5)
			test_info6 = DataFrame(ks_result = pvalue_test6, pT = DoF6)
			Test_Info6 = vcat(Test_Info5,test_info6)
			test_info7 = DataFrame(ks_result = pvalue_test7, pT = DoF7)
			Test_Info7 = vcat(Test_Info6,test_info7)
			println(Test_Info7)
			#@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and pT on ", file, Dates.today(),".csv"), Test_Info7)
			
			result8 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test8 = pvalue(result8)
			println(pvalue_test8)
			result9 = MannWhitneyUTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test9 = pvalue(result9)
			println(pvalue_test9)
			result10 = MannWhitneyUTest(grouped_data[3][!,1],grouped_data[4][!,1])
			pvalue_test10 = pvalue(result10)
			println(pvalue_test10)
			result11 = MannWhitneyUTest(grouped_data[4][!,1],grouped_data[5][!,1])
			pvalue_test11 = pvalue(result11)
			println(pvalue_test11)
			result12 = MannWhitneyUTest(grouped_data[5][!,1],grouped_data[6][!,1])
			pvalue_test12 = pvalue(result12)
			println(pvalue_test12)
			result13 = MannWhitneyUTest(grouped_data[6][!,1],grouped_data[7][!,1])
			pvalue_test13 = pvalue(result13)
			println(pvalue_test13)
			pvalue_test14 = NaN
			try
				result14 = MannWhitneyUTest(grouped_data[7][!,1],grouped_data[8][!,1])
				pvalue_test14 = pvalue(result14)
				println(pvalue_test14)
			catch err
				pvalue_test14 = NaN
				println("Error occured; pvalue = NaN")
			end
			testinfo8 = DataFrame(ks_result = pvalue_test8, pT = ~)
			TestInfo8 = vcat(grouped_data_DataFrame, testinfo8)
			testinfo9 = DataFrame(ks_result = pvalue_test9, pT = ~)
			TestInfo9 = vcat(TestInfo8, testinfo9)
			testinfo10 = DataFrame(ks_result = pvalue_test10, pT = ~)
			TestInfo10 = vcat(TestInfo9, testinfo10)
			testinfo11 = DataFrame(ks_result = pvalue_test11, pT = ~)
			TestInfo11 = vcat(TestInfo10, testinfo11)
			testinfo12 = DataFrame(ks_result = pvalue_test12, pT = ~)
			TestInfo12 = vcat(TestInfo11, testinfo12)
			testinfo13 = DataFrame(ks_result = pvalue_test13, pT = ~)
			TestInfo13 = vcat(TestInfo12, testinfo13)
			testinfo14 = DataFrame(ks_result = pvalue_test14, pT = ~)
			TestInfo14 = vcat(TestInfo13, testinfo14)
			println(TestInfo14)
			#@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of KS Test and pT on ", file, Dates.today(),".csv"), TestInfo14)
	 	end
	 end

	# 26 - Unequal Variance Test and Mann Whitney of Chebyshev and pT (pairwise)
	for file in ReadingFiles3
		if occursin(".csv", file)
			filepath = joinpath(directory, file)
			query_file = CSV.File(filepath) |> DataFrame
			column_names = strip.(names(query_file))
			println(column_names)
				
			selected_columns = select(query_file, :cbs_result, :pT)
			#println(selected_columns)
			grouped_data = groupby(selected_columns, :pT)
			grouped_data_DataFrame = DataFrame(grouped_data)
			println(grouped_data_DataFrame)
				
			result1 = UnequalVarianceTTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test1 = NaN 
			try
				pvalue_test1 = pvalue(result1)
				println(pvalue_test1)
			catch err
				pvalue_test1 = NaN
				println("Error occurred; pvalue = NaN")
			end
			n1=length(grouped_data[1][!,1])
			n2 = length(grouped_data[2][!,1])
			DoF1 = n1+n2-2
			println(DoF1)
			result2 = UnequalVarianceTTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test2 = pvalue(result2)
			println(pvalue_test2)
			n1=length(grouped_data[2][!,1])
			n2 = length(grouped_data[3][!,1])
			DoF2 = n1+n2-2
			println(DoF2)
			result3 = UnequalVarianceTTest(grouped_data[3][!,1],grouped_data[4][!,1])
			pvalue_test3 = pvalue(result3)
			println(pvalue_test3)
			n1=length(grouped_data[3][!,1])
			n2 = length(grouped_data[4][!,1])
			DoF3 = n1+n2-2
			println(DoF3)
			result4 = UnequalVarianceTTest(grouped_data[4][!,1],grouped_data[5][!,1])
			pvalue_test4 = NaN
			try
				pvalue_test4 = pvalue(result4)
				println(pvalue_test4)
			catch err
				pvalue_test4 = NaN
				println("Error occurred; pvalue = NaN")
			end
			n1=length(grouped_data[4][!,1])
			n2 = length(grouped_data[5][!,1])
			DoF4 = n1+n2-2
			println(DoF4)
			result5 = UnequalVarianceTTest(grouped_data[5][!,1],grouped_data[6][!,1])
			pvalue_test5 = NaN
			try
				pvalue_test5 = pvalue(result5)
				println(pvalue_test5)
			catch err
				pvalue_test5 = NaN
				println("Error occurred; pvalue = NaN")
			end
			n1=length(grouped_data[5][!,1])
			n2 = length(grouped_data[6][!,1])
			DoF5 = n1+n2-2
			println(DoF5)
			result6 = UnequalVarianceTTest(grouped_data[6][!,1],grouped_data[7][!,1])
			pvalue_test6 = NaN
			try
				pvalue_test6 = pvalue(result6)
				println(pvalue_test6)
			catch err
				pvalue_test6 = NaN
				println("Error occured; pvalue = NaN")
			end
			n1=length(grouped_data[6][!,1])
			n2 = length(grouped_data[7][!,1])
			DoF6 = n1+n2-2
			println(DoF6)
			pvalue_test7 = NaN
			DoF7 = NaN
			try
				result7 = UnequalVarianceTTest(grouped_data[7][!,1],grouped_data[8][!,1])
				pvalue_test7 = pvalue(result7)
				println(pvalue_test7)
				n1=length(grouped_data[7][!,1])
				n2 = length(grouped_data[8][!,1])
				DoF7 = n1+n2-2
				println(DoF7)
			catch err
				pvalue_test7 = NaN
				println("Error occurred; PValue = NaN")
			end
			test_info1 = DataFrame(cbs_result = pvalue_test1, pT = DoF1)
			Test_Info = vcat(grouped_data_DataFrame, test_info1)
			test_info2 = DataFrame(cbs_result = pvalue_test2, pT = DoF2)
			Test_Info2 = vcat(Test_Info, test_info2)
			test_info3 = DataFrame(cbs_result = pvalue_test3, pT = DoF3)
			Test_Info3 = vcat(Test_Info2,test_info3)
			test_info4 = DataFrame(cbs_result = pvalue_test4, pT = DoF4)
			Test_Info4 = vcat(Test_Info3,test_info4)
			test_info5 = DataFrame(cbs_result = pvalue_test5, pT = DoF5)
			Test_Info5 = vcat(Test_Info4,test_info5)
			test_info6 = DataFrame(cbs_result = pvalue_test6, pT = DoF6)
			Test_Info6 = vcat(Test_Info5,test_info6)
			test_info7 = DataFrame(cbs_result = pvalue_test7, pT = DoF7)
			Test_Info7 = vcat(Test_Info6,test_info7)
			println(Test_Info7)
			#@time CSV.write(string(JuliaStatsDir,"\\","Unequal Variance of KS Test and pT on ", file, Dates.today(),".csv"), Test_Info7)
			
			result8 = MannWhitneyUTest(grouped_data[1][!,1],grouped_data[2][!,1])
			pvalue_test8 = pvalue(result8)
			println(pvalue_test8)
			result9 = MannWhitneyUTest(grouped_data[2][!,1],grouped_data[3][!,1])
			pvalue_test9 = pvalue(result9)
			println(pvalue_test9)
			result10 = MannWhitneyUTest(grouped_data[3][!,1],grouped_data[4][!,1])
			pvalue_test10 = pvalue(result10)
			println(pvalue_test10)
			result11 = MannWhitneyUTest(grouped_data[4][!,1],grouped_data[5][!,1])
			pvalue_test11 = pvalue(result11)
			println(pvalue_test11)
			result12 = MannWhitneyUTest(grouped_data[5][!,1],grouped_data[6][!,1])
			pvalue_test12 = pvalue(result12)
			println(pvalue_test12)
			result13 = MannWhitneyUTest(grouped_data[6][!,1],grouped_data[7][!,1])
			pvalue_test13 = pvalue(result13)
			println(pvalue_test13)
			result14 = MannWhitneyUTest(grouped_data[7][!,1],grouped_data[8][!,1])
			pvalue_test14 = pvalue(result14)
			println(pvalue_test14)
			testinfo8 = DataFrame(cbs_result = pvalue_test8, pT = ~)
			TestInfo8 = vcat(grouped_data_DataFrame, testinfo8)
			testinfo9 = DataFrame(cbs_result = pvalue_test9, pT = ~)
			TestInfo9 = vcat(TestInfo8, testinfo9)
			testinfo10 = DataFrame(cbs_result = pvalue_test10, pT = ~)
			TestInfo10 = vcat(TestInfo9, testinfo10)
			testinfo11 = DataFrame(cbs_result = pvalue_test11, pT = ~)
			TestInfo11 = vcat(TestInfo10, testinfo11)
			testinfo12 = DataFrame(cbs_result = pvalue_test12, pT = ~)
			TestInfo12 = vcat(TestInfo11, testinfo12)
			testinfo13 = DataFrame(cbs_result = pvalue_test13, pT = ~)
			TestInfo13 = vcat(TestInfo12, testinfo13)
			testinfo14 = DataFrame(cbs_result = pvalue_test14, pT = ~)
			TestInfo14 = vcat(TestInfo13, testinfo14)
			println(TestInfo14)
			#@time CSV.write(string(JuliaStatsDir,"\\","MannWhitney of Chebyshev Test and pT on ", file, Dates.today(),".csv"), TestInfo14)
		end
	end
end

directory8 = "C:\\Users\\camara.casson\\Dropbox (UFL)\\research-share\\Camara\\ccRCC-TIME-analysis\\Results-and-Analysis\\Panel2-Tumor-Cutoff5\\interdist"
RunStatistics_Inter(directory8)


# #plotting method?
# # histogram([data for data in dist if length(data)>50],normalize=:pdf,bins=:scott)

# # plot([data for data in intradist if length(data)>50],normalize=:pdf,bins=:sturges,layout=(3,1),label=nothing,seriestype=:barhist,xrotation=45,xtickfontsize=6,ytickfontsize=6,seriescolor=[:green :red :blue])