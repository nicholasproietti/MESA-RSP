# Code that creates color-magnitude diagrams using data from Modules for Experiments in Stellar Astrophysics

using DataFrames
using CSV
using Pkg
using Plots
# for single line comments
#=
for multi-line comments
=#

cwd = pwd()
sets = ['A','B','C','D']

UsedFiles = 0

AbsV = []
AbsI = []
AbsVMinusI = []

for j ∈  1:length(sets)
	println("Creating Color-magnitude diagram for Set " * sets[j])
	CMDName = "CMD_Set"*sets[j]*".png"

	ParentDir = joinpath(cwd,"LOGS_"*sets[j])
	dirs = readdir(ParentDir, join=true, sort=true)

	global AbsV = []
	global AbsI = []
	global AbsVMinusI = []

	global UsedFiles = 0

	println("Loading files...")
	for i ∈  1:length(dirs)
		LINAFile = joinpath(dirs[i], "LINA_period_growth.data")
		CurrentFile = joinpath(dirs[i],"history.data")

		if !isfile(LINAFile)
			continue
		else
		end

		LINAData = CSV.read(LINAFile, DataFrame, skipto = 2, header = ["","P(days)", "growth"], delim=' ', ignorerepeated=true)
		GrowthRates = LINAData[:, :3]

		if (GrowthRates[1] <= 0)
			continue
		else
		end

		if !isfile(CurrentFile)
			continue
		else
		end

		global UsedFiles += 1

		data = DataFrame(CSV.File(CurrentFile, skipto = 7, header = 6, delim=' ', ignorerepeated=true))

		TempV = log10.(data[:, :42])
		TempI =	log10.(data[:, :43])

		TempVMinusI = TempV - TempI

		global AbsV = push!(AbsV, TempV)
		global AbsI = push!(AbsI, TempI)
		global AbsVMinusI = push!(AbsVMinusI,TempVMinusI)

		println("Reading file: " * string(i) * " of " * string(length(dirs)))
	end

	MyPlot = scatter(AbsVMinusI, AbsV, title = "Set " * sets[j] * " Color-magnitude diagram", xlabel = "V - I (mag)", ylabel = "V (mag)", yflip = true, seriescolor = RGBA(255/255, 117/255, 255/255, 1.0), markershape= :+, legend = false)
	scatter!(MyPlot, AbsVMinusI, AbsV, title = "Set " * sets[j] * " Color-magnitude diagram", xlabel = "V - I (mag)", ylabel = "V (mag)", yflip = true, seriescolor = RGBA(255/255, 117/255, 255/255, 1.0), markershape= :+, legend = false)
	savefig(CMDName)

	println("Color magnitude diagram saved as " * CMDName * " in " * cwd)
	println("Total model directories: " * string(length(dirs)))
	println("Total positive FU growth rate model files: " * string(UsedFiles))

end
