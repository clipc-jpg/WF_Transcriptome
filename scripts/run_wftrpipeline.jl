
using JSON

#TODO put into configuration
nattempts = 1
isempty(ARGS) && throw("Script requires one argument: the path to a configuration file")
config_path = ARGS[1]
config = open(JSON.parse, config_path, "r")
#workdir = config["working_directory"]
workdir = realpath("./")
logfile = joinpath(workdir, "logs", "runlog.txt") #config["logging_location"]
outputdir = joinpath(workdir, "output")
finaloutput_dir = joinpath(outputdir, "de_analysis")
expected_outputs = ["coldata.tsv", "dtu_plots.pdf", 
                   "results_dexseq.tsv", "results_dge.pdf", 
                   "results_dge.tsv", "results_dtu.pdf", 
                   "results_dtu_gene.tsv", "results_dtu_stageR.tsv", 
                   "results_dtu_transcript.tsv"]

for k in 1:nattempts
    println("Processing attempt: ", k)
    try
        include(joinpath(@__DIR__, "wftrpipeline.jl"))
        opmode = isfile(logfile) ? "a" : "w"
        open(logfile,opmode) do io
            println(io, "\nPipeline done!")
        end
    catch e
        println("Error: ", e)
    end
    if isfile(logfile) && 
                logfile |> eachline |> !isempty &&
                logfile |> eachline |> last |> contains("Pipeline done!") && 
                isdir(finaloutput_dir) &&
                all(in(readdir(finaloutput_dir)),expected_outputs)
        break
    else
        k0 = maximum(readdir(workdir;join=true);init=0) do dir
            if isdir(dir)
                mtch = match(r"_(\d+)\z",dir)
                isnothing(mtch) ? 0 : parse(Int64,first(mtch))
            else
                0
            end
        end
        isdir(outputdir) && mv(outputdir, join((outputdir, "_run_", k0+1)))
    end
end





