




using Random

################################################################################
## creating FIFO that will be used to mask /dev/random and /dev/urandom inside the container
## this hopefully makes random algorithms deterministic (not completely deterministic, due to out-of-order execution)
################################################################################

println("Preparing FIFO")

fifo_id = rand('a':'z',4) |> join
fifo_name = "singularity_random_fifo_" * fifo_id

cd("/WFTDeps") do
    isfifo(fifo_name) && rm(fifo_name)
    ispath(fifo_name) && throw("Some object already occupies \"/WFTDeps/$(fifo_name)\"")
    run(`mkfifo $(fifo_name)`)
    # filling the FIFO must be done in a separate process
    task = Threads.@spawn run(`julia -q /WFTDeps/fill_urandom.jl /WFTDeps/$(fifo_name) 42`)
end

cd(raw"/mnt/c/Dokumente/TRR156_local/Projekte/SingularityPipelines/Genomics/wf_transcriptome0/") do
    singcomm = `singularity exec --bind /mnt:/mnt  --bind /WFTDeps/$(fifo_name):/dev/random --bind /WFTDeps/$(fifo_name):/dev/urandom wf_deps4.sif julia "/mnt/c/Dokumente/TRR156_local/Projekte/SingularityPipelines/Genomics/wf_transcriptome/repeat_wfpipeline.jl"`

    run(singcomm)
end










