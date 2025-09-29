
#singularity exec --bind /mnt:/mnt wf_deps4.sif julia "/mnt/c/Dokumente/TRR156_local/Projekte/SingularityPipelines/Genomics/wf_transcriptome/repeat_wfpipeline.jl"

using JSON, Dates

datadir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/epi2medaten"
codedir = joinpath(@__DIR__, "scripts")
rundir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/epi2me_test_runs"
logdir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/epi2me_test_runs/epi2metestleer"
ref_genome = "hg38_chr20.fa"
ref_transcriptome = "ref_transcriptome.fasta"
ref_annotation = "gencode.v22.annotation.chr20.gtf"
sample_sheet= "sample_sheet.csv"


datadir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/B03/scRNAseq/run1 - Kopie"
codedir = joinpath(@__DIR__, "scripts")
rundir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/b03_e1_runs"
logdir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/b03_e1_runs/epi2metestleer"
ref_genome = "Mus_musculus.GRCm39.dna.toplevel.fa"
ref_transcriptome = nothing
ref_annotation = "Mus_musculus.GRCm39.110.gtf"
sample_sheet= "sample_sheet.csv"

datadir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/B03/scRNAseq/run1 - Kopie"
codedir = joinpath(@__DIR__, "scripts")
rundir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/b03_e1_runs"
logdir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/b03_e1_runs/epi2metestleer"
ref_genome = nothing
ref_transcriptome = "Mus_musculus.GRCm39.cdna.all.fa"
ref_annotation = "Mus_musculus.GRCm39.111.gff3"
sample_sheet= "sample_sheet.csv"

datadir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/B03/scRNAseq/run1 - Kopie"
codedir = joinpath(@__DIR__, "scripts")
rundir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/b03_e1_runs"
logdir = raw"/mnt/c/Dokumente/TRR156_local/Testdaten/b03_e1_runs/epi2metestleer"
ref_genome = nothing
ref_transcriptome = "Mus_musculus.GRCm39.cdna.all.fa"
ref_annotation = "Mus_musculus.GRCm39.111.gff3"
sample_sheet= "sample_sheet.csv"

function maybe_joinpath(path, maybe_path)
    return isnothing(maybe_path) ? nothing : joinpath(path, maybe_path)
end

function frmt_date(dt) 
    return replace(string(dt), r"[^\w\d]" => '_')
end

function copydir(src,dst)
    println("Start copying...")
    for (r,dirs,files) in walkdir(src)
        for d in dirs
            println("Copying ", basename(d), "...")
            cp(joinpath(r,d), joinpath(dst,d))
        end
        for f in files
            println("Copying ", basename(f), "...")
            cp(joinpath(r,f), joinpath(dst,f))
        end
    end
end

function run_experiments(Ns)
    for k in Ns
        println("Preparing Run ", k, "...")
        outdir = joinpath(rundir, "epi2metest$k")
        if isdir(outdir)
            #println("Target directory ", outdir, " already exists. Skipping Run.")
            println("Target directory ", outdir, " already exists.")
            #continue
        else
            mkpath(outdir)
        end
        inputdir = joinpath(outdir, "input")
        if !isdir(inputdir)
            try
                cp(datadir, inputdir)
            catch e
                println("Error: ", e)
            end
        end
        fastq_directory =  "differential_expression_fastq"
        println("Start processing Run ", k, "...")
        try
            now = frmt_date(Dates.now())
            isdir(logdir) || mkpath(logdir)
            logfile = joinpath(logdir,string("run_", k, '_', now, ".txt"))
            write(logfile, "")
            config = Dict("working_directory" => outdir,
                          #"input_directory" => joinpath(inputdir,
                          "fastq_directory" => joinpath(inputdir,fastq_directory),
                          "output_directory" => joinpath(outdir,"output"),
                          "ref_genome" => maybe_joinpath(inputdir,ref_genome),
                          "ref_transcriptome" => maybe_joinpath(inputdir,ref_transcriptome),
                          "ref_annotation" => joinpath(inputdir,ref_annotation),
                          "mapping" => Dict(
                                "indexing" => "-k14",
                                "alignment" => "-uf",
                                "min_mapping_quality" => 40,
                                "transcriptome_assembly" => "--conservative"
                          ),
                          "differential_expression_analysis" => Dict(
                                "sample_sheet" => maybe_joinpath(inputdir,sample_sheet),
                                "minimum_gene_expression" => 10,
                                "minimum_gene_expression_samples" => 3,
                                "minimum_feature_expression" => 3,
                                "minimum_feature_expression_samples" => 1,
                          ),
            )
            workflow_config_path = "workflow_configuration.json"
            open(workflow_config_path, "w") do io
                JSON.print(io, config)
            end
            script = joinpath(codedir, "run_wftrpipeline.jl")
            #run(pipeline(`julia -t 8 $(script)`; stdout="WFTR_RUNLOG$(k).txt", stderr="WFTR_RUNERRS$(k).txt");wait=true)
            run(pipeline(`julia -t 8 $(script) $(workflow_config_path)`);wait=true)
        catch e
            println(e)
        end
        
        now = frmt_date(Dates.now())
        isfile("runlog.txt") && mv("runlog.txt", "runlog_run_$(k)_$(now).txt")
        
    end
end


println("Starting...")
run_experiments(301)











