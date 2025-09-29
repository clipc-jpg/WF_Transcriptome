
using CSV, DataFrames, JSON
include(joinpath(@__DIR__, "JSONLogger.jl"))
include(joinpath(@__DIR__, "pipeline_steps.jl"))

config_path = isempty(ARGS) ? "workflow_configuration.json" : ARGS[1]




################################################################################
## Input validation
################################################################################

wfinfo("Starting Input validation"; operation="Input validation",
       inputs=(),
       outputs=(),
       attempt_num=-1)

config = open(JSON.parse, config_path, "r")



wf_assert(config, "Configuration specifies necessary inputs") do cnf
    all(in(keys(cnf)), ("fastq_directory_path", "ref_annotation_path", "ref_genome_path", "ref_transcriptome_path", "mapping", "differential_expression_analysis"))
end

wf_assert(config, "Configuration specifies fastq_directory_path and it is a directory") do cnf
    obj = cnf["fastq_directory_path"]
    isdir(obj)
end

wf_assert(config, "Configuration specifies ref_annotation_path and it is a file") do cnf
    obj = cnf["ref_annotation_path"]
    isfile(obj)
end

wf_assert(config, "Configuration specifies either ref_genome or ref_transcriptome, and the specified files exist") do cnf
    rgp = cnf["ref_genome_path"]
    rtp = cnf["ref_transcriptome_path"]
    if !isnothing(rgp) && !isnothing(rtp)
        isfile(rgp) && isfile(rtp)
    elseif !isnothing(rgp)
        isfile(rgp)
    elseif !isnothing(rtp)
        isfile(rtp)
    else
        false
    end
end


mappingconfig = config["mapping"]

wf_assert(mappingconfig, "Configuration specifies all necessary mapping settings") do cnf
    all(in(keys(cnf)), ("indexing", "alignment", "min_mapping_quality", "transcriptome_assembly",))
end



wf_assert((config["ref_genome_path"],mappingconfig), "Minimum mapping quality is a whole number when a reference genome is used") do (ref_g, mcnf)
    mmq = mcnf["min_mapping_quality"]
    if isnothing(ref_g)
        true
    elseif mmq isa Integer
        true
    elseif mmq isa String
        !isnothing(tryparse(Int64, mmq))
    else
        false
    end
end


dea_config = config["differential_expression_analysis"]

if !isnothing(dea_config)

    wf_assert(dea_config, "DEA Configuration specifies necessary inputs") do cnf
        all(in(keys(cnf)), ("sample_sheet_path", "minimum_gene_expression", "minimum_gene_expression_samples", "minimum_feature_expression", "minimum_feature_expression_samples"))
    end


    wf_assert(dea_config, "DEA Configuration sample sheet exists") do cnf
        isfile(cnf["sample_sheet_path"])
    end

    for k in ("minimum_gene_expression", "minimum_gene_expression_samples", "minimum_feature_expression", "minimum_feature_expression_samples")
        wf_assert(dea_config, string("DEA configuration input '", k , "' is a whole number")) do cnf
            obj = cnf[k]
            if obj isa Integer
                true
            elseif obj isa String
                !isnothing(tryparse(Int64, obj))
            else
                false
            end
        end
    end
end


###  Example configuration
#
# {
#     "fastq_directory_path": "/mnt/c/Dokumente/TRR156_local/Testdaten/B03/Exp25_06_30/colony_analysis/input/full_data_pass_basecalling_filtered",
#     "ref_annotation_path": "/mnt/c/Dokumente/TRR156_local/Testdaten/B03/Exp25_06_30/colony_analysis/input/Mus_musculus_balbcj.BALB_cJ_v1.113.gtf",
#     "ref_genome_path": null,
#     "ref_transcriptome_path": "/mnt/c/Dokumente/TRR156_local/Testdaten/B03/Exp25_06_30/colony_analysis/input/Mus_musculus_balbcj.BALB_cJ_v1.cdna.all.fa",
#     "mapping": {
#         "indexing": "-k14",
#         "alignment": "-uf",
#         "min_mapping_quality": 40,
#         "transcriptome_assembly": "--conservative"
#     },
#     "differential_expression_analysis": {
#         "sample_sheet_path": "/mnt/c/Dokumente/TRR156_local/Testdaten/B03/Exp25_06_30/colony_analysis/input/sample_sheet.csv",
#         "minimum_gene_expression": 10,
#         "minimum_gene_expression_samples": 3,
#         "minimum_feature_expression": 3,
#         "minimum_feature_expression_samples": 1
#     }
# }




################################################################################
## Checking and creating output directory
################################################################################

working_directory = realpath("./")

logdir = joinpath(working_directory, "logs")
isdir(logdir) || mkpath(logdir)

logpath = begin
        logfiles = readdir(logdir)
        lognum = 1 + maximum(logfiles; init=0) do lfile
            mtch = match(r"run_(\d+)\.txt\z",lfile)
            if isnothing(mtch) 
                0
            else
                parse(Int64, first(mtch))
            end
        end
        joinpath(logdir, string("runlog_", lognum, ".txt"))
end

loggerio = open(logpath, "a")
logger = JSONLogger(;io=loggerio,ioerr=loggerio,flush_size=-1)

Logging.global_logger(logger)

wfinfo("Starting Pipeline"; operation=basename(@__FILE__), 
                            inputs=(),
                            outputs=(),
                            attempt_num=-1)

outputdir = created_directory(joinpath(working_directory, "output"))


srcdir = config["fastq_directory_path"]

scriptdir = checked_directory(@__DIR__)

subdirs = ("caching", "unzippedfiles", "concatenateddata", "trimmeddata",
            "minimapindex", "alignments", "wflogs", "gff_files",
            "transcriptome_output", "query_annotation_outputdir")

cachingdir, unzipdir, concatdir, trimdir, minimapindexdir, 
    alignmentdir, logdir, gffdir, transcriptome_outputdir, 
    query_annotation_outputdir = created_subdirectory.(outputdir, subdirs)


################################################################################
## Checking existence of reference files and optionally creating some
################################################################################

emptyfile = joinpath(working_directory, "emptyfile.txt")
open(io-> print(io,""), emptyfile, "w")

################################################################################
## Stocktaking sample fastq files
## fastq files are organized inside one directory for each sample source
## their relation to comparison groups is inside the sample_sheet.csvm if that exists
################################################################################

srcfilegroups = map(walkdir(srcdir)) do (r,ds,fs)
    filter!(contains(".fastq"),fs)
    sampleid = replace(basename(r), "barcode" => "sample")
    sampleid => joinpath.(r,fs)
end



# srcfilegroups = if isnothing(dea_config)
#         map(walkdir(srcdir)) do (r,ds,fs)
#             filter!(contains(".fastq"),fs)
#             sampleid = replace(basename(r), "barcode" => "sample")
#             sampleid => joinpath.(r,fs)
#         end
#     else
#         sample_sheet = checked_file(dea_config["sample_sheet_path"])
#
#         sample_df = CSV.read(sample_sheet, DataFrame)
#
#         wf_assert(names(sample_df), "Samplesheet contains columns 'barcode', 'alias' and 'condition'") do nms
#             all(in(nms), ("barcode", "alias", "condition"))
#         end
#
#         wf_assert(allunique, sample_df.barcode, "Entries in Samplesheet column 'barcode' are unique")
#
#         sample_map = Dict(zip(sample_df.barcode, sample_df.alias))
#
#
#         srcs = Pair{String, String}}[]
#         for (r,ds,fs) in walkdir(srcdir)
#             filter!(contains(".fastq"),fs)
#
#             barcode = basename(r)
#             if haskey(sample_map, barcode)
#                 sampleid = sample_map[barcode]
#                 entry = sampleid => joinpath.(r,fs)
#                 push!(srcs, entry)
#             end
#         end
#
#         unused_keys = setdiff(keys(sample_map), first.(srcs))
#
#         if !isempty(unused_keys)
#             for barcode in unused_keys
#                 wferror(string("Failed assertion: Barcode '", barcode, "' specified in sample_sheet, and a directory of that name was found");
#                         operation="assertion",
#                         inputs=inputs,
#                         outputs=(),
#                         attempt_num=attempt_num)
#             end
#         end
#
#         wf_assert(isempty, unused_keys, string("All specified barcode subdirectories could be found in the given fastq directory"))
# end



wf_assert(srcfilegroups, "Fastq sample groups have been found") do sfgs
    !isempty(sfgs)
end

filter!(!isempty∘last,srcfilegroups)

wf_assert(srcfilegroups, "No fastq sample group is empty") do sfgs
    !isempty(sfgs) && !any(isempty∘last, sfgs)
end

wf_assert(srcfilegroups, "No fastq file is empty") do sfgs
    all(sfgs) do (sample_id, srcfiles)
        !any(==(0)∘filesize, srcfiles)
    end
end

#wfinfo(string("Found these fastq source files: ", srcfilegroups))

################################################################################
## Processing each sample source individually, each with possibly many fastq files
################################################################################

foreach(srcfilegroups) do (sample_id, srcfiles)

    wfinfo(string("Processing ", sample_id); operation="Start first Loop")
    
    ################################################################################
    ## Decompressing files
    ################################################################################
    
    wfinfo("Step a.1")
    
    foreach(srcfiles) do file
        decompress_file(file,joinpath(unzipdir, sample_id))
    end
    

    ################################################################################
    ## Concatenating all files for a sample source to one large file
    ################################################################################
    wfinfo("Step a.2")
    unzippedfiles = readdir(joinpath(unzipdir,sample_id), join=true)
    while length(unzippedfiles) < length(srcfiles)
        global unzippedfiles
        sleep(0.5)
        unzippedfiles = readdir(joinpath(unzipdir,sample_id), join=true)
    end
	
    concatfile = concatenate_files(unzippedfiles, joinpath(concatdir, "all_data_$(sample_id).fastq"))

    ################################################################################
    ## Pychopper must be used if and only if cDNA is being used
    ## Pychopper finds out which direction the cDNA points towards
    ## and then trims away primer fragments left over from their PCR generation
    ################################################################################
    wfinfo("Step a.3")
    trimfile = joinpath(trimdir, "all_data_$(sample_id)_trimmed.fastq")
    pychopper_pdf = joinpath(trimdir, "pychopper.pdf")
    pychopper_tsv = joinpath(trimdir, "pychopper.tsv")
    #TODO: parameterize number of threads
    cd(trimdir) do
        runprogr("pychopper", "-t", "4", concatfile, trimfile;
                    outputs=(trimfile,pychopper_pdf,pychopper_tsv),
                    cmdout=IndirectWFLogStump(wfinfo,joinpath(logdir, "pychopper_templog_out.txt")),
                    cmderr=IndirectWFLogStump(wfinfo,joinpath(logdir, "pychopper_templog_err.txt")))
    end
end

if !isnothing(config["ref_genome_path"])
    wfinfo("Step a.0")
    
    referencefile = checked_file(config["ref_genome_path"])
    minimapindexfile = joinpath(minimapindexdir, "referenceindex.mmi")

    runprogr("minimap2", "-d", minimapindexfile, referencefile; 
                outputs=(minimapindexfile,),
                cmdout=IndirectWFLogStump(wfinfo,joinpath(logdir,"templog1.txt")),
                cmderr=IndirectWFLogStump(wfinfo,joinpath(logdir,"templog2.txt")))
    
    foreach(srcfilegroups) do (sample_id, srcfiles)
        
        trimfile = joinpath(trimdir, "all_data_$(sample_id)_trimmed.fastq")
        pychopper_pdf = joinpath(trimdir, "pychopper.pdf")
        pychopper_tsv = joinpath(trimdir, "pychopper.tsv")
        
        ################################################################################
        ## Minimap2 (in default mode) takes fastq files containing RNA/cDNA-reads and finds where in the 
        ## reference genome those fragments supposedly originate from (Alignment)
        ################################################################################
        wfinfo("Step a.4")
        alignmentfile = joinpath(alignmentdir, "alignment_$(sample_id).sam")
        runprogr("minimap2", "-a", referencefile, trimfile;
                    outputs=(alignmentfile,),
                    cmdout=alignmentfile,
                    cmderr=IndirectWFLogStump(wfinfo, joinpath(logdir,"templog5.txt")))

        seqkitlogfile = joinpath(logdir, "internal_priming_fail.tsv")

        sortedreadsfile = joinpath(alignmentdir, "reads_$(sample_id)_aln_sorted.bam")
        sortingreadslog = joinpath(logdir, "read_aln_$(sample_id)_stats.tsv")

        cachingfiles = [joinpath(cachingdir, "file$(k)_$(sample_id).throwaway") for k in 1:5]

        wfinfo("Step 1")
        runprogr("minimap2", "-t", "4", "-ax", "splice", minimapindexfile, trimfile;
                    outputs=(cachingfiles[1],),
                    cmdout=cachingfiles[1],
                    cmderr=IndirectWFLogStump(wfinfo, joinpath(logdir,"templog6.txt")))
        wfinfo("Step 2")
        ################################################################################
        ## "samtools view" takes a bam file and decompresses it into 
        ## a human-readable sam file
        ################################################################################
        runprogr("samtools", "view", "-q", mappingconfig["min_mapping_quality"], "-F", "2304", "-Sb";
                    cmdin=cachingfiles[1],
                    outputs=(cachingfiles[2],),
                    cmdout=cachingfiles[2])
        
        ################################################################################
        ## "seqkit bam" creates some summary statistics on a bam file
        ## like: how many alignments have failed,
        ## or how many led to ambiguous results
        ################################################################################
        
        seqkityamlfile = create_seqkityamlfile(joinpath(logdir,"seqkitToolbox.yml"), referencefile, seqkitlogfile)
        
        wfinfo("Step 3")
        #TODO: parameterize number of threads
        runprogr("seqkit", "bam", "-j", "4", "-x", "-T", "{Yaml: $(seqkityamlfile)}";
                    cmdin=cachingfiles[2],
                    outputs=(cachingfiles[3],),
                    cmdout=cachingfiles[3])
        
        ################################################################################
        ## sort alignments by leftmost coordinate in reference genome
        ################################################################################
        wfinfo("Step 4")
		
        #TODO: parameterize number of threads
        runprogr("samtools", "sort", "-@", "4", "-o", sortedreadsfile, "-T", "{Yaml: $(seqkityamlfile)}";
                    cmdin=cachingfiles[3],
                    outputs=(cachingfiles[4],sortedreadsfile),
                    cmdout=cachingfiles[4])
        
        wfinfo("Step 5")
		
        runprogr("seqkit", "bam", "-s", "-j", "4", "-";
                    cmdin=sortedreadsfile,
                    outputs=(cachingfiles[5],),
                    cmdout=cachingfiles[5],
                    cmderr=cachingfiles[5])
        
        wfinfo("Step 6")
        #isfile(sortingreadslog) || run(pipeline(`tee $(sortingreadslog)`)) # stdout=stdout
        

        
        
        ################################################################################
        ## stringtie takes positional mappings from the alignment step
        ## and estimates which fragments belong to the same transcript
        ## TODO: set G_FLAG to reference transcriptome if one is used
        ## TODO: split bam file and call once for each bundle with their prefix
        ################################################################################
        wfinfo("Step 7")
        # TODO: set G_FLAG to reference transcriptome if one is used
        G_FLAG = ""
        prefix = "CompleteBAM_$(sample_id)"
        complBAMout = joinpath(outputdir, prefix * ".gff")
        #run(pipeline(`stringtie --rf -L -v -p 1 -o $(complBAMout) -l $(prefix) $(sortedreadsfile)`;stderr=devnull); wait=true)
        runprogr("stringtie", "--rf", "-L", "-v", "-p", "1", mappingconfig["transcriptome_assembly"], "-o", complBAMout,
                               "-l", prefix, sortedreadsfile;
                    outputs=(complBAMout,),
                    cmderr=devnull)
        #stringtie --rf ${G_FLAG} -L -v -p ${task.cpus} ${params.stringtie_opts} \
        #    -o  ${prefix}.gff -l ${prefix} ${bam} 2>/dev/null
        
        
        ################################################################################
        ## this step by hand appends outputs from stringtie
        ## that have been made for the same sample
        ## currently useless, since the bam file had not been split into
        ## smaller, more manageable chunks
        ################################################################################
        wfinfo("Step 8")
        gffmerged = merge_gff_files((complBAMout,), joinpath(gffdir, "merged_$(sample_id).gff"))

        ################################################################################
        ## GffCompare provides classification and reference annotation mapping and
        ## matching statistics for RNA-Seq assemblies (transfrags) or other generic
        ## GFF/GTF files.
        ## GffCompare also clusters and tracks transcripts across multiple GFF/GTF
        ## files (samples), writing matching transcripts (identical intron chains) into
        ## <outprefix>.tracking, and a GTF file <outprefix>.combined.gtf which
        ## contains a nonredundant set of transcripts across all input files (with
        ## a single representative transfrag chosen for each clique of matching transfrags
        ## across samples).
        ################################################################################
        wfinfo("Step 9")
        #if !isnothing(config["ref_transcriptome_path"])
            #ref_annotation = checked_file(joinpath(inputdir,"ref_transcriptome.fasta"))
            ref_annotation = isnothing(config["ref_transcriptome_path"]) ? emptyfile : checked_file(config["ref_transcriptome_path"])
            query_annotation = gffmerged

            gffcoutp = joinpath(query_annotation_outputdir, "str_merged_$(sample_id)")
            gffcoutpfile = joinpath(query_annotation_outputdir, "str_merged_$(sample_id).combined.gtf")
            #isfile(gffcoutpfile) || run(`gffcompare -o $(gffcoutp) -r $(ref_annotation) $(query_annotation)`; wait=true)
            runprogr("gffcompare", "-o", gffcoutp, "-r", ref_annotation, query_annotation;
                        outputs=(gffcoutpfile,))
            
            #gffcompare -o ${out_dir}/str_merged -r ${ref_annotation} \
            #            ${params.gffcompare_opts} ${query_annotation}


            ################################################################################
            ## here, gffread extracts spliced exons for each transcript
            ################################################################################
            wfinfo("Step 10")
            transcriptomeoutp = joinpath(transcriptome_outputdir, "transcriptome_$(sample_id).fas")
            transcripts_gff = gffcoutpfile

            #isfile(transcriptomeoutp) || run(`gffread -g $(referencefile) -w $(transcriptomeoutp) $(transcripts_gff)`; wait=true)
            runprogr("gffread", "-g", referencefile, "-w", transcriptomeoutp, transcripts_gff;
                        outputs=(transcriptomeoutp,))
            
            #=
            get_transcriptome(
                            merge_gff_bundles.out.gff
                            .join(run_gffcompare.out.gffcmp_dir)
                            .join(seq_for_transcriptome_build))

            tuple val(sample_id), path(transcripts_gff), path(gffcmp_dir), path(reference_seq)
            gffread -g ${reference_seq} -w ${transcriptome} ${transcripts_gff}
                    if  [ "\$(ls -A $gffcmp_dir)" ];
                        then
                            gffread -F -g ${reference_seq} -w ${merged_transcriptome} $gffcmp_dir/str_merged.annotated.gtf
                    fi
            =#
    #    else
    #        query_annotation_file = joinpath(query_annotation_outputdir, basename(gffmerged))
    #        isfile(query_annotation_file) || cp(gffmerged, query_annotation_file)
        #end
    end

    ################################################################################
    ## Merge all alignment results from each sample
    ## into one transcriptome for later pseudo-alignment (i.e. alignment of reads with the help of a pre-assembled transcriptome)
    ##
    ## 1. merge all transcripts into one file
    ## 2. remove everything not labelled transcript
    ## 3. remove/merge duplicate transcripts
    ## 4. remove empty sequences with sed regex
    ## 5. ??? probably remove empty lines
    ## 6. clean up intermediate files
    ################################################################################
    wfinfo("Step 11")
    stringtieoutp = joinpath(outputdir, "stringtie.gtf")
    ref_annotation = checked_file(config["ref_annotation_path"])
    final_transcriptome = joinpath(outputdir, "final_non_redundant_transcriptome.fasta")
    query_annotations = readdir(query_annotation_outputdir; join=true)
    println("Query annotations: ", query_annotations)
    filter!(endswith(r".g[tf]f3?"), query_annotations)

    temp_transcriptome = joinpath(outputdir, "temp_transcriptome.fasta")
    temp_del_repeats   = joinpath(outputdir, "temp_del_repeats.fasta")
    temp_rm_empty_seq  = joinpath(outputdir, "temp_rm_empty_seq.fasta")

    wfinfo("Step 11.1")

        runprogr("stringtie", "--merge", "-G", ref_annotation, "-p", "1", "-o", stringtieoutp, query_annotations;
                    outputs=(stringtieoutp,))

    wfinfo("Step 11.2")
    # TODO: currently does not get skipped when already completed
    runprogr("seqkit", "subseq", "--feature", "transcript", "--gtf-tag", "transcript_id", "--gtf", stringtieoutp, referencefile;
                outputs=(temp_transcriptome,),
                cmdout=temp_transcriptome)
    wfinfo("Step 11.3")
    runprogr("seqkit", "rmdup", "-s";
                cmdin=temp_transcriptome,
                outputs=(temp_del_repeats,),
                cmdout=temp_del_repeats)
                
                
    wfinfo("Step 11.4")

    clean_temp_del_repeats(temp_del_repeats, temp_rm_empty_seq)

    wfinfo("Step 11.5")
    runprogr("awk", "BEGIN {RS = \">\" ; FS = \"\\n\" ; ORS = \"\"} \$2 {print \">\"\$0}",  temp_rm_empty_seq;
                outputs=(final_transcriptome,),
                cmdout=final_transcriptome)

    wfinfo("Step 11.6")
    delete_file(temp_transcriptome)
    delete_file(temp_del_repeats)
    delete_file(temp_rm_empty_seq)

else
    cp(config["ref_transcriptome_path"], joinpath(outputdir, "final_non_redundant_transcriptome.fasta");force=true)
    #TODO: in reference_genome path and here and in the following, rename stringtie.gtf to annotation.gtf
    cp(config["ref_annotation_path"], joinpath(outputdir, "stringtie.gtf"); force=true)
end

################################################################################
## 
## 
## at this point, reference-guided and precomputed workflows converge
## from here on, they differ in the source of "final_transcriptome"
## 
## 
################################################################################





################################################################################
## checking input files for second  phase: pseudo-alignment
## TODO: move section towards beginning for earlier error management
################################################################################
#sample_sheet = checked_file(joinpath(inputdir, "sample_sheet.csv")) # that one is universal
#sample_sheet = checked_file(joinpath(inputdir, config["sample_sheet_path"])) # that one is universal
#ref_annotation = checked_file(joinpath(inputdir, "gencode.v22.annotation.chr20.gtf")) #referencetranscriptome
ref_annotation = checked_file(config["ref_annotation_path"]) #referencetranscriptome
final_transcriptome = checked_file(joinpath(outputdir, "final_non_redundant_transcriptome.fasta"))
    

foreach(srcfilegroups) do (sample_id, srcfiles)

    ################################################################################
    ## finding input reads again
    ## TODO: merge this redundant step into an earlier one
    ################################################################################
    wfinfo("Step 12")
    fastqreads = map(srcfiles) do file
        joinpath(unzipdir, sample_id, replace(basename(file),".gz" => ""))
    end
    wfinfo(string("Found fastqreads: ", fastqreads))
    ################################################################################
    ## index combined transcriptome with minimap2 for faster search
    ################################################################################
    mmioutp = joinpath(outputdir, "genome_index.mmi")
    minimapoutp = joinpath(outputdir, "output_$(sample_id).bam")
    samtoolsoutp = joinpath(outputdir, "reads_aln_sorted_$(sample_id).bam")
    # TODO: implement
    #run(pipeline(`checkSampleSheetCondition(sample_sheet) `); wait=true; wait=true) # refers to workflowglue
    #isfile(mmioutp) || run(pipeline(`minimap2 -t 1 -I 1000G -d $(mmioutp) "$(final_transcriptome)"`); wait=true)
    wfinfo("Step 12.1")
    runprogr("minimap2", "-t", "1", "-I", "1000G", mappingconfig["indexing"], "-d", mmioutp, final_transcriptome;
                outputs=(mmioutp,),
                cmderr=IndirectWFLogStump(wfinfo, joinpath(logdir,"templog7.txt")))
    #run(pipeline(pipeline(`minimap2 -t 1 -ax splice -uf -p 1.0 "${index}" "${fastq_reads}"`,
    #isfile(minimapoutp) || run(pipeline(pipeline(`minimap2 -t 1 -ax splice -uf -p 1.0 $(mmioutp) $(fastqreads)`,
    #                       `samtools view -Sb`); stdout=minimapoutp); wait=true)
    wfinfo("Step 12.2")
    cachefile = joinpath(cachingdir, "minimap2samtoolscache_$(sample_id).throwaway")
    runprogr("minimap2", "-t", "1", "-ax", "splice", mappingconfig["alignment"], "-p", "1.0", mmioutp, fastqreads;
                outputs=(cachefile,),
                cmdout=cachefile,
                cmderr=IndirectWFLogStump(wfinfo, joinpath(logdir,"templog8.txt")))
    wfinfo("Step 12.3")
    runprogr("samtools", "view", "-Sb";
                cmdin=cachefile,
                outputs=(minimapoutp,),
                cmdout=minimapoutp)
    #isfile(samtoolsoutp) || run(pipeline(`samtools sort -@ 1 $(minimapoutp) -o $(samtoolsoutp)`); wait=true)
    wfinfo("Step 12.4")
    runprogr("samtools", "sort", "-@", "1", minimapoutp, "-o", samtoolsoutp;
                outputs=(samtoolsoutp,))
    #=
    minimap2 -t ${task.cpus} -ax splice -uf -p 1.0 "${index}" "${fastq_reads}" \
    | samtools view -Sb > "output.bam"
    samtools sort -@ ${task.cpus} "output.bam" -o "${meta.alias}_reads_aln_sorted.bam"
    =#

    ################################################################################
    ## salmon pseudo-aligns reads with a reference transcriptome
    ## the output is a statistically useful table with read counts
    ################################################################################
    wfinfo("Step 13")
    salmonoutpdir = joinpath(outputdir, "counts")
    salmonoutp = joinpath(outputdir, "transcript_counts_$(sample_id).tsv")
    #isfile(salmonoutp) || run(pipeline(`salmon quant --noErrorModel -p 1 -t $(final_transcriptome) -l SF -a $(samtoolsoutp) -o $(salmonoutpdir)`); wait=true)
    #isfile(salmonoutp) || mv(joinpath(salmonoutpdir, "quant.sf"), salmonoutp)
    
    # TODO: currently does not get skipped when already completed
    # => combine steps into one OR set output to the moved file
    runprogr("salmon", "quant", "--noErrorModel", "-p", "1", "-t", final_transcriptome, 
                            "-l", "SF", "-a", samtoolsoutp, "-o", salmonoutpdir;
                outputs=(salmonoutp,))
                #outputs=(salmonoutpdir,),
                #skip_predicate=Returns(false))
    if !isfile(salmonoutp)
        salmonoutp_quant = checked_file(joinpath(salmonoutpdir, "quant.sf"))
        copy_file(salmonoutp_quant, salmonoutp)
    end
    
    ################################################################################
    ## "seqkit bam" creates some summary statistics on the new bam file
    ## like: how many alignments have failed,
    ## or how many led to ambiguous results
    ################################################################################
    seqkitstats = joinpath(outputdir, "seqkit.stats")
    seqkitstats_output = joinpath(outputdir, "sample_$(sample_id)_seqkit.stats")
    #run(pipeline(`seqkit bam $(samtoolsoutp)`;stderr=seqkitstats); wait=true)
    runprogr("seqkit", "bam", samtoolsoutp;
            outputs=(seqkitstats_output,),
            cmderr=seqkitstats)
    copy_file(seqkitstats, seqkitstats_output)
end
wfinfo("Done with per-sample processing")


################################################################################
## find all tables with read counts by name
################################################################################
wfinfo("Step 14")
countfiles = filter!(readdir(outputdir;join=true)) do file
    isfile(file) && contains(basename(file), r"transcript_counts_sample\d+.tsv")
end

wf_assert(countfiles, "Countfiles have been created") do files
    !isempty(files)
end

################################################################################
## merge all such count tables
##   - "Reference" is the gene identifying column
##   - filter out all genes for which there was no observation
##   - due to outerjoin, genes will reappear if there had been an observation
##     inside another sample file, but with "missing/NULL/NA" instead of 0
##   - rename all "Count"-column to their respective sample identifier to be able to compare them
################################################################################

#countdf = mapreduce((a,b)->outerjoin(a,b,on="Reference"),countfiles) do file
#    df = CSV.read(file, DataFrame)
#    rename!(df, "Name" => "Reference")
#    df.Count = df.NumReads
#    filter!(:Count => >(0),df)
#    sort!(df,:Count;rev=true)
#    name = match(r"(sample\d+)\.[^\.]+\z",file) |> first
#    rename!(df, "Count" => name)
#    df[!,["Reference", name]]
#end
#
#################################################################################
### during join, zero counts were removed but may have reappeared as missing
### => replace "missing" with 0
#################################################################################
#filter!(:Reference=>!ismissing, countdf)
#foreach(names(countdf)) do col
#    replace!(countdf[!,col], missing => 0.0)
#end
#
#de_transcript_counts = joinpath(outputdir, "de_transcript_counts.tsv")
#CSV.write(de_transcript_counts, countdf; delim = '\t')

de_transcript_counts = joinpath(outputdir, "de_transcript_counts.tsv")
merge_count_tables(countfiles, de_transcript_counts)

#run(pipeline(`workflow-glue merge_count_tsvs -z -o de_transcript_counts.tsv -tsvs ${counts}`); wait=true)

################################################################################
## repeat the upper steps with TPM-counts
## TODO: step is identical to non-TPM treatment
## TODO: find out where data came from
## (will not be used anymore => TODO: what are TPM counts?)
################################################################################
wfinfo("Step 15")

#TPMdf = mapreduce((a,b)->outerjoin(a,b,on="Reference"),countfiles) do file
#    df = CSV.read(file, DataFrame)
#    rename!(df, "Name" => "Reference")
#    df.Count = df.TPM
#    filter!(:Count => >(0),df)
#    sort!(df,:Count;rev=true)
#    name = match(r"(sample\d+)\.[^\.]+\z",file) |> first
#    rename!(df, "Count" => name)
#    df[!,["Reference", name]]
#end
#de_tpm_counts = joinpath(outputdir,"de_tpm_transcript_counts.tsv")
#CSV.write(de_tpm_counts, TPMdf; delim = '\t')

de_tpm_counts = joinpath(outputdir, "de_tpm_transcript_counts.tsv")
merge_count_tables(countfiles, de_tpm_counts)

#run(pipeline(`workflow-glue merge_count_tsvs -o de_tpm_transcript_counts.tsv -z -tpm True -tsvs $counts`); wait=true)

################################################################################
## TODO?: move variable definition upwards
################################################################################
merged_tsv = de_transcript_counts #it is not "de_tpm_transcript_counts.tsv"
annotation_gtf = checked_file(joinpath(outputdir, "stringtie.gtf"))


################################################################################
## Sometimes, direction of the transcripts cannot be inferred
## this will not be accepted by downstream tools
## also, removing these transcripts will not be tolerated
## ergo, the NULL-symbol "." is arbitrarily being replaced by "+"
## this fix is obviously wrong
## TODO: find out how to filter these cases out instead
## TODO: find out, whether those exceptions should occur in the first place
## TODO: find out, whether tool is stochastic and rerunning
##       some previous instructions can improve/save the annotation
################################################################################
annotstr = read(annotation_gtf, String)
annotstr = replace(annotstr, r"(\d+\t\d+\t\d+\t)\.(\t)" => s"\1+\2")
open(annotation_gtf, "w") do io
    print(io, annotstr)
end

annotation_gtf = repair_annotation_gtf(annotation_gtf, annotation_gtf)

dea_config = config["differential_expression_analysis"]


################################################################################
## In this section differential gene expression analysis is done in R
## Subsequently, another R script creates visualizations
################################################################################

if !isnothing(dea_config)

    wfinfo("Starting differential analysis in R...")

    de_analysisdir = created_subdirectory(outputdir,"de_analysis")
    #sample_sheet = checked_file(joinpath(inputdir, "sample_sheet.csv")) # that one is universal
    sample_sheet = checked_file(dea_config["sample_sheet_path"]) # that one is universal

    cp(merged_tsv, joinpath(de_analysisdir, "all_counts.tsv"); force=true)
    cp(merged_tsv, joinpath(outputdir, "all_counts.tsv"); force=true)
    #cp(sample_sheet, joinpath(de_analysisdir, "coldata.tsv"); force=true)
    cp(sample_sheet, joinpath(de_analysisdir, "sample_sheet.tsv"); force=true)
    cp(sample_sheet, joinpath(outputdir, "sample_sheet.csv"); force=true)


    min_gene_expr = dea_config["minimum_gene_expression"]
    min_samps_gene_expr = dea_config["minimum_gene_expression_samples"]
    min_feature_expr = dea_config["minimum_feature_expression"]
    min_samps_feature_expr = dea_config["minimum_feature_expression_samples"]

    cd(outputdir) do
#    cd(de_analysisdir) do
        de_analysis_script = checked_file(joinpath(scriptdir, "de_analysis.R"))
        #run(pipeline(`Rscript $(de_analysis_script) $(annotation_gtf) $(min_samps_gene_expr) $(min_samps_feature_expr) $(min_gene_expr) $(min_feature_expr) $(annotation_type) $(strip_version)`); wait=true)
        mkdir("merged")
        runprogr("Rscript", de_analysis_script, annotation_gtf, min_samps_gene_expr, min_samps_feature_expr, min_gene_expr, min_feature_expr,
                    outputs=("cts.csv","txdf.csv", "all_counts_before_filtering.tsv",
                             "merged/filtered_transcript_counts_with_genes.tsv", "merged/all_gene_counts.tsv",
                             "de_analysis/results_dge.pdf", "de_analysis/results_dge.tsv",
                             "de_analysis/results_dtu.pdf", "de_analysis/results_dtu_gene.tsv",
                             "de_analysis/results_dtu_transcript.tsv",
                             "de_analysis/results_dexseq.tsv",
                             "de_analysis/results_dtu_stageR.tsv"))
        cp("merged/filtered_transcript_counts_with_genes.tsv", "de_analysis/filtered_transcript_counts_with_genes.tsv"; force=true)
    end
    
    
    cd(de_analysisdir) do
        wfinfo("Starting to visualize results in R...")
        result_visualization_script = checked_file(joinpath(scriptdir, "plot_dtu_results.R"))
        #run(pipeline(`Rscript $(result_visualization_script)`); wait=true)
        runprogr("Rscript", result_visualization_script;
                    inputs=(result_visualization_script, 
                            "de_analysis/coldata.tsv", "de_analysis/results_dtu_stageR.tsv",
                            "merged/all_counts_filtered.tsv"),
                    outputs=("de_analysis/dtu_plots.pdf", ))
    end
end

println("WF-Tr-Pipeline completed.")







