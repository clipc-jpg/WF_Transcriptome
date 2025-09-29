
using CSV, DataFrames
using UUIDs

include(joinpath(@__DIR__, "JSONLogger.jl"))

################################################################################
## appinfo and global state
################################################################################

function app_info()
    return "0.1.0"
end

const OPCODES = Dict()

################################################################################
## creating/checking directories
################################################################################

function checked_directory(dir::AbstractString, error_message="Directory does not exist: $(dir)"; attempt_num=-1)
    if !isdir(dir)
        wferror(error_message; operation="checked_directory",
                              inputs=(dir,),
                              outputs=(),
                              attempt_num=attempt_num)
        throw(error_message)
    end
    
    return dir
end

# TODO: how to find out order of operations, if for one operation,
#       input and output file are the same??
function created_directory(dir::AbstractString, error_message="Directory does not exist: $(dir)"; attempt_num=-1)
    new_dir = try
        isdir(dir) || mkpath(dir)
        dir
    catch e
        wferror(string(e); operation="created_directory",
                              inputs=(),
                              outputs=(dir,),
                              attempt_num=attempt_num)
        throw(e)
    end
    
    return new_dir
end

function checked_subdirectory(parent_dir, subdir, error_message="Directory does not exist: $(joinpath(parent_dir, subdir))"; attempt_num=-1)
    new_dir = try
        new_dir = joinpath(parent_dir, subdir)
        isdir(new_dir) || throw(string(error_message, new_dir))
        new_dir
    catch e
        wferror(string(e); operation="checked_subdirectory",
                              inputs=(parent_dir,),
                              outputs=(subdir,),
                              attempt_num=attempt_num)
        throw(e)
    end
    
    return new_dir
end

function created_subdirectory(parent_dir, subdir, error_message="Directory could not be created: $(joinpath(parent_dir, subdir))"; attempt_num=-1)
    new_dir = try
        new_dir = joinpath(parent_dir, subdir)
        isdir(new_dir) || mkpath(new_dir)
        new_dir
    catch e
        wferror(string(e); operation="created_subdirectory",
                              inputs=(parent_dir,),
                              outputs=(subdir,),
                              attempt_num=attempt_num)
        throw(e)
    end
    
    return new_dir
end

################################################################################
## checking, moving, deleting files
################################################################################

function checked_file(predicate::Function, filepath::AbstractString, error_message="Failing assertion for file: $(filepath)"; attempt_num=-1)
    
    if !predicate(filepath)
        wferror(error_message; operation="checked_file",
                              inputs=(filepath,),
                              outputs=(),
                              attempt_num=attempt_num)
        throw(error_message)
    end
    
    return filepath
end

function checked_file(filepath::AbstractString, error_message="File does not exist: $(filepath)", attempt_num=-1)
    return checked_file(isfile, filepath, error_message; attempt_num=attempt_num)
end


function copy_file(file, target; attempt_num=-1, force=true)
    !force && isfile(target) && return target
    try
        cp(file, target; force)
    catch e
        wferror(string(e); operation="move_file",
                              inputs=(file,),
                              outputs=(),
                              attempt_num=attempt_num)
        throw(e)
    end
    
    return target
end


function move_file(file, target; attempt_num=-1, force=true)
    !force && isfile(target) && return target
    try
        mv(file, target; force)
    catch e
        wferror(string(e); operation="move_file",
                              inputs=(file,),
                              outputs=(),
                              attempt_num=attempt_num)
        throw(e)
    end
    
    return target
end

function delete_file(file, required=false)
    try
        enforced = required || isfile(file)
        enforced && rm(file)
    catch e
        wferror(string(e); operation="delete_file",
                              inputs=(file,),
                              outputs=(),
                              attempt_num=attempt_num)
        throw(e)
    end
    
    return nothing
end

################################################################################
## general assertions
################################################################################

function wf_assert(predicate::Function, inputs, assertion_description; attempt_num=-1)
    
    condition_fulfilled = try
            predicate(inputs)
        catch e
            false
    end

    if !condition_fulfilled
        wferror("Failed assertion: "*assertion_description; operation="assertion",
                              inputs=inputs,
                              outputs=(),
                              attempt_num=attempt_num)
        throw("Failed assertion: "*assertion_description)
    else
        wfinfo("Passed assertion: "*assertion_description; operation="assertion",
                              inputs=inputs,
                              outputs=(),
                              attempt_num=attempt_num)
    end
    
    return nothing
end


################################################################################
## running a command line step
################################################################################


function runprogr(cmdstring::AbstractString, 
                   argstrings...; 
                                log_message="Running pipeline step $(cmdstring)", error_message="Failed pipeline step $(cmdstring)", skip_message="All outputs already created. Skipping $(cmdstring).", 
                                inputs=argstrings, outputs=(), attempt_num=-1, required=false, 
                                skip_predicate=(outs->!isempty(outs) && all(ispath,outs)), 
                                wait=true, cmdin=stdin, cmdout=stdout, cmderr=stderr)
    args = String[]
    for arg in argstrings
        if arg isa String
            push!(args, arg)
        elseif arg isa Array
            append!(args, string.(arg))
        else
            push!(args, string(arg))
        end
    end
    cmd = isempty(argstrings) ? `$(cmdstring)` : `$(cmdstring) $(args)`
    # often (e.g with minimap2) outputs cannot be captured without error
    # TODO move code chunks into separate block
    pipein = cmdin
    pipeout = if cmdout isa WFLogStump
            WFLogPipe(cmdout, operation=cmdstring, inputs=inputs, outputs=outputs, attempt_num=attempt_num)
        elseif cmdout isa IndirectWFLogStump
            cmdout.dump_file
        else
            cmdout
    end
    pipeerr = if cmderr isa WFLogStump
            WFLogPipe(cmderr, operation=cmdstring, inputs=inputs, outputs=outputs, attempt_num=attempt_num)
        elseif cmderr isa IndirectWFLogStump
            cmderr.dump_file
        else
            cmderr
    end
    p = pipeline(cmd; stdin=pipein, stdout=pipeout, stderr=pipeerr)
    wfinfo(log_message; operation=cmdstring, inputs=inputs, outputs=outputs, attempt_num=attempt_num)
    if skip_predicate(outputs)
        wfinfo(skip_message; operation=cmdstring, inputs=inputs, outputs=outputs, attempt_num=attempt_num)
    else
        try
            result = run(p; wait=wait)
            if cmdout isa IndirectWFLogStump
                msg = read(cmdout.dump_file, String)
                cmdout.log_function(msg; operation=cmdstring, inputs=inputs, outputs=outputs, exitcode=result.exitcode, attempt_num=attempt_num)
            end
            if cmderr isa IndirectWFLogStump
                msg = read(cmderr.dump_file, String)
                cmderr.log_function(msg; operation=cmdstring, inputs=inputs, outputs=outputs, exitcode=result.exitcode, attempt_num=attempt_num)
            end
        catch e
            # there is always exactly one process involved
            exitcode = e isa ProcessFailedException ? e.procs[1].exitcode : -24
            wferror(string(error_message, e); operation=cmdstring, inputs=inputs, outputs=outputs, exitcode=exitcode, attempt_num=attempt_num )
            required && throw(e)
        end
    end
end

################################################################################
## specific, encapsulated steps in the pipeline
################################################################################

function collect_fastq_inputs(srcdir)
    # is that the most sensible way? might someone organize their files in a way such that 
    # this block will behave unexpectedly?
    wfinfo("Starting to collect fastq inputs..."; operation="collect_fastq_inputs")
    srcfilegroups = map(walkdir(srcdir)) do (r,ds,fs)
        filter!(contains(".fastq"),fs)
        # here might be some different naming and nesting scheme at work
        # e.g. /barcode1/organ1/fastq
        sampleid = replace(basename(r), "barcode" => "sample")
        sampleid => joinpath.(r,fs)
    end

    wf_assert(srcfilegroups, "Some fastq sample groups have been found") do sfgs
        !isempty(sfgs)
    end

    filter!(!isempty∘last,srcfilegroups)
    
    wf_assert(srcfilegroups, "No fastq sample group is empty") do sfgs
        !isempty(sfgs) && !any(isempty∘last, sfgs)
    end

    wf_assert(srcfilegroups, "No decompressed fastq file is empty") do sfgs
        all(sfgs) do (sample_id, srcfiles)
            !any(==(0)∘filesize, srcfiles)
        end
    end

    wfinfo("Collecting fastq inputs"; operation="collect_fastq_inputs", inputs=(srcdir,), outputs=last.(srcfilegroups), attempt_num=-1)
    
    return srcfilegroups
end

function decompress_file(file, targetdir)
    tardir = created_directory(targetdir)
    tarfile = joinpath(tardir, replace(basename(file),r"\.fastq\..*" => ".fastq"))
    isfile(tarfile) && return tarfile
    if endswith(file, ".fastq")
        mv(file, tarfile)
    elseif endswith(file, ".gz")
        #p = pipeline(`gzip -dk $(file)`)
        #run(p; wait=true)
        runprogr("gzip", "-dk", file, cmdout=stdout, cmderr=stderr)
        decompressed = replace(file, r"\.fastq\..*" => ".fastq")
        for _ in 1:5
            isfile(decompressed) && continue
            sleep(0.5)
        end
        mv(decompressed,tarfile)
    else
        wf_assert(Returns(false),(file,),"Unexpected file type: "*file)
    end
    
    return tarfile
end

function concatenate_files(files, targetfile; overwrite=false)
    !overwrite && isfile(targetfile) && return targetfile
    open(targetfile,"w") do tario
        for src in files
            open(src,"r") do srcio
                write(tario,srcio)
            end
        end
    end
    return targetfile
end

function create_seqkityamlfile(targetfile, referencefile, seqkitlogfile; overwrite=true)
    #!overwrite && isfile(targetfile) && return targetfile
    open(targetfile,"w") do io
        println(io,"AlnContext:")
        println(io,"  Ref: \"$(referencefile)\"")
        println(io,"  LeftShift: -10")
        println(io,"  RightShift: 10")
        println(io,"  RegexEnd: \"[Aa]{5,}\"")
        println(io,"  Stranded: True")
        println(io,"  Invert: True")
        println(io,"  Tsv: \"$(seqkitlogfile)\"")
    end
    return targetfile
end

function merge_gff_files(gff_files, targetfile; overwrite=true)
    #!overwrite && isfile(targetfile) && return targetfile
    open(targetfile,"w") do tario
        println(tario, "##gff-version 2")
        println(tario, "#pipeline-nanopore-isoforms: stringtie")
        for gff_file in gff_files
            for line in eachline(gff_file)
                contains(line,'#') || println(tario, line)
            end
        end
    end
    return targetfile
end


function merge_count_tables(countfiles, targetfile, join_column="Reference")
    countdf = mapreduce((a,b)->outerjoin(a,b,on=join_column),countfiles) do file
        df = CSV.read(file, DataFrame)
        rename!(df, "Name" => "Reference")
        df.Count = df.NumReads
        filter!(:Count => >(0),df)
        sort!(df,:Count;rev=true)
        name = match(r"(sample\d+)\.[^\.]+\z",file) |> first
        rename!(df, "Count" => name)
        df[!,["Reference", name]]
    end

    ################################################################################
    ## during join, zero counts were removed but may have reappeared as missing
    ## => replace "missing" with 0
    ################################################################################
    filter!(:Reference=>!ismissing, countdf)
    foreach(names(countdf)) do col
        replace!(countdf[!,col], missing => 0.0)
    end

    CSV.write(targetfile, countdf; delim = '\t')
    
    return targetfile
end

function repair_annotation_gtf(annotation_gtf, targetfile)
    annotstr = read(annotation_gtf, String)
    # replaces the fourth column in a tab-delimited file
    # from '.' to '+', which is obviously wrong
    annotstr = replace(annotstr, r"(\d+\t\d+\t\d+\t)\.(\t)" => s"\1+\2")
    open(targetfile, "w") do io
        print(io, annotstr)
    end
    
    return targetfile
end


#        pipeline(`cat $(temp_del_repeats)`, 
#                    `sed 's/>.* />/g'`, 
#                    `sed -e 's/_[0-9]* \[/ \[/g'`)
#run(pipeline(; stdout=temp_rm_empty_seq); wait=true)

function clean_temp_del_repeats(temp_del_repeats, targetfile)
    str = read(temp_del_repeats, String)
    str = replace(str, r">.* " => '>', r"_[0-9]* \[" => " [")
    write(targetfile, str)
    return targetfile
end









