


# abhängigkeiten inside singularity container
#singularity exec --bind /mnt:/mnt wf_deps4.sif julia "/mnt/c/Dokumente/TRR156_local/Projekte/SingularityPipelines/Genomics/wf_transcriptome/repeat_wfpipeline.jl"

using Dates

datadir = "/mnt/sds-hd/sd23k007/repeat_wftr/input"
superdir = "/mnt/sds-hd/sd23k007/repeat_wftr"
workdirs = superdir .* "repetition_c_" .* string.(1:6)

container_pth = "/mnt/sds-hd/sd23k007/repeat_simpleaf/simpleaf_chromium/alevin_fry_apps.sif"


slurm_template = """#!/bin/bash
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3:00:05
#SBATCH --mem=24gb

    # you need to change the "SV" to yours or else the script will fail
    # then run this file with sbatch
    # the printing outputs will be written into a slurm-*out file, at ~
    # the cellranger results are written to /mnt/sds-hd/\${SV}/cellranger_workdir/\${projectname}

    projectname={{project_name}}
	project_dir={{project_dir}}
	julia_script={{julia_script}}
	container_path={{container_path}}
	workflow_configuration={{workflow_configuration}}
	project_workdir={{project_workdir}}
	workdir={{workdir}}
	
    SV="sd23k007"

    echo "Project \${projectname} started"

    #load needed modules
        module load math/julia
		module load system/singularity
    echo "Modules loaded"
	
	# create a temporary job workspace of the day for 2 days
	#(workspaces may not live longer than 30 days)
	# creating a workspace for each job is an intended use, see 'man ws_allocate'

	echo "Operating inside workspace \${project_dir}"

	# singularity may not look into your storage,
	# at most, it can see your home directory
	# copy all data from sds
		# already done outside
		#{{copy_commands}}

	# run the actual application
        julia \${julia_script} \${container_path} \${workflow_configuration}
	
    echo "Completed running job \${projectname}"
	
	# copy all data from workspace back into sds
	
	cp -r \${project_workdir} \${workdir}
"""

for workdir in workdirs
	# files will be copied into here
	isdir(workdir) || mkpath(workdir)
	
	project_name = "simpleaf_" * basename(workdir)
	
	today = Dates.now() |> Date
	println("Allocating workspace: $(project_name)_$(today)")
	project_dir = read(`ws_allocate "$(project_name)_$(today)" 2`, String) |> strip
	#project_dir = match(r"(.*)\n?", project_dir) |> first
	project_workdir = joinpath(project_dir, basename(workdir))
	project_datadir = joinpath(project_dir, basename(datadir))
	project_container = joinpath(project_workdir, basename(container_pth))
	
	println("Starting to copy workdir")
	cp(workdir, project_workdir; force=true)
	println("Starting to copy datadir")
	cp(datadir, project_datadir; force=true)
	println("Starting to copy container")
	cp(container_pth, project_container; force=true)
	println("Copy operations complete")
	
	copy_operations = (workdir => project_workdir, 
					   datadir => project_datadir, 
					   container_pth => project_container)
	copy_commands = join( map(copy_operations) do (src, tar)
			"cp -r $(src) $(dirname(tar))"
	end, '\n')
	
	slurm_file_path = joinpath(project_workdir, "julia_autogen.slurm")
	config_pth = joinpath(project_workdir, "workflow_config.json")
	
	config = Dict()
	config["working_directory"] = project_workdir
	config["data_directory"] = project_datadir
	config["config_template"] = joinpath(project_datadir, "10x-feature-barcode-antibody_template.jsonnet")
	
	config_pth = joinpath(project_workdir, "workflow_config.json")
	open(config_pth, "w") do io
		JSON3.pretty(io, config)
	end
	
	batch_file = replace(slurm_template,
				"{{project_name}}" => project_name,
				"{{project_dir}}" => project_dir,
				"{{julia_script}}" => joinpath(@__DIR__, "run_simpleaf.jl"),
				"{{container_path}}" => project_container,
				"{{workflow_configuration}}" => config_pth,
				"{{project_workdir}}" => project_workdir,
				"{{workdir}}" => workdir,
				"{{copy_commands}}" => ""
	)
	
	batch_file_pth = joinpath(project_workdir,"sbatch_simpleaf.slurm")
	open(batch_file_pth, "w") do io
		print(io, batch_file)
	end
	
	slurm_output_pth = joinpath(project_workdir, "slurm-$(project_name).out")
	
	
	run(`sbatch --job-name $(project_name) -o $(slurm_output_pth) $(batch_file_pth)`)
	
end










