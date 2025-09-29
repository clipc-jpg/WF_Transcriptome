
# WF_Transcriptome
In this repository the epi2me workflow wf_transcriptome has been refactored into a simple script and containerized into singularity suitable for use with the Colony Launcher Application.

## Build Steps
The build steps have been straamlined into a custom cargo subcommand as defined in the .cargo directory and implemented in the build_tasks module.
You first need to install rust, which comes with the cargo package manager. Install the dioxus bundler with "cargo install dioxus-cli".
This project was created with Dioxus Version 6.3. The dioxus-cli version will likely need to match.
Then you should be able to create a fully functional container with "cargo containerize".

The full cycle consists of:
	1. Creating a sif file with all dependencies installed (which takes the majority of build time)
	   You do this in the container directory by calling "sudo singularity build wf_deps.sif wf_deps.def"
	2. Compiling the selfconfigurator. At the top-level directory run "dx build -p self_configurator --release".
	3. Creating the fully functional container. In the container directory call "sudo singularity build wf_apps.sif wf_apps.def"



    
    
