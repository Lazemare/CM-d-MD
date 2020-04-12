""" 
    inputparser(configfile::String)

The input file parser of CMdMD.
"""
function inputparser(configfile::String)

    system = MDsystem()
    # Initialize the system
    system.namdconfigurefile      = ""
    system.iter                   = 1
    system.steps                  = 1
    system.dcdfreq                = 1
    system.nproc                  = 1
	system.gpu                    = -1
	system.orderparam             = ""
    system.selec                  = "all"
	system.alignref               = ""
	system.refframe               = 0
    system.rmsdprecision          = "atoms"
    system.rmsdregionselection    = "all"
    system.rmsdstructureselection = "all"
    system.disprecision           = "atoms"
    system.disregionselection     = "all"
	system.disstructureselection  = "all"
	system.posprecision           = "atoms"
    system.posregionselection     = "all"
    system.posstructureselection  = "all"
    system.stretchprecision       = "atoms"
    system.stretchfactor          = 1.10
    system.freq                   = 1
    system.commutemap             = true
    system.lag                    = 1

    # Read cmdmd configure file.
    config = readdlm(configfile)
    for i = 1:size(config,1)
        if config[i,1] |> lowercase == "namdconfigurefile"
            system.namdconfigurefile = config[i,2] |> String
        elseif config[i,1] |> lowercase == "iterations"
            system.iter = config[i,2]
        elseif config[i,1] |> lowercase == "steps"
            system.steps = config[i,2]
        elseif config[i,1] |> lowercase == "dcdfreq"
            system.dcdfreq = config[i,2]
        elseif config[i,1] |> lowercase == "frequency"
            system.freq = config[i,2]
        elseif config[i,1] |> lowercase == "nproc"
            system.nproc = config[i,2]
        elseif config[i,1] |> lowercase == "gpu"
            system.gpu = config[i,2]
        elseif config[i,1] |> lowercase == "selection"
			system.selec = config[i,2] |> String |> lowercase
		elseif config[i,1] |> lowercase == "orderparameter"
			system.orderparam = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "rmsdprecision"
            system.rmsdprecision = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "rmsdregionselection"
            system.rmsdregionselection = config[i,2] |> String
        elseif config[i,1] |> lowercase == "rmsdstructureselection"
            system.rmsdstructureselection = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "referenceframe"
            system.refframe = config[i,2]
        elseif config[i,1] |> lowercase == "disprecision"
            system.disprecision = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "disregionselection"
            system.disregionselection = config[i,2] |> String
        elseif config[i,1] |> lowercase == "disstructureselection"
			system.disstructureselection = config[i,2] |> String |> lowercase
		elseif config[i,1] |> lowercase == "posprecision"
            system.posprecision = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "posregionselection"
            system.posregionselection = config[i,2] |> String
        elseif config[i,1] |> lowercase == "posstructureselection"
            system.posstructureselection = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "alignreference"
            system.alignref = config[i,2] |> String
        elseif config[i,1] |> lowercase == "stretchprecision"
            system.stretchprecision = config[i,2] |> String |> lowercase
        elseif config[i,1] |> lowercase == "stretchfactor"
            system.stretchfactor = config[i,2]
		elseif config[i,1] |> lowercase == "commutemap"
			if (typeof(config[i,2]) == String && lowercase(config[i,2]) == "true")
				system.commutemap = true
			else
				system.commutemap = config[i,2]
			end
        elseif config[i,1] |> lowercase == "lagtime"
            system.lag = config[i,2]
        end
	end
	
    # Read and check initial NAMD configure file.
    system.initcoord, system.topo = read_conf_namd(system.namdconfigurefile)
    # Make some checks.
    if system.namdconfigurefile == ""
        error("You must assign the NAMD configure file name!")
	end
	if (system.orderparam != "rmsd" && system.orderparam != "distance" &&
		system.orderparam != "position")
		error("You must assign which order parameter to use!")
	end
	if (system.rmsdprecision != "atoms" && system.rmsdprecision != "residues")
		error("Acceptable values of `RMSDPrecision` are: `atoms` and `residues`.")
	end
	if (system.disprecision != "atoms" && system.disprecision != "residues")
		error("Acceptable values of `DisPrecision` are: `atoms` and `residues`.")
	end
	if (system.posprecision != "atoms" && system.posprecision != "residues")
		error("Acceptable values of `PosPrecision` are: `atoms` and `residues`.")
	end
	if system.iter == 1
		@warn "Only 1 iteration of CMdMD will be performed, which make no sense.
		Did you just forget to set the `Iterations` value?"
	end
	if system.steps == 1
		@warn "Only 1 step of MD simulation will be performed, which make no sense.
		Did you just forget to set the `Steps` value?"
	end
	if system.lag == 1
		@warn "The lag time has been set to 1, which is a very small value.
		Did you just forget to set the `LagTime` value?"
	end
	if ((system.orderparam == "rmsd" || system.orderparam == "position") &&
		system.alignref == "")
		@warn "You haven't assigned `AlignReference`. This value has been set to `all`."
		system.alignref = "all"
	end
	if (system.orderparam == "distance" && system.disprecision == "atoms")
		@warn "You are using distance between atoms in your selection as the order parameter! 
		This could be VERY slow!"
    end
    if (system.stretchprecision != "atoms" && system.stretchprecision != "residues")
		error("Acceptable values of `StretchPrecision` are: `atoms` and `residues`.")
	end
    if system.stretchfactor <= 1.0
        @warn "Setting `StretchFactor` smaller than 1.0 will lead no acceleration!"
    end
    if system.stretchfactor > 1.3
        @warn "Setting `StretchFactor` larger than 1.3 could lead bad dynamics."
    end
    
return system
end
