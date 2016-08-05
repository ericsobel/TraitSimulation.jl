using ArgParse

function main(args)

    # get command line
    parsed_args = get_command_line(args)
    
    println(parsed_args)
end

# parse command line arguments
function get_command_line(args)
    
    # description of the GWAS simulation program
    description = ArgParseSettings("Simulate GWAS phenotypes")

    # command line arguments
    @add_arg_table description begin
        
        # input genotype file in PLINK format
        "--bfile"
            help = "an option"     # used by the help screen
        "--opt2"
            action = :store_true   # this makes it a flag
            help = "a flag"
    end

    parsed_args = parse_args(args, description)

    return parsed_args

end

main(ARGS)
