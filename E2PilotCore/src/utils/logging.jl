using Logging, LoggingExtras
using Dates

function clean_debug_log_file(name::String = "e2pilot")
	logger_type = [Logging.Debug, Logging.Info][1]
	log_path = joinpath(k_root_dir, "log",  "$logger_type-$name.log")
    rm(log_path; force=true)
end

function get_logger(name::String; console_level=Logging.Info, dir_path::String=k_root_dir)
	timestamp_logger(logger) = TransformerLogger(logger) do log
		time_now = now()
	  	Base.merge(log, (; message = "$(time_now) $(log.message)"))
	end

	loggers = []
	logger_types = [Logging.Debug, Logging.Info]
    log_dir = joinpath(dir_path, "log", "$(today())")
    if !isdir(log_dir)
        try 
            mkdir(log_dir)
        catch e
            @warn "failed to create log dir $log_dir with error $e. But we will try to continue." 
        end
    end
	for type in logger_types
		log_path = joinpath(log_dir, "$type-$name.log")
        @info "getting logger for $name with path $log_path"
		logger_tmp = FileLogger(log_path; always_flush=true, append=true)
		logger_tmp = MinLevelLogger(logger_tmp, type)
		logger_tmp = timestamp_logger(logger_tmp)
		push!(loggers, logger_tmp)
	end

	logger1 = ConsoleLogger(stdout, console_level) |> timestamp_logger
	push!(loggers, logger1)
	logger = TeeLogger(loggers...)
	return logger
end

function log_debug()
    @debug "test" 
end
