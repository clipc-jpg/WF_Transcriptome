
import Logging
import JSON
import Dates

mutable struct JSONLogger{IO1<:IO,IO2<:IO, 
                   TLL<:Union{Base.CoreLogging.LogLevel,Int64}}<:Logging.AbstractLogger
    io::IO1
    ioerr::IO2
    min_loglevel::TLL
    flush_size::Int64
    flush_counter::Int64
    lock::ReentrantLock
    print_to_screen::Bool
    function JSONLogger(;io=stdout,ioerr=stderr,min_loglevel=Logging.Info,flush_size=1,reelock=ReentrantLock(),print_to_screen=true)
        IO1 = typeof(io)
        IO2 = typeof(ioerr)
        TLL = typeof(min_loglevel)
        return new{IO1,IO2,TLL}(io, ioerr, min_loglevel, flush_size, zero(flush_size), reelock,print_to_screen)
    end
end

function Logging.handle_message(jsl::JSONLogger, level, message, 
                                _module, group, id, 
                                file, line; 
                                kwargs...)
    data = Dict(:timestamp => Dates.now(), :loglevel => level, 
                    :message => message, kwargs...)
    lock(jsl.lock) do
        jsl.print_to_screen && println(stdout, message)
        if level < Logging.Error
            JSON.Writer.print(jsl.io, data)
            println(jsl.io)
        else
            JSON.Writer.print(jsl.ioerr, data)
            println(jsl.ioerr)
        end
        jsl.flush_counter += 1
        if jsl.flush_counter >= jsl.flush_size
            flush(jsl.io)
            jsl.io != jsl.ioerr && flush(jsl.ioerr)
            jsl.flush_counter = zero(jsl.flush_counter)
        end
    end
end

function Logging.shouldlog(jsl::JSONLogger, level, _module, group, id)
    return level >= jsl.min_loglevel
end

function Logging.min_enabled_level(jsl::JSONLogger)
    return jsl.min_loglevel
end


function wfdebug(event, message=event; operation=event, inputs=(), outputs=(), exitcode=0, attempt_num=-1)
    @debug message event=event operation=operation inputs=inputs outputs=outputs exitcode=exitcode attempt_num=attempt_num
end

function wfinfo(event, message=event; operation=event, inputs=(), outputs=(), exitcode=0, attempt_num=-1)
    @info message event=event operation=operation inputs=inputs outputs=outputs exitcode=exitcode attempt_num=attempt_num
end

function wfwarn(event, message=event; operation=event, inputs=(), outputs=(), exitcode=0, attempt_num=-1)
    @warn message event=event operation=operation inputs=inputs outputs=outputs exitcode=exitcode attempt_num=attempt_num
end

function wferror(event, message=event; operation=event, inputs=(), outputs=(), exitcode=1, attempt_num=-1)
    @error message event=event operation=operation inputs=inputs outputs=outputs exitcode=exitcode attempt_num=attempt_num
end


struct WFLogPipe{F,TI,TO}<:IO
    log_function::F
    operation::String
    inputs::TI
    outputs::TO
    attempt_num::Int64
    
    function WFLogPipe(log_function; operation, inputs=(), outputs=(), attempt_num=-1)
        return new{typeof(log_function), typeof(inputs), typeof(outputs)}(log_function, operation, inputs, outputs, attempt_num)
    end
end

function Base.write(wflp::WFLogPipe, message)
    (; operation, inputs, outputs, attempt_num) = wflp
    wflp.log_function(message; operation, inputs, outputs, attempt_num)
end

function Base.write(wflp::WFLogPipe, message_io::IO)
    (; operation, inputs, outputs, attempt_num) = wflp
    message = read(message_io, String)
    wflp.log_function(message; operation, inputs, outputs, attempt_num)
end

struct WFLogStump{F}
    log_function::F
end

function WFLogPipe(wfls::WFLogStump, operation, inputs, outputs, attempt_num)
    return WFLogPipe(wfls.log_function; operation, inputs, outputs, attempt_num)
end

struct IndirectWFLogStump{F}
    log_function::F
    dump_file::String
    function IndirectWFLogStump(log_function::F; parentdir=pwd(), cleanup=true) where {F}
        dump_file = tempname(parentdir; cleanup)
        return new{F}(log_function, dump_file)
    end
    function IndirectWFLogStump(log_function::F, dump_file; cleanup=false) where {F}
        cleanup && Base.Filesystem.temp_cleanup_later(dump_file)
        return new{F}(log_function, dump_file)
    end
end

function WFLogPipe(iwfls::IndirectWFLogStump, operation, inputs, outputs, attempt_num)
    return WFLogPipe(wfls.log_function; operation, inputs, outputs, attempt_num)
end






