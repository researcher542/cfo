

"""
re-implement the map function with multi-threads

gc_flag: If we want to perform garbage collection before the tasks.

verbose: If we want to print the progress.

chunk_size: we divide the whole arg_vec into multiple chunks. Useful when there is limited number of threads.
"""
function e2map(f::Base.Callable, args::AbstractArray, gc_flag::Bool = true, verbose::Bool = false; prefix::String = "", nthread::Int = -1)
    # Manually performance garbage collection, so there will be less overhead for tasks with small allocation
    gc_flag && GC.gc(false)
    if Threads.nthreads() == 1 || length(args) <= 1 || nthread == 1
        return map(f, args)
    end

    if nthread == -1
        nthread = Threads.nthreads()
    end

    ntask = length(args)
    res_vec = Any[nothing for i in 1:ntask]

    # Threads.@threads :greedy for i in eachindex(args)
    #     res_vec[i] = f(args[i])
    # end
    
    ####

     chnl = Channel{Int}(nthread, spawn=true) do ch
         for iarg in 1:ntask
             put!(ch, iarg)
         end
         ## put the end signal
         for ithread in 1:nthread
             put!(ch, -1)
         end
     end;

     Threads.@threads :greedy for ithread in 1:nthread
         while true
             itask = take!(chnl)
             if itask == -1
                 break
             end
             res_vec[itask] = f(args[itask])
             if verbose
                 msg = "$(prefix) progress of e2map: $(itask)/$ntask"
                 @info msg
             end
             # GC.safepoint()
         end
     end

    ###########

    # @show chnl

    # jj = Threads.Atomic{Int}(0)
    # l = Threads.SpinLock()
    # chunk_vec = collect(Iterators.partition(1:ntask, chunk_size))
    # Threads.@threads for i_chunk in 1:length(chunk_vec)
    # Threads.@threads for i in 1:ntask
    #     itask = take!(chnl)
    #     res_vec[i] = f(args[itask])
    #     if verbose
    #         Threads.lock(l)
    #         Threads.atomic_add!(jj, 1)
    #         msg = "$(prefix) progress of e2map: $(jj[])/$ntask"
    #         @info msg
    #         Threads.unlock(l)
    #     end
    #     GC.safepoint()
    #     # chunk_i = chunk_vec[i_chunk]
    #     # for i in chunk_i
    #     #     res_vec[i] = f(args[i])
    #     #     if verbose
    #     #         Threads.lock(l)
    #     #         Threads.atomic_add!(jj, 1)
    #     #         msg = "$(prefix) progress of e2map: $(jj[])/$ntask"
    #     #         @info msg
    #     #         Threads.unlock(l)
    #     #     end
    #     #     GC.safepoint()
    #     # end
    # end
    return [res for res in res_vec]
end

"""
map the function f to args, with self-defined multiple thread 

    verbose: if true: print the progress.
"""
function e2map_old(f::Function, args, verbose::Bool = false)
    # Manully performance garbage collection, so there will be less overhead for tesks with small allocation
    GC.gc(false)
    if Threads.nthreads() == 1 || length(args) == 1
        return map(f, args)
    end

    task_list = [@task f(arg) for arg in args]

    nthread = Threads.nthreads()
    ntask = length(task_list)
    ntask_each_thread = zeros(nthread)
    task_thread_id = zeros(Int, ntask)
    ntask_each_thread[1] = Inf

    pbar = ProgressMeter.Progress(ntask; dt=0.01, showspeed=true)

    for (itask,t) in enumerate(task_list)
        # t may now get run on any thread
        # task.stikcy means we statically set its 
        t.sticky = true
        # Dynamically schedule the task.
        for i in 2:nthread
            ntask_each_thread[i] = 0
            for id in task_thread_id
                if id != 0
                    ntask_each_thread[id] += 1
                end
            end
        end
        # @show ntask_each_thread
        if all(ntask_each_thread .>= 4) || nthread == 1
            cur_id = Threads.threadid()
            set_task_tid(t, cur_id)
            schedule(t)
            # @debug "yielding with" 
            # @show ntask_each_thread
            task_thread_id[itask] = cur_id
            yield()
            # yieldto(t)
        else
            # find the thread with minimum number of running task.
            idx = argmin(ntask_each_thread)
            task_thread_id[itask] = idx
            set_task_tid(t, idx)
            schedule(t)
        end
        # if the number of running task is larger than cpu cores, yield the current core and wait for available cores.
        # if task is done, we set its thread id to zero.
        task_thread_id[istaskdone.(task_list)] .= 0
        ntask_done = sum(istaskdone.(task_list))
        if verbose
            ProgressMeter.update!(pbar, ntask_done)
            # force to print the current progress
            flush(pbar.output)
        end
    end
    if verbose
        ProgressMeter.finish!(pbar)
        println()
    end

    res_vec = fetch.(task_list)
    return res_vec
end