using SnoopPrecompile

@precompile_setup begin
    @precompile_all_calls begin
        init()
        # net = readoldmap("oldmap")
    end
end