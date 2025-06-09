using Test
using E2PilotCore: subregion_poly, isinHK, isinchina

@testset "geometry" begin
	g_hk_poly = subregion_poly(:china, "香港特别行政区")
	@test isinHK(22.3367, 114.1724) # city u
	@test isinHK(22.4196, 114.2068) # CUHK
	@test !isinHK(22.5429, 114.0596) # shenzheng
	@test isinchina(22.5429, 114.0596) # shenzheng
    latlon = (36.7783, -119.4179)  # California
	@test !isinchina(latlon...) # California
    state = ep.latlon2subregion(latlon..., :us)
    @test state == "California"
end