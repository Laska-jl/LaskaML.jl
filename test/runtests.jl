using LaskaML
using Test

@testset "LaskaML.jl" begin
    # Setup
    alpha_n(v) = (0.02 * (v - 25.0)) / (1.0 - exp((-1.0 * (v - 25.0)) / 9.0))
    beta_n(v) = (-0.002 * (v - 25.0)) / (1.0 - exp((v - 25.0) / 9.0))

    alpha_m(v) = (0.182 * (v + 35.0)) / (1.0 - exp((-1.0 * (v + 35.0)) / 9.0))
    beta_m(v) = (-0.124 * (v + 35.0)) / (1.0 - exp((v + 35.0) / 9.0))

    alpha_h(v) = 0.25 * exp((-1.0 * (v + 90.0)) / 12.0)
    beta_h(v) = (0.25 * exp((v + 62.0) / 6.0)) / exp((v + 90.0) / 12.0)

    # Examnple potassium channel
    chanK = LaskaML.hh.hhchannel(
        "K",
        35.0,
        -77.0,
        alpha_n,
        beta_n,
        4,
        alpha_h,
        beta_h,
        0
    )
    # Example sodium channel
    chanNa = LaskaML.hh.hhchannel(
        "Na",
        40.0,
        55.0,
        alpha_m,
        beta_m,
        3,
        alpha_h,
        beta_h,
        1
    )
    # Example leaky channel
    chanL = LaskaML.hh.hhchannel(
        "L",
        0.3,
        -65.0
    )

    @test chanK isa LaskaML.hh.ChannelType
    @test chanNA isa LaskaML.hh.ChannelType
    @test chanL isa LaskaML.hh.ChannelType

    LaskaML.hh.parsechannel(chan1)
end
