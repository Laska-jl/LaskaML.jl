using LaskaML
using Test

@testset "HH channels and model" begin
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
        alpha_m=alpha_n,
        beta_m=beta_n,
        p=4,
        alpha_h=alpha_h,
        beta_h=beta_h,
        q=0
    )
    # Example sodium channel
    chanNa = LaskaML.hh.hhchannel(
        "Na",
        40.0,
        55.0,
        alpha_m=alpha_m,
        beta_m=beta_m,
        p=3,
        alpha_h=alpha_h,
        beta_h=beta_h,
        q=1
    )
    # Example leaky channel
    chanL = LaskaML.hh.hhchannel(
        "L",
        0.3,
        -65.0
    )

    @test LaskaML.hh.beta_m(chanK) == beta_n
    @test LaskaML.hh.alpha_m(chanK) == alpha_n
    @test LaskaML.hh.alpha_h(chanK) == alpha_h
    @test LaskaML.hh.beta_h(chanK) == beta_h

    @test LaskaML.hh.alpha_m(chanL)() == 1
    @test LaskaML.hh.beta_m(chanL)() == 1
    @test LaskaML.hh.alpha_h(chanL)() == 1
    @test LaskaML.hh.beta_h(chanL)() == 1

    @test LaskaML.hh.parsechannel(chanK) == "(g_K*(m_K^4)*(v-E_K))"
    @test LaskaML.hh.parsechannel(chanNa) == "(g_Na*(m_Na^3)*(h_Na^1)*(v-E_Na))"
    @test LaskaML.hh.parsechannel(chanL) == "(g_L*(v-E_L))"

    model = LaskaML.hh.HHModel(
        [chanK]
    )

    model += chanNA

end
