

@testset "utility functions" begin

    N = 5
    c1 = FerriteIGA.diagonal_beo(N)
    c2 = FerriteIGA.diagonal_beo(N)

    _diagonalmatrix(N) = [ i==j ? 1.0 : 0.0 for i in 1:N, j in 1:N]

    @test FerriteIGA.beo2matrix(c1) == _diagonalmatrix(N)

    c3 = FerriteIGA.combine_beo(c1, c2)
    @test FerriteIGA.beo2matrix(c3) == _diagonalmatrix(2N)

end