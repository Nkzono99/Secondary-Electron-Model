module m_random_distributions
    use iso_fortran_env, only: dp => real64
    implicit none

    private
    public :: set_rng
    public :: sample_uniform
    public :: sample_normal
    public :: sample_binomial
    public :: sample_hemisphere

    ! 手続きインターフェースの定義
    abstract interface
        subroutine random_uniform_interface(u)
            import dp
            real(dp), intent(out) :: u
        end subroutine
    end interface

    ! 手続きポインタ変数（初期値は intrinsic の random_number）
    procedure(random_uniform_interface), pointer :: rng => default_rng

    ! 正規乱数キャッシュ用
    logical :: have_saved_normal = .false.
    real(dp) :: saved_normal

contains

    subroutine default_rng(u)
        real(dp), intent(out) :: u

        call random_number(u)
    end subroutine

    ! RNG を差し替えるためのサブルーチン
    subroutine set_rng(proc)
        procedure(random_uniform_interface), pointer, intent(in) :: proc
        rng => proc
    end subroutine

    subroutine sample_uniform(u)
        real(dp), intent(out) :: u

        call rng(u)
    end subroutine

    ! μ, σ 指定可能な正規乱数サンプラー
    subroutine sample_normal(mu, sigma, x)
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: sigma
        real(dp), intent(out) :: x

        real(dp) :: u, v
        real(dp) :: s
        real(dp) :: factor
        real(dp) :: z0, z1

        ! すでに保存分があればそれを返す
        if (have_saved_normal) then
            x = mu + sigma*saved_normal
            have_saved_normal = .false.
            return
        end if

        ! Marsaglia の極座標法 本体
        do
            call rng(u); u = 2.0_dp*u - 1.0_dp
            call rng(v); v = 2.0_dp*v - 1.0_dp

            s = u*u + v*v

            if (s > 0.0_dp .and. s < 1.0_dp) exit
        end do

        factor = sqrt(-2.0_dp*log(s)/s)
        z0 = u*factor
        z1 = v*factor

        ! 片方を返し、もう片方を保存
        x = mu + sigma*z0

        saved_normal = z1
        have_saved_normal = .true.
    end subroutine

    !> Binomial(n, p) の乱数サンプルを返す関数
    !! n: 試行回数
    !! p: 成功確率 (0<=p<=1)
    !! 戻り値: 成功回数 (0..n)
    function sample_binomial(n, p) result(k)
        integer, intent(in)    :: n
        real(dp), intent(in)   :: p
        integer                :: k, i
        real(dp)               :: u

        k = 0
        do i = 1, n
            call rng(u)
            if (u < p) k = k + 1
        end do
    end function

    subroutine sample_hemisphere(dir)
        real(dp), intent(out) :: dir(3)

        real(dp) :: phi, theta
        real(dp) :: u1, u2
        real(dp) :: mu
        real(dp) :: sin_t

        call sample_uniform(u1)
        phi = 2.0_dp*acos(-1.0_dp)*u1

        call sample_uniform(u2)
        mu = sqrt(u2)

        theta = acos(mu)

        sin_t = sin(theta)
        dir(1) = sin_t*cos(phi)
        dir(2) = sin_t*sin(phi)
        dir(3) = mu
    end subroutine

end module
