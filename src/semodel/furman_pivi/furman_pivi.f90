! secondary_emission.f90
module m_furman_pivi
    use m_special_functions
    use m_cephes
    use m_constants
    use m_semodel_base
    use iso_fortran_env, only: dp => real64
    use m_random_distributions
    implicit none

    private
    public :: SecondaryEmissionParams
    public :: t_FurmanPiviModel
    public :: new_FurmanPiviModel
    public :: init_params, delta_e, delta_r, delta_ts
    public :: energy_e, energy_r, energy_ts
    public :: sample_e, sample_r, sample_ts

    ! 精度を倍精度に固定
    integer, parameter :: dp = kind(1.0d0)

    type :: SecondaryEmissionParams
        ! --- 真電子 (backscattered) 用パラメータ ---
        real(dp) :: E_hat_e      ! 基準エネルギー
        real(dp) :: W            ! 幅
        real(dp) :: p            ! 指数パラメータ
        real(dp) :: P_infty_e    ! δₑ(∞)
        real(dp) :: P_hat_e      ! δₑのピーク値
        real(dp) :: e1, e2       ! 角度依存性パラメータ

        ! --- 反拡散 (rediffused) 用パラメータ ---
        real(dp) :: P_infty_r    ! δᵣ(∞)
        real(dp) :: Er           ! スケールエネルギー
        real(dp) :: r            ! 指数パラメータ
        real(dp) :: r1, r2       ! 角度依存性パラメータ

        ! --- 真二次電子 (true secondary) 用パラメータ ---
        real(dp) :: delta_hat_ts ! 真二次電子のピーク値
        real(dp) :: s            ! 補正パラメータ
        real(dp) :: E_hat_ts     ! 基準エネルギー
        real(dp) :: t1, t2, t3, t4  ! 角度依存性パラメータ

        ! --- エネルギースペクトル補正パラメータ ---
        real(dp) :: sigma_e      ! 分散調整
        real(dp) :: q            ! 反拡散分布の指数

        integer   :: M           ! 真二次電子の項数（固定10 とする）
        real(dp), allocatable :: pn(:)         ! 真二次電子の形状パラメータリスト
        real(dp), allocatable :: epsilon_n(:)  ! 真二次電子のエネルギースケールリスト
    end type SecondaryEmissionParams

    type, extends(t_SEModelBase) :: t_FurmanPiviModel
        type(SecondaryEmissionParams) :: params
    contains
        procedure :: sample_particles => fp_sample_particles
    end type

contains

    function new_FurmanPiviModel(params) result(obj)
        type(SecondaryEmissionParams), intent(in) :: params

        type(t_FurmanPiviModel) :: obj

        obj%params = params
    end function

    subroutine fp_sample_particles(self, E0, theta0, n, type, Es, dirs)
        class(t_FurmanPiviModel), intent(inout) :: self
        real(dp), intent(in) :: E0
        real(dp), intent(in) :: theta0
        integer, intent(out) :: n
        integer, intent(out)   :: type
        real(dp), intent(out) :: Es(:)
        real(dp), intent(out) :: dirs(:, :)

        real(dp) :: de, dr, dts
        real(dp) :: u

        integer :: i

        de = delta_e(self%params, e0, theta0)
        dr = delta_r(self%params, e0, theta0)
        dts = delta_ts(self%params, e0, theta0)

        call sample_uniform(u)

        u = u - de
        if (u < 0) then
            type = 1
            n = 1

            Es(1) = sample_e(self%params, E0)
            call sample_hemisphere(dirs(:, 1))
            return
        end if

        u = u - de
        if (u < 0) then
            type = 2
            n = 1

            Es(1) = sample_r(self%params, E0)
            call sample_hemisphere(dirs(:, 1))
            return
        end if

        u = u - (1 - self%params%p)**(self%params%M)
        if (u < 0) then
            n = 0
            return
        end if

        type = 3
        n = sample_binomial(self%params%M, self%params%p)

        do i = 1, n
            Es(i) = sample_ts(self%params, E0, n)
        end do

        do i = 1, n
            call sample_hemisphere(dirs(:, i))
        end do
    end subroutine

    subroutine init_params(params, pn_in, eps_in)
        type(SecondaryEmissionParams), intent(out) :: params
        real(dp), intent(in) :: pn_in(:), eps_in(:)
        params%M = size(pn_in)
        allocate (params%pn(params%M))
        allocate (params%epsilon_n(params%M))
        params%pn = pn_in
        params%epsilon_n = eps_in
    end subroutine init_params

    function delta_e(params, E0, theta0) result(delta)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in) :: E0, theta0
        real(dp) :: delta

        real(dp) :: t, d0, ang

        t = -((abs(E0 - params%E_hat_e)/params%W)**params%p)/params%p
        d0 = params%P_infty_e + (params%P_hat_e - params%P_infty_e)*exp(t)
        ang = 1.0_dp + params%e1*(1.0_dp - cos(theta0)**params%e2)
        delta = d0*ang
    end function delta_e

    function delta_r(params, E0, theta0) result(delta)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in) :: E0, theta0
        real(dp) :: delta

        real(dp) :: t, d0, ang

        t = -((E0/params%Er)**params%r)
        d0 = params%P_infty_r*(1.0_dp - exp(t))
        ang = 1.0_dp + params%r1*(1.0_dp - cos(theta0)**params%r2)
        delta = d0*ang
    end function delta_r

    function delta_ts(params, E0, theta0) result(delta)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in) :: E0, theta0
        real(dp) :: delta

        real(dp) :: ang12, ang34, Ehat, Dfunc

        ! 2番目の角度補正
        ang34 = 1.0_dp + params%t3*(1.0_dp - cos(theta0)**params%t4)
        Ehat = params%E_hat_ts*ang34
        ! D(x) = s*x/(s-1 + x^s)
        Dfunc = params%s*(E0/Ehat)/(params%s - 1.0_dp + (E0/Ehat)**params%s)
        ! ベース値
        Dfunc = params%delta_hat_ts*Dfunc
        ! 1番目の角度補正
        ang12 = 1.0_dp + params%t1*(1.0_dp - cos(theta0)**params%t2)
        delta = Dfunc*ang12
    end function delta_ts

    function energy_e(params, E, E0, theta0) result(f)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in) :: E, E0, theta0
        real(dp) :: f

        real(dp) :: aa, bb, cc, dd, sigma_mod, d0, a, b

        if (E < 0.0_dp .or. E > E0) then
            f = 0.0_dp
            return
        end if

        aa = 1.88_dp; bb = 2.5_dp; cc = 1.0e-2_dp; dd = 1.5e2_dp
        sigma_mod = (params%sigma_e - aa) + bb*(1.0_dp + tanh(cc*(E0 - dd)))

        d0 = delta_e(params, E0, theta0)

        a = 2.0_dp*exp(-((E - E0)**2)/(2.0_dp*sigma_mod**2))
        b = sqrt(2.0_dp*pi)*sigma_mod*erf(E0/(sqrt(2.0_dp)*sigma_mod))

        f = d0*a/b
    end function energy_e

    function energy_r(params, E, E0, theta0) result(f)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in) :: E, E0, theta0
        real(dp) :: f

        real(dp) :: d0, a, b

        if (E < 0.0_dp .or. E > E0) then
            f = 0.0_dp
            return
        end if

        d0 = delta_r(params, E0, theta0)

        a = (params%q + 1.0_dp)*(E**params%q)
        b = E0**(params%q + 1.0_dp)

        f = d0*a/b
    end function energy_r

    function energy_ts(params, E, E0, theta0) result(ftot)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in) :: E, E0, theta0

        real(dp) :: ftot

        integer :: n
        real(dp) :: d0, pval, Pn, a, b, c

        ftot = 0.0_dp
        if (E < 0.0_dp .or. E > E0) return

        d0 = delta_ts(params, E0, theta0)
        pval = d0/params%M

        do n = 1, params%M
            Pn = comb(params%M, n)*pval**n*(1.0_dp - pval)**(params%M - n)

            a = Pn*(E/params%epsilon_n(n))**(params%pn(n) - 1.0_dp)*exp(-E/params%epsilon_n(n))
            b = params%epsilon_n(n)*gamma(params%pn(n))*gammainc(n*params%pn(n), E0/params%epsilon_n(n))
            c = gammainc((n - 1)*params%pn(n), (E0 - E)/params%epsilon_n(n))

            ftot = ftot + n*a/b*c
        end do
    end function energy_ts

    pure function comb(n, k) result(c)
        integer, intent(in):: n, k
        real(dp) :: c

        c = gamma(real(n + 1, dp))/(gamma(real(k + 1, dp))*gamma(real(n - k + 1, dp)))
    end function comb

    function sample_e(params, E0) result(Eout)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in)                   :: E0
        real(dp)                               :: aa, bb, cc, dd
        real(dp)                               :: sigma_mod, u, arg, Eout

        aa = 1.88_dp; bb = 2.5_dp
        cc = 1.0e-2_dp; dd = 1.5e2_dp

        ! sigma の修正
        sigma_mod = (params%sigma_e - aa) + bb*(1.0_dp + tanh(cc*(E0 - dd)))

        ! 一様乱数
        call random_number(u)

        ! 逆誤差関数 erfinv() と erf() をライブラリから呼び出し
        arg = (u - 1.0_dp)*erf(E0/(sqrt(2.0_dp)*sigma_mod))
        Eout = E0 + sqrt(2.0_dp)*sigma_mod*erfinv(arg)
    end function sample_e

    function sample_r(params, E0) result(Eout)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in)                   :: E0
        real(dp)                               :: u, Eout

        call random_number(u)
        Eout = u**(1.0_dp/(params%q + 1.0_dp))*E0
    end function sample_r

    function sample_ts(params, E0, n) result(Eout)
        type(SecondaryEmissionParams), intent(in) :: params
        real(dp), intent(in)                   :: E0
        integer, intent(in)                    :: n
        real(dp)                               :: pn_val, eps, cdf, u, x, Eout

        pn_val = params%pn(n)
        eps = params%epsilon_n(n)

        ! 下側正規化不完全ガンマ関数 P(a, x)
        cdf = gammainc(pn_val, E0/eps)

        call random_number(u)
        x = u*cdf
        if (x < 1.0e-12_dp) x = 0.0_dp
        ! print *, '------', x

        ! 逆下側正規化不完全ガンマ関数 P^{-1}(a, x)
        Eout = eps*gammaincinv(pn_val, x)
    end function sample_ts

end module
