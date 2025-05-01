module m_secondary_electron_parameters
    use m_furman_pivi
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    private

    public set_copper
    public set_stainless_steel

contains

    !> Copperパラメータで初期化する
    subroutine set_copper(params)
        type(SecondaryEmissionParams), intent(out) :: params

        integer, parameter :: M_CU = 10

        real(dp), dimension(M_CU), parameter :: pn_vals = (/ &
                                                2.5_dp, 3.3_dp, 2.5_dp, 2.5_dp, 2.8_dp, &
                                                1.3_dp, 1.5_dp, 1.5_dp, 1.5_dp, 1.5_dp/)

        real(dp), dimension(M_CU), parameter :: eps_vals = (/ &
                                                1.5_dp, 1.75_dp, 1.0_dp, 3.75_dp, 8.5_dp, &
                                                11.5_dp, 2.5_dp, 3.0_dp, 2.5_dp, 3.0_dp/)

        ! 真二次電子リストの確保とコピー
        params%M = M_CU
        allocate (params%pn(M_CU))
        allocate (params%epsilon_n(M_CU))
        params%pn = pn_vals
        params%epsilon_n = eps_vals

        ! Backscattered electrons
        params%P_infty_e = 0.02_dp
        params%P_hat_e = 0.496_dp
        params%E_hat_e = 0.0_dp
        params%W = 60.86_dp
        params%p = 1.0_dp
        params%e1 = 0.26_dp
        params%e2 = 2.0_dp
        params%sigma_e = 2.0_dp

        ! Rediffused electrons
        params%P_infty_r = 0.2_dp
        params%Er = 0.041_dp
        params%r = 0.104_dp
        params%r1 = 0.26_dp
        params%r2 = 2.0_dp
        params%q = 0.5_dp

        ! True-secondary electrons
        params%delta_hat_ts = 1.8848_dp
        params%E_hat_ts = 276.8_dp
        params%s = 1.54_dp
        params%t1 = 0.66_dp
        params%t2 = 0.8_dp
        params%t3 = 0.7_dp
        params%t4 = 1.0_dp
    end subroutine

    !> Stainless steelパラメータで初期化する
    subroutine set_stainless_steel(params)
        type(SecondaryEmissionParams), intent(out) :: params
        integer, parameter :: M_SS = 10
        real(dp), dimension(M_SS), parameter :: pn_vals = (/ &
                                                1.6_dp, 2.0_dp, 1.8_dp, 4.7_dp, 1.8_dp, &
                                                2.4_dp, 1.8_dp, 1.8_dp, 2.3_dp, 1.8_dp/)
        real(dp), dimension(M_SS), parameter :: eps_vals = (/ &
                                                3.9_dp, 6.2_dp, 13.0_dp, 8.8_dp, 6.25_dp, &
                                                2.25_dp, 9.2_dp, 5.3_dp, 17.8_dp, 10.0_dp/)

        params%M = M_SS
        allocate (params%pn(M_SS))
        allocate (params%epsilon_n(M_SS))
        params%pn = pn_vals
        params%epsilon_n = eps_vals

        ! Backscattered electrons
        params%P_infty_e = 0.07_dp
        params%P_hat_e = 0.5_dp
        params%E_hat_e = 0.0_dp
        params%W = 100.0_dp
        params%p = 0.9_dp
        params%e1 = 0.26_dp
        params%e2 = 2.0_dp
        params%sigma_e = 1.9_dp

        ! Rediffused electrons
        params%P_infty_r = 0.74_dp
        params%Er = 40.0_dp
        params%r = 1.0_dp
        params%r1 = 0.26_dp
        params%r2 = 2.0_dp
        params%q = 0.4_dp

        ! True-secondary electrons
        params%delta_hat_ts = 1.22_dp
        params%E_hat_ts = 310.0_dp
        params%s = 1.813_dp
        params%t1 = 0.66_dp
        params%t2 = 0.8_dp
        params%t3 = 0.7_dp
        params%t4 = 1.0_dp
    end subroutine

end module
