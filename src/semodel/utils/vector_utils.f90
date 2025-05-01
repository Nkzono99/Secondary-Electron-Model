module m_vector_utils
    implicit none
    integer, parameter :: dp = kind(1.0d0)

contains

    !---------------------------------------------------
    ! 3成分ベクトルの長さを返す関数
    !---------------------------------------------------
    pure function norm3(x) result(n)
        real(dp), intent(in) :: x(3)
        real(dp)             :: n
        n = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
    end function norm3

    !---------------------------------------------------
    ! 3成分ベクトルの外積
    !---------------------------------------------------
    pure subroutine cross3(a, b, c)
        real(dp), intent(in)  :: a(3), b(3)
        real(dp), intent(out) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine cross3

    !---------------------------------------------------
    ! エネルギー E[eV]，偏角 theta(1)=polar, theta(2)=azimuth，
    ! 法線ベクトル n(3) から速度ベクトル v(3) を返すサブルーチン
    !
    ! E は電子の運動エネルギー（eV），
    ! 速度大きさは v = sqrt(2 E e / m_e) で計算
    !---------------------------------------------------
    subroutine compute_velocity(E, dir, n, qm, v)
        real(dp), intent(in)  :: E         ! eV
        real(dp), intent(in)  :: dir(3)  ! [polar, azimuth] in rad
        real(dp), intent(in)  :: n(3)      ! 法線ベクトル
        real(dp), intent(in)  :: qm
        real(dp), intent(out) :: v(3)      ! 出力速度ベクトル
        !--- 物理定数 ---
        real(dp), parameter :: m_e = 9.10938356e-31_dp  ! 電子質量 [kg]
        real(dp), parameter :: e_c = 1.602176634e-19_dp ! 電子電荷 [C]

        real(dp) :: vmag       ! 速度の大きさ [m/s]
        real(dp) :: un(3)      ! 単位法線ベクトル
        real(dp) :: r(3), t1(3), t2(3) ! 法線に直交する基底ベクトル
        real(dp) :: len

        !---- 速度の大きさ計算 ----
        vmag = sqrt(2.0_dp*E*qm)

        !---- 法線ベクトルを単位化 ----
        len = norm3(n)
        if (len <= 0.0_dp) then
            stop "Error: zero-length normal vector"
        end if
        un = n/len

        !---- 法線に直交する任意のベクトル r を選択 ----
        if (abs(un(1)) < abs(un(2))) then
            r = [1.0_dp, 0.0_dp, 0.0_dp]
        else
            r = [0.0_dp, 1.0_dp, 0.0_dp]
        end if

        !---- t1 = r × un を単位化 ----
        call cross3(r, un, t1)
        t1 = t1/norm3(t1)

        !---- t2 = un × t1 ----
        call cross3(un, t1, t2)

        !---- 方向ベクトル dir を作成 ----
        ! dir = sinθ cosφ t1 + sinθ sinφ t2 + cosθ un
        v = vmag*(dir(1)*t1 + dir(2)*t2 + dir(3)*un)
    end subroutine compute_velocity

end module
