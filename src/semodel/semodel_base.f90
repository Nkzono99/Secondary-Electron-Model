module m_semodel_base
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    private
    public :: t_SEModelBase

    type, abstract :: t_SEModelBase
    contains
        procedure(sample_particles_interface), deferred :: sample_particles
    end type t_SEModelBase

    interface
        subroutine sample_particles_interface(self, E0, theta0, n, type, Es, dirs)
            import :: t_SEModelBase, dp
            class(t_SEModelBase), intent(inout) :: self
            real(dp), intent(in) :: E0       ! 入射エネルギー
            real(dp), intent(in) :: theta0   ! 入射角度
            integer, intent(out) :: n        ! 生成粒子数
            integer, intent(out) :: type ! 粒子タイプ
            real(dp), intent(out) :: Es(:) ! 速度ベクトル (vx,vy,vz)×n
            real(dp), intent(out) :: dirs(:, :) ! 速度ベクトル (vx,vy,vz)×n
        end subroutine
    end interface

end module m_semodel_base
