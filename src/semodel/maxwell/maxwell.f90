module m_maxwell
    use m_semodel_base
    use iso_fortran_env, only: dp => real64
    use m_random_distributions
    implicit none

    private

    public t_MaxwellModel
    public new_MaxwellModel

    type, extends(t_SEModelBase) :: t_MaxwellModel
        real(dp) :: yield
        real(dp) :: mu
        real(dp) :: sigma
    contains
        procedure :: sample_particles => mx_sample_particles
    end type

contains

    function new_MaxwellModel(yield, mu, sigma) result(obj)
        real(dp), intent(in) :: yield
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: sigma
        type(t_MaxwellModel) :: obj

        obj%yield = yield
        obj%mu = mu
        obj%sigma = sigma
    end function

    subroutine mx_sample_particles(self, E0, theta0, n, type, Es, dirs)
        class(t_MaxwellModel), intent(inout) :: self
        real(dp), intent(in) :: E0
        real(dp), intent(in) :: theta0
        integer, intent(out) :: n
        integer, intent(out)   :: type
        real(dp), intent(out) :: Es(:)
        real(dp), intent(out) :: dirs(:, :)

        real(dp) :: u

        call sample_uniform(u)

        if (u > self%yield) then
            n = 0
            return
        end if

        n = 1
        type = 1
        call sample_normal(self%mu, self%sigma, Es(1))
        call sample_hemisphere(dirs(:, 1))
    end subroutine

end module
