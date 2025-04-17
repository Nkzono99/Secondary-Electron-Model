module m_constants
    use iso_fortran_env, only: dp => real64
    implicit none

    private
    public pi

    real(dp), parameter :: pi = acos(-1.0_dp)
end module