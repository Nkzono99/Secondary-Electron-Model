module m_cephes
    use iso_c_binding, only: c_double
    implicit none
    private

    public :: gammainc
    public :: gammaincinv

    interface
        function c_igam(a, x) bind(C, name="igam")
            import :: c_double
            real(c_double), value :: a, x
            real(c_double)         :: c_igam
        end function
        function c_igami(a, y) bind(C, name="igami")
            import :: c_double
            real(c_double), value :: a, y
            real(c_double)         :: c_igami
        end function
    end interface

contains

    function gammainc(a, x) result(y)
        real(c_double), intent(in) :: a, x
        real(c_double)             :: y
        y = c_igam(a, x)
    end function

    function gammaincinv(a, y) result(x)
        real(c_double), intent(in) :: a, y
        real(c_double)             :: x
        x = c_igami(a, 1 - y)
    end function

end module
