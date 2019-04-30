module functions
    implicit none
    contains
    function analitico (t_i, t_l,a, l, x) result(phi)
        implicit none
        double precision::t_i, t_l, a, l, phi, x
        phi=a*((((t_l-t_i)/l)*x)+t_i)
    end function analitico
    function dif_ter(x,t) result(d_k)
        implicit none
        double precision::d_k, x, t
        d_k=1.0+0*t+0*x
    end function dif_ter
    function dirichlet(t_i, t_l, w, i,n) result (g_c)
        implicit none
        double precision:: g_c, t_i, t_l
        integer:: i,n
        double precision, allocatable, dimension(:)::w
        if (i<n) then
            g_c=2*t_i-w(1-i)
        end if
        if (i>=n) then
            g_c=2*t_l-w(-i+2*n+1)
        end if
    end function
end module functions