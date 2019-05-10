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
    subroutine dirichlet(n_g, w, k, t_i, t_l,n)
        implicit none
        double precision:: t_i, t_l
        integer::i, n_g, n
        double precision, allocatable, dimension(:)::w, k
        do i=1-n_g,0
            w(i)=(2*t_i)-k((-i)+1)
        end do
        
        do i=n+1, n+n_g
            w(i)=(2*t_l)-k((-i)+(2*n))
        end do
    end subroutine 
end module