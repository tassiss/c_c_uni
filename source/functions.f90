module functions
    implicit none
    contains
    function analitico (t_i, t_l, l, x) result(phi)
        implicit none
        double precision::t_i, t_l, l, phi, x
        phi=((((t_l-t_i)/l)*x)+t_i)
    end function analitico
    function dif_ter(x,t,k) result(d_k)
        implicit none
        double precision::d_k, x, t, k
        d_k=k+0*t+0*x
    end function dif_ter
    subroutine dirichlet(n_g, w, b, t_i, t_l,n)
        implicit none
        double precision:: t_i, t_l
        integer::i, n_g, n
        double precision, allocatable, dimension(:)::w, b
        do i=1-n_g,0
            w(i)=(2*t_i)-b((-i)+1)
        end do
        
        do i=n+1, n+n_g
            w(i)=(2*t_l)-b((-i)+(2*n))
        end do
   
    end subroutine 
    subroutine g_s (n,a,b,c, w, z)
        implicit none
        integer:: n, i !número de elementos, contador
        double precision:: a,b,c !constantes embutidos na matriz
        double precision, allocatable, dimension (:):: w, z !vetores de calculo
        write(*,*) n, size(z), size(w)
        do i=2,n-1
            w(i)=(1.0d0/b)*(z(i)-(a*w(i-1)+c*z(i+1)))
        end do
    end subroutine
end module
