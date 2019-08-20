subroutine g_s (n,a,b,c, w, z)
    implicit none
    integer:: n, i !número de elementos, contador
    double precision:: a,b,c !constantes embutidos na matriz
    double precision, allocatable, dimension (:):: w, z !vetores de calculo
    write(*,*) n
    do i=3,n-2
        w(i)=(1.0d0/b)*(z(i)-(a*w(i-1)+c*z(i+1)))
    end do
end subroutine  