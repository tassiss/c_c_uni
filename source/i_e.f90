module i_e
    implicit none
    contains
    subroutine implicito (d_k)
        implicit none
        integer:: i, i_e, n_g, n, show
        double precision:: l, t_i, t_l, k , tol, dt, dx, time,rho, cp
        double precision, allocatable, dimension()::
        call entrada (n,show,l,t_i,t_l,k,tol,dt,dx, time, n_g, rho, cp, i_e)
        do i=1-n_g,n+n_g !a temperatura nas pontas da barra variam
            if (i>=1 .and. i<=n) then
                dkp=(d_k(i)+d_k(i+1))/2.0d0
                dkn=(d_k(i-1)+d_k(i))/2.0d0
            end if
        
        
            if(i<1 .or. i>n)then
                a(i,i)=1.0d0
            else
                z=(dt/((dx**2.0d0)*rho*cp))
                a(i,i)=(z*(dkp+dkn))+1
                a(i,i+1)=-z*dkp
                a(i,i-1)=-z*dkn
            end if
            write(9,*) a(i,:)
            call lapack_lin_sys(n+(2*n_g),a,b,w)
            if (all(b==w)) then
                write(*,*)'Erro'
                stop
            end if
        end do
