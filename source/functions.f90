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
    subroutine dirichlet(n_g, w, b, t_i, t_l,n, i_e)
        implicit none
        double precision:: t_i, t_l
        integer::i, n_g, n, i_e
        double precision, allocatable, dimension(:)::w, b
        do i=1-n_g,0
            if (i_e==2) then
                w(i)=(2*t_i)-b((-i)+1)
            else
                w(i)=(t_i)-b((-i)+1)
            end if
        end do
        
        do i=n+1, n+n_g
            if (i_e==2) then
                w(i)=(2*t_l)-b((-i)+(2*n))
            else
                w(i)=(t_l)-b((-i)+(2*n))
            end if
        end do
   
    end subroutine 
    subroutine g_s (n,a,b,c, w, z, n_g, tol, n_i)
        implicit none
        integer:: n, i,j, n_g, n_i !número de elementos, contador
        double precision:: a,b,c, m_n1, m_n2, tol!constantes embutidos na matriz !m_n será o maior número obtido na iteraçao
        double precision, allocatable, dimension (:):: w, z !vetores de calculo
        !write(*,*) n, size(z), size(w)
        w(1-n_g:0)=z(1-n_g:0)               !celulas ghosts
        w((n+1)-(2*n_g):n-n_g)=z((n+1)-(2*n_g):n-n_g)   !células ghosts
        m_n2=1
        m_n1=0
        write(*,*) '       iteracao              ','tol'
        j=1
        do while (abs(m_n2-m_n1)>tol)
            
            m_n2=m_n1
            do i=1,n-(2*n_g)
                w(i)=(1.0d0/b)*(z(i)-(a*w(i-1)+c*w(i+1)))   
            end do
            m_n1=maxval(abs(w(1:n-2*n_g)))
            write(*,*)j, '       ',abs(m_n2-m_n1)
            j=j+1
            if(j>n_i) then
                exit
            end if
        end do
        
    end subroutine
end module
