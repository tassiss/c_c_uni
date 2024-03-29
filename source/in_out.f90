module in_out
    implicit none
    contains
    !#############################################################################################
    subroutine entrada (n,show,l,t_i,t_l,k,tol,dt,dx, time, n_g, rho, cp, i_e, n_i)
        implicit none
        integer::n, show, input_id=22, status=0, n_g, i_e, n_i
        double precision:: l, t_i,t_l,k, tol, dt, dx, time, cfl, rho, cp
        open(unit=input_id, action='read', status='old', iostat=status, file='../input.dat')!chamada dos dados de entrada
        read(input_id, *)l !lê cada linha do arquivo e associa com uma variavel de entrada
        read(input_id, *)n  !
        read(input_id, *)t_i    
        read(input_id, *)t_l    !
        read(input_id, *)k  ! difusividade térmica
        read(input_id, *)tol    !
        read(input_id, *)cfl !
        read(input_id, *)time  !
        read(input_id, *)show   !
        read(input_id,*) n_g
        read(input_id,*) rho
        read(input_id,*) cp
        read(input_id,*) i_e
        read(input_id,*) n_i
        close(input_id)
        dx=l/(real(n-1)) !distância entre o nós
        dt=cfl*(dx**2.0d0/(k/(rho*cp)))
        write(*,*) dx, dt 
    end subroutine entrada
    !################################################################################################
    subroutine saida(t,dt,x,w,cont, erro, n_g)
        implicit none
        integer, intent(in)::n_g
        integer::status=0, cont, j,  fid=20,n, show, i_e, n_i
        double precision::x,t_l,t_i, dx, k, tol, dt, l, time, rho, cp
        double precision, intent(in)::t, erro
        double precision,dimension(:)::w
        character*2048::out,name,file
        call entrada(n,show,l,t_i,t_l,k,tol,dt,dx,time, n_g, rho, cp, i_e, n_i)
        out='../output/'
        
        call system ('mkdir '//out)
        
        write(file,'(I0.10)') cont !tranforma cont em character
        name=trim(out)//trim(adjustl(file))//'.dat' !cria o nome do arquivo
        open(unit=fid, file=Trim((name)), action='write', status='unknown', iostat=status) !abre o arquivo
        x=0.0d0
        j=n_g+1
        do while (j<=n+n_g)!varre todo o vetor, inclusive as pontas
            write(fid,*)t, x, w(j), cont, erro !escreve tempo, posição do nó e sua respectiva temperatura
            j=j+1
            x=x+dx
        end do
        close(fid, iostat=status)
    end subroutine saida

end module in_out
