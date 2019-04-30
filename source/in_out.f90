module in_out
    implicit none
    contains
    !#############################################################################################
    subroutine entrada (n,show,l,t_i,t_l,a,tol,dt,dx, time, n_g)
        implicit none
        integer::n, show, input_id=22, status=0, n_g
        double precision:: l, t_i,t_l,a, tol, dt, dx, time, cfl
        open(unit=input_id, action='read', status='old', iostat=status, file='../input.dat')!chamada dos dados de entrada
        read(input_id, *)l !lê cada linha do arquivo e associa com uma variavel de entrada
        read(input_id, *)n  !
        read(input_id, *)t_i    !
        read(input_id, *)t_l    !
        read(input_id, *)a  !
        read(input_id, *)tol    !
        read(input_id, *)cfl !
        read(input_id, *)time  !
        read(input_id, *)show   !
        read(input_id,*) n_g
        close(input_id)
        dx=l/(real(n-1)) !distância entre o nós
        dt=cfl*(dx**2/a**2)
    end subroutine entrada
    !################################################################################################
    subroutine saida(t,dt,x,w,cont, erro)
        implicit none
        integer::status=0, cont, j=1,  fid=20,n, show, n_g
        double precision::x,t_l,t_i, dx, a, tol, dt, l, time
        double precision, intent(in)::t, erro
        double precision,dimension(:)::w
        character*2048::out,name,file
        call entrada(n,show,l,t_i,t_l,a,tol,dt,dx,time, n_g)
        out='../output/'
        write(file,'(I0.10)') cont !tranforma cont em character
        name=trim(out)//trim(adjustl(file))//'.dat' !cria o nome do arquivo
        open(unit=fid, file=Trim((name)), action='write', status='unknown', iostat=status) !abre o arquivo
        x=0.00000
        do while (j-1/=n)!varre todo o vetor, inclusive as pontas
            write(fid,*)t, x, w(j), cont, erro !escreve tempo, posição do nó e sua respectiva temperatura
            w(n)=t_l
            w(1)=t_i
            j=j+1
            x=x+dx
        end do
        j=1
        close(fid, iostat=status)
    end subroutine saida

end module in_out
