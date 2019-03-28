module conducao
implicit none

contains
subroutine ent_cond(l,n,t_i,t_l,a,tol,dt,t, k, w, dx, show)
    implicit none
    integer::n, status=0, input_id=21, show
    double precision, allocatable, dimension(:):: k, w
        double precision:: l, t_i, t_l, t, dt, dx, tol, a
    open(unit=input_id, action='read', status='old', iostat=status, file='../input.dat')!chamada dos dados de entrada
    read(input_id, *)l !lê cada linha do arquivo e associa com uma variavel de entrada
    read(input_id, *)n  !
    read(input_id, *)t_i    !
    read(input_id, *)t_l    !
    read(input_id, *)a  !
    read(input_id, *)tol    !
    read(input_id, *)dt !
    read(input_id, *)t  !
    read(input_id, *)show   !
    close(input_id) !fecha o arquivo para economizar memória
    allocate(k(n),stat=status) !aloca o vetor k utilizado para calcular a temperatura
    allocate(w(n), stat=status) !aloca o vetor w que é o vetor temperatura
    k=t_l !todos os pontos do vetor são equivalente a t_l
    k(1)=t_i ! o primeiro nó recebe a temperatura t_i
    dx=l/(real(n-1)) !distância entre o nós
    w=k
end subroutine 
subroutine calc_cond (n,t_i,t_l,a,tol,dt,t, k, w, dx, fid, d, i, show)
implicit none
integer::n, fid,i, show, j=1, status=0, cont=1
double precision, allocatable,dimension(:)::w,k,dif
double precision:: d,x,dx, tol, t, dt, a, t_l, t_i
character*2048::file, name,out
    !########################################################################################
    !método para criação de arquivo output
    out='../output/'
    write(file,'(I0.8)') cont !tranforma cont em character
    name=trim(out)//trim(adjustl(file))//'.dat' !cria o nome do arquivo
    open(unit=22, file=Trim((name)), action='write', status='unknown', iostat=status) !abre o arquivo
    x=0.00000
    do while (j<=n)!varre todo o vetor, inclusive as pontas
    write(22,*)t, x, w(j) !escreve tempo, posição do nó e sua respectiva temperatura
    w(j)=t_l
    w(1)=t_i
    j=j+1
    x=x+dx
    end do
    j=1
    close(22, iostat=status)
    !fim
    !#########################################################################################
    do
        do while(i/=n) !enquanto i for menor que n:
            w(i)=((dt*a*(k(i+1)-(2*(k(i)))+k(i-1)))/dx**2)+k(i) ! o vetor w é igua a isso tudo que envolve o vetor k
            i=i+1             
        end do
        i=2
        dif=w-k !vetor diferença
        d=(sum(dif))/n !soma dos elementos do vetor diferença
        k=w ! k igual w diz que na próxima entrada os valores de k serão iguais os valores calculados para w nessa saída
        write(fid,*)w
        if(mod(cont,show)==0 .or. abs(d)<tol)then !se o resto da divisão de cont/show for exata, se a diferença entre linhas for menor que a tolerância estabecida & o primeiro calculo
            write(file,'(I0.8)') cont
            name=trim(out)//trim(adjustl(file))//'.dat'
            open(unit=22, file=(name), action='write', status='unknown', iostat=status)
            x=0.00000
            do while (j-1/=n)
                write(22,*)t, x, w(j)
                j=j+1       
                x=x+dx
            end do
            j=1
            close(22, iostat=status)
        end if
        x=dx
        if(abs(d)<tol)then !se a soma foir igual a zero, quero dizer que a peóxima linha i é igual a linha i-1
            exit !encerra o loop e mostra o tempo necessario para estabilizar
        end if
        cont=cont+1
        t=t+dt
    end do
end subroutine
function analitico (dx, n, t_i, t_l, a, l)  result(phi)
    integer::n, j, status=0
    double precision:: x
    double precision, intent(in):: t_i, t_l, a, l, dx
    double precision, allocatable,dimension(:):: phi
    x=dx
    allocate(phi(n), stat=status)
    j=2
    phi=t_l
    phi(1)=t_i
    do while (j/=n)
        phi(j)=a*((((t_l-t_i)/l)*x)+t_i)
        j=j+1
        x=x+dx
    end do
end function
end module conducao
