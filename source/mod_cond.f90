module conducao
implicit none

contains
subroutine ent_cond(l,n,t_i,t_l,a,tol,dt,t, k, w, phi, dx)
    implicit none
    integer::n, status=0, input_id=21, cont=1
    double precision, allocatable, dimension(:):: k, w, phi
    double precision:: l, t_i, t_l, t, dt, dx, tol, a, in
    open(unit=input_id, action='read', status='old', iostat=status, file='../input.dat')!chamada dos dados de entrada
    do
        read(input_id, *,iostat=status)in
        if (status/=0)then
            exit
        end if
        if(cont==1)then         !linha 1 
            l=in                !tamanho da barra
        else if (cont==2)then   !linha 2
            n=nint(in)                !número de nós
        else if (cont==3) then  
            t_i=in              !temperatura no inicio da barra
        else if(cont==4) then   
            t_l=in              !temperatura ao longo da barra
        else if (cont==5)then
            a=in                !difusidade térmica
        else if (cont==6) then
            tol=in              !tolerância
        else if (cont==7)then
            dt=in               !passo de tempo
        else if(cont==8) then
            t=in                !tempo final
        end if
        cont=cont+1            !linha do arquivo
    end do
    close(input_id)
    allocate(k(n),stat=status)
    allocate(w(n), stat=status)
    allocate(phi(n), stat=status)
    k=t_l
    k(1)=t_i
    dx=l/(n-1)
    w=k
    phi=t_l
    phi(1)=t_i
end subroutine 
subroutine calc_cond (l,n,t_i,t_l,a,tol,dt,t, k, w, phi, dx, fid, d, i)
implicit none
integer::n, fid,i
double precision, allocatable,dimension(:)::w,k,dif,phi
double precision:: d,x,dx, tol, t, dt, l, a, t_l, t_i
x=dx
    do
        do while(i/=n) !enquanto i for menor que n:
            w(i)=((dt*a*(k(i+1)-2*(k(i))+k(i-1)))/dx**2)+k(i) ! o vetor w é igua a isso tudo que envolve o vetor k
            i=i+1             
        end do
        i=2
        dif=w-k !vetor diferença
        d=sum(dif) !soma dos elementos do vetor diferença
        if(d<tol)then !se a soma foir igual a zero, quero dizer que a peóxima linha i é igual a linha i-1
            exit !encerra o loop e mostra o tempo necessario para estabilizar
        end if
        k=w ! k igual w diz que na próxima entrada os valores de k serão iguais os valores calculados para w nessa saída
        write(fid,*)w
        t=t+dt
    end do
    do while (i/=n)
        phi(i)=a*((((t_l-t_i)/l)*x)+t_i) !resolução do modelo matemático
        i=i+1
        x=x+dx
    end do
    dif=phi-w !diferença entre o modelo matemático e númerico
    write(*,*) dif
    deallocate(dif)
end subroutine
end module conducao