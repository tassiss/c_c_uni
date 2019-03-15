program cond
    implicit none
    integer::status=0, n=10, i=2, fid=20
    double precision:: dx, dt, l, t=0.0, a, d, t_i, t_l, x
    double precision, allocatable, dimension (:)::k, w, dif, phi
    allocate(k(n), stat=status)
    allocate(w(n), stat=status)
    allocate(phi(n), stat=status)
    open(unit=fid, action='write', status='unknown', iostat=status,file='temp.dat')
    t_i=200 !temperatura no inicio da barra
    t_l=100 !temperatura no final da barra
    a=1 !difusidade térmica
    l=2.0 !comprimento da barra
    k=100.0 ! temperatura da barra
    k(1)=t_i ! temperatura no primeiro no
    k(n)=t_l !temperatura no último nó
    dx=l/(n-1) !tamanho do dx
    dt=0.001 !tamnaho do dt
    w=k
    x=dx
    phi(1)=t_i
    phi(n)=t_l
    write(fid,*) w
    do
        do while(i/=n) !enquanto i for menor que n:
            w(i)=((dt*a*(k(i+1)-2*(k(i))+k(i-1)))/dx**2)+k(i) ! o vetor w é igua a isso tudo que envolve o vetor k
            i=i+1             
        end do
        i=2
        dif=w-k !vetor diferença
        d=sum(dif) !soma dos elemento do vetor diferença
        if(d<1E-14)then !se a soma foir igual a zero, quero dizer que a peóxima linha i é igual a linha i-1
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
    deallocate(k)
    deallocate(w)
    deallocate(dif)
    deallocate(phi)
    close(fid)
    write(*,*) t
end program