program cond
    implicit none
    integer::status=0, n=10, i=2
    double precision:: dx, dt, l, t=0.0, a
    double precision, allocatable, dimension (:)::k, w
    allocate(k(n), stat=status)
    allocate(w(n), stat=status)
    a=1.0 !difusidade térmica
    l=2.0 !comprimento da barra
    k=100.0 ! temperatura da barra
    k(1)=200.0 ! temperatura no primeiro no
    k(n)=100.0 !temperatura no último nó
    dx=l/n !tamanho do dx
    dt=0.001 !tamnaho do dt
    w=k
    write(20,*) w
    do while (t<5) !enquanto t for menor que 5:
        do while(i/=n) !enquanto i for menor que n:
            w(i)=((dt*a*(k(i+1)-2*(k(i))+k(i-1)))/dx**2)+k(i) ! o vetor w é igua a isso tudo que envolve o vetor k
            i=i+1 
        end do
        i=2
        k=w ! k igual w diz que na próxima entrada os valores de k serão iguais os valores calculados para w nessa saída
        write(20,*)w
        t=t+dt 
    end do
end program