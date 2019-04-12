program cond
    use in_out
    implicit none      
    integer::i=2, n, status=0, show, cont=1 !i=posição no vetor, n=número de nós
    double precision:: a,dt,dx, d, tol, l, t_i, t_f, t, x
    double precision, allocatable, dimension(:)::w, k, dif !w é o vetor temperatura na barra, k é um vetor utilizado para calcular w
    allocate(k(n),stat=status) !aloca o vetor k utilizado para calcular a temperatura
    allocate(w(n), stat=status) !aloca o vetor w que é o vetor temperatura
    allocate(dif(n), stat=status)
    call entrada(n,show,l,t_i,t_f,a,tol,dt,t,dx)
    w(n)=t_f
    w(1)=t_i
    k=w
    call saida(t,x,w, cont)
    do
        do while (i/=n) !a temperatura nas pontas da barra não variam
            w(i)=((dt*a*(k(i+1)-(2*(k(i)))+k(i-1)))/dx**2)+k(i) !equação da tempratura na barra
            i=i+1
        end do
        i=2

        dif=w-k !vetor diferença entre os loops do calculo
        d=sum(dif)/n !d é média entre as diferenças entre os vetores w e k
        k=w !! k igual w diz que na próxima entrada os valores de k serão iguais os valores calculados para w nessa saída
        t=t+dt
        call saida(t,x,w, cont)
        if(abs(d)<=tol) then !se a diferença entre os vetores for menor ou igual a tlerância estabelcida
            exit !o calculo para
        end if
        cont=cont+1
    end do
    deallocate(w,stat=status)
    deallocate(k,stat=status)
    deallocate(dif,stat=status)
end program