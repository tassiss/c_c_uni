program cond
    use in_out
    use functions
    use lapack_parcer
    implicit none      
    integer::i=2, n, status=0, show, cont=1, n_g!i=posição no vetor, n=número de nós
    double precision:: a,dt,dx, d, tol, l, t_i, t_l, t=0.0d0, x, erro=0, rm, time, rho, cp, dkp, dkn, z
    double precision, allocatable, dimension(:)::w, k, dif, r_a, e_m,d_k !r_a é "resultado analitico" !w é o vetor temperatura na barra, k é um vetor utilizado para calcular w
    double precision, allocatable, dimension(:,:):: m !matriz 'A'

    call entrada(n,show,l,t_i,t_l,a,tol,dt,dx, time, n_g, rho, cp)
    allocate(m(1-n_g:n+n_g,1-n_g:n+n_g), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação m'
        stop
    end if
    write(8,*) m
    allocate(d_k(1-n_g:n+n_g), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação d_k'
        stop
    end if
    allocate(k(1-n_g:n+n_g),stat=status) !aloca o vetor k utilizado para calcular a temperatura
    if (status/=0) then
        write(*,*) 'erro de alocação k'
        stop
    end if
    allocate(w(1-n_g:n+n_g), stat=status) !aloca o vetor w que é o vetor temperatura
    if (status/=0) then
        write(*,*) 'erro de alocação w'
        stop
    end if
    allocate(dif(n), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação dif'
        stop
    end if
    allocate(r_a(n), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação r_a'
        stop
    end if
    
    allocate(e_m(n), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação e_m'
        stop
    end if
    w(1:n)=t_l
    call dirichlet(n_g, w, k, t_i, t_l, n)
    k=w
    
    call saida(t,dt,x,w, cont, erro, n_g)
    r_a=w(1:n)
    !######CÁLCULO DO MÉTODO ANALITICO#######
    do i=1,n
        x=dx/2.d0 +(i-1)*dx
        r_a(i)=analitico(t_i, t_l, l, x)
    end do
    !#######################################
    
    
    m= 0.0d0
    
    
    do
    !#####CÁLCULO DA DIFUSIVIDADE###########
        do i=1-n_g,n+n_g
            x=dx/2.d0 +(i-1)*dx
            d_k(i)=dif_ter(x,t,a)
        end do

    !#######################################
        

    !#####CÁLCULO DO MÉTODO NUMÉRICO########   
        do i=1-n_g,n+n_g !a temperatura nas pontas da barra variam
            if (i>=1 .and. i<=n) then
                dkp=(d_k(i)+d_k(i+1))/2.0d0
                dkn=(d_k(i-1)+d_k(i))/2.0d0
            end if
            if(i<1 .or. i>n)then
                 m(i,i)=1.0d0
            else
                z=(dt/((dx**2.0d0)*rho*cp))
                m(i,i)=(z*(dkp+dkn))+1
                m(i,i+1)=-z*dkp
                m(i,i-1)=-z*dkn
            end if
            write(9,*) m(i,:)
        end do
        
        call lapack_lin_sys(n+(2*n_g),m,k,w)
        if (all(k==w)) then
            write(*,*)'Erro'
            stop
        end if
    !#######################################
        dif=w-k !vetor diferença entre os loops do calculo
        e_m=w(1:n)-r_a !vetor diferença em relação ao analitico
        rm=abs(sum(e_m)/n) !média do erro em relação analitico
        erro=abs(sum(dif)/n) !d é média entre as diferenças entre os vetores w e k
        t=t+dt
        cont=cont+1
        write(*,*) 'erro ao analitico', rm
        write(*,*) 'diferença', erro
        write(*,*) 'tempo[s]:      ', t
        if(mod(cont,show)==0)then
            call saida(t,dt,x,w,cont, erro, n_g)
        end if

        if((abs(erro)<=tol .and. abs(rm)<=tol) .and. cont>100 .or.time<=t) then !se a diferença entre os vetores for menor ou igual a tlerância estabelcida
            call saida(t,dt,x,w, cont, erro, n_g)
            write(*,*) t
            exit !o calculo para
        end if 
        call dirichlet(n_g, w, k, t_i, t_l, n)
        k=w
        d=erro
    end do 
    deallocate(w,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
    deallocate(k,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
    deallocate(dif,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
    
    deallocate (r_a, stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
    deallocate(e_m, stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
    deallocate(d_k,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if    
    deallocate(m,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
end program
