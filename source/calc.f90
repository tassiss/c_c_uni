program cond
    use in_out
    use functions
    !use g_s
    implicit none      
    integer::i=2, n, status=0, show, cont=1, n_g, i_e!i=posição no vetor, n=número de nós
    double precision:: k,dt,dx, d, tol, l, t_i, t_l, t=0.0d0, x, erro=0, rm, time, rho, cp, dkp, dkn, z
    double precision, allocatable, dimension(:)::w, b, dif, r_a, e_m,d_k !r_a é "resultado analitico" !w é o vetor temperatura na barra, k é um vetor utilizado para calcular w
    double precision, allocatable, dimension(:,:):: a !matriz 'A'

    call entrada(n,show,l,t_i,t_l,k,tol,dt,dx, time, n_g, rho, cp, i_e)
    allocate(a(1-n_g:n+n_g,1-n_g:n+n_g), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação a'
        stop
    end if
    !write(8,*) m
    allocate(d_k(1-n_g:n+n_g), stat=status)
    if (status/=0) then
        write(*,*) 'erro de alocação d_k'
        stop
    end if
    allocate(b(1-n_g:n+n_g),stat=status) !aloca o vetor k utilizado para calcular a temperatura
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
    call dirichlet(n_g, w, b, t_i, t_l, n)
    b=w
    
    call saida(t,dt,x,w, cont, erro, n_g)
    r_a=w(1:n)
    !######CÁLCULO DO MÉTODO ANALITICO#######
    do i=1,n
        x=dx/2.d0 +(i-1)*dx
        r_a(i)=analitico(t_i, t_l, l, x)
    end do
    !#######################################
    
    
    a= 0.0d0
    
    
    do
    !#####CÁLCULO DA DIFUSIVIDADE###########
        do i=1-n_g,n+n_g
            x=dx/2.d0 +(i-1)*dx
            d_k(i)=dif_ter(x,t,k)
        end do

    !#######################################
        

    !#####CÁLCULO DO MÉTODO NUMÉRICO########   
        if (i_e==1) then  !seleção para resolver méto implicito
            do i=1-n_g,n+n_g !a temperatura nas pontas da barra variam
                if (i>=1 .and. i<=n) then
                    dkp=(d_k(i)+d_k(i+1))/2.0d0
                    dkn=(d_k(i-1)+d_k(i))/2.0d0
                end if
            
            
                if(i<1 .or. i>n)then
                    a(i,i)=1.0d0
                else
                    z=(dt/((dx**2.0d0)*rho*cp))
                    a(i,i)=(z*(dkp+dkn))+1.0d0 !B
                    a(i,i+1)=-z*dkp            !C
                    a(i,i-1)=-z*dkn            !A
                end if
                write(9,*) a(i,:)
            end do
            call g_s (n+(2*n_g),-z*dkn,(z*(dkp+dkn))+1.0d0,-z*dkp,w,b)
            !call lapack_lin_sys(n+(2*n_g),a,b,w)
            if (all(b==w)) then
                write(*,*)'Erro'
                stop
            end if
        end if
        
        if (i_e==2)then !seleção para resolver méto explicito
            do i=1,n !a temperatura nas pontas da barra variam
                !#método explicito###########################
                dkp=(d_k(i)+d_k(i+1))/2.0d0
                dkn=(d_k(i-1)+d_k(i))/2.0d0
                w(i)= (((dt/(4.0d0*(dx**2.0d0)))*((dkp*(b(i+1)-b(i)))-(dkn*(b(i)-b(i-1)))))+b(i))/(rho*cp)!equação da tempratura na barra
            end do
        end if
        
        
       !###########################################
         
    !#######################################
        dif=w-b !vetor diferença entre os loops do calculo
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
        b=w
        call dirichlet(n_g, w, b, t_i, t_l, n)
        d=erro
    end do 
    deallocate(w,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
    deallocate(b,stat=status)
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
    deallocate(a,stat=status)
    if (status/=0) then
        write(*,*) 'erro de dealocação'
        stop 
    end if
end program
