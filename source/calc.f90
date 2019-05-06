program cond
    use in_out
    use functions
    implicit none      
    integer::i=2, n, status=0, show, cont=1, n_g !i=posição no vetor, n=número de nós
    double precision:: a,dt,dx, d, tol, l, t_i, t_l, t=0.000, x, erro=0, rm, time, rho, cp, dkp, dkn
    double precision, allocatable, dimension(:)::w, k, dif, r_a, e_m,d_k,g_c !r_a é "resultado analitico" !w é o vetor temperatura na barra, k é um vetor utilizado para calcular w
    call entrada(n,show,l,t_i,t_l,a,tol,dt,dx, time, n_g)
    allocate(g_c(1-n_g:n+n_g))
    allocate(d_k(1-n_g:n+n_g), stat=status)
    allocate(k(n),stat=status) !aloca o vetor k utilizado para calcular a temperatura
    allocate(w(n), stat=status) !aloca o vetor w que é o vetor temperatura
    allocate(dif(n), stat=status)
    allocate(r_a(n), stat=status)
    allocate(e_m(n), stat=status)
    w=t_l
    k=w
    rho=1.0
    cp=1.0
    call saida(t,dt,x,w, cont, erro)
    r_a=w
    !######CÁLCULO DO MÉTODO ANALITICO#######
    do i=2,(n-1)
        r_a(i)=analitico(t_i, t_l,a, l, x)
        x=(i)*dx
    end do
    !#######################################
    g_c(1:n)=w
    !######CÁLCULO DAS CELULAS GHOSTS#######
    do i=1-n_g,0
        g_c(i)=dirichlet(t_i, t_l, w, i, n)

    end do
    do i=n+1,n+n_g
        g_c(i)=dirichlet(t_i, t_l,w,i, n)

    end do
    !########################################
    do
    !#####CÁLCULO DA DIFUSIVIDADE###########
        do i=1-n_g,n+n_g
            d_k(i)=dif_ter(x,t)
            x=(i-1)*dx
        end do

    !#######################################
        
    !######CÁLCULO DAS CELULAS GHOSTS#######
        do i=1-n_g,0
            g_c(i)=dirichlet(t_i, t_l, w, i, n)
    
        end do
        do i=n+1,n+n_g
            g_c(i)=dirichlet(t_i, t_l,w,i, n)
    
        end do
        !#########################################

    !#####CÁLCULO DO MÉTODO NUMÉRICO########   
        do i=1,n !a temperatura nas pontas da barra variam
            dkp=(d_k(i)+d_k(i+1))/2.0
            dkn=(d_k(i-1)+d_k(i))/2.0
            w(i)= (((dt/(4*(dx**2)))*((dkp*(g_c(i+1)-g_c(i)))-(dkn*(g_c(i)-g_c(i-1)))))+g_c(i))/(rho*cp)!equação da tempratura na barra
        end do
    !#######################################
        dif=w-k !vetor diferença entre os loops do calculo
        e_m=w-r_a !vetor diferença em relação ao analitico
        rm=abs(sum(e_m)/n) !média do erro em relação analitico
        erro=abs(sum(dif)/n) !d é média entre as diferenças entre os vetores w e k
        t=t+dt
        cont=cont+1
        write(*,*) 'erro ao analitico', rm
        write(*,*) 'erro', erro
        write(*,*) 'tempo[s]:      ', t
        if(mod(cont,show)==0)then
            call saida(t,dt,x,w,cont, erro)
        end if
        ! if(abs(erro)>abs(5*d) .and. cont>100) then
        !     call saida(t,dt,x,w,cont, erro)
        !     write(*,*) 'o programa apresenta erro // alta divergência'
        !     exit
        ! end if 
        if((abs(erro)<=tol .and. abs(rm)<=tol) .and. cont>100 .or.time<=t) then !se a diferença entre os vetores for menor ou igual a tlerância estabelcida
            call saida(t,dt,x,w, cont, erro)
            write(*,*) t
            exit !o calculo para
        end if
        !g_c(1:n)=w
         do i=1,n
             g_c(i)=(2*g_c(i)+g_c(i+1)+g_c(i-1))/4
         end do !! k igual w diz que na próxima entrada os valores de k serão iguais os valores calculados para w nessa saída
        write(*,*) g_c(0)
         d=erro      
    end do
    !write(*,*) erro
    deallocate(w,stat=status)
    deallocate(k,stat=status)
    deallocate(dif,stat=status)
    deallocate (r_a, stat=status)
    deallocate(e_m, stat=status)
    deallocate(d_k,stat=status)
    deallocate(g_c)
end program
