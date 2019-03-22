program cond
    use conducao
    implicit none
    integer::status=0, n, i=2, fid=20, show
    double precision:: dx, dt, l,d, t, a, t_i, t_l, tol
    double precision, allocatable, dimension (:)::k, w
    double precision,allocatable,dimension(:):: phi
    open(unit=fid, action='write', status='unknown', iostat=status,file='temp.dat')
    !##########################################################################
    call ent_cond(l,n,t_i,t_l,a,tol,dt,t, k, w, dx, show) !chamada de rotina entrada de dados
    !#######################################################################################
    call calc_cond (n,t_i,t_l,a,tol,dt,t, k, w, dx, fid, d, i, show)
    write(fid,*) w
    phi= analitico(dx, n, t_i, t_l, a, l)
    write(*,*) phi
    deallocate(k)
    deallocate(w)
    !deallocate(phi)
    close(fid)
    write(*,*) t
end program