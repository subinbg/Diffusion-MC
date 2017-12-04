! Diffuion Quantum Monte Carlo Solver
! For Ground State of Hydrogen Molecule
! Subin Bang, Seoul National University @ 2017.

#define Hartree 27.211385
#define tpi 6.28318530718
!compile with -cpp(gfortran) or -fpp(ifort) flag.

module iso_fortran_env

  ! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
  ! See Subclause 13.8.2 of the Fortran 2003 standard. 

  !implicit NONE 
  !public 

  integer, parameter :: Character_Storage_Size = 8 
  integer, parameter :: Error_Unit = 0 
  integer, parameter :: File_Storage_Size = 8 
  integer, parameter :: Input_Unit = 5 
  integer, parameter :: IOSTAT_END = -1 
  integer, parameter :: IOSTAT_EOR = -2 
  integer, parameter :: Numeric_Storage_Size = 32 
  integer, parameter :: Output_Unit = 6 

end module iso_fortran_env

module functions
  USE IFPORT
  implicit none
  !private
  
  public init_random_seed

contains

    subroutine init_random_seed()
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid, t(2), s
        integer(8) :: count, tms

        call random_seed(size = n)
        allocate(seed(n))
        open(newunit=un, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            call system_clock(count)
        if (count /= 0) then
            t = transfer(count, t)
        else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))

        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
                seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
        else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
        end if
        call random_seed(put=seed)
    end subroutine init_random_seed

    function potential(x, nelec, ndim, R)
        ! Calculating potential operator of H2 for given x.
        real(8) :: potential
        real(8), intent(in) :: x(:,:), R
        integer, intent(in) :: nelec, ndim

        integer :: i, j
        real(8) :: dist_proton1, dist_proton2, dist_elec

        potential = 0.0
        !if(ndim .ne. 3) stop 'Potential Only for 3-dim.'
        do i=1,nelec
            dist_proton1 = sqrt((x(i,1)-R/2)**2+(x(i,2))**2+(x(i,3))**2)
            dist_proton2 = sqrt((x(i,1)+R/2)**2+(x(i,2))**2+(x(i,3))**2)
            potential = potential - 1.0/dist_proton1 - 1.0/dist_proton2
        enddo

        do i=1,nelec
            do j=i+1,nelec
                dist_elec = sqrt((x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+(x(i,3)-x(j,3))**2)
                potential = potential + 1.0/dist_elec
            enddo
        enddo
            
    end function potential

    function potential_h(x, nelec, ndim, R)
        ! Calculating potential operator of H for given x.
        real(8) :: potential_h
        real(8), intent(in) :: x(:,:), R
        integer, intent(in) :: nelec, ndim

        integer :: i, j
        real(8) :: dist_proton1, dist_proton2, dist_elec

        potential_h = 0.0
        !if(ndim .ne. 3) stop 'Potential Only for 3-dim.'
        do i=1,nelec
            dist_proton1 = sqrt((x(i,1))**2+(x(i,2))**2+(x(i,3))**2)
        !    dist_proton2 = sqrt((x(i,1)+R/2)**2+(x(i,2))**2+(x(i,3))**2)
            potential_h = potential_h - 1.0/dist_proton1 !- 1.0/dist_proton2
        enddo

        !do i=1,nelec
        !    do j=i+1,nelec
        !        dist_elec = sqrt((x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+(x(i,3)-x(j,3))**2)
        !        potential = potential + 1.0/dist_elec
        !    enddo
        !enddo
    end function potential_h

    subroutine life_count(life, n_life)
        integer, intent(in) :: life(:)
        real(8), intent(inout) :: n_life

        integer :: i, test
        n_life = 0.0
        test = size(life)
        do i=1,size(life)
            if(life(i).ge.0) n_life = n_life + 1.0
        enddo

    end subroutine life_count

    subroutine box_muller(ranseed)
        real(8), intent(inout) :: ranseed

        real(8) :: U, V
        
        U = 0; V = 0
        do while(.true.)
            call random_number(U)
            call random_number(V)
            if(U>epsilon(V)) exit
        enddo
        ranseed = sqrt(-2*log(U))*cos(tpi*V)
    end subroutine box_muller

    subroutine result(x, N, nelec, ndim, life, E)
        real(8), intent(in) :: x(:,:,:), E
        integer, intent(in) :: N, nelec, ndim, life(:)
        integer :: un, i, j
        open(un, file='coord.dat', status='new')

        !write(un,*) E
        do i=1,N
            if(life(i)<0) cycle
            do j=1,nelec
                write(un,*) x(i,j,1), x(i,j,2), x(i,j,3)
            enddo
        enddo

        close(un,status='keep')
    end subroutine result

    subroutine calc_probd(temp, lattice, x, prob_x, niter, avg_start, r_max)
        real(8), intent(in) ::  lattice(:), x(:)
        real(8), intent(inout) :: temp(:), prob_x(:)
        integer, intent(in) :: niter, avg_start
        real(8), intent(in) :: r_max

        real(8) :: normed
        integer :: ix

        temp(:) = lattice(:)
        normed = norm2(x(1:2))
        if(normed .le. r_max) then
            temp(:) = normed - temp(:)
            temp = pack([(ix,ix=1,size(temp))], temp<0)
            prob_x(temp(1)-1) = prob_x(temp(1)-1) + 1.0/(niter-avg_start)
        endif

    end subroutine calc_probd

    subroutine calc_probd_2(temp, tempp, lattice, x, prob_x2, niter, avg_start, r_max)
        real(8), intent(in) ::  lattice(:), x(:)
        real(8), intent(inout) :: temp(:), prob_x2(:,:), tempp(:)
        integer, intent(in) :: niter, avg_start
        real(8), intent(in) :: r_max

        real(8) :: x1,x2
        integer :: ix

        temp(:) = lattice(:); tempp(:) = lattice(:)
        x1 = abs(x(1))
        x2 = abs(x(2))
        if((x1.lt.r_max/2) .and. (x2.lt.r_max/2)) then
            temp(:) = x(1) - temp(:)
            tempp(:) = x(2) - tempp(:)
            temp = pack([(ix,ix=1,size(temp))], temp<0)
            tempp = pack([(ix,ix=1,size(tempp))], tempp<0)
            prob_x2(temp(1)-1, tempp(1)-1) = prob_x2(temp(1)-1, tempp(1)-1) + 1.0/(niter-avg_start)
        endif
    end subroutine calc_probd_2

    subroutine calc_probd_h(temp, lattice, x, prob_x, niter, avg_start, r_max)
        real(8), intent(in) ::  lattice(:), x(:)
        real(8), intent(inout) :: temp(:), prob_x(:)
        integer, intent(in) :: niter, avg_start
        real(8), intent(in) :: r_max

        real(8) :: normed
        integer :: ix

        temp(:) = lattice(:)
        normed = norm2(x(:))
        if(normed .le. r_max) then
            temp(:) = normed - temp(:)
            temp = pack([(ix,ix=1,size(temp))], temp<0)
            prob_x(temp(1)-1) = prob_x(temp(1)-1) + 1.0/(niter-avg_start)
        endif

    end subroutine calc_probd_h

end module functions

program dqmc_H2
    ! Diffuion Quantum Monte Carlo Solver
    ! For Ground State of Hydrogen Molecule
    ! Subin Bang, Seoul National University @ 2017.
    use functions
    use iso_fortran_env

    ! N = number of replica
    ! nelec = number of electrons
    ! niter = number of iteration
    ! replica_test = 
    ! array_size = 
    integer :: N, nelec, niter, replica_test, array_size, avg_start, avg_step

    ! R = equilibrium distance between two protons
    ! It is already known that R = 0.74A, so R = 1.4 Bohr Radii.
    ! dt = time step defined in propagator
    ! potential_avg = 
    ! W = 
    ! n_life = 
    ! E = 
    real(8) :: R, dt, potential_avg, W, n_life_before, n_life_now, E

    ! x  = array storing positions of 
    !      2N*nelec*3(spatial degree of freedom) pseudoparticles.
    real(8), allocatable :: x(:,:,:)

    ! life = array storing information about
    !        whether i-th replica is dead or alive with 2N dimension.
    !        If alive, 1. If dead, -1.
    ! life_index = array storing index of elements of life, with value +-1.
    integer, allocatable :: life(:), life_index(:), test(:)

    integer :: i_iter, i, j, k, count, choice, ix
    real(8) :: ranseed, pot_value

    real(8) :: temp1, temp2 ! debugging variables

    real(8), allocatable :: prob_x(:), prob_x2(:,:), lattice(:), temp(:), tempp(:)
    integer :: r_step
    real(8) :: r_max!,r1,r2,normed

    abstract interface
      function func (x, nelec, ndim, R)
        real(8) :: func
        real(8), intent(in) :: x(:,:), R
        integer, intent(in) :: nelec, ndim
      end function func

      subroutine probd(temp, lattice, x, prob_x, niter, avg_start, r_max)
        real(8), intent(in) ::  lattice(:), x(:)
        real(8), intent(inout) :: temp(:), prob_x(:)
        integer, intent(in) :: niter, avg_start
        real(8), intent(in) :: r_max
      end subroutine probd
    end interface
    procedure (func), pointer :: f_ptr => null ()
    procedure (probd), pointer :: probd_ptr => null ()

    call init_random_seed

    open(10,file='energy.dat',status='new')
    open(11,file='probd.dat',status='new')
    open(12,file='input',status='old')

    read(12,*) N != 10000
    read(12,*) nelec != 1 !!!!
    read(12,*) niter != 100000
    read(12,*) dt != 0.01
    read(12,*) R != 1.4
    array_size = 10*N
    read(12,*) avg_start != 1000
    read(12,*) avg_step
    read(12,*) r_step != 1000
    read(12,*) r_max != 10.0
    read(12,*) choice

    allocate(x(1:array_size,1:nelec,1:3), life(1:array_size), life_index(1:array_size), test(1:array_size))
    life(1:N) = 1; life(N+1:array_size) = -1

    if(choice .eq. 1) then
        f_ptr => potential_h
        probd_ptr => calc_probd_h
        allocate(prob_x(1:r_step),lattice(1:r_step+1),temp(1:r_step+1))
        x(:,:,:) = 1.0
        do i=1,r_step+1
            lattice(i) = r_max/r_step*(i-1)
        enddo
    else
        f_ptr => potential
        !probd_ptr => calc_probd
        allocate(prob_x2(1:r_step,1:r_step), lattice(1:r_step+1),temp(1:r_step+1), tempp(1:r_step+1))
        x(:,:,:) = 0.0
        do i=1,r_step+1
            lattice(i) = r_max/r_step*(i-1-r_step/2)
        enddo
    endif

    life_index(:) = 0
    test(:) = 0
    prob_x(:) = 0

    E = 0.0
    ! Main loop
    do i_iter = 1,niter
        ! Propose new positions
        do i = 1,array_size
            if (life(i) < 0) cycle
            do j = 1,nelec
                do k = 1,3
                    call box_muller(ranseed)
                    x(i,j,k) = x(i,j,k) + sqrt(dt) * ranseed
                enddo
            enddo
        enddo

        if (i_iter .eq. 1) then
            potential_avg = 0.0
            do i = 1,N !array_size
                !if (life(i) < 0) cycle
                potential_avg = potential_avg + f_ptr(x(i,:,:), nelec, 3, R) !!!!
            enddo
            potential_avg = potential_avg/N
            n_life_before = N
            n_life_now = N
        else
            call life_count(life, n_life_now)
            !potential_avg = potential_avg + (1.0-n_life_now/n_life_before)/dt
            potential_avg = 0.0
            do i = 1,array_size
                if (life(i) < 0) cycle
                potential_avg = potential_avg + f_ptr(x(i,:,:), nelec, 3, R) !!!!
            enddo
            potential_avg = potential_avg/n_life_now
            potential_avg = potential_avg !+ (1.0-n_life_now/n_life_before)/dt
            n_life_before = n_life_now
        endif

        test(:) = 0
        count = 1
        do i=1,array_size
            if (life(i) < 0) cycle
            if ( ANY( test(1:count)==i ) ) cycle
            W = exp(-(f_ptr(x(i,:,:), nelec, 3, R)-potential_avg)*dt) !!!!

            !write(Output_Unit,*) i_iter, i, f_ptr(x(i,:,:), nelec, 3, R), potential_avg

            call random_number(ranseed)
            replica_test = min(int(W+ranseed), 3)

            if(replica_test .eq. 0) then
                life(i) = -1
            elseif(replica_test .eq. 2) then
                life_index = pack([(ix,ix=1,size(life))], life<0)
                x(life_index(1),:,:) = x(i,:,:)
                life(life_index(1)) = 1

                test(count) = life_index(1)
                count = count + 1
            elseif(replica_test .eq. 3) then
                life_index = pack([(ix,ix=1,size(life))], life<0)
                x(life_index(1),:,:) = x(i,:,:)
                x(life_index(2),:,:) = x(i,:,:)
                life(life_index(1)) = 1
                life(life_index(2)) = 1

                test(count) = life_index(1)
                test(count+1) = life_index(2)
                count = count + 2
            endif
        enddo

        if(n_life_now.le.1) stop

        if (i_iter>avg_start) then
            pot_value = 0.0
            if(choice .eq. 1) then
                do i = 1,array_size
                    if (life(i) < 0) cycle
                    pot_value = pot_value + f_ptr(x(i,:,:), nelec, 3, R) !!!!
                    do j=1,nelec
                        ! do k=1,r_step
                        !     r1 = r_max/real(r_step)*(k)
                        !     r2 = r_max/real(r_step)*(k-1)
                        !     normed = norm2(x(i,j,:))
                        !     if(r2.le.normed) then
                        !         if(r1.gt.normed) prob_x(k) = prob_x(k) + 1.0/(niter-avg_start)
                        !     endif
                        ! enddo
                        ! temp(:) = lattice(:)
                        ! normed = norm2(x(i,j,:))
                        ! temp(:) = normed - temp(:)
                        ! temp = pack([(ix,ix=1,size(temp))], temp<0)
                        ! prob_x(temp(1)-1) = prob_x(temp(1)-1) + 1.0/(niter-avg_start)
                        call probd_ptr(temp, lattice, x(i,j,:), prob_x, niter, avg_start, r_max)
                    enddo 
                enddo
            else
                do i = 1,array_size
                    if (life(i) < 0) cycle
                    pot_value = pot_value + f_ptr(x(i,:,:), nelec, 3, R) !!!!
                    do j=1,nelec
                        ! do k=1,r_step
                        !     r1 = r_max/real(r_step)*(k)
                        !     r2 = r_max/real(r_step)*(k-1)
                        !     normed = norm2(x(i,j,:))
                        !     if(r2.le.normed) then
                        !         if(r1.gt.normed) prob_x(k) = prob_x(k) + 1.0/(niter-avg_start)
                        !     endif
                        ! enddo
                        ! temp(:) = lattice(:)
                        ! normed = norm2(x(i,j,:))
                        ! temp(:) = normed - temp(:)
                        ! temp = pack([(ix,ix=1,size(temp))], temp<0)
                        ! prob_x(temp(1)-1) = prob_x(temp(1)-1) + 1.0/(niter-avg_start)
                        call calc_probd_2(temp, tempp, lattice, x(i,j,:), prob_x2, niter, avg_start, r_max)
                    enddo 
                enddo
            endif
            
            pot_value = pot_value/n_life_now
            E = E + pot_value * Hartree
            if(choice .ne. 1) E = E + Hartree / R
            if(mod(i_iter,avg_step).eq.0) write(10,*) E/(i_iter-avg_start)
            write(Output_Unit, *) E/(i_iter-avg_start)
            write(Output_Unit, *) n_life_before
            write(Output_Unit, *) n_life_now
            write(Output_Unit, *) potential_avg
            write(Output_Unit, *) replica_test
            write(Output_Unit, *) potential_h(x(i,:,:), nelec, 3, R) !!!!
            write(Output_Unit, *) i_iter
            write(Output_Unit, *) '================='
        endif
        
    enddo
    E = E/(niter-avg_start)

    if(choice .eq. 1) then
        do i=1,r_step
            r2 = r_max/real(r_step)*i
            write(11,*) r2, prob_x(i), (r2**2)*exp(-r2)
        enddo
    else
        do i=1,r_step
            r1 = r_max/real(r_step)*(i-1-r_step/2)
            do j=1,r_step
                r2 = r_max/real(r_step)*(j-1-r_step/2)
                write(11,*) r1, r2, prob_x2(i,j)
            enddo
        enddo
    endif
    

    call result(x, array_size, nelec, ndim, life, E)

    deallocate(x, life, life_index)
    if(choice .eq. 1) deallocate(test, prob_x, lattice)
    close(10,status='keep')
    close(11,status='keep')

end program dqmc_H2
