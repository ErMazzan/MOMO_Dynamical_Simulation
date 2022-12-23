c	Dynamical Simulation proram
c	Author: David Mazzanti Tarancón

c	The program is divided in two parts. The first one solves the first exercise
c	of the final project and the second part solves the second exercise. Once the
c	code is compiled, if the user executes the program it will be asked to solve
c	one of the exercises. The input must be 1 or 2.
	
	implicit none
	integer n, i, j, Ntimesteps, seed, exer, msd_flag
	parameter (n = 125)
	real*8 rho, l, cutoff, upot, dt, ekin, etotal, T, time, sigma_T
	real*8 nu, sigma, t_ins, p, eps, m, kb, nv, press, msd_val
	real*8, dimension(n,3) :: r, auxr, oldr, r_eul, v_eul, f, v, r_ini
	character(8) fmt,ext
	
	seed=784378
	call srand(seed)
	
	write(*,*) "Enter the number exercise to solve:"
	read(*,*) exer
	
c	EXERCISE 1	
	if (exer == 1) then
	
c	Here we define the initial conditions for all the simulation

	rho = 0.7
	T = 100.0
	
	Ntimesteps = 200000
	dt = 0.001
	time = 0.d0
	
	r = 0.d0
	v = 0.d0
	f = 0.d0
	upot = 0.d0
	ekin = 0.d0

c	Initially, the velocities are distributed following a bimodal distribution with
c	mean = 0 and sigma = 1. Particles are set on a square lattice	
	call bimodal_v(v,T,n)
	call square_lattice(rho,n,r,l)
	cutoff = l/3.d0
	
	r_eul = r
	v_eul = v
	open(10, file="exercise_1/dt_10-3/verlet_thermodynamics.dat")
	open(11, file="exercise_1/dt_10-3/euler_thermodynamics.dat")
	open(12, file="exercise_1/dt_10-3/initial_distribution.dat")
	open(13, file="exercise_1/dt_10-3/final_distribution.dat")

c	Writting of the initial velocities	
	do i = 1,n
		do j = 1,3
		write(12,*) v(i,j)
		enddo
	enddo
	
	close(12)

c	Main core of the simulation, loop that perfoms N steps.

	do i = 1, Ntimesteps
		call verl_vel(r,v,l,cutoff,dt,n)
		call pbc2(r,l,n)
		
		call euler_step(n,r_eul,v_eul,l,cutoff,upot,f,dt)
		call pbc2(r_eul,l,n)
		
		if (mod(i,1000)==0) then
			print*,i
			call kin_energy(n,v,ekin)
			call leonard_jones_force(n,r,l,cutoff,upot,f)
			call total_mom(ekin,p)
			etotal = ekin + upot
			write(10,*) time, ekin, upot, etotal, p
			
			call kin_energy(n,v_eul,ekin)
			call leonard_jones_force(n,r_eul,l,cutoff,upot,f)
			call total_mom(ekin,p)
			etotal = ekin + upot
			write(11,*) time, ekin, upot, etotal, p
		endif
		
		time = time + dt
	enddo
	
	close(10)
	close(11)
	
	do i = 1,n
		do j = 1,3
		write(13,*) v(i,j)
		enddo
	enddo
	
	close(13)
	
	endif

c	EXERCISE 2	
	if (exer == 2) then

c	Here we define the initial conditions for all the simulation	
	Ntimesteps = 500000
	dt = 0.0001
	time = 0.d0
	
	r = 0.d0
	v = 0.d0
	f = 0.d0
	upot = 0.d0
	ekin = 0.d0
	
	eps = 0.998
	sigma = 3.4
	m = 40.0
	nv = 6.02214d23
	kb = 1.380649d-23
	
	rho = 0.8
	
c	MSD flag is a variable that activates the calculus of the MSD (set 1 to compute it
c	set 0 to not do it)
	msd_flag = 1
	
	call square_lattice(rho,n,r,l)
	cutoff = l/3.d0
	
	fmt='(f5.3)'
	write(ext,fmt) rho
	open(10, file="exercise_2/thermodynamics_"//trim(ext)//".dat")
	
	if (msd_flag == 1) then
		open(11, file="exercise_2/msd_"//trim(ext)//".dat")
		dt = 0.001
	endif

c	Nu parameter for Andersen thermostat	
	nu = 0.1

c	Save the initial configuration (for MSD calculus)	
	do i = 1, n
		do j = 1, 3
		r_ini(i,j) = r(i,j)
		enddo
	enddo

c	Main core of the simulation	
	do i = 1, Ntimesteps
		if (i.lt.(10000)) then
			T = 100.0
		else
			T = 1.2
		endif
		sigma_T = dsqrt(T)
		
		call verl_vel(r,v,l,cutoff,dt,n)
		
c	In order to compute the MSD we do not apply PBC. Otherwise, the sistem would
c	not be open and the diffusion coefficient would be incorrect.

		if (msd_flag == 0) then
			call pbc2(r,l,n)
		endif
		
		call therm_Andersen(n,v,nu,sigma_T)
		
		if (mod(i,1000)==0) then
			print*, i
			call kin_energy(n,v,ekin)
			call leonard_jones_force(n,r,l,cutoff,upot,f)
			call temp_inst(ekin,n,t_ins)
			call p_inst(n, rho, l, r, f, t_ins, press)

			etotal = ekin + upot
			
			t_ins = 1d3*eps/(kb*nv)*t_ins
			press = 1d3*eps/(nv*(sigma*1d-10)**3d0)*press
			time = dsqrt(m*1d-3/(eps*1d3))*sigma*1d-10*dt*i
			
			ekin = ekin*eps
			upot = upot*eps
			etotal = etotal*eps
			
			write(10,*) time, ekin, upot, etotal, t_ins, press
			
			if (msd_flag == 1) then
				call msd(n, r_ini, r, msd_val)
				msd_val = msd_val*sigma*sigma
				write(11,*) time, msd_val
				print*, msd_val
			endif
		endif
	
	enddo
	
	close(10)
	
	if ((rho == 0.05).or.(rho == 0.8)) then
		close(11)
	endif
	
	endif
	
	end

C	Subroutine that sets initially the system on a square lattice
C	INPUT: rho (density of the system), n (# of particles)
C	OUTPUT: r (positions of the particles), l (size of the box)
	subroutine square_lattice(rho,n,r,l)
	implicit none
	integer m, n, nx, ny, nz, i
	real*8 l, a, rho
	real*8 r(n,3)
	
	l = (n/rho)**(1.d0/3.d0)
	m = nint(n**(1.d0/3.d0))
	a = l/m
	
	i = 1
	
	do nx = 1, m
		do ny = 1, m
			do nz = 1, m
			r(i, 1) = (nx-1)*a
			r(i, 2) = (ny-1)*a
			r(i, 3) = (nz-1)*a
			i = i+1
			enddo
		enddo
	enddo
	
	end

C	Subroutine that sets the initial velocities of the particles
C	INPUT: T (temperature of the system) n (# of particles)
C 	OUTPUT: v (velocities of the particles)

	subroutine bimodal_v(v,T,n)
	implicit none
	integer n, i, j
	real*8 vel, T, x
	real*8, dimension(n,3) :: v
	
	vel = dsqrt(T)
	do i = 1, n
		do j = 1, 3
			x = rand()
			if (x.gt.0.5) then
			v(i,j) = vel
			else
			v(i,j) = vel*(-1.d0)
			endif
		enddo
	enddo
	end

C	Subroutine that computes the momenta of the system
C	INPUT: ekin (kinetic energy of the system)
C	OUTPUT: p_tot (total momenta of the system)
	
	subroutine total_mom(ekin,p_tot)
	implicit none
	real*8 ekin, p_tot
	
	p_tot = 0.d0
	p_tot = dsqrt(dsqrt(2*ekin))
	end

C	Subroutine that computes the lennard-jones interactions (potential & force)
C	INPUT: n (# of particles), r (position of the particles), l (size of the system)
C		   cutoff (parameter that sets the max distance to compute the interaction)
C	OUPUT: upot (potential energy of the system), f (force interacting on each particle)
	
	subroutine leonard_jones_force(n,r,l,cutoff,upot,f)
	implicit none
	integer i, j, n
	real*8 cutoff, upot, d2, d4, d6, d8, d12, d14
	real*8 cutoff2, cf4, cf6, cf12, l
	real*8 r(n,3), rij(3), f(n,3)
	
	cutoff2 = cutoff*cutoff
	cf4 = cutoff2*cutoff2
	cf6 = cutoff2*cf4
	cf12 = cf6*cf6
	upot = 0.d0
	f = 0.d0
	
	do i = 1, n
		do j = i+1, n
			rij = r(i,:) - r(j,:)
			call pbc(rij,l) 
			d2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
			if (d2.lt.cutoff2) then
				d4 = d2*d2
				d6 = d4*d2
				d8 = d6*d2
				d12 = d6*d6
				d14 = d8*d6
				upot = upot+4.d0*(1.d0/d12-1.d0/d6) - 4.d0*(1.d0/cf12-1.d0/cf6)
				f(i,:) = f(i,:) + (48.d0/d14 - 24.d0/d8)*rij
				f(j,:) = f(j,:) - (48.d0/d14 - 24.d0/d8)*rij
			endif
		enddo
	enddo
	end

C	Subroutine that computes the Periodic Boundary Conditions (for Lennard-Jones
C	subroutine)
C	INPUT: r (positions of the particles), l (size of the system)
C	OUTPUT: r (positions after performing the boundaries)

	subroutine pbc(r,l)
	integer i
	real*8 l
	real*8 r(3)
	
	do i = 1,3
		if (r(i).gt.(l/2.d0)) then 
			r(i) = r(i) - l 
		endif
		
		if (r(i).lt.(-l/2.d0)) then 
			r(i) = r(i) + l
		endif
	enddo
	end

C	Subroutine that computes the Periodic Boundary Conditions (for main program)
C	INPUT: r (positions of the particles), l (size of the system), n (# of particles)
C	OUTPUT: r (positions after performing the boundaries)	

	subroutine pbc2(r,l,n)
	implicit none

	real*8 l 
	integer n
	real*8 r(n,3)
	integer i, j


	do i=1,n
		do j=1,3

		if (r(i,j).gt.l/2.d0) then
			r(i,j)=r(i,j)-l
		endif

		if (r(i,j).lt.(-l/2.d0)) then
			r(i,j)=r(i,j)+l
		endif

		enddo
	enddo 
	end

C	Subroutine that do a step of the Velocity Verlet integrator
C	INPUT: r (positions of the particles), v (velocities of the particles), l (size of the 
C	system), n (# of particles), dt (time step), cutoff (parameter that sets the max
C	distance to compute the interaction)
C	OUTPUT: r (positions of the particles), v (velocities of the particles)
	
	subroutine verl_vel(r,v,l,cutoff,dt,n)
	implicit none
	
	integer n
	real*8 r(n,3), v(n,3), v_new(n,3), f(n,3)
	real*8 l, cutoff, dt, nu, mu, sigma, upot
	
	call leonard_jones_force(n, r, l, cutoff, upot, f)
	r = r + v*dt + 0.5*f*dt*dt
	v = v + 0.5*f*dt
	call leonard_jones_force(n, r, l, cutoff, upot, f)
	v = v + 0.5*f*dt
	end
	
C	Subroutine that do a step of the Euler integrator
C	INPUT: r (positions of the particles), v (velocities of the particles), l (size of the 
C	system), n (# of particles), dt (time step), cutoff (parameter that sets the max
C	distance to compute the interaction)
C	OUTPUT: r (positions of the particles), v (velocities of the particles), upot 
C	(potential energy of the system), f (force acting to each particle)	
	
	subroutine euler_step(n,r,v,l,cutoff,upot,f,dt)
	implicit none
	
	integer n
	real*8 l, cutoff, upot,dt
	real*8 r(n,3), v(n,3), f(n,3)
	
	call leonard_jones_force(n, r, l, cutoff, upot, f)
	r = r + v*dt + 0.5*f*dt*dt
	v = v + 0.5*f*dt
	
	end
	
C	Subroutine that computes the kinetic energy of the system
C	INPUT: n (# of particles), v (velocities of the particles)
C	OUTPUT: ekin (kinetic energy of the system)
	
	subroutine kin_energy(n,v,ekin)
	implicit none
	integer n
	real*8 ekin
	real*8 v(n,3), v2(n,3)
	
	v2 = 0.5*v*v
	ekin = sum(v2)
	
	end
	
C	Subroutine that performs a Box-Muller random numbers
C	INPUTS: sigma (variance), x1, x2 (two random uniform numbers)
C	OUPUTS: xout1, xou2 (two random normal numbers)

	subroutine box_muller(sigma, x1, x2, xout1, xout2)
	implicit none
	
	real*8 sigma,x1,x2,xout1,xout2,pi

    	pi=4.d0*atan(1.d0)

    	xout1=sigma*sqrt(-2.d0*(log(1.d0-x1)))*cos(2.d0*pi*x2)
    	xout2=sigma*sqrt(-2.d0*(log(1.d0-x1)))*sin(2.d0*pi*x2)
	end
	
C	Subroutine that applies a Andersen Thermostat to the system
C	INPUT: n (# of particles), v (velocities of the particles), nu (probability of 
C	changing the velocity), sigma (variance)
C	OUTPUT: v (velocities of the particles)

	subroutine therm_Andersen(n,v,nu,sigma)
	implicit none
	
    	real*8 nu,sigma,x1,x2,xout1,xout2,x3,xout3,xout4,x4,x5
    	real*8 v(n,3)
    	integer i, n 

    	do i=1,n
        	x1=rand()
        	x2=rand()
		x3=rand()
		x4=rand()
		x5=rand()

		if (x5.lt.nu) then
			call box_muller(sigma,x1,x2,xout1,xout2)
			call box_muller(sigma,x3,x4,xout3,xout4)
		    	v(i,:)=(/xout1,xout2,xout3/)
		endif
	    enddo
	end
	
C	Subroutine that computes the instantaneous temperature of the system
C	INPUT: ekin (kinetic energy), n (# of particles)
C	OUTPUT: temp (temperature of the system)
 
	subroutine temp_inst(ekin,n,temp)
    	implicit none
    	integer n, i
    	real*8 temp,ekin

   	temp=2.d0/3.d0*ekin/dble(n)

	end
	
C	Subroutine that computes the instantaneous pressure of the system
C	INPUT: n (# of particles), rho (density), l (size of the system), r (positions of
C	the particles), f (force acting on the particles), t_inst (temperature acting on
C	the system
C	OUTPUT: P (pressure of the system)

	subroutine p_inst(n, rho, l, r, f, t_inst, P)
	implicit none
	integer n
	real*8 r(n,3), f(n,3)
	real*8 t_inst, P, l, rho
	
	P = rho*t_inst + (sum(abs(f*r))/dble(n))/(3.d0*l*l*l)
	end
	
C	Subroutine that computes the MSD
C	INPUT: n (# of particles), r_ini (positions for t = 0), r_fin (positions a time t)
C	OUTPUT: msd_val (value of the MSD)
	subroutine msd(n, r_ini, r_fin, msd_val)
	implicit none
	integer n, i
	real*8 r_ini(n,3), r_fin(n,3), dr(n)
	real*8 drx, dry, drz, drtotal, msd_val
	
	drtotal = 0.0
	do i = 1, n
		drx = r_fin(i,1) - r_ini(i,1)
		dry = r_fin(i,2) - r_ini(i,2)
		drz = r_fin(i,3) - r_ini(i,3)
		dr(i) = drx*drx + dry*dry + drz*drz
		drtotal = drtotal + dr(i)
	enddo
	msd_val = 1.d0/dble(n)*(drtotal)
	end
	
