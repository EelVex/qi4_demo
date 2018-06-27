program demo
   use QI4
   implicit none

   ! declarations for master equation solution
   real(dp)::t
   real(dp),parameter::t_fin = 1e-5
   real(dp),parameter::dt = 0.1e-5
   
   complex(dp)::rho(4,4)=reshape(&
    [1,0,0,0,&
     0,1,0,0,&
     0,0,1,0,&
     0,0,0,1],shape(rho))
   complex(dp)::sigma(4,4) = 0
   complex(dp)::A(2,2)=reshape(&
    [33,24,&
     48,57],shape(A))
   complex(dp)::B(3,3)=reshape(&
    [6, 9, 13, &
     12,17,30, &
     7, 9, 21],shape(B))
   complex(dp)::psip(2,2)=reshape(&
    [0.5_dp, 0.5_dp, &
     0.5_dp, 0.5_dp],shape(psip))
   complex(dp)::psim(2,2)=reshape(&
    [0.5_dp, 0.0_dp, &
     0.0_dp, 0.5_dp],shape(psip))

   rho = ket(bell1) .o. bra(bell1)
   sigma(1,1) = 0.75_dp; sigma(2,2) = 0.10_dp ;  sigma(3,3) = 0.10_dp ; sigma(4,4) = 0.05_dp
   
   print *, "======= some constants from module qi_constants: ============"
   print *, "pi is ", PI_c
   print *, "e is ", E_c
   print *, "q (charge) is ", Q_c
   print *, "pauli x matrix: "; call Mdisplay(PAULIX)
   print *, "pauli y matrix: "; call Mdisplay(PAULIY)
   print *; print *

   print *, "======= some QI functions from module qi_functions: ============"
   print *, "concurrence(rho) = ", concurrence(rho)
   print *, "concurrence(sigma) = ", concurrence(sigma)
   print *, "fidelity(rho, sigma) = ", fidelity(rho, sigma)
   print *, "entropy of a pure state psip", entropy_von_neumann(psip)
   print *, "entropy of a mixture psim", entropy_von_neumann(psim)
   print *, "where the matrices are ";
   print *, "rho = "; call Mdisplay(rho); print *, "sigma = "; call Mdisplay(sigma)
   print *; print *

   print *, "======= some mathematical functions from module qi_functions: ============"
   print *, "let A = "; call Mdisplay(A)
   print *, "let B = "; call Mdisplay(B)
   print *, "eigenvalues (A) = ",  real(eigenvalues(A))
   print *, "eigenvectors (A) = ";  call Mdisplay(eigenvectors(A))

   print *, "square root of A = "; call Mdisplay(Msqrt(A))
   print *, "square root of B = "; call Mdisplay(Msqrt(B))
   print *, "square root of PAULIX = "; call Mdisplay(Msqrt(PAULIX))
   print *, "log of PAULIX = "; call Mdisplay(Mlog(PAULIX))
   print *, "operator '.o.' for kronecker product C = A .o. B"; call Mdisplay(A .o. B)

   print *, "partial trace of C over A"
   call Mdisplay(ptrace(A .o. B, 2, 2))

   stop
   
   print *; print *

   print *, "======= some tools from module qi_tools: ============"
   print *, "*All above matrices are printed through 'call Mdisplay(matrix)'"
   print *
   print *, "Ground state rho for a 3x3 system: rho_ground(3) = "; call Mdisplay(rho_ground(3))
   print *, "kets, bras and tensor product: ket(bell1) .o. bra(bell1) = "; call Mdisplay(ket(bell1)  .o. (bra(bell1)))
   print *; print *

   print *, "======= solving a master equation with module integrators: ============"
   ! initialize rho to ground state
   rho = rho_ground(4); t = 0
   ! replace with something like:
   ! initiliaze() ! ? set default variables?
   ! solve(t_0=0,t_fin=t_fin,dt=dt,steps=400,rhodot=eqs,method='rk4')
   do while (t < t_fin)
        t = t + dt
        rho = RKH4_step(t, rho, rhodot, dt=dt)
        print *, t
   enddo

   print *, "at t = ", t
   call Mdisplay(rho)

   contains
           function rhodot(t, rho)
                      real(dp),intent(in)::t
                      complex(dp),dimension(:,:),intent(in)::rho
                      complex(dp),dimension(size(rho,1),size(rho,2))::rhodot

                      rhodot = exp(-hbar_c*t*(0,1));

           end function


end program demo
