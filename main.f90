!-------------------------------
 program main
!-------------------------------
 use precision, only: pr => dp
 use mtmod,only: grnd
 use modulo_ti
 
 implicit none
 real(pr)              :: alpha, tol,F, gm
 integer               :: nmax, kmax,i,k,n_h, n_j, m, idum
 real(pr), allocatable :: J(:), H(:,:), J0(:)
 character(len=15)     :: filename
 character(len=15)     :: filename2
 character (len=20)     :: aux

 
 
 open (7,file='infid.dat',status='unknown',form='formatted')
 open (13,file='eff.dat',status='unknown',form='formatted')

 !--------------------------------------
 !Parameter definition
 !--------------------------------------
  alpha = 0.05
  tol   = 1E-4_pr
  nmax  = 40
  kmax  = 1000000
  allocate (J0(3))

 !-------------------------------------
 !Loop over system's dimensions
 !-------------------------------------
 
 
  do n_j = 3, nmax
        
      n_h = n_j+1
      
      write(*,*) n_j
      
      allocate(J(1:n_j), H(1:n_h,1:n_h))
     
      call inicializar(J,J0) !initialize couplings

     
      !write(*,*) J, J0
      
      do k = 1, kmax
        
          m = int(grnd()*n_j)+1 !choose random coupling to modify
                
          call grad(m,H,J,gm) !calculates associated gradient
                
          J(m) = J(m)-alpha*gm !modifies chosen coupling
        
          call fidelidad(J,H,n_h,F) !recalculates fidelity
          
          !------------------------------------------------------
          ! generates fidelity evolution file 
          !------------------------------------------------------
           write(aux,'(I15)') n_h
           filename2 = 'evol_' // trim(adjustl(aux))//'.dat'
           open (10,file= filename2,status='unknown',form='formatted')
           write(10,'(I15,E16.8)') k,1-F
          !---------------------------------------------------------
          
          
        if ((1._pr-F) < tol) exit !  tolerance check
        
      end do  
      
      
      write(7,*) n_h, 1-F 
           
      write(13,'(I15,2E16.8)') n_h, maxval(J), acos(-1._pr)*n_h*4  !efficiency file

      !----------------------------------------------------------------
      ! Final couplings file
      !----------------------------------------------------------------
      
       write(aux,'(I15)') n_h
       filename = 'J_' // trim(adjustl(aux))//'.dat'
       open (9,file= filename,status='unknown',form='formatted')
             
       do i = 1,n_j
            write(9,*) i,J(i)
       end do
      !------------------------------------------------------------- ---------------
      !saves previous couplings to initialize next dimension
      !--------------------------------------------------------------------------------
      
       deallocate(J0) 
     
       allocate (J0(n_j))
    
       J0 = J 
     
       deallocate(J,H)
       
  end do
  
  end program
