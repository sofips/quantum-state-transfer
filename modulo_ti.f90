!-----------------------------------------
 module modulo_ti
 use precision , pr => dp
 implicit none
 
 contains
 
!---------------------------------------------------
 subroutine inicializar(J,J_0)
!---------------------------------------------------
!Inicializa los acoplamientos entre spines,
!como J=1 si n=3, si n>3 en funcion de los anteriores
!--------------------------------------------------
 real(pr), dimension(:), intent(inout) :: J, J_0
 integer                               :: n_j,i
 
  n_j = size(J)
  
  if (n_j==3) then 
    J = 1._pr
  else
  
    do i= 1,n_j/2 !int redondea para abajo
         J(i)=J_0(i)
    end do
  
    do i= (n_j/2+1),n_j
         J(i)=J_0(i-1)
    end do  
   
  end if
   
 end subroutine

!------------------------------------------------ 
 subroutine grad(m,H,J,gm)
!------------------------------------------------
 real(pr), dimension(:), intent(in)      :: J
 real(pr), dimension(:,:), intent(inout) :: H
 integer, intent(in)                     :: m
 real(pr), intent(out)                   :: gm
 real(pr), dimension(size(J))            :: Jp, Jl
 real(pr)                                :: beta, lambda, Lp, Ll, Fp, Fl
 integer                                 :: n_h
 !defino los parametros--------
  
  beta   = 0.5_pr
  lambda = 1E-6_pr
  n_h    = size(H,1)
 !------------------------------- 
 !calculo la variacion
 !-------------------------------
   Jp= J
  
   Jp(m) = J(m) + beta  
 
   Jl= J
 
   Jl(m) = J(m) - beta  
 !------------------------------------------------------------
 !calculo la fidelidad y la funcion de perdida con variaciones
 !-------------------------------------------------------------
  call fidelidad(Jp,H,n_h,Fp)
  
   Lp = 1-Fp +lambda*maxval(Jp) 
   
  call fidelidad(Jl,H,n_h,Fl)
   
   Ll = 1-Fl+lambda*maxval(Jl)
 !----------------------------------------
 !calculo el gradiente de la fn de perdida
 !----------------------------------------
   gm = (Lp - Ll)/(2._pr*beta)
   
 end subroutine grad 
 
!--------------------------------------------
 subroutine fidelidad(J,H,n_h,F)
!--------------------------------------------
 real(pr), dimension(:), intent(in)      :: J
 real(pr), dimension(:,:), intent(inout) :: H
 integer, intent(in)                     :: n_h
 real(pr), intent(out)                   :: F
 real(pr), dimension(n_h)                :: w, c1cn
 real(pr), allocatable                   :: work(:)
 integer                                 :: lwork, lwmax, INFO, i
 character                               :: JOBZ, UPLO
 real(pr)                                :: Fr, Fi
 
  !-----------------------------------------------------
  !construyo el hamiltoniano con los acoplamientos dados
  !------------------------------------------------------
   call hamiltoniano(J,H)
 
  !Para calcular F, necesito autovalores y autovectores, uso dsyev de las lib. Lapack
  
  !-----------------------------------------------------------------------
  !'V' calcula los autovalores, 'U' se guarda solo la triangular superior
  !-----------------------------------------------------------------------
   JOBZ = 'V'
   UPLO = 'U'
  !-----------------------------------
  !genera un espacio de trabajo
  !-----------------------------------
   lwmax = 10000
   lwork = 5
   allocate(work(lwmax))
  
  !----------------------------------
  !optimizar el espacio de trabajo
  !------------------------------------
   lwork = -1
 
   call dsyev( JOBZ,  UPLO, n_h, H, n_h, w, work, lwork, INFO )
 
   lwork = work(1) 
  !-----------------------------------------
  !diagonalizar la matriz
  !-----------------------------------------
   call dsyev( JOBZ,  UPLO, n_h, H, n_h, w, WORK, lwork, INFO )
 
    do i = 1,n_h
         c1cn(i) =H(1,i)*H(n_h,i)
    end do
    
    Fr = 0._pr
    Fi = 0._pr

    do i = 1,n_h
       Fr = Fr + cos(w(i))*c1cn(i)
       
       Fi = Fi + sin(w(i))*c1cn(i) 
    end do
    
    F = Fr*Fr + Fi*Fi
    
   
  
 end subroutine fidelidad
  
!---------------------------------------------------
 subroutine hamiltoniano(J,H)
!---------------------------------------------------
  real(pr), dimension(:), intent(in)      :: J
  real(pr), dimension(:,:), intent(inout) :: H
  integer                                 :: n_j,i
   
   n_j =size(J)
   H =0._pr

   do i = 1,n_j
     H(i,i+1) = J(i)*(-0.5_pr)
     H(i+1,i) = H(i,i+1)
   end do
   
  end subroutine hamiltoniano 
 !------------------------------------------------ 
 end module
