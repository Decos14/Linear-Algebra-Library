module linalg
    USE, INTRINSIC :: IEEE_ARITHMETIC
    implicit none


contains

    function Transpose (A) result(AT) 
        implicit none
        real, intent(in) :: A(:,:)
        real::AT(size(A,2),size(A,1)) 
        integer :: i,j,n,m
        n = size(A,1)
        m = size(A,2)
        do i = 1,n
            do j = 1,m 
                AT(j,i) = A(i,j)
            enddo
        enddo
    end function Transpose

    function MatMul(A,B) result(C)
        implicit none
        real, intent(in) :: A(:,:), B(:,:)
        real:: C(size(A,1), size(B,2))
        integer :: i,j,k
        
        C(:,:) = 0
        do i = 1,size(A,1)
            do j = 1,size(B,2)
                do k = 1,size(A,2)
                    C(i,j) = C(i,j) + (A(i,k) * b(k,j))
                enddo
            enddo
        enddo
    end function MatMul

    subroutine LU_Factor(A, L, U)
        implicit none
        real, intent(in) :: A(:,:)
        real:: L(:,:),U(:,:)
        real :: B(size(A,1),size(A,1))
        integer::i,j,k
        B = A           
        L(:,:) = 0
        U(:,:) = 0
        do k = 1, size(B,1)-1
            do i = k+1, size(B,1)
                if (B(k,k) == 0) then
                    B(i,k) = 0
                else
                    B(i,k) = B(i,k)/B(k,k)
                endif
                do j = k+1,size(B,1)
                    B(i,j) = B(i,j) - B(i,k)*B(k,j)
                enddo
            enddo
        enddo
        
        do i = 1, size(U,2)
            do j = i,size(U,2)
                L(j,i) = B(j,i)
                U(i,j) = B(i,j)
                if (i == j) then
                    L(i,j) = 1
                endif
            enddo
        enddo
    
        
    end subroutine LU_Factor

    function Determinant(A) result(d)
        implicit none
        real, intent(in) :: A(:,:)
        real:: d
        integer::i
        real::L(size(A,1), size(A,1)), U(size(A,1), size(A,1))
        call LU_Factor(A,L,U)

        d = 1
        do i = 1,size(U,1)
            d = d*U(i,i)
        enddo

    end function Determinant

    ! Produces an nxn Identity matrix
    function Identity(n) result(I)
        implicit none
        integer, intent(in) :: n
        real:: I(n,n)
        integer::j
        I(:,:) = 0
        do j = 1,n
            I(j,j) = 1
        enddo
    end function Identity

    ! takes a matrix A to an integer power
    function Mat_Pwr(A,n) result(B)
        implicit none
        integer, intent(in) :: n
        real, intent(in) :: A(:,:)
        real:: B(size(A,1), size(A,2))
        integer::i
        if(n == 0) then
            B = Identity(size(A,1))
        else
            B = A
            do i = 1,n-1
                B = MatMul(A,B)
            enddo
        endif
    end function Mat_Pwr

    function Trace(A) result(n)
        implicit none
        real, intent(in) :: A(:,:)
        real:: n
        integer::i
        n = 0

        do i = 1,size(A,1)
            n = n + A(i,i)
        enddo
    end function Trace

    function Is_Square(A) result(b)
        real, intent(in) :: A(:,:)
        logical:: b
        b = (size(A,1) == size(A,2))
    end function Is_Square

    function Is_Symmetric(A) result(b)
        real, intent(in) :: A(:,:)
        logical:: b
        integer:: i,j
        b = .False.
        
        if (Is_Square(A)) then
            b = .True.
            do i = 1,size(A,1)
                do j = 1,size(A,2)
                    b = (b .and. (A(i,j) == A(j,i)))
                enddo
            enddo
        endif
    end function Is_Symmetric

    !Solves upper tringular systems
    function U_solve(U,b) result(x)
        real, intent(in) :: U(:,:), b(:,:)
        real:: x(size(b),1)
        real::temp
        integer::i,j
        x(:,1) = 0

        do i = size(U,1), 1, -1
            temp = b(i,1)
            do j = size(U,1), i, -1
                temp = temp - x(j,1)*U(i,j)
            enddo
            x(i,1) = temp/U(i,i)
        enddo
    end function U_solve

    ! Solves lower triangular systems
    function L_solve(L,b) result(x)
        real, intent(in) :: L(:,:), b(:,:)
        real:: x(size(b),1)
        real::temp
        integer::i,j
        x(:,1) = 0

        do i = 1, size(L,1)
            temp = b(i,1)
            do j = 1, i
                temp = temp - x(j,1)*L(i,j)
            enddo
            x(i,1) = temp/L(i,i)
        enddo
    end function L_solve

    function LU_solve(A,b) result(x)
        real, intent(in) :: A(:,:), b(:,:)
        real:: x(size(A,2),1), y(size(A,1),1)
        real::L(size(A,1), size(A,1)), U(size(A,1),size(A,2))

        call LU_Factor(A,L,U)
        
        y = L_solve(L,b)
        x = U_solve(U,y)
    end function LU_solve

    function Inverse(A) result(AT)
        real, intent(in) :: A(:,:)
        real:: AT(size(A,1), size(A,1))
        real::temp(size(A,1), 1), ei(size(A,1), 1)
        integer::i,j

        
        do i = 1, size(A,1)
            ei(:,:) = 0
            ei(i,1) = 1
            temp = LU_solve(A,ei)
            do j = 1, size(A,2)
                AT(j,i) = temp(j,1)
            enddo
        enddo

    end function Inverse

    ! Finds the innter product of two matricies defined by matrix A
    function Inner_Product(u,v,A) result(c)
        real, intent(in) :: A(:,:), u(:,:), v(:,:)
        real::temp(1,1)
        real::c

        temp = MatMul(MatMul(u,A),v)
        c = temp(1,1)
    end function Inner_Product

    function Inf_Norm_Vect(v) result(c)
        real, intent(in) :: v(:,:)
        real:: c
        integer::i
        c = v(1,1)

        do i = 1,size(v)
            if(v(i,1) > c) then
                c = v(i,1)
            endif
        enddo   
    end function Inf_Norm_Vect

    ! Computes the p-norm of a vector for any integer p
    function P_Norm_Vect(v,p) result(c)
        real, intent(in) :: v(:,:)
        integer, intent(in) :: p
        real:: c
        integer::i
        c = 0

        do i = 1,size(v)
            c = c + (abs(v(i,1))**p)
        enddo

        c = (c**(1.0/p))
    end function P_Norm_Vect

    function One_Norm_Mat(A) result(c)
        real, intent(in) :: A(:,:)
        real:: c
        real::v(size(A,1),1)
        real::s
        integer::i,j
        c = 0

        do i = 1,size(A,1)
            s = 0
            do j = 1, size(A,2)
                s = s + A(i,j)
            enddo
            v(i,1) = abs(s)
        enddo

        c = Inf_Norm_Vect(v)
    end function One_Norm_Mat

    function Inf_Norm_Mat(A) result(c)
        real, intent(in) :: A(:,:)
        real:: c
        real::v(size(A,1),1)
        real::s
        integer::i,j
        c = 0

        do i = 1,size(A,2)
            s = 0
            do j = 1, size(A,1)
                s = s + A(i,j)
            enddo
            v(i,1) = abs(s)
        enddo

        c = Inf_Norm_Vect(v)
    end function Inf_Norm_Mat

    ! Computes the Frobenius norm of any given Matrix
    function F_Norm_Mat(A) result(c)
        real, intent(in) :: A(:,:)
        real:: c
        integer::i,j

        c = 0
        do i = 1, size(A,1)
            do j = 2, size(A,2)
                c = c + (A(i,j)**2)
            enddo
        enddo

        c = c**(1.0/2.0)
    end function F_Norm_Mat

    ! Finds the least squares solution to a system 
    function Least_Squares(A,b) result(x)
        real, intent(in) :: A(:,:), b(:,:)
        real:: x(size(A,2), 1)
        
        x = LU_solve(MatMul(Transpose(A),A), MatMul(Transpose(A),b))
    end function Least_Squares

    subroutine QR_Factor(A,Q,R)
        real, intent(in) :: A(:,:)
        real, intent(inout) :: Q(:,:), R(:,:)
        integer::i,j,k
        Q = A

        do j = 1,size(A,2)
            R(j,j) = 0
            do i = 1,size(A,1)
                R(j,j) = R(j,j) + (Q(i,j)**2)
            enddo
            R(j,j) = sqrt(R(j,j))
            do i = 1,size(A,1)
                Q(i,j) = Q(i,j)/R(j,j)
            enddo
            do k = j+1,size(A,2)
                R(j,k) = 0
                do i = 1,size(A,1)
                    R(j,k) = R(j,k) + (Q(i,j)*Q(i,k))
                enddo
                do i = 1,size(A,1)
                    Q(i,k) = Q(i,k) - Q(i,j)*R(j,k)
                enddo
            enddo
        enddo
    end subroutine QR_Factor

    function Eigenvalues(A) result(v)
        real, intent(in) :: A(:,:)
        real::v(size(A,1),1)
        real::Q(size(A,1), size(A,1)),R(size(A,1), size(A,1)), Ak(size(A,1), size(A,1))
        real::v_prev(size(A,1),1),eps,err
        integer::i
        eps = 0.00000000001
        err = 1
        Ak = A

        do while(err > eps)
            v_prev = v
            call QR_Factor(Ak,Q,R)
            Ak = MatMul(R,Q)
            do i = 1,size(A,1)
                v(i,1) = Ak(i,i)
            enddo
            err = P_Norm_Vect(v - v_prev,2)
        enddo
    end function Eigenvalues

    !Uses the shifted inverse power method to find the nearest eigenvector to a given eigenvalue guess
    function SI_Power_Method(A, e) result(v)
        real, intent(in) :: A(:,:), e
        real:: v(size(A,1),1)
        real::v_prev(size(A,1),1),err,eps, coeff(size(A,1), size(A,1))
        integer::i
        eps = 0.0000001
        err = 1
        v(:,1) = 2
        v(1,1) = 1
        coeff = A

        do i = 1,size(A,1)
            coeff(i,i) = coeff(i,i) - e
        enddo
        coeff = Inverse(coeff)

        do while(err > eps)
            v_prev = v
            v = MatMul(coeff,v_prev)
            v = v/P_Norm_Vect(v,2)
            err = P_Norm_Vect(v - v_prev,2)
        enddo
    end function SI_Power_Method

    function Eigenvectors(A) result(E)
        real, intent(in) :: A(:,:)
        real:: E(size(A,1), size(A,1))
        real::v(size(A,1),1), eig(size(A,1),1)
        integer::i,j

        v = Eigenvalues(A)

        do i = 1, size(v)
            eig = SI_Power_Method(A,v(i,1)-0.1)
            do j = 1, size(E,1)
                E(j,i) = eig(j,1)
            enddo
        enddo
    end function Eigenvectors

    ! Simultanuesly finds both eigenvalues and eigenvectors for posative definite matricies
    subroutine Pos_QR_eig(K, e_val, e_vect)
        real, intent(in) :: K(:,:)
        real, intent(inout) :: e_vect(:,:), e_val(:,:)
        real::Q(size(K,1), size(K,1)),R(size(K,1), size(K,1)), Ki(size(K,1),size(K,1)) 
        real::e_prev(size(e_val),1),eps,err
        integer::i
        Ki = K
        eps = 0.00000000001
        err = 1
        e_vect = Identity(size(e_vect,1))

        do while(err > eps)
            e_prev = e_val
            call QR_Factor(Ki,Q,R)
            e_vect = MatMul(e_vect,Q)
            Ki = MatMul(R,Q)
            do i = 1,size(Ki,1)
                e_val(i,1) = Ki(i,i)
            enddo
            err = P_Norm_Vect(e_val - e_prev,2)
        enddo
    end subroutine Pos_QR_eig
    
    !Inverse of a diagonal matrix (not neccesarily square)
    function Diag_Inv(A) result(A_inv)
        real, intent(in) :: A(:,:)
        real::A_inv(size(A,2), size(A,1))
        integer::i
        A_inv(:,:) = 0

        do i = 1, size(A,2)
            A_inv(i,i) = 1.0/A(i,i)
        enddo
    end function Diag_Inv

    function Singular_Values(A) result(s)
        real, intent(in) :: A(:,:)
        real::s(size(A,2),1)
        real::K(size(A,2), size(A,2))
        integer::i

        K = MatMul(Transpose(A),A)
        s = Eigenvalues(K)

        do i = 1, size(s)
            s(i,1) = sqrt(s(i,1))
        enddo
    end function Singular_Values

    !Singular value decomposition of a matrix
    subroutine SVD(A,U,Sigma,V)
        real, intent(in) :: A(:,:)
        real, intent(inout) :: U(:,:),Sigma(:,:),V(:,:)
        real:: K(size(A,2), size(A,2)), e_val(size(K,1), 1), ui(size(U,1),1), vi(size(V,1),1)
        integer::i,j

        K = MatMul(Transpose(A), A)
        call Pos_QR_eig(K,e_val,V)

        Sigma(:,:) = 0
        do i = 1, size(Sigma,1)
            Sigma(i,i) = sqrt(e_val(i,1))
        enddo

        do i = 1,size(U,2)
            do j = 1, size(V,1)
                vi(j,1) = V(j,i)
            enddo
            ui = MatMul(A,vi)/Sigma(i,i)
            do j = 1, size(U,1)
                U(j,i) = ui(j,1)
            enddo
        enddo
    end subroutine SVD

    ! Moore-Penrose Psuedoinverse of a matrix
    function MP_Inv(A) result(A_inv)
        real, intent(in) :: A(:,:)
        real::A_inv(size(A,2), size(A,1))
        real::U(size(A,1),size(A,2)), E(size(A,2),size(A,2)), V(size(A,2),size(A,2))

        call SVD(A,U,E,V)

        A_inv = MatMul(V,MatMul(Diag_Inv(E),Transpose(U)))
    end function MP_Inv

    function factorial(n) result(f)
        real, intent(in) :: n
        real::f, i
        f = 1
        i = n
        do while(i > 0)
            f = f * i
            i = i - 1
        enddo
    end function factorial

    ! Computes the matrix e^(tA)
    function exp(A, t) result(E)
        real, intent(in) :: A(:,:), t
        real::E(size(A,1), size(A,1)), pwr(size(A,1), size(A,2)), eps, err, f_val, f_prev, n
        eps = 0.0001
        err = 1
        n = 1
        E = Identity(size(E,1))

        do while(err > eps)
            pwr = Mat_Pwr(A,int(n))
            if(pwr(1,1) /= pwr(1,1) .or. ieee_is_finite(pwr(1,1))) then
                E = E + (((t**n)/factorial(n))*pwr)
                f_val = F_Norm_Mat(E)
                err = abs(f_val - f_prev)
            else
                err = 0
            endif
            n = n + 1
        enddo
    end function exp

    subroutine print_matrix(A)
        real, intent(in) :: A(:,:)
        integer :: i

        do i = 1, size(A,1)
            print *, A(i,:)
        end do

    end subroutine print_matrix

    subroutine populate(A)
        real, intent(inout) :: A(:,:)
        integer::i,j

        do i = 1,size(A,1)
            do j = 1,size(A,2)
                A(i,j) = size(A,2)*(i-1) + j
            end do
        end do
    end subroutine populate


end module 

program use_mod
    use linalg
    implicit none

    real::A(3,3), A_exp(3,3)
    integer :: i

    call populate(A)

    ! errors in populate matrix going to inf????
    ! errors in incorrect solution by orders of magnitude????
    A_exp = exp(A,1.0)

    call print_matrix(A)
    Write(*,*)
    call print_matrix(A_exp)
    Write(*,*)


end program use_mod