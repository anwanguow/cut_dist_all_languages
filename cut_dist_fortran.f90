module cutnorm_mod
  implicit none
  private
  public :: cut_distance

contains

  function cut_distance(G, H) result(sdp)
    implicit none
    real(8), intent(in) :: G(:,:), H(:,:)
    real(8) :: sdp
    integer :: n
    real(8), allocatable :: diff(:,:)
    n = size(G,1)
    allocate(diff(n,n))
    diff = (G - H) / dble(n*n)
    sdp = cut_norm(diff)
    deallocate(diff)
  end function cut_distance

  function cut_norm(A) result(sdp_result)
    implicit none
    real(8), intent(in) :: A(:,:)
    real(8) :: sdp_result
    integer :: n1, N, n2, p, temp_int, i, j
    real(8) :: temp
    real(8), allocatable :: A_col_sum(:), A_row_sum(:)
    real(8), allocatable :: A_ext(:,:), A_final(:,:)
    real(8), allocatable :: x0(:,:), x(:,:)
    real(8), allocatable :: nrmx0(:)
    real(8), allocatable :: U(:,:), V(:,:), M(:,:)

    n1 = size(A,1)
    allocate(A_col_sum(n1), A_row_sum(n1))
    do j = 1, n1
      A_col_sum(j) = sum(A(:,j))
    end do
    do i = 1, n1
      A_row_sum(i) = sum(A(i,:))
    end do
    temp = sum(A_col_sum)

    allocate(A_ext(n1, n1+1))
    A_ext(:,1:n1) = A
    do i = 1, n1
      A_ext(i, n1+1) = -A_row_sum(i)
    end do

    N = n1 + 1
    allocate(A_final(N,N))
    A_final(1:n1,:) = A_ext
    do j = 1, n1
      A_final(N,j) = -A_col_sum(j)
    end do
    A_final(N,N) = temp

    temp = sqrt(2.0d0*dble(n1)) / 2.0d0
    temp_int = nint(temp)
    if (temp_int > 100) then
      p = 100
    else if (temp_int < 1) then
      p = 1
    else
      p = temp_int
    end if

    n2 = 2 * N
    allocate(x0(p,n2))
    call random_normal(x0)

    allocate(nrmx0(n2))
    do j = 1, n2
      nrmx0(j) = sqrt(sum(x0(:,j)**2))
      if (nrmx0(j) > 1.0d-12) x0(:,j) = x0(:,j) / nrmx0(j)
    end do
    deallocate(nrmx0)

    call optimize(x0, A_final, x)

    allocate(U(p,N), V(p,N))
    U = x(:,1:N)
    V = x(:,N+1:n2)
    allocate(M(N,N))
    M = matmul(transpose(U), V)
    sdp_result = abs(sum(A_final * M)) / 4.0d0

    deallocate(A_col_sum, A_row_sum, A_ext, A_final, x0, x, U, V, M)
  end function cut_norm

  subroutine cut_norm_quad(V, A, f, g)
    implicit none
    real(8), intent(in) :: V(:,:), A(:,:)
    real(8), intent(out) :: f
    real(8), allocatable, intent(out) :: g(:,:)
    integer :: p, n
    real(8), allocatable :: Us(:,:), Vs(:,:)

    p = size(V,1)
    n = size(A,1)
    allocate(Us(p,n), Vs(p,n))
    Us = V(:,1:n)
    Vs = V(:,n+1:2*n)

    allocate(g(p,2*n))
    g(:,1:n)   = 2.0d0 * matmul(Vs,    transpose(A))
    g(:,n+1:2*n) = 2.0d0 * matmul(Us,    A)
    f = ( sum(g(:,1:n) * Us) + sum(g(:,n+1:2*n) * Vs) ) / 2.0d0

    deallocate(Us, Vs)
  end subroutine cut_norm_quad

  subroutine optimize(x, A, x_out)
    implicit none
    real(8), intent(inout) :: x(:,:)
    real(8), intent(in)    :: A(:,:)
    real(8), allocatable, intent(out) :: x_out(:,:)

    integer :: p, n2, itr, mxitr, nt, n_iter_range, ls_iter, i
    real(8) :: xtol, ftol, gtol, rho, eta, gamma, tau, tau_orig
    real(8), allocatable :: nrmx(:)
    real(8) :: f, fp, deriv, XDiff, FDiff, sy
    real(8) :: Q, Qp, Cval, temp_tau
    real(8), allocatable :: g(:,:), gp(:,:), dtX(:,:), dtXP(:,:), s(:,:), y(:,:)
    real(8), allocatable :: xtg(:), gg(:), xx(:), xxgg(:)
    real(8), allocatable :: a1(:), a2(:), xp(:,:)
    real(8), allocatable :: crit(:,:), mcrit(:)
    real(8) :: nrmG

    mxitr = 600
    xtol   = 1.0d-8;  ftol = 1.0d-10; gtol = 1.0d-8
    rho    = 1.0d-4;  eta  = 0.1d0; gamma = 0.85d0
    tau    = 1.0d-3;  tau_orig = tau
    nt     = 5

    p = size(x,1)
    n2 = size(x,2)

    allocate(nrmx(n2))
    do ls_iter = 1, n2
      nrmx(ls_iter) = sqrt(sum(x(:,ls_iter)**2))
      if (nrmx(ls_iter) > 1.0d-12) x(:,ls_iter) = x(:,ls_iter) / nrmx(ls_iter)
    end do
    deallocate(nrmx)

    call cut_norm_quad(x, A, f, g)

    allocate(xtg(n2), gg(n2), xx(n2), xxgg(n2))
    do ls_iter = 1, n2
      xtg(ls_iter) = sum(x(:,ls_iter) * g(:,ls_iter))
      gg(ls_iter)   = sum(g(:,ls_iter)**2)
      xx(ls_iter)   = sum(x(:,ls_iter)**2)
      xxgg(ls_iter) = xx(ls_iter) * gg(ls_iter)
    end do

    allocate(dtX(p,n2))
    do ls_iter = 1, n2
      dtX(:,ls_iter) = x(:,ls_iter)*xtg(ls_iter) - g(:,ls_iter)
    end do
    nrmG = sqrt(sum(dtX**2))

    allocate(crit(mxitr,3))
    crit = 0.0d0
    Q    = 1.0d0
    Cval = f

    ! 主迭代
    do itr = 1, mxitr
      dtXP = dtX; xp = x; gp = g; fp = f
      deriv = rho * nrmG**2

      ls_iter = 0
      do
        ls_iter = ls_iter + 1
        temp_tau = tau / 2.0d0
        allocate(a1(n2), a2(n2))
        do i = 1, n2
          a1(i) = ((1.0d0 + temp_tau*xtg(i))**2 - temp_tau**2*xxgg(i)) &
                  / (1.0d0 + temp_tau**2*(-xtg(i)**2 + xxgg(i)))
          a2(i) = -tau * xx(i) / (1.0d0 + temp_tau**2*(-xtg(i)**2 + xxgg(i)))
        end do
        do i = 1, n2
          x(:,i) = xp(:,i)*a1(i) + gp(:,i)*a2(i)
        end do
        deallocate(a1, a2)
        call cut_norm_quad(x, A, f, g)
        if (f <= Cval - tau*deriv .or. ls_iter >= 5) exit
        tau = eta * tau
      end do

      ! 更新梯度与度量
      do i = 1, n2
        xtg(i) = sum(x(:,i) * g(:,i))
        gg(i)   = sum(g(:,i)**2)
        xx(i)   = sum(x(:,i)**2)
        xxgg(i) = xx(i) * gg(i)
      end do
      do i = 1, n2
        dtX(:,i) = x(:,i)*xtg(i) - g(:,i)
      end do
      nrmG = sqrt(sum(dtX**2))
      s     = x - xp
      XDiff = sqrt(sum(s**2)) / sqrt(dble(p))
      FDiff = abs(fp - f) / (abs(fp) + 1.0d0)
      crit(itr,1) = nrmG
      crit(itr,2) = XDiff
      crit(itr,3) = FDiff
      allocate(mcrit(3))
      mcrit = 0.0d0
      n_iter_range = min(nt, itr)
      do i = itr-n_iter_range+1, itr
        mcrit = mcrit + crit(i,:)
      end do
      mcrit = mcrit / dble(n_iter_range)
      if ((XDiff < xtol .and. FDiff < ftol) .or. &
          (nrmG < gtol) .or. &
          (mcrit(2) < 10.0d0*xtol .and. mcrit(3) < 10.0d0*ftol)) then
        deallocate(mcrit)
        exit
      end if
      deallocate(mcrit)

      y  = dtX - dtXP
      sy = abs(sum(s * y))
      tau = tau_orig
      if (sy > 0.0d0) then
        if (mod(itr,2) == 0) then
          tau = sum(s**2) / sy
        else
          tau = sy / sum(y**2)
        end if
        tau = max(min(tau,1.0d20),1.0d-20)
      end if

      Qp    = Q
      Q     = gamma * Qp + 1.0d0
      Cval  = (gamma*Qp*Cval + f) / Q
    end do

    deallocate(xtg, gg, xx, xxgg, dtX, crit, g, gp)
    x_out = x
  end subroutine optimize

  subroutine random_normal(mat)
    implicit none
    real(8), intent(out) :: mat(:,:)
    integer :: m, n, total, k, row, col
    real(8) :: u1, u2, z, pi, eps
    m = size(mat,1)
    n = size(mat,2)
    total = m*n
    pi = 4.0d0 * atan(1.0d0)
    eps = 1.0d-12
    k = 1
    do while (k <= total)
      call random_number(u1)
      call random_number(u2)
      if (u1 < eps) u1 = eps
      if (k < total) then
        z = sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
        row = mod(k-1,m)+1; col = ((k-1)/m)+1
        mat(row,col) = z; k = k+1
        z = sqrt(-2.0d0*log(u1)) * sin(2.0d0*pi*u2)
        row = mod(k-1,m)+1; col = ((k-1)/m)+1
        mat(row,col) = z; k = k+1
      else
        z = sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
        row = mod(k-1,m)+1; col = ((k-1)/m)+1
        mat(row,col) = z; k = k+1
      end if
    end do
  end subroutine random_normal

end module cutnorm_mod

program main
  use cutnorm_mod
  implicit none
  integer, parameter :: n = 400
  real(8), dimension(n,n) :: A, B
  integer :: i, j
  real(8) :: s, r, p1, p2

  call random_seed()

  p1 = 0.2d0
  p2 = 0.5d0
  A = 0.0d0
  do i = 1, n
    do j = i+1, n
      call random_number(r)
      if (r < p1) then
        A(i,j) = 1.0d0; A(j,i) = 1.0d0
      end if
    end do
  end do

  B = 0.0d0
  do i = 1, n
    do j = i+1, n
      call random_number(r)
      if (r < p2) then
        B(i,j) = 1.0d0; B(j,i) = 1.0d0
      end if
    end do
  end do

  s = cut_distance(A, B)
  print *, "The cut distance between A and B: ", s
end program main
