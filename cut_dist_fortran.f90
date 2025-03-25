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
    diff = (G - H) / (dble(n*n))
    sdp = cut_norm(diff)
    deallocate(diff)
  end function cut_distance

  function cut_norm(A) result(sdp_result)
    implicit none
    real(8), intent(in) :: A(:,:)
    real(8) :: sdp_result
    integer :: n1, N, n2, p, i, j, temp_int
    real(8) :: temp
    real(8), allocatable :: A_col_sum(:), A_row_sum(:)
    real(8), allocatable :: A_ext(:,:), A_final(:,:)
    real(8), allocatable :: x0(:,:), x(:,:)
    real(8), allocatable :: nrmx0(:)
    real(8), allocatable :: U(:,:), V(:,:), M(:,:)
    n1 = size(A,1)
    allocate(A_col_sum(n1))
    allocate(A_row_sum(n1))
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
    allocate(A_final(N, N))
    A_final(1:n1, :) = A_ext
    do j = 1, n1
       A_final(N, j) = -A_col_sum(j)
    end do
    A_final(N, N) = temp
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
    allocate(x0(p, n2))
    call random_normal(x0)
    allocate(nrmx0(n2))
    do j = 1, n2
       nrmx0(j) = sqrt(sum(x0(:,j)**2))
       if (nrmx0(j) > 1.0d-12) then
          x0(:,j) = x0(:,j) / nrmx0(j)
       end if
    end do
    deallocate(nrmx0)
    call optimize(x0, A_final, x)
    allocate(U(p, N))
    allocate(V(p, N))
    U = x(:, 1:N)
    V = x(:, N+1:n2)
    allocate(M(N, N))
    M = matmul(transpose(U), V)
    sdp_result = abs(sum(A_final * M)) / 4.0d0
    deallocate(A_col_sum, A_row_sum, A_ext, A_final, x0, x, U, V, M)
  end function cut_norm

  subroutine cut_norm_quad(V, A, f, g)
    implicit none
    real(8), intent(in) :: A(:,:)
    real(8), intent(in) :: V(:,:)
    real(8), intent(out) :: f
    real(8), allocatable, intent(out) :: g(:,:)
    integer :: p, n, i, j
    real(8), allocatable :: Us(:,:), Vs(:,:)
    n = size(A,1)
    p = size(V,1)
    allocate(Us(p, n))
    allocate(Vs(p, n))
    Us = V(:, 1:n)
    Vs = V(:, n+1:2*n)
    allocate(g(p, 2*n))
    g(:,1:n)   = 2.0d0 * matmul(Vs, transpose(A))
    g(:,n+1:2*n) = 2.0d0 * matmul(Us, A)
    f = ( sum(g(:,1:n) * Us) + sum(g(:,n+1:2*n) * Vs) ) / 2.0d0
    deallocate(Us, Vs)
  end subroutine cut_norm_quad

  subroutine optimize(x, A, x_out)
    implicit none
    real(8), intent(inout) :: x(:,:)
    real(8), intent(in) :: A(:,:)
    real(8), allocatable, intent(out) :: x_out(:,:)
    integer :: p, n2, i, j, itr, mxitr, nt, n_iter_range, nls
    real(8) :: xtol, ftol, gtol, rho, eta, gamma, tau, tau_orig
    real(8), allocatable :: nrmx(:)
    real(8) :: f, fp, deriv, XDiff, FDiff, sy, temp_tau, Q, Qp, Cval
    real(8), allocatable :: g(:,:), gp(:,:), dtX(:,:), dtXP(:,:), s(:,:)
    real(8), allocatable :: xtg(:), gg(:), xx(:), xxgg(:)
    real(8), allocatable :: a1(:), a2(:)
    real(8) :: nrmG
    real(8), allocatable :: xp(:,:), y(:,:)
    real(8), allocatable :: crit(:,:), mcrit(:)
    p = size(x,1)
    n2 = size(x,2)
    xtol = 1.0d-8
    ftol = 1.0d-10
    gtol = 1.0d-8
    rho  = 1.0d-4
    eta  = 0.1d0
    gamma = 0.85d0
    tau = 1.0d-3
    tau_orig = tau
    nt = 5
    mxitr = 600
    allocate(nrmx(n2))
    do j = 1, n2
       nrmx(j) = sqrt(sum(x(:,j)**2))
       if (nrmx(j) > 1.0d-12) then
          x(:,j) = x(:,j) / nrmx(j)
       end if
    end do
    deallocate(nrmx)
    call cut_norm_quad(x, A, f, g)
    allocate(xtg(n2), gg(n2), xx(n2), xxgg(n2))
    do j = 1, n2
       xtg(j) = sum(x(:,j) * g(:,j))
       gg(j)   = sum(g(:,j)**2)
       xx(j)   = sum(x(:,j)**2)
       xxgg(j) = xx(j) * gg(j)
    end do
    allocate(dtX(p, n2))
    do j = 1, n2
       dtX(:,j) = x(:,j) * xtg(j) - g(:,j)
    end do
    nrmG = sqrt(sum(dtX**2))
    Q = 1.0d0
    Cval = f
    allocate(crit(mxitr,3))
    crit = 0.0d0
    do itr = 1, mxitr
       dtXP = dtX
       allocate(xp(p, n2))
       xp = x
       fp = f
       allocate(gp(size(g,1), n2))
       gp = g
       nls = 1
       deriv = rho * nrmG**2
       do
          temp_tau = tau / 2.0d0
          allocate(a1(n2), a2(n2))
          do j = 1, n2
             a1(j) = ((1.0d0 + temp_tau*xtg(j))**2 - temp_tau**2 * xxgg(j)) / (1.0d0 + temp_tau**2 * (-xtg(j)**2 + xxgg(j)))
             a2(j) = -tau * xx(j) / (1.0d0 + temp_tau**2 * (-xtg(j)**2 + xxgg(j)))
          end do
          do j = 1, n2
             x(:,j) = xp(:,j)*a1(j) + gp(:,j)*a2(j)
          end do
          deallocate(a1, a2)
          call cut_norm_quad(x, A, f, g)
          if ( f <= Cval - tau * deriv .or. nls >= 5 ) exit
          tau = eta * tau
          nls = nls + 1
       end do
       do j = 1, n2
          xtg(j) = sum(x(:,j) * g(:,j))
          gg(j)   = sum(g(:,j)**2)
          xx(j)   = sum(x(:,j)**2)
          xxgg(j) = xx(j) * gg(j)
       end do
       do j = 1, n2
          dtX(:,j) = x(:,j)*xtg(j) - g(:,j)
       end do
       nrmG = sqrt(sum(dtX**2))
       allocate(s(p, n2))
       s = x - xp
       XDiff = sqrt(sum(s**2)) / sqrt(dble(p))
       FDiff = abs(fp - f) / (abs(fp) + 1.0d0)
       crit(itr,1) = nrmG
       crit(itr,2) = XDiff
       crit(itr,3) = FDiff
       n_iter_range = min(nt, itr)
       allocate(mcrit(3))
       mcrit = 0.0d0
       do j = itr - n_iter_range + 1, itr
          mcrit = mcrit + crit(j, :)
       end do
       mcrit = mcrit / dble(n_iter_range)
       if ((XDiff < xtol .and. FDiff < ftol) .or. (nrmG < gtol) .or. (mcrit(2) < 10.0d0*xtol .and. mcrit(3) < 10.0d0*ftol)) then
          deallocate(mcrit)
          deallocate(xp, gp, s)
          exit
       end if
       deallocate(mcrit)
       deallocate(xp, gp, s)
       allocate(y(p, n2))
       y = dtX - dtXP
       sy = abs(sum(x*y))
       tau = tau_orig
       if (sy > 0.0d0) then
          if (mod(itr,2) == 0) then
             tau = sum(x**2) / sy
          else
             tau = sy / sum(y**2)
          end if
          if (tau > 1.0d20) tau = 1.0d20
          if (tau < 1.0d-20) tau = 1.0d-20
       end if
       Qp = Q
       Q = gamma * Qp + 1.0d0
       Cval = (gamma * Qp * Cval + f) / Q
       deallocate(y)
    end do
    x_out = x
    deallocate(xtg, gg, xx, xxgg, dtX, crit, g)
  end subroutine optimize

  subroutine random_normal(mat)
    implicit none
    real(8), intent(out) :: mat(:,:)
    integer :: m, n, total, k, row, col
    real(8) :: u1, u2, z, pi, eps
    m = size(mat,1)
    n = size(mat,2)
    total = m * n
    pi = 4.0d0 * atan(1.0d0)
    eps = 1.0d-12
    k = 1
    do while (k <= total)
       call random_number(u1)
       call random_number(u2)
       if (u1 < eps) then
          u1 = eps
       end if
       if (k < total) then
          z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
          row = mod(k-1, m) + 1
          col = ((k-1) / m) + 1
          mat(row, col) = z
          k = k + 1
          z = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * pi * u2)
          row = mod(k-1, m) + 1
          col = ((k-1) / m) + 1
          mat(row, col) = z
          k = k + 1
       else
          z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
          row = mod(k-1, m) + 1
          col = ((k-1) / m) + 1
          mat(row, col) = z
          k = k + 1
       end if
    end do
  end subroutine random_normal

end module cutnorm_mod

program main
  use cutnorm_mod
  implicit none
  integer, parameter :: n = 30
  real(8), dimension(n,n) :: A, B
  real(8) :: s, r
  real(8) :: p1, p2
  integer :: i, j
  call random_seed()
  p1 = 0.2d0
  p2 = 0.5d0
  A = 0.0d0
  do i = 1, n
     do j = i+1, n
        call random_number(r)
        if (r <  p1) then
           A(i,j) = 1.0d0
           A(j,i) = 1.0d0
        end if
     end do
  end do
  B = 0.0d0
  do i = 1, n
     do j = i+1, n
        call random_number(r)
        if (r < p2) then
           B(i,j) = 1.0d0
           B(j,i) = 1.0d0
        end if
     end do
  end do
  s = cut_distance(A, B)
  print *, "The cut distance between A and B: ", s
end program main
