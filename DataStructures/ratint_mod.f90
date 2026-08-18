module ratint_mod
    
  use limb_class
  use numPrecision
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic

  implicit none 

  type ratint_t
    type(limb_t) :: p, q
  end type ratint_t

  interface signed
    module procedure signed_longInt
    module procedure signed_rational
    module procedure signed_shortInt
  end interface signed

  interface convert_ieee 
    module procedure convert_ieee64
    module procedure convert_ieee64Vector
  end interface convert_ieee

  interface isZero 
    module procedure isZero_flat 
    module procedure isZero_Vector
  end interface isZero

  interface operator (+)
    module procedure addpure 
    module procedure addmixedL
    module procedure addmixedR
  end interface operator (+)

  interface operator (-)
    module procedure subtractpure
    module procedure subtractmixedL
    module procedure subtractmixedR
  end interface operator (-)

  interface operator (*)
    module procedure multiplypure
    module procedure multiplymixedL
    module procedure multiplymixedR
    module procedure multiplyVectorL
    module procedure multiplyVectorR
  end interface operator (*)

  interface operator (/)
    module procedure dividepure
    module procedure dividemixedL
    module procedure dividemixedR
  end interface operator (/)

  interface assignment (=)
    module procedure assignpure
  end interface assignment (=)

  interface operator (==)
    module procedure equality
  end interface operator (==)

  interface operator (>)
    module procedure gneq
  end interface operator (>)

  interface operator (>=)
    module procedure geq
  end interface operator (>=)

  interface dot_product 
    module procedure dot_product_ratint
  end interface dot_product

  interface swapSign 
    module procedure swapSign_flat 
    module procedure swapSign_vector 
  end interface swapSign
  
  contains
  !! Function 'absoluteValue'
  !!
  !! Description:
  !!   Returns the absolute value of a rational number.
  !!
  !! Arguments:
  !!   r [in] -> Rational number.
  !!
  !! Result:
  !!   res -> Absolute value of r.
  !!
  elemental function absoluteValue(r) result(res)
    type(ratint_t), intent(in) :: r
    type(ratint_t)             :: res

    ! Copy number then ensure signs are positive.
    res = r
    res % p % sign = 1
    res % q % sign = 1

  end function absoluteValue

  !! Function 'signed_longInt'
  !!
  !! Description:
  !!   Returns |n| with the sign of r.
  !!
  !! Arguments:
  !!   n [in] -> Integer.
  !!   r [in] -> Rational number.
  !!
  !! Result:
  !!   signedInt -> |n| with sign of r.
  !!
  !! Note:
  !!   Assumes that r1 % q % sign = 1, does not check.
  !!
  elemental function signed_longInt(n, r) result(signedInt)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in)   :: r
    integer(longInt)             :: signedInt

    signedInt = abs(n) * int(r % p % sign, longInt)

  end function signed_longInt

  !! Function 'signed_rational'
  !!
  !! Description:
  !!   Returns |r1| with the sign of r2.
  !!
  !! Arguments:
  !!   r1 [in] -> Rational number.
  !!   r2 [in] -> Rational number.
  !!
  !! Result:
  !!   res -> r1 with sign of r2.
  !!
  !! Note:
  !!   Assumes that r1 % q % sign = 1 and r2 % q % sign = 1, does not check.
  !!
  elemental function signed_rational(r1, r2) result(res)
    type(ratint_t), intent(in) :: r1, r2
    type(ratint_t)             :: res

    res = absoluteValue(r1)
    res % p % sign = r2 % p % sign

  end function signed_rational

  !! Function 'signed_shortInt'
  !!
  !! Description:
  !!   Returns |n| with the sign of r.
  !!
  !! Arguments:
  !!   n [in] -> Integer.
  !!   r [in] -> Rational number.
  !!
  !! Result:
  !!   signedInt -> n with sign of r.
  !!
  !! Note:
  !!   Assumes that r % q % sign = 1, does not check.
  !!
  elemental function signed_shortInt(n, r) result(signedInt)
    integer(shortInt), intent(in) :: n
    type(ratint_t), intent(in)    :: r
    integer(shortInt)             :: signedInt

    signedInt = abs(n) * r % p % sign

  end function signed_shortInt

  !! Enforces the denominator to always be positive by moving any negative
  !! sign into the numerator. All comparison operators cross-multiply by q 
  !! and are only valid while this holds.
  pure subroutine normaliseSign(r)
    type(ratint_t), intent(inout) :: r

    if(checkInvalidRatint(r)) return
    if(r % q % sign == -1) then
      r % q % sign = 1
      if(.not. limbiszero(r % p)) r % p % sign = -r % p % sign

    end if

  end subroutine normaliseSign

  pure function def_ratint(n, d, s) result(r)
      integer(longInt), intent(in) :: n, d, s
      type(ratint_t) :: r 

      r%p = initlimb(n) * s
      r%q = initlimb(d)

      !r = simplify(r)

  end function def_ratint

  pure function def_ratint_large() result (r)
      type(ratint_t) :: r 

      r%p = initlimb() 
      r%p%limbs(1) = 0 
      r%p%limbs(2) = 0
      r%p%limbs(3) = 0
      r%p%limbs(4) = 1 
      r%p%front = 4 
      r%q = initlimb(1_8)

  end function def_ratint_large

  pure function initratint_vector() result(r)
      type(ratint_t), dimension(3) :: r 
      integer(shortInt) :: i 

      do i =1, 3 
          r(i) = convert_int(0_8)
      end do


  end function initratint_vector


  pure function convert_int(n) result(r)
      integer(longInt), intent(in) :: n 
      type(ratint_t) :: r 

      r%q = initlimb(1_8) 
      r%p = initlimb(n)

  end function convert_int

  pure subroutine swapSign_flat(a)
      type(ratint_t), intent(inout) :: a 

      a%p%sign = a%p%sign * (-1)
      
  end subroutine swapSign_flat


  pure subroutine swapSign_vector(a)
      type(ratint_t), dimension(:), intent(inout) :: a 
      integer(shortInt) :: i

      do i =1, size(a)
          call swapSign_flat(a(i))
      end do
      
  end subroutine swapSign_vector

  !! Function 'twoToThePower'
  !!
  !! Basic description:
  !!   Returns 2 ** n as a limb number, for 0 <= n. Never computes 2 ** n directly. The base of an integer literal is a default
  !!   integer, so the result overflows for 31 <= n whatever the kind of n. Here 2 ** n = 2 ** (31 * j + s) is built as 2 ** s 
  !!   shifted up by j whole limbs, with s < 31 so the remaining power is safe.
  !!
  pure function twoToThePower(n) result(l)
    integer(shortInt), intent(in) :: n
    type(limb_t)                  :: l

    l = shiftbyn(initlimb(2_longInt ** mod(n, 31)), n / 31)

  end function twoToThePower

  !! Function 'convert_ieee64'
  !!
  !! Basic description:
  !!   Converts a double precision number into an exact rational.
  !!
  !! Detailed description:
  !!   Every double is exactly mantissa * 2 ** exponentOfTwo with mantissa an integer of
  !!   at most 53 bits, so the conversion is lossless and needs no search over bits.
  !!
  pure function convert_ieee64(n) result(res)
    real(defReal), intent(in) :: n
    integer(longInt)          :: mantissa
    integer(shortInt)         :: exponentOfTwo
    type(ratint_t)            :: res

    ! Handle zero.
    if(n == ZERO) then
      call setZero(res)
      return

    end if

    ! Decompose into an integer mantissa and a binary exponent.
    mantissa = int(scale(fraction(abs(n)), 53), longInt)
    exponentOfTwo = exponent(n) - 53

    ! Strip trailing zero bits to keep the denominator as small as possible.
    do while(iand(mantissa, 1_longInt) == 0_longInt)
      mantissa = shiftr(mantissa, 1)
      exponentOfTwo = exponentOfTwo + 1

    end do

    ! Assemble the rational.
    if(exponentOfTwo < 0) then
      res % p = initlimb(mantissa)
      res % q = twoToThePower(-exponentOfTwo)

    else
      res % p = initlimb(mantissa) * twoToThePower(exponentOfTwo)
      res % q = initlimb(1_longInt)

    end if

    ! Apply the sign.
    if(n < ZERO) res % p % sign = -1

  end function convert_ieee64

  pure function convert_ieee64Vector(n) result(r)
      real(defReal), dimension(:), intent(in) :: n 
      type(ratint_t), allocatable:: r(:)
      integer(shortInt) :: i, height 

      height = size(n)
      allocate(r(height))

      do i=1, height 
          r(i) = convert_ieee64(n(i))
      end do


  end function convert_ieee64Vector

  pure function checkInvalidRatint(r) result(x)
      type(ratint_t), intent(in) :: r 
      logical :: x 

      x = .false.

      if (checkInvalid(r%p) .or. checkInvalid(r%q)) then 
          x = .true. 
          return
      end if 

  end function checkInvalidRatint


  pure subroutine setInvalidRatint(r)
      type(ratint_t), intent(inout) :: r 

      call setInvalid(r%p)

      call setINvalid(r%q)

  end subroutine setInvalidRatint

  elemental function isZero_flat(r) result(isIt)
    type(ratint_t), intent(in) :: r 
    logical(defBool)           :: isIt

    isIt = limbiszero(r % p)

  end function isZero_flat


  pure function isZero_Vector(r) result(x)
      type(ratint_t), dimension(:), intent(in) :: r 
      logical :: x
      integer(shortInt) :: i, height

      x = .true.

      height = size(r)

      do i=1, height
          if (.not. isZero_flat(r(i))) then 
              x = .false. 
              return 
          end if 
      end do 

  end function isZero_Vector



  pure function get_numerator(r) result(n)
      type(ratint_t), intent(in) :: r 
      type(limb_t) :: n 

      n = r%p 
  end function get_numerator

  pure function get_denominator(r) result(d)
      type(ratint_t), intent(in) :: r 
      type(limb_t) :: d

      d = r%q
  end function get_denominator

  pure function evaluate(r) result(v)
    type(ratint_t), intent(in) :: r 
    real(defReal) :: v 

    v = r%p / r%q

  end function evaluate

  elemental function addpure(r1, r2) result(res)
    type(ratint_t), intent(in) :: r1, r2
    integer(shortInt)          :: n
    type(ratint_t)             :: res

    if(checkInvalidRatint(r1) .or. checkInvalidRatint(r2)) then 
      call setInvalidRatint(res)
      return

    end if

    ! Check if both r1 and r2 have the same denominator.
    if(r1 % q == r2 % q) then
      res % p = r1 % p + r2 % p
      res % q = r1 % q

    else
      res % p = r1 % p * r2 % q + r2 % p * r1 % q
      res % q = r1 % q * r2 % q

      ! Strip common factors of two.
      n = min(trailingZeroBits(res % p), trailingZeroBits(res % q))
      if(0 < n) then
        call shiftRightBits(res % p, n)
        call shiftRightBits(res % q, n)

      end if

    end if

    ! Normalise signs.
    call normaliseSign(res)

  end function addpure

  ! Allows addition between : int + ratint
  elemental type(ratint_t) function addmixedL(n, r1)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in) :: r1
    type(ratint_t) :: rn 

    rn = convert_int(n)
    addmixedL = addpure(rn, r1)

  end function addmixedL


  ! Allows addition between : ratint + int
  elemental type(ratint_t) function addmixedR(r1, n)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in) :: r1
    type(ratint_t) :: rn 

    rn = convert_int(n)
    addmixedR = addpure(r1, rn)

  end function



  ! negates the second value and adds the results
  elemental type(ratint_t) function subtractpure(r1, r2)
    type(ratint_t), intent(in) :: r1, r2 
    type(ratint_t) :: r2t

    r2t%p = r2%p
    r2t%p%sign = r2%p%sign * (-1)
    r2t%q = r2%q

    subtractpure = addpure(r1, r2t)
  
  end function subtractpure


  ! Allows subtraction between : int - ratint
  elemental type(ratint_t) function subtractmixedL(n, r1)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in) :: r1
    type(ratint_t) :: rn 

    rn = convert_int(n)

    subtractmixedL = subtractpure(rn, r1)

  end function subtractmixedL

  
  ! Allows subtraction between : ratint - int
  elemental type(ratint_t) function subtractmixedR(r1, n)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn =convert_int(n)

      subtractmixedR = subtractpure(r1, rn)
  end function subtractmixedR

  
  !! multiplies numerator and denominator then simplifies the fraction
  elemental function multiplypure(r1, r2) result(res)
    type(ratint_t), intent(in) :: r1, r2
    integer(shortInt)          :: n
    type(ratint_t)             :: res

    if(checkInvalidRatint(r1) .or. checkInvalidRatint(r2)) then 
      call setInvalidRatint(res)
      return

    end if 

    res % p = r1 % p * r2 % p 
    res % q = r1 % q * r2 % q

    ! Strip common factors of two.
    n = min(trailingZeroBits(res % p), trailingZeroBits(res % q))
    if(0 < n) then
      call shiftRightBits(res % p, n)
      call shiftRightBits(res % q, n)

    end if

    ! Normalise signs.
    call normaliseSign(res)

  end function multiplypure



  ! Allows multiplication between : int * ratint
  elemental type(ratint_t) function multiplymixedL(n, r1)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn = convert_int(n)

      multiplymixedL = multiplypure(rn, r1)

  end function multiplymixedL

  

  ! Allows multiplication between : ratint * int
  elemental type(ratint_t) function multiplymixedR(r1, n)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn =convert_int(n)

      multiplymixedR = multiplypure(r1, rn)
  end function multiplymixedR


  pure function multiplyVectorL(n, r1) result(rn)
      integer(longInt), intent(in) :: n
      type(ratint_t), dimension(:), intent(in) :: r1
      type(ratint_t), dimension(size(r1)) :: rn 
      integer :: i 

      do i=1, size(r1)
          rn(i) = r1(i) * n 

      end do


  end function multiplyVectorL

  pure function multiplyVectorR(r1, n) result(rn)
      integer(longInt), intent(in) :: n
      type(ratint_t), dimension(:), intent(in) :: r1
      type(ratint_t), dimension(size(r1)) :: rn 
      integer :: i 

      do i=1, size(r1)
          rn(i) = r1(i) * n 

      end do


  end function multiplyVectorR


  

  ! Follows keep, change, flip rule, then applies multiplication
  ! NOTE: division by 0 causes NaN via modulo() call in gcd

  elemental type(ratint_t) function dividepure(r1, r2)
    type(ratint_t), intent(in) :: r1, r2
    type(ratint_t)             :: r3

    ! Handle division by 0.
    if(isZero(r2)) then
      call setInvalidRatint(dividepure)
      return

    end if

    r3 % p = r2 % q 
    r3 % q = r2 % p 

    dividepure = multiplypure(r1, r3)

    ! Normalise signs.
    call normaliseSign(dividepure)

  end function dividepure


  ! Allows division between : int / ratint
  elemental type(ratint_t) function dividemixedL(n, r1)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn = convert_int(n)

      dividemixedL = dividepure(rn, r1)

  end function dividemixedL

  
  ! Allows division between : ratint / int
  elemental type(ratint_t) function dividemixedR(r1, n)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn = convert_int(n)

      dividemixedR = dividepure(r1, rn)
  end function dividemixedR


  ! Copies over the values from rin (r input) into rout (r output)
  pure subroutine assignpure(rout, rin)
      type(ratint_t), intent(out) :: rout 
      type(ratint_t), intent(in) :: rin

      if (checkInvalidRatint(rin)) then 
          call setInvalidRatint(rout)
          return 
      end if

      rout%p = rin%p 
      rout%q = rin%q

  end subroutine assignpure

  pure function equality(a, b) result(r)
      type(ratint_t), intent(in) :: a,b 
      type(limb_t) :: prod1, prod2
      logical :: r 

      r = .false.

      if (checkInvalidRatint(a) .or. checkInvalidRatint(b)) then 
          return 
      end if

      if (a%q == b%q) then 
          if (a%p == b%p) then 
              r = .true. 
              return 
          end if 
      else 
          prod1 = a%p * b%q 
          prod2 = b%p * a%q 
          if (prod1 == prod2) then 
              r = .true.
              return 
          end if 
      end if

  
  end function equality


  pure function gneq(a, b) result(r)
      type(ratint_t), intent(in) :: a,b 
      type(limb_t) :: prod1, prod2
      logical :: r 

      r = .false.

      if (checkInvalidRatint(a) .or. checkInvalidRatint(b)) then 
          return 
      end if


      if (a%q == b%q) then 
          if (a%p > b%p) then 
              r = .true. 
              return 
          end if 
      else 
          prod1 = a%p * b%q 
          prod2 = b%p * a%q 
          if (prod1 > prod2) then 
              r = .true.
              return 
          end if 
      end if
  end function gneq 




  pure function geq(a,b) result(r)
      type(ratint_t), intent(in) :: a,b 
      type(limb_t) :: prod1, prod2
      logical :: r 

      r = .false.

      if (checkInvalidRatint(a) .or. checkInvalidRatint(b)) then 
          return 
      end if


      if (a%q == b%q) then 
          if (a%p >= b%p) then 
              r = .true. 
              return 
          end if 
      else 
          prod1 = a%p * b%q 
          prod2 = b%p * a%q 
          if (prod1 >= prod2) then 
              r = .true.
              return 
          end if 
      end if

  end function geq

  pure subroutine setZero(r)
      type(ratint_t), intent(inout) :: r 
      r%p = initlimb(0_8)
      r%q = initlimb(1_8)
  end subroutine setZero



  subroutine printRatInt(a)
      type(ratint_t), intent(in) :: a 
      print *, 'Numerator'
      call printlimb(a%p)
      print *, '/////'
      print *, 'Denominator'
      call printlimb(a%q)
  end subroutine printRatInt



  pure function dot_product_ratint(a,b) result(r)
      type(ratint_t), intent(in) :: a(:), b(:) 
      type(ratint_t) :: r
      integer(shortInt) :: height, i

      call setZero(r)

      if (size(a) /= size(b)) then 
          call setInvalidRatint(r)
          return 
      end if 

      height = size(a)

      do i = 1, height
          r = r + (a(i) * b(i))
      end do

  end function dot_product_ratint

end module ratint_mod



! program test 
!     use limb_class
!     use ratint
!     use, intrinsic :: iso_fortran_env
!     use, intrinsic :: ieee_arithmetic

!     implicit none 

!     real(defReal) :: v1, v2, v3, eval
!     type(ratint_t) :: ratint1, ratint2, vres
!     logical :: result


!     v1 = 1.0_real64 / 75.0_real64
!     print *, '----'
!     print *, v1

!     ratint1 = convert_ieee64(v1)
!     call printRatInt(ratint1)
!     eval = evaluate(ratint1)
!     print *, 'evaluated:'
!     print *, eval 

    
    
!     result = v1 == eval



! end program test