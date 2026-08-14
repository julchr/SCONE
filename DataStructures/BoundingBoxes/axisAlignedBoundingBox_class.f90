module axisAlignedBoundingBox_class

  use genericProcedures,  only : anyAreEqual, areEqual, swap
  use numPrecision
  use publicObjects,      only : intersectionTestPayload, intersectionTestResult, &
                                 newRationalIntersectionTestPayload, rationalIntersectionTestPayload
  use ratint_mod
  use universalVariables, only : INF, NUDGE, ONE, SURF_TOL, ZERO

  implicit none
  private

  !!
  !!
  !!
  type, public :: axisAlignedBoundingBox
    private
    real(defReal), dimension(3)     :: centre = ZERO, halfwidths = ZERO
    real(defReal), dimension(3, 2)  :: bounds = ZERO
    type(ratint_t), dimension(3, 2) :: rationalBounds
  contains
    ! Build procedure.
    procedure          :: init
    procedure          :: kill
    ! Runtime procedures.
    generic            :: computeBounds => computeBoundsFromCoords, computeBoundsFromBoundingBoxes
    procedure, private :: computeBoundsFromBoundingBoxes
    procedure, private :: computeBoundsFromCoords
    generic            :: contains => containsCoords, containsBoundingBox
    procedure, private :: containsBoundingBox
    procedure, private :: containsCoords
    generic            :: distanceSquared => distanceSquared_Vertex
    procedure, private :: distanceSquared_Vertex
    procedure          :: getBounds
    procedure          :: getCentre
    procedure          :: getHalfwidths
    generic            :: intersects => intersects_BoundingBox, intersects_Ray, intersects_Ray_rational
    procedure, private :: intersects_BoundingBox
    procedure, private :: intersects_Ray
    procedure, private :: intersects_Ray_rational
    procedure          :: pushFromBoundary
  end type axisAlignedBoundingBox

contains
  !! Subroutine 'computeBoundingBox'
  !!
  !! Basic description:
  !!   Computes the bounding box of the node along a specified dimension.
  !!
  !! Detailed description:
  !!   First sorts coordinates in the node in increasing order along the specified dimension. The
  !!   bounding box of the node along this dimension is then simply given by the minimum and maximum
  !!   coordinates along the dimension.
  !!
  !! Arguments:
  !!   nVertices [in]   -> Number of vertices in the node.
  !!   boundingBox [in] -> Bounding boxes of each coordinate in the node.
  !!
  pure subroutine computeBoundsFromBoundingBoxes(self, boundingBoxes)
    class(axisAlignedBoundingBox), intent(inout)           :: self
    type(axisAlignedBoundingBox), dimension(:), intent(in) :: boundingBoxes
    integer(shortInt)                                      :: i, nBoundingBoxes
    real(defReal), dimension(3)                            :: minCoords, maxCoords
    real(defReal), dimension(3, 2)                         :: bounds

    ! Check against zero-sized arrays.
    nBoundingBoxes = size(boundingBoxes)
    if (nBoundingBoxes == 0) return

    bounds = boundingBoxes(1) % getBounds()
    minCoords = bounds(:, 1)
    maxCoords = bounds(:, 2)

    do i = 2, nBoundingBoxes
      bounds = boundingBoxes(i) % getBounds()
      minCoords = min(minCoords, bounds(:, 1))
      maxCoords = max(maxCoords, bounds(:, 2))

    end do
    call self % init(reshape([minCoords, maxCoords], shape=[3, 2]))

  end subroutine computeBoundsFromBoundingBoxes

  !!
  !!
  !!
  pure subroutine computeBoundsFromCoords(self, primitives)
    class(axisAlignedBoundingBox), intent(inout) :: self
    real(defReal), dimension(:, :), intent(in)   :: primitives

    ! Check against zero-sized arrays.
    if (size(primitives, 2) == 0) return
    call self % init([minval(primitives, dim = 2), maxval(primitives, dim = 2)])

  end subroutine computeBoundsFromCoords

  !!
  !!
  !!
  elemental function containsBoundingBox(self, boundingBox) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    type(axisAlignedBoundingBox), intent(in)  :: boundingBox
    logical(defBool)                          :: doesIt

    doesIt = all(self % bounds(:, 1) <= boundingBox % bounds(:, 1)) .and. all(boundingBox % bounds(:, 2) <= self % bounds(:, 2))

  end function containsBoundingBox

  !!
  !!
  !!
  pure function containsCoords(self, r, u) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3), intent(in)   :: r, u
    integer(shortInt)                         :: i
    logical(defBool)                          :: doesIt

    do i = 1, 3
      if (areEqual(self % bounds(i, 1), r(i))) then
        doesIt = ZERO < u(i)

      elseif (areEqual(self % bounds(i, 2), r(i))) then
        doesIt = u(i) < ZERO

      else
        doesIt = self % bounds(i, 1) < r(i) .and. r(i) < self % bounds(i, 2)

      end if
      if (.not. doesIt) return

    end do

  end function containsCoords

  !!
  !!
  !!
  pure function distanceSquared_Vertex(self, r) result(dSquared)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3), intent(in)   :: r
    real(defReal)                             :: dSquared
    real(defReal), dimension(3)               :: diff
    integer(shortInt)                         :: i

    ! Initialise diff = ZERO then loop over all dimensions.
    diff = ZERO
    do i = 1, 3
      if (r(i) < self % bounds(i, 1)) then
        diff(i) = self % bounds(i, 1) - r(i)

      elseif (r(i) > self % bounds(i, 2)) then
        diff(i) = r(i) - self % bounds(i, 2)

      end if

    end do
    dSquared = dot_product(diff, diff)

  end function distanceSquared_Vertex

  !!
  !!
  !!
  pure function getBounds(self) result(bounds)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3, 2)            :: bounds

    bounds = self % bounds

  end function getBounds

  !!
  !!
  !!
  pure function getCentre(self) result(centre)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3)               :: centre

    centre = self % centre

  end function getCentre
  
  !!
  !!
  !!
  pure function getHalfwidths(self) result(halfwidths)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3)               :: halfwidths

    halfwidths = self % halfwidths

  end function getHalfwidths 

  !!
  !!
  !!
  pure subroutine init(self, bounds)
    class(axisAlignedBoundingBox), intent(inout) :: self
    real(defReal), dimension(3, 2), intent(in)   :: bounds

    self % bounds = bounds
    self % centre = HALF * (bounds(:, 1) + bounds(:, 2))
    self % halfwidths = HALF * (bounds(:, 2) - bounds(:, 1))

    ! Pre-compute rational bounds.
    self % rationalBounds(:, 1) = convert_ieee(self % bounds(:, 1))
    self % rationalBounds(:, 2) = convert_ieee(self % bounds(:, 2))

  end subroutine init

  !!
  !!
  !!
  elemental function intersects_BoundingBox(self, boundingBox) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    type(axisAlignedBoundingBox), intent(in)  :: boundingBox
    logical(defBool)                          :: doesIt

    doesIt = all(self % bounds(:, 1) <= boundingBox % bounds(:, 2)) .and. all(boundingBox % bounds(:, 1) <= self % bounds(:, 2))

  end function intersects_BoundingBox

  !!
  !!
  !!
  pure function intersects_Ray(self, payload) result(res)
    class(axisAlignedBoundingBox), intent(in)  :: self
    class(intersectionTestPayload), intent(in) :: payload
    integer(shortInt)                          :: i
    real(defReal)                              :: inverseU, tFar, tNear, t1, t2
    type(intersectionTestResult)               :: res

    ! Initialise tNear = -INF, tFar = INF, then loop over all halfwidths.
    tNear = -INF
    tFar = INF

    do i = 1, 3
      if(areEqual(payload % u(i), ZERO)) then
        if(areEqual(payload % r(i), self % bounds(i, 1)) .or. &
           areEqual(payload % r(i), self % bounds(i, 2))) then
          ! If ray is parallel to current dimension and is within tolerance of either plane of the
          ! bounding box along the current dimension, launch exact computation and return.
          res = self % intersects_ray_rational(newRationalIntersectionTestPayload(convert_ieee(payload % r), &
                                                                                  convert_ieee(payload % u), &
                                                                                  convert_ieee(payload % dMax)))
          return

        end if
        ! If ray is parallel and is definitely outside the box, there can be no intersection so return early.
        ! Else, cycle to the next dimension.
        if(payload % r(i) < self % bounds(i, 1) .or. self % bounds(i, 2) < payload % r(i)) return
        cycle

      else
        inverseU = ONE / payload % u(i)
        t1 = (self % bounds(i, 1) - payload % r(i)) * inverseU
        t2 = (self % bounds(i, 2) - payload % r(i)) * inverseU
        if(t2 < t1) call swap(t1, t2)

        tNear = max(tNear, t1)
        tFar = min(tFar, t2)

        ! Return early if intersection is impossible (far intersection is definitely greater than near intersection, or far
        ! intersection is definitely negative).
        if((tFar < tNear .and. .not. areEqual(tNear, tFar)) .or. &
           (tFar < ZERO .and. .not. areEqual(tFar, ZERO))) return

      end if

    end do

    ! If results are ambiguous, launch an exact computation.
    if(areEqual(tNear, ZERO) .or. areEqual(tFar, ZERO) .or. areEqual(tNear, tFar)) then
      res = self % intersects_ray_rational(newRationalIntersectionTestPayload(convert_ieee(payload % r), &
                                                                              convert_ieee(payload % u), &
                                                                              convert_ieee(payload % dMax)))
      ! TODO: change this; needs exact computed distance. Do this when reworking acceleration structures.
      res % d = merge(tFar, tNear, tNear < ZERO)

    else
      res % d = merge(tFar, tNear, tNear < ZERO)
      res % intersects = .true.

    end if

  end function intersects_Ray

  !!
  !!
  !!
  pure function intersects_ray_rational(self, payload) result(res)
    class(axisAlignedBoundingBox), intent(in)         :: self
    type(rationalIntersectionTestPayload), intent(in) :: payload
    integer(shortInt)                                 :: i
    logical(defBool)                                  :: areDistancesInvalid
    type(intersectionTestResult)                      :: res
    type(ratint_t)                                    :: temp, tFar, tNear, t1, t2, ZERO_rational

    ! Pre-compute ZERO_rational.
    ZERO_rational = convert_int(0_longInt)

    ! Initialise areDistancesInvalid = .true. then loop over all dimensions.
    areDistancesInvalid = .true.
    do i = 1, 3
      ! Check if ray is parallel to the current dimension.
      if(isZero(payload % u(i))) then
        ! If ray starts outside the box simply return, else cycle to the next dimension.
        if(self % rationalBounds(i, 1) > payload % r(i) .or. &
           payload % r(i) > self % rationalBounds(i, 2)) return
        cycle

      end if

      ! Compute t1 and t2.
      t1 = (self % rationalBounds(i, 1) - payload % r(i)) / payload % u(i)
      t2 = (self % rationalBounds(i, 2) - payload % r(i)) / payload % u(i)

      ! Swap if necessary.
      if(t1 > t2) then
        temp = t1
        t1 = t2
        t2 = temp

      end if

      ! Update values.
      if(areDistancesInvalid) then
        tNear = t1
        tFar = t2

        ! Update flag.
        areDistancesInvalid = .false.

      else
        if(t1 > tNear) tNear = t1
        if(tFar > t2) tFar = t2

      end if

      ! Return early if crossing is impossible.
      if(tNear > tFar .or. ZERO_rational > tFar) return

    end do

    ! If reached, the ray intersects the bounding box.
    res % intersects = .true.

  end function intersects_ray_rational

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(axisAlignedBoundingBox), intent(inout) :: self

    self % bounds = ZERO
    self % centre = ZERO
    self % halfwidths = ZERO

  end subroutine kill

  !!
  !!
  !! Note:
  !!   Assumes that the coordinates are already inside the bounding box to begin with.
  !!
  subroutine pushFromBoundary(self, u, r, inside)
    class(axisAlignedBoundingBox), intent(in)  :: self
    real(defReal), dimension(3), intent(in)    :: u
    real(defReal), dimension(3), intent(inout) :: r
    logical(defBool), intent(out)              :: inside
    real(defReal), dimension(3)                :: nudgeDirection
    logical(defBool)                           :: isOnBoundary
    integer(shortInt)                          :: i

    ! Initialise hasEscaped = .false. and check if coordinates are on the boundary.
    inside = .true.
    isOnBoundary = anyAreEqual(self % bounds(:, 1), r) .or. anyAreEqual(self % bounds(:, 2), r)

    ! If the coordinates are not on the boundary simply return.
    if (.not. isOnBoundary) return

    ! Nudge the coordinates until they are not on any boundaries anymore.
    do while (isOnBoundary)
      nudgeDirection = ZERO
      do i = 1, 3
        if (areEqual(self % bounds(i, 1), r(i))) then
          if (areEqual(u(i), ZERO)) nudgeDirection(i) = ONE

        elseif (areEqual(self % bounds(i, 2), r(i))) then
          if (areEqual(u(i), ZERO)) nudgeDirection(i) = -ONE

        end if

      end do

      if (any(nudgeDirection /= ZERO)) then
        nudgeDirection = nudgeDirection / norm2(nudgeDirection)

      else
        nudgeDirection = u

      end if
      r = r + nudgeDirection * NUDGE
      isOnBoundary = anyAreEqual(self % bounds(:, 1), r) .or. anyAreEqual(self % bounds(:, 2), r)

    end do

    inside = self % contains(r, u)

  end subroutine pushFromBoundary

end module axisAlignedBoundingBox_class