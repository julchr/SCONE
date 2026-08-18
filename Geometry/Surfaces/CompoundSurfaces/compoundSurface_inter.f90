module compoundSurface_inter

  use errors_mod,         only : fatalError
  use numPrecision
  use genericProcedures,  only : areEqual, numToChar
  use ratint_mod
  use surface_inter,      only : kill_super => kill, surface
  use universalVariables, only : INF, MISS_TOL, PERIODIC_BC, REFLECTIVE_BC, VACUUM_BC

  implicit none
  private

  !!
  !! Abstract interface to group all surfaces containing halfwidths. These surfaces are called
  !! compound surfaces. A compound surface contains one (eg, truncated quadratic surfaces), two
  !! (eg, square cylinders) or three (eg, boxes) halfwidths.
  !!
  !! Private members:
  !!   halfwidths          -> Halfwidth(s) of the surface. It is an allocatable array as different 
  !!                          compound surfaces can have different numbers of halfwidths.
  !!   halfwidthsOrigin    -> Origin or offset of the surface's halfwidth(s). Also allocatable.
  !!
  !! Interface:
  !!   closestDist         -> Returns the minimum and maximum distances to the intersected surface 
  !!                          halwidth.
  !!   isHalfspacePositive -> Returns .true. if a particle is in the positive halfspace with respect
  !!                          to the surface's halfwidth(s). Note that it only checks halfspace sign
  !!                          with respect to the components of a surface having halfwidths.
  !!   evaluateHalfwidths  -> Evaluates the expression F(r) = c with respect to the halfwidths of a
  !!                          surface.
  !!   killHalfwidths      -> Returns the halfwidths component to an unitialised state.
  !!   getHalfwidths       -> Returns the halfwidth(s) of the surface.
  !!   setHalfwidths       -> Sets the halfwidth(s) of the surface and their origin.
  !!
  type, public, abstract, extends(surface)       :: compoundSurface
    private
    integer(shortInt)                            :: nHalfwidths = 0
    integer(shortInt), dimension(:), allocatable :: BCs
    real(defReal), dimension(:), allocatable     :: halfwidths, halfwidthsOrigin
    real(defReal), dimension(:, :), allocatable  :: bounds
    type(ratint_t), dimension(:), allocatable    :: rationalHalfwidths, rationalHalfwidthsOrigin
    type(ratint_t), dimension(:, :), allocatable :: rationalBounds
  contains
    procedure          :: distancesCompound
    generic            :: isHalfspacePositive => isHalfspacePositive_defReal, isHalfSpacePositive_rational
    procedure, private :: isHalfspacePositive_defReal
    procedure, private :: isHalfspacePositive_rational
    procedure          :: evaluateCompound
    procedure          :: setCompoundBCs
    procedure          :: killCompound
    procedure          :: getBounds
    procedure          :: getHalfwidths
    procedure          :: getOrigins
    procedure          :: getRationalBounds
    procedure          :: setHalfwidths
    procedure          :: explicitCompoundBCs
    procedure          :: transformCompoundBCs
  end type compoundSurface

contains

  !! Function 'distancesCompound'
  !!
  !! Basic description:
  !!   Computes the distance to the nearest intersection between a ray of origin r and direction
  !!   u and a set of halfwidths.
  !!
  !! Detailed description:
  !!   Intersection test is based on the slab method. In essence, an axis-aligned bounding box can
  !!   be represented by two triples (l_x, l_y, l_z) and (h_x, h_y, h_z) which are
  !!   the low and high bounds of the box along each dimension. Each pair of parallel planes (eg,
  !!   x_min and x_max) defines a slab.
  !!
  !!   A point p(d) = [p_x(d), p_y(d), p_z(d)] along the ray can be represented in parametric form
  !!   as:
  !!
  !!   p(d) = r + d * u
  !!
  !!   Assuming that none of the direction components are ZERO, then we can solve for t as:
  !!
  !!   d = (p - r) / u
  !!
  !!   Therefore, the intersection of the ray with the two slabs orthogonal to the i-th dimension
  !!   are given by:
  !!
  !!   d_l_i = (l_i - r_i) / u
  !!   d_h_i = (h_i - r_i) / u
  !!
  !!   The minimum and maximum extrema of the segment within the i-th slab are then:
  !!
  !!   d_min_i = min(d_l_i, d_h_i)
  !!   d_max_i = max(d_l_i, d_h_i)
  !!
  !!   And the intersection of all these segments therefore is:
  !!
  !!   d_min = max(d_min_x, d_min_y, d_min_z)
  !!   d_max = min(d_max_x, d_max_y, d_max_z)
  !!
  !! Arguments:
  !!   offsetCoords [in]     -> Origin-centred coordinates components of the ray along the set of 
  !!                            halfwidths.
  !!   u [in]                -> Direction of the ray.
  !!   dMin [out]            -> Minimum distance of intersection to the halfwidth.
  !!   dMax [out]            -> Maximum distance of intersection to the halfwidth.
  !!   cannotIntersect [out] -> .true. if the ray cannot intersect the halfwidth(s).
  !!
  !! Notes:
  !!   Special handling is required when one or more component(s) of the ray's direction is (are) ZERO.
  !!   In this case, the ray is parallel to one or more pair(s) of orthogonal planes (eg, the ray only
  !!   moves in the y-direction).
  !!
  pure function distancesCompound(self, rComponents, uComponents, dMinIn, dMaxIn, surfTolCondition) result(d)
    class(compoundSurface), intent(in)                       :: self
    real(defReal), dimension(self % nHalfwidths), intent(in) :: rComponents, uComponents
    real(defReal), intent(in)                                :: dMinIn, dMaxIn
    logical(defBool), intent(in)                             :: surfTolCondition
    integer(shortInt)                                        :: i
    real(defReal)                                            :: d, dMin, dMax, halfwidth, offsetCoord, &
                                                                uComponent, inverseU

    ! Initialise d = INF then loop over all halfwidths.
    d = INF
    dMin = dMinIn
    dMax = dMaxIn
    do i = 1, self % nHalfwidths
      halfwidth = self % halfwidths(i)
      offsetCoord = rComponents(i) - self % halfwidthsOrigin(i)
      uComponent = uComponents(i)

      ! Perform early check to see if the particle is outside the slab and moving away from it along
      ! the current dimension. If yes the particle cannot intersect the slab and we can return early.
      if ((offsetCoord <= -halfwidth .and. uComponent <= ZERO) .or. &
      (offsetCoord >= halfwidth .and. uComponent >= ZERO)) return

      ! Cycle to next dimension if the current direction component is ZERO, since in this case there
      ! can be no intersection along the current dimension.
      if (areEqual(uComponent, ZERO)) cycle

      ! Pre-compute inverse direction component and update dMin and dMax.
      inverseU = 1 / uComponent
      dMin = max(dMin, (-sign(halfwidth, uComponent) - offsetCoord) * inverseU)
      dMax = min(dMax, (sign(halfwidth, uComponent) - offsetCoord) * inverseU)

    end do

    ! If particle intersects very close to a corner return. TODO: change this in a future update.
    if (dMax <= dMin * MISS_TOL) return

    ! If reached here, ray intersects the halfwidth(s). Update d = dMin. Check if particle is on 
    ! surface or already inside the surface and update d = dMax if yes.
    d = dMin
    if ((surfTolCondition .and. abs(dMax) >= abs(dMin)) .or. (.not. surfTolCondition .and. dMin <= ZERO)) d = dMax
    
    ! Cap distance to INF if d <= ZERO or d > INF.
    if (d <= ZERO .or. d > INF) d = INF

  end function distancesCompound

  !! Function 'isHalfspacePositive'
  !!
  !! Basic description:
  !!   Returns .true. if particle lies in the positive halfspace with respect to a set of halfwidths.
  !!
  !! Detailed description:
  !!   Performs a check along all dimensions along which the surface contains a halfwidth (we need to do
  !!   this to properly account for particles at corners). If the component of the direction vector along
  !!   a given dimension is ZERO, halfspace allocation is decided based on particle's location. If not, it
  !!   is determined by the sign of the projection of the particle's direction component on the halfwidth's
  !!   normal vector.
  !!
  !! Arguments:
  !!   offsetCoords [in] -> Origin-centred coordinate components of the particle along the dimensions 
  !!                        containing halfwidths.
  !!   uComponents [in]  -> Components of the particle's direction vector along the dimensions containing
  !!                        halfwidths.
  !!
  !! Result:
  !!   isIt              -> .true. if particle is in positive halfspace.
  !!
  pure function isHalfspacePositive_defReal(self, rComponents, uComponents) result(isIt)
    class(compoundSurface), intent(in)                       :: self
    real(defReal), dimension(self % nHalfwidths), intent(in) :: rComponents, uComponents
    logical(defBool)                                         :: isIt
    real(defReal)                                            :: surfTol, halfwidth, offsetCoord, &
                                                                dist, uComponent
    integer(shortInt)                                        :: i

    ! Initialise isIt = .true., retrieve surface tolerance then loop over all compound dimensions (note, 
    ! we need to do this to properly account for particles close to corners).
    isIt = .true.
    surfTol = self % getSurfTol()
    do i = 1, self % nHalfwidths
      ! Retrieve offset coordinates and halfwidth components for the current dimension. Compute
      ! distance to halfwidth along current dimension and cycle to the next dimension if
      ! abs(dist) >= surface tolerance.
      halfwidth = self % halfwidths(i)
      offsetCoord = rComponents(i) - self % halfwidthsOrigin(i)
      dist = abs(offsetCoord) - halfwidth
      if (abs(dist) >= surfTol) cycle
      
      ! Retrieve direction component for the current dimension.
      uComponent = uComponents(i)
      if (areEqual(uComponent, ZERO)) then
        ! If ray is approximately parallel to current dimension, launch exact computation and return.
        isIt = self % isHalfspacePositive(convert_ieee(rComponents), convert_ieee(uComponents))
        return

      end if

      ! If reached here, halfspace is given by the sign of the projection of uComponent onto
      ! halfwidth's normal. Return early if projection > ZERO.
      if (uComponent * sign(ONE, offsetCoord) > ZERO) return

    end do

    ! If reached here, the particle is in the negative halfspace. Update isIt = .false.
    isIt = .false.

  end function isHalfspacePositive_defReal

  !!
  !!
  !!
  pure function isHalfspacePositive_rational(self, rComponents, uComponents) result(isIt)
    class(compoundSurface), intent(in)                        :: self
    type(ratint_t), dimension(self % nHalfwidths), intent(in) :: rComponents, uComponents
    integer(shortInt)                                         :: i
    logical(defBool)                                          :: isIt
    logical(defBool), dimension(self % nHalfwidths)           :: isOnSurfaces
    type(ratint_t)                                            :: d, ZERO_rational
    type(ratint_t), dimension(self % nHalfwidths)             :: offsetCoords

    ! Pre-compute ZERO_rational.
    ZERO_rational = convert_int(0_longInt)

    ! First check that the point is actually on the surface.
    isIt = .true.
    isOnSurfaces = .false.
    do i = 1, self % nHalfwidths
      offsetCoords(i) = rComponents(i) - self % rationalHalfwidthsOrigin(i)
      d = absoluteValue(offsetCoords(i)) - self % rationalHalfwidths(i)

      ! If the distance is positive along any dimension, the point cannot be on the surface so return.
      if(d > ZERO_rational) return

      ! Check if point is on surface along the current dimension.
      if(isZero(d)) isOnSurfaces(i) = .true.

    end do

    ! If the point is not on any surfaces here it is in the negative halfspace so return.
    if(.not. any(isOnSurfaces)) then
      isIt = .false.
      return

    end if

    ! Loop over all compound dimensions (note: we need to do this to properly account for particles close to corners).
    do i = 1, self % nHalfwidths
      ! If point is not on surface along the current dimension cycle.
      if(.not. isOnSurfaces(i)) cycle
      
      ! If ray is parallel to current dimension the point is considered in the positive halfspace so return.
      if(isZero(uComponents(i))) return

      ! If reached here, halfspace is given by the sign of the projection of uComponent onto
      ! halfwidth's normal. Return early if projection > ZERO_rational.
      if(uComponents(i) * signed(1_longInt, offsetCoords(i)) > ZERO_rational) return

    end do

    ! If reached here, the particle is in the negative halfspace. Update isIt = .false.
    isIt = .false.

  end function isHalfspacePositive_rational

  !! Function 'evaluateCompound'
  !!
  !! Basic description:
  !!   Evaluates the expression F(r) = c for all dimensions containing halfwidths.
  !!
  !! Arguments:
  !!   offsetCoords -> Origin-centred coordinate components of the particle along the dimensions 
  !!                   containing halfwidths.
  !! Result:
  !!   c            -> Evaluated expression F(r).
  !!
  pure function evaluateCompound(self, rComponents) result(c)
    class(compoundSurface), intent(in)                       :: self
    real(defReal), dimension(self % nHalfwidths), intent(in) :: rComponents
    real(defReal)                                            :: c

    c = maxval(abs(rComponents - self % halfwidthsOrigin) - self % halfwidths)

  end function evaluateCompound

  !! Subroutine 'setCompoundBCs'
  !!
  !! Basic description:
  !!   Sets boundary conditions of the compound surface according to BCs. If BCs is absent, initialises
  !!   the boundary conditions to VACUUM_BC.
  !!
  !! Arguments:
  !!   nRequired [in]      -> Number of required boundary conditions for the compound surface.
  !!   BCs [in] [optional] -> Integer array containing the boundary condition flags. If absent boundary
  !!                          conditions are set to VACUUM_BC by default.
  !!
  !! Errors:
  !!   - fatalError if number of required boundary conditions does not match number of supplied boundary
  !!     conditions.
  !!   - fatalError if any boundary condition flags are unrecognised.
  !!   - fatalError if periodic boundary conditions are not applied on opposite surfaces.
  !!
  subroutine setCompoundBCs(self, nRequired, BCs)
    class(compoundSurface), intent(inout)                 :: self
    integer(shortInt), intent(in)                         :: nRequired
    integer(shortInt), dimension(:), intent(in), optional :: BCs
    integer(shortInt)                                     :: nBCs, i
    character(*), parameter                               :: here = 'setCompoundBCs (compoundSurface_inter.f90)'

    if (.not. present(BCs)) then
      allocate(self % BCs(nRequired))
      self % BCs = VACUUM_BC
      return

    end if

    nBCs = size(BCs)
    if (nBCs < nRequired) call fatalError(Here, &
    'Boundary conditions must have size '//numToChar(nRequired)//'. Has: '//numToChar(nBCs)//'.')

    ! Verify BC flags.
    do i = 1, nBCs
      select case(BCs(i))
        case (VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC)
          ! Do nothing, pass
        case default
          call fatalError(here, 'Unrecognised BC: '//numToChar(BCs(i))//' in position: '//numToChar(i)//'.')

      end select

    end do

    ! Verify periodic BCs.
    if (.not. all((BCs([(i, i = 1, nBCs, 2)]) == PERIODIC_BC) .eqv. (BCs([(i, i = 2, nBCs, 2)]) == PERIODIC_BC))) &
    call fatalError(here, 'Periodic BCs need to be applied on opposite surfaces.')

    ! Load BC flags.
    self % BCs = BCs(1:nBCs)

  end subroutine setCompoundBCs

  !! Subroutine 'killCompound'
  !!
  !! Basic description:
  !!   Returns compound surface to an unitialised state.
  !!
  elemental subroutine killCompound(self)
    class(compoundSurface), intent(inout) :: self

    self % nHalfwidths = 0
    if(allocated(self % BCs)) deallocate(self % BCs)
    if(allocated(self % halfwidths)) deallocate(self % halfwidths)
    if(allocated(self % halfwidthsOrigin)) deallocate(self % halfwidthsOrigin)
    if(allocated(self % bounds)) deallocate(self % bounds)
    if(allocated(self % rationalHalfwidths)) deallocate(self % rationalHalfwidths)
    if(allocated(self % rationalHalfwidthsOrigin)) deallocate(self % rationalHalfwidthsOrigin)
    if(allocated(self % rationalBounds)) deallocate(self % rationalBounds)

  end subroutine killCompound

  !!
  !!
  !!
  pure function getBounds(self) result(bounds)
    class(compoundSurface), intent(in)          :: self
    real(defReal), dimension(:, :), allocatable :: bounds

    if(allocated(self % bounds)) then
      bounds = self % bounds

    else
      allocate(bounds(0, 0))

    end if

  end function getBounds

  !! Function 'getHalfwidths'
  !!
  !! Basic description:
  !!   Returns the halfwidths of the surface.
  !!
  !! Result:
  !!   halfwidths -> Halfwidths of the surface.
  !!
  pure function getHalfwidths(self) result(halfwidths)
    class(compoundSurface), intent(in)       :: self
    real(defReal), dimension(:), allocatable :: halfwidths

    halfwidths = self % halfwidths

  end function getHalfwidths

  !!
  !!
  !!
  pure function getOrigins(self) result(origins)
    class(compoundSurface), intent(in)       :: self
    real(defReal), dimension(:), allocatable :: origins

    origins = self % halfwidthsOrigin

  end function getOrigins

  !!
  !!
  !!
  pure function getRationalBounds(self) result(rationalBounds)
    class(compoundSurface), intent(in)           :: self
    type(ratint_t), dimension(:, :), allocatable :: rationalBounds

    if(allocated(self % rationalBounds)) then
      rationalBounds = self % rationalBounds

    else
      allocate(rationalBounds(0, 0))

    end if

  end function getRationalBounds

  !! Subroutine 'setHalfwidths'
  !!
  !! Basic description:
  !!   Sets the halfwidths of the surface.
  !!
  !! Arguments:
  !!   halfwidths [in] -> Array of halfwidths.
  !!   nRequired [in]  -> Number of required halfwidths for the surface.
  !!
  !! Errors:
  !!   - fatalError if the number of supplied halfwidths does not match the number of required
  !!     halfwidths.
  !!   - fatalError if any of the supplied halfwidths are negative.
  !!
  subroutine setHalfwidths(self, halfwidths, origin, nRequired)
    class(compoundSurface), intent(inout)   :: self
    real(defReal), dimension(:), intent(in) :: halfwidths, origin
    integer(shortInt), intent(in)           :: nRequired
    integer(shortInt)                       :: nHalfwidths, nOrigin
    character(*), parameter                 :: here = 'setHalfwidths (compoundSurface_inter.f90)'

    nHalfwidths = size(halfwidths)
    if (nHalfwidths /= nRequired) call fatalError(here, &
    'Halfwidths must have size: '//numToChar(nRequired)//'. Has: '//numToChar(nHalfwidths)//'.')
    if (any(halfwidths <= ZERO)) call fatalError(here, 'Halfwidths cannot contain negative entries.')

    nOrigin = size(origin)
    if (nOrigin /= nRequired) call fatalError(here, &
    'Halfwidths origin must have size: '//numToChar(nRequired)//'. Has: '//numToChar(nOrigin)//'.')

    self % halfwidths = halfwidths
    self % halfwidthsOrigin = origin
    self % nHalfwidths = nHalfwidths

    ! Set rational components.
    self % rationalHalfwidths = convert_ieee(halfwidths)
    self % rationalHalfwidthsOrigin = convert_ieee(origin)

    ! Allocate bounds and compute them.
    allocate(self % bounds(nRequired, 2), self % rationalBounds(nRequired, 2))
    self % bounds(:, 1) = self % halfwidthsOrigin - self % halfwidths
    self % bounds(:, 2) = self % halfwidthsOrigin + self % halfwidths

    self % rationalBounds(:, 1) = self % rationalHalfwidthsOrigin - self % rationalHalfwidths
    self % rationalBounds(:, 2) = self % rationalHalfwidthsOrigin + self % rationalHalfwidths

  end subroutine setHalfwidths

  !! Subroutine 'explicitCompoundBCs'
  !!
  !! Basic description:
  !!   Explicitly applies boundary conditions on particle position and direction for each
  !!   halfwidth of the compound surface along specified dimensions.
  !!
  !! Arguments:
  !!   dimensions [in] -> Dimensions along which to apply boundary conditions.
  !!   r [inout]       -> Particle position.
  !!   u [inout]       -> Particle direction.
  !!
  pure subroutine explicitCompoundBCs(self, dimensions, r, u)
    class(compoundSurface), intent(in)                           :: self
    integer(shortInt), dimension(self % nHalfwidths), intent(in) :: dimensions
    real(defReal), dimension(3), intent(inout)                   :: r, u
    integer(shortInt)                                            :: i, dim, bc
    real(defReal)                                                :: surfTol, offsetCoord, halfwidth

    ! Retrieve surface tolerance then loop over all surface halfwidths.
    surfTol = self % getSurfTol()
    do i = 1, self % nHalfwidths
      ! Retrieve halfwidth and offset coordinate with respect to box centre
      ! for the current dimension.
      dim = dimensions(i)
      halfwidth = self % halfwidths(i)
      offsetCoord = r(dim) - self % halfwidthsOrigin(i)

      ! Skip if particle is well inside the domain.
      if (abs(offsetCoord) <= halfwidth - surfTol) cycle

      ! Retrieve appropriate BC based on whether offsetCoord < ZERO and apply it.
      bc = self % BCs(merge(2 * i - 1, 2 * i, offsetCoord < ZERO))
      if (bc == REFLECTIVE_BC) then
        u(dim) = -u(dim)
        cycle

      end if
      if (bc == PERIODIC_BC) r(dim) = r(dim) - TWO * sign(halfwidth, offsetCoord)

    end do

  end subroutine explicitCompoundBCs

  !! Subroutine 'transformCompoundBCs'
  !!
  !! Basic description:
  !!   Transforms particle position and direction using boundary conditions.
  !!
  !! Arguments:
  !!   dimensions [in] -> Dimensions along which to apply boundary conditions.
  !!   r [inout]       -> Particle position.
  !!   u [inout]       -> Particle direction.
  !!
  pure subroutine transformCompoundBCs(self, dimensions, r, u)
    class(compoundSurface), intent(in)                           :: self
    integer(shortInt), dimension(self % nHalfwidths), intent(in) :: dimensions
    real(defReal), dimension(3), intent(inout)                   :: r, u
    real(defReal), dimension(self % nHalfwidths)                 :: halfwidths, halfwidthsOrigin
    integer(shortInt), dimension(self % nHalfwidths)             :: nTransforms
    integer(shortInt)                                            :: i, dim, j, bc
    real(defReal)                                                :: offsetCoord, originComponent, halfwidth

    ! Calculate number of transforms in each direction. Reduce halfwidths by surface
    ! tolerance to catch particles at the surface.
    halfwidths = self % halfwidths
    halfwidthsOrigin = self % halfwidthsOrigin
    nTransforms = ceiling(abs(r(dimensions) - halfwidthsOrigin) / (halfwidths - self % getSurfTol())) / 2

    ! Loop over all dimensions.
    do i = 1, self % nHalfwidths
      ! Retrieve origin and halfwidth and apply number of transforms for the current dimension.
      dim = dimensions(i)
      originComponent = halfwidthsOrigin(i)
      halfwidth = halfwidths(i)
      
      do j = 1, nTransforms(i)
        ! Offset coordinates with respect to box centre.
        offsetCoord = r(dim) - originComponent

        ! Retrieve appropriate BC based on whether offsetCoord < ZERO and apply it.
        bc = self % BCs(merge(2 * i - 1, 2 * i, offsetCoord < ZERO))
        if (bc == REFLECTIVE_BC) then
          ! Perform reflection based on distance to the closest plane to offsetCoord.
          r(dim) = -r(dim) + TWO * (sign(halfwidth, offsetCoord) + originComponent)
          u(dim) = -u(dim)
          cycle

        end if

        ! If periodic BC perform translation.
        if (bc == PERIODIC_BC) r(dim) = r(dim) - TWO * sign(halfwidth, offsetCoord)

      end do

    end do

  end subroutine transformCompoundBCs
  
end module compoundSurface_inter