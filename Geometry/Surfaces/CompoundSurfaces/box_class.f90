module box_class

  use compoundSurface_inter, only : compoundSurface
  use dictionary_class,      only : dictionary
  use errors_mod,            only : fatalError
  use genericProcedures,     only : append, areEqual, numToChar, swap
  use numPrecision
  use ratint_mod,            evaluate_ratint => evaluate
  use surface_inter,         only : kill_super => kill
  use universalVariables,    only : INF

  implicit none
  private

  !!
  !! Axis Aligned Box
  !!
  !! F(r) = maxval(abs(r - o) - a)
  !!
  !! Where: a -> halfwidth vector, o-> origin position
  !!        maxval -> maximum element (L_inf norm)
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   aab { type box; id 92; origin (0.0 0.0 9.0); halfwidth(1.0 2.0 0.3);}
  !!
  !! Boundary Conditions:
  !!   BC order: x_min, x_max, y_min, y_max, z_min, z_max
  !!   Each face can have diffrent BC. Any combination is supported with coordinates transform.
  !!
  !! Private Members:
  !!   origin     -> 3-D coordinates of the middle of the box
  !!   halfwidths -> Half-lengths of the box in each direction (must be > ZERO)
  !!   BCs        -> Boundary conditions flags for each face (x_min, x_max, y_min, y_max, z_min, z_max)
  !!
  !! Interface:
  !!   compoundSurface interface
  !!
  type, public, extends(compoundSurface) :: box
    private
    integer(shortInt) :: nBCs = 6
  contains
    ! Superclass procedures.
    procedure          :: cropsBoundingBox
    procedure          :: init
    procedure          :: evaluate
    procedure          :: distance
    procedure, private :: distance_rational
    procedure          :: entersPositiveHalfspace
    procedure          :: kill
    procedure          :: setBCs
    procedure          :: explicitBC
    procedure          :: transformBC
  end type box

contains
  !! Subroutine 'cropsBoundingBox'
  !!
  !! Basic description:
  !!   Checks whether the box crops a given bounding box.
  !!
  !! Arguments:
  !!   boundingBox -> A defReal array containing the x_min, y_min, z_min, x_max, y_max and z_max
  !!                  coordinates of the bounding box.
  !!
  !! Result:
  !!   doesIt      -> .true. if the box crops the bounding box.
  !!
  pure function cropsBoundingBox(self, boundingBox) result(doesIt)
    class(box), intent(in)                  :: self
    real(defReal), dimension(6), intent(in) :: boundingBox
    logical(defBool)                        :: doesIt
    real(defReal), dimension(3)             :: halfwidths
    integer(shortInt)                       :: i

    ! Initialise doesIt = .true., retrieve the box's halfwidths then loop through all dimensions.
    doesIt = .true.
    halfwidths = self % getHalfwidths()
    do i = 1, 3
      if (any(abs(boundingBox([i, 3 + i])) >= halfwidths(i))) return

    end do

    ! If reached here, update doesIt = .false.
    doesIt = .false.

  end function cropsBoundingBox

  !!
  !! Initialise box from a dictionary
  !!
  !! See surface_inter for more details
  !!
  subroutine init(self, dict)
    class(box), intent(inout)                :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: origin, halfwidths
    real(defReal), dimension(6)              :: boundingBox

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)

    ! Load origin.
    call dict % get(origin, 'origin')
    call self % setOrigin(origin, 3)

    ! Load halfwidths.
    call dict % get(halfwidths, 'halfwidth')
    call self % setHalfwidths(halfwidths, origin, 3)

    ! Set bounding box.
    boundingBox(1:3) = origin - halfwidths
    boundingBox(4:6) = origin + halfwidths
    call self % setBoundingBox(boundingBox)
    call self % setType('box')

    ! Initialise BCs.
    call self % setCompoundBCs(self % nBCs)

  end subroutine init
  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    ! Compute c from origin-centred coordinates and box halfwidths.
    c = self % evaluateCompound(r)

  end function evaluate

  !!
  !!
  !!
  pure function distance(self, r, u) result(d)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt)                       :: i
    logical(defBool)                        :: surfTolCondition
    real(defReal)                           :: d, inverseU, tFar, tNear, t1, t2
    real(defReal), dimension(3, 2)          :: bounds
    
    ! Initialise d = INF, tNear = -INF, tFar = INF.
    d = INF
    tNear = -INF
    tFar = INF

    surfTolCondition = abs(self % evaluate(r)) < self % getSurfTol()

    ! Retrieve bounds then loop over all dimensions.
    bounds = self % getBounds()

    do i = 1, 3
      if((areEqual(r(i), bounds(i, 1)) .or. areEqual(r(i), bounds(i, 2))) .and. areEqual(u(i), ZERO)) then
        d = self % distance_rational(surfTolCondition, convert_ieee(r), convert_ieee(u))
        return

      end if

      ! Perform early check to see if the particle is outside the slab and moving away from it along
      ! the current dimension. If yes the particle cannot intersect the slab and we can return early.
      if((r(i) <= bounds(i, 1) .and. u(i) <= ZERO) .or. (r(i) >= bounds(i, 2) .and. u(i) >= ZERO)) return

      if(areEqual(u(i), ZERO)) cycle
      inverseU = ONE / u(i)
      t1 = (bounds(i, 1) - r(i)) * inverseU
      t2 = (bounds(i, 2) - r(i)) * inverseU
      if(t2 < t1) call swap(t1, t2)

      tNear = max(tNear, t1)
      tFar = min(tFar, t2)

      ! Return early if intersection is impossible (far intersection is definitely greater than near intersection, or far
      ! intersection is definitely negative).
      if((tFar < tNear .and. .not. areEqual(tNear, tFar)) .or. &
          (tFar < ZERO .and. .not. areEqual(tFar, ZERO))) return

    end do

    ! If results are ambiguous, launch an exact computation.
    if(areEqual(tNear, tFar) .or. &
       (.not. surfTolCondition .and. (areEqual(tNear, ZERO) .or. areEqual(tFar, ZERO)))) then
      d = self % distance_rational(surfTolCondition, convert_ieee(r), convert_ieee(u))

    else
      ! Take the far intersection if the particle is on the surface or already inside it.
      if((surfTolCondition .and. abs(tFar) >= abs(tNear)) .or. &
         (.not. surfTolCondition .and. tNear <= ZERO)) then
        d = tFar

      else
        d = tNear

      end if

    end if

    ! Cap distance to INF if d <= ZERO or d > INF.
    if(d <= ZERO .or. d > INF) d = INF
  
  end function distance

  !!
  !!
  !!
  pure function distance_rational(self, surfTolCondition, r, u) result(d)
    class(box), intent(in)                   :: self
    logical(defBool), intent(in)             :: surfTolCondition
    type(ratint_t), dimension(3), intent(in) :: r, u
    integer(shortInt)                        :: i
    logical(defBool)                         :: areDistancesInvalid
    real(defReal)                            :: d
    type(ratint_t)                           :: temp, tFar, tNear, t1, t2, ZERO_rational
    type(ratint_t), dimension(3, 2)          :: bounds

    ! Pre-compute ZERO_rational.
    ZERO_rational = convert_int(0_longInt)
    
    ! Initialise areDistancesInvalid = .true. and d = INF.
    areDistancesInvalid = .true.
    d = INF

    ! Retrieve bounds then loop over all dimensions.
    bounds = self % getRationalBounds()
    do i = 1, 3
      if((r(i) == bounds(i, 1) .or. r(i) == bounds(i, 2)) .and. isZero(u(i))) return

      ! Perform early check to see if the particle is outside the slab and moving away from it along
      ! the current dimension. If yes the particle cannot intersect the slab and we can return early.
      if((bounds(i, 1) >= r(i) .and. ZERO_rational >= u(i)) .or. (r(i) >= bounds(i, 2) .and. u(i) >= ZERO_rational)) return

      if(isZero(u(i))) cycle
      t1 = (bounds(i, 1) - r(i)) / u(i)
      t2 = (bounds(i, 2) - r(i)) / u(i)

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

    ! Return if distance are still invalid (should never happen for a unit vector).
    if(areDistancesInvalid) return

    ! Take the far intersection if the particle is on the surface or already inside it.
    if((surfTolCondition .and. absoluteValue(tFar) >= absoluteValue(tNear)) .or. &
       (.not. surfTolCondition .and. ZERO_rational >= tNear)) then
      d = evaluate_ratint(tFar)

    else
      d = evaluate_ratint(tNear)

    end if

    ! Cap distance to INF if d <= ZERO or d > INF.
    if(d <= ZERO .or. d > INF) d = INF
  
  end function distance_rational

  !! Function 'entersPositiveHalfspace'
  !!
  !! Basic description:
  !!   Returns .true. if the particle is going into positive halfspace.
  !!
  !! Detailed description:
  !!   Works similarly to the 'distance' function, except that here we are only interested
  !!   in the value of tMax (renamed t). The algorithm loops over all dimensions and checks
  !!   whether the particle is within surface tolerance of a halfwidth for the current 
  !!   dimension. If the direction component of a dimension for which the particle is within 
  !!   surface tolerance to a halfwidth is ZERO (ie, the ray is parallel to the halfwidth), 
  !!   then halfspace is determined based on which side of the halfwidth the particle lies. 
  !!   Else, t is updated for the current dimension. 
  !!
  !!   Then, the particle will cross into a positive halfspace if t < maxDist, where maxDist 
  !!   is the maximum possible distance that the particle can travel before either crossing 
  !!   one (or more) halfwidths or being outside surface tolerance of all three halfwidths. 
  !!   This distance is obtained by considering a particle lying within surface tolerance to 
  !!   all three halfwidths (ie, very close to a corner of the box) and having a direction 
  !!   vector u = [sqrt(3) / 3, sqrt(3) / 3, sqrt(3) / 3]; in this case applying Pythagoras 
  !!   theorem in 3-D gives maxDist = sqrt(3) * (surface tolerance). 
  !!
  !! Arguments:
  !!   r [in] -> Coordinates of the ray's origin.
  !!   u [in] -> Direction of the ray.
  !!
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: isHalfspacePositive

    isHalfspacePositive = self % isHalfspacePositive(r, u)

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(box), intent(inout) :: self

    ! Superclass.
     call kill_super(self)

    ! Compound.
    call self % killCompound()

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBCs(self, BCs)
    class(box), intent(inout)                   :: self
    integer(shortInt), dimension(:), intent(in) :: BCs

    call self % setCompoundBCs(self % nBCs, BCs)
    
  end subroutine setBCs

  !!
  !! Apply explicit boundary conditions.
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!
  pure subroutine explicitBC(self, r, u)
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r, u

    call self % explicitCompoundBCs([1, 2, 3], r, u)

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BC
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Order of transformations does not matter
  !!   - Calculate distance (in # of transformations) for each direction and apply them
  !!
  pure subroutine transformBC(self, r, u)
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r, u

    call self % transformCompoundBCs([1, 2, 3], r, u)

  end subroutine transformBC

end module box_class