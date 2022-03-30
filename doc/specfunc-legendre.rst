.. index::
   single: Legendre polynomials
   single: Legendre functions
   single: spherical harmonics
   single: conical functions
   single: hyperbolic space

The Legendre Functions and Legendre Polynomials are described in
Abramowitz & Stegun, Chapter 8.  These functions are declared in 
the header file :file:`gsl_sf_legendre.h`.

Legendre Polynomials
--------------------

.. function:: double gsl_sf_legendre_P1 (double x)
              double gsl_sf_legendre_P2 (double x)
              double gsl_sf_legendre_P3 (double x)
              int gsl_sf_legendre_P1_e (double x, gsl_sf_result * result)
              int gsl_sf_legendre_P2_e (double x, gsl_sf_result * result)
              int gsl_sf_legendre_P3_e (double x, gsl_sf_result * result)

   These functions evaluate the Legendre polynomials
   :math:`P_l(x)` using explicit representations for :math:`l = 1, 2, 3`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_legendre_Pl (int l, double x)
              int gsl_sf_legendre_Pl_e (int l, double x, gsl_sf_result * result)

   These functions evaluate the Legendre polynomial :math:`P_l(x)`
   for a specific value of :data:`l`, :data:`x` subject to :math:`l \ge 0` and
   :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_Pl_array (int lmax, double x, double result_array[])
              int gsl_sf_legendre_Pl_deriv_array (int lmax, double x, double result_array[], double result_deriv_array[])

   These functions compute arrays of Legendre polynomials
   :math:`P_l(x)` and derivatives :math:`dP_l(x)/dx`
   for :math:`l = 0, \dots, lmax` and :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Q0 (double x)
              int gsl_sf_legendre_Q0_e (double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_0(x)` for
   :math:`x > -1` and :math:`x \ne 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Q1 (double x)
              int gsl_sf_legendre_Q1_e (double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_1(x)` for
   :math:`x > -1` and :math:`x \ne 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Ql (int l, double x)
              int gsl_sf_legendre_Ql_e (int l, double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_l(x)` for
   :math:`x > -1`, :math:`x \ne 1` and :math:`l \ge 0`.
.. Exceptional Return Values: GSL_EDOM

Associated Legendre Polynomials and Spherical Harmonics
-------------------------------------------------------

The following functions compute the associated Legendre polynomials
:math:`P_l^m(x)` which are solutions of the differential equation

.. only:: not texinfo

   .. math:: (1 - x^2) {d^2 \over dx^2} P_l^m(x) - 2x {d \over dx} P_l^m(x) +
             \left( l(l+1) - {m^2 \over 1 - x^2} \right) P_l^m(x) = 0

.. only:: texinfo

   ::

      (1 - x^2) d^2 P_l^m(x) / dx^2 P_l^m(x) - 2x d/dx P_l^m(x) +
      ( l(l+1) - m^2 / (1 - x^2) ) P_l^m(x) = 0

where the degree :math:`l` and order :math:`m` satisfy :math:`0 \le l` and
:math:`0 \le m \le l`.
The functions :math:`P_l^m(x)` grow combinatorially with
:math:`l` and can overflow for :math:`l` larger than about 150.
Alternatively, one may calculate normalized associated Legendre
polynomials. There are a number of different normalization conventions,
and these
functions can be stably computed up to degree and order 2700. The
following normalizations are provided:

* Schmidt semi-normalization

  Schmidt semi-normalized associated Legendre polynomials are often
  used in the geophysical community and are defined as (Winch et al, 2005)

  .. only:: not texinfo

     .. math::

        S_l^0(x) &= P_l^0(x) \\
        S_l^m(x) &= (-1)^m \sqrt{2 {(l-m)! \over (l+m)!}} P_l^m(x), m > 0 

  .. only:: texinfo

     ::

        S_l^0(x) = P_l^0(x)
        S_l^m(x) = (-1)^m \sqrt((2(l-m)! / (l+m)!)) P_l^m(x), m > 0 

  The factor of :math:`(-1)^m` is called the Condon-Shortley phase
  factor and can be included if desired by setting the flag
  :macro:`GSL_SF_LEGENDRE_FLG_CSPHASE` in the function :func:`gsl_sf_legendre_precompute`
  below. These functions satisfy the normalization condition,

  .. only:: not texinfo

     .. math::

        \int_{-1}^1 S_k^m(x) S_l^m(x) dx =
        \left\{
          \begin{array}{ll}
            \frac{2}{2l+1} \delta_{kl}, & m = 0 \\
            \frac{4}{2l+1} \delta_{kl}, & m > 0
          \end{array}
        \right.

  .. only:: texinfo

     ::

        \int_{-1}^1 S_k^m(x) S_l^m(x) dx = { 2/(2l+1) \delta_{kl}, m = 0
                                           { 4/(2l+1) \delta_{kl}, m > 0

* Spherical Harmonic Normalization

  The associated Legendre polynomials suitable for calculating spherical
  harmonics are defined as

  .. only:: not texinfo

     .. math:: Y_l^m(x) = (-1)^m \sqrt{{2l + 1 \over 4 \pi} {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo

     ::

        Y_l^m(x) = (-1)^m \sqrt((2l + 1) * (l-m)! / (4 \pi) / (l+m)!) P_l^m(x)

  where again the phase factor :math:`(-1)^m` can be included or excluded
  if desired. These functions satisfy the normalization condition,

  .. only:: not texinfo

     .. math::

        \int_{-1}^1 Y_k^m(x) Y_l^m(x) dx = \frac{\delta_{kl}}{2\pi}

  .. only:: texinfo

     ::

        \int_{-1}^1 Y_k^m(x) Y_l^m(x) dx = \delta_{kl} / (2 \pi)

  Note that these functions, when coupled with the factor
  :math:`e^{i m \phi}` produce the orthonormalized complex spherical
  harmonics.

* Full Normalization

  The fully normalized associated Legendre polynomials are defined as

  .. only:: not texinfo

     .. math:: N_l^m(x) = (-1)^m \sqrt{(l + {1 \over 2}) {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo
  
     ::
     
        N_l^m(x) = (-1)^m \sqrt((l + 1/2) (l-m)! / (l+m)!) P_l^m(x)

  and satisfy the normalization condition,

  .. math:: \int_{-1}^1 N_k^m(x) N_l^m(x) dx = \delta_{kl}

* :math:`4 \pi` Normalization

  The :math:`4 \pi` normalized associated Legendre polynomials are often used in geodesy and are defined as

  .. only:: not texinfo

     .. math:: R_l^m(x) = (-1)^m \sqrt{(2 - \delta_{m0}) (2 l + 1) {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo
  
     ::
     
        R_l^m(x) = (-1)^m \sqrt((2 - \delta_{m0}) (2 l + 1) (l-m)! / (l+m)!) P_l^m(x)

  These functions satisfy the normalization condition,

  .. only:: not texinfo

     .. math:: \int_{-1}^1 R_k^m(x) R_l^m(x) dx = 2 \left( 2 - \delta_{m0} \right) \delta_{kl}

  .. only:: texinfo

     ::

        \int_{-1}^1 R_k^m(x) R_l^m(x) dx = 2 (2 - \delta_{m0}) \delta_{kl}

  When used in the definition of real spherical harmonics, they satisfy a
  :math:`4\pi` normalization condition when integrated over the unit sphere.
  More information on these functions can be found in Hofmann-Wellenhof and Moritz, 2006.

The normalized associated Legendre routines below use a recurrence
relation which is stable up to a degree and order of about 2700.
Beyond this, the computed functions could suffer from underflow
leading to incorrect results. Routines are provided to compute
first and second derivatives
:math:`dP_l^m(x)/dx` and :math:`d^2 P_l^m(x)/dx^2` as well as their alternate
versions :math:`d P_l^m(\cos{\theta})/d\theta` and
:math:`d^2 P_l^m(\cos{\theta})/d\theta^2`. While there is a simple
scaling relationship between the two forms, the derivatives
involving :math:`\theta` are heavily used in spherical harmonic
expansions, and also do not suffer from singularities at the poles,
:math:`x = \pm 1`, and so these routines are also provided.

In the functions below, a parameter of type :type:`gsl_sf_legendre_t`
specifies the type of normalization to use. The possible values are

.. type:: gsl_sf_legendre_t

   ================================== ===============================================================================
   Value                              Description
   ================================== ===============================================================================
   :code:`GSL_SF_LEGENDRE_NONE`       The unnormalized associated Legendre polynomials :math:`P_l^m(x)`
   :code:`GSL_SF_LEGENDRE_SCHMIDT`    The Schmidt semi-normalized associated Legendre polynomials :math:`S_l^m(x)`
   :code:`GSL_SF_LEGENDRE_SPHARM`     The spherical harmonic associated Legendre polynomials :math:`Y_l^m(x)`
   :code:`GSL_SF_LEGENDRE_FULL`       The fully normalized associated Legendre polynomials :math:`N_l^m(x)`
   :code:`GSL_SF_LEGENDRE_FOURPI`     The :math:`4\pi` normalized associated Legendre polynomials :math:`R_l^m(x)`
   ================================== ===============================================================================

In the routines below which return an array of ALFs, there are two possible
indexing schemes which can be used.

* M-major indexing

  This scheme uses the following indexing function to locate an ALF of degree :math:`l` and order
  :math:`m`,

  .. math:: \mathcal{I}_m(l,m,L) = m L - \frac{m(m-1)}{2} + l

  where :math:`L` is the maximum degree, :code:`lmax`. This corresponds to the following memory layout,

  .. math:: l \quad \overbrace{0 \; 1 \; 2 \; \cdots \; L}^{m = 0} \quad \overbrace{1 \; 2 \; \cdots \; L}^{m = 1} \quad \overbrace{2 \; 3 \; \cdots \; L}^{m = 2} \quad \cdots \quad \overbrace{L}^{m = L}

  This is the default method and is recommended since it corresponds to the order that ALFs
  are computed in their recurrence relations, and is therefore cache-efficient. The following
  code demonstrates how to access the array elements in order,

  .. code::

     idx = 0;
     for (m = 0; m <= lmax; ++m) {
       for (l = m; l <= lmax; ++l) {
         double value = Plm[idx]; /* (l,m) element */
         ++idx;
       }
     }

* L-major indexing

  This scheme uses the following indexing function to locate an ALF of degree :math:`l` and order
  :math:`m`,

  .. math:: \mathcal{I}_l(l,m) = \frac{l(l+1)}{2} + m

  This corresponds to the following memory layout,

  .. math:: m \quad \overbrace{0}^{l = 0} \quad \overbrace{0 \; 1}^{l = 1} \quad \overbrace{0 \; 1 \; 2}^{l = 2} \quad \cdots \quad \overbrace{0 \; 1 \; 2 \; \cdots \; L}^{l = L}

  The following code demonstrates how to access array elements in order,

  .. code::

     idx = 0;
     for (l = 0; l <= lmax; ++l) {
       for (m = 0; m <= l; ++m) {
         double value = Plm[idx]; /* (l,m) element */
         ++idx;
       }
     }

  .. important::

     The L-major indexing scheme was the default in GSL versions 2.7 and earlier. However, it
     is not as cache efficient as the M-major indexing scheme, and so GSL v2.8 and later
     use M-major indexing by default. It is recommended to use M-major indexing to get
     maximum performance when computing ALFs.

.. function:: int gsl_sf_legendre_precompute (const gsl_sf_legendre_t norm, const size_t lmax, const size_t flags, double result_array[])

   This function precomputes the multiplicative factors needed for the
   associated Legendre recurrence relations. The input :data:`norm`
   specifies the ALF normalization. The input :data:`lmax` specifies
   the maximum ALF degree. The input :data:`flags` is a bitmask which
   specifies how the ALFs are computed and stored. It can contain the
   following values,

   .. macro:: GSL_SF_LEGENDRE_FLG_CSPHASE

      This flag will include the Condon-Shortley phase factor in the calculation
      of the ALFs

   .. macro:: GSL_SF_LEGENDRE_FLG_INDEXL

      This flag will store the output arrays using L-major indexing,
      :math:`\mathcal{I}_l(l,m)`. If this flag is not set, the output
      arrays will use M-major indexing, :math:`\mathcal{I}_m(l,m,L)`.

   The output array :data:`result_array` should have a length
   as returned by the function :func:`gsl_sf_legendre_array_n`. The computed
   recurrence factors are stored at the end of :data:`result_array`, leaving
   room at the front for the calculation of the ALFs.

   This routine must be called prior to calling ALF functions with the suffix
   :code:`arrayx`.

.. function:: int gsl_sf_legendre_arrayx (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[])

   This function calculates all associated Legendre polynomials
   for :math:`0 \le l \le lmax` and :math:`0 \le m \le l` for :math:`|x| \le 1`.
   The :data:`norm` parameter specifies which normalization is used.
   The normalized :math:`P_l^m(x)` values are stored in :data:`result_array`, whose
   minimum size can be obtained from calling :func:`gsl_sf_legendre_array_n`.
   The array index of :math:`P_l^m(x)` is obtained from calling
   :code:`gsl_sf_legendre_array_index(l, m)`.

   The function :func:`gsl_sf_legendre_precompute` must be called first
   using the same :data:`norm` and :data:`lmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_legendre_deriv_alt_arrayx (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[])

   This function calculates all associated Legendre
   functions and their (alternate) first derivatives for :math:`0 \leq l \leq lmax`
   and :math:`0 \leq m \leq l` for :math:`|x| \leq 1`.
   The :data:`norm` parameter specifies which normalization
   is used. The :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(\cos{\theta})/d\theta` are stored in :data:`result_array` and
   :data:`result_deriv_array` respectively.

   The function :func:`gsl_sf_legendre_precompute` must be called first
   using the same :data:`norm` and :data:`lmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_legendre_deriv_arrayx (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[])

   This function calculates all associated Legendre
   functions and their first derivatives for :math:`0 \leq l \leq lmax`
   and :math:`0 \leq m \leq l` for :math:`|x| < 1`.
   The :data:`norm` parameter specifies which
   normalization is used. The :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(x)/dx` are stored in :data:`result_array` and
   :data:`result_deriv_array` respectively.

   Note that for some orders :math:`m`, the derivatives :math:`dP_l^m(x)/dx` have
   singularities at the end points :math:`x = \pm 1`, and so this function only
   accepts interior points as input, :math:`x \in (-1,1)`.

   The function :func:`gsl_sf_legendre_precompute` must be called first
   using the same :data:`norm` and :data:`lmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_legendre_deriv2_alt_arrayx (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   This function calculates all associated Legendre
   functions and their (alternate) first and second derivatives for
   :math:`0 \leq l \leq lmax` and :math:`0 \leq m \leq l` for :math:`|x| \leq 1`.
   The :data:`norm` parameter specifies which normalization is used. The normalized :math:`P_l^m(x)`
   values and their derivatives :math:`dP_l^m(\cos{\theta})/d\theta` and
   :math:`d^2 P_l^m(\cos{\theta})/d\theta^2` are stored in :data:`result_array`,
   :data:`result_deriv_array`, and :data:`result_deriv2_array` respectively.

   The function :func:`gsl_sf_legendre_precompute` must be called first
   using the same :data:`norm` and :data:`lmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_legendre_deriv2_arrayx (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   This function calculates all associated Legendre
   functions and their first and second derivatives for
   :math:`0 \leq l \leq lmax` and :math:`0 \leq m \leq l` for :math:`|x| < 1`.
   The :data:`norm` parameter specifies which
   normalization is used. The :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(x)/dx` and :math:`d^2 P_l^m(x)/dx^2` are stored in
   :data:`result_array`, :data:`result_deriv_array`, and
   :data:`result_deriv2_array` respectively.

   Note that for some orders :math:`m`, the derivatives :math:`dP_l^m(x)/dx`
   and :math:`d^2 P_l^m(x)/dx^2` have
   singularities at the end points :math:`x = \pm 1`, and so this function only
   accepts interior points as input, :math:`x \in (-1,1)`.

   The function :func:`gsl_sf_legendre_precompute` must be called first
   using the same :data:`norm` and :data:`lmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_legendre_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[])
              int gsl_sf_legendre_deriv_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[])
              int gsl_sf_legendre_deriv_alt_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[])
              int gsl_sf_legendre_deriv2_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])
              int gsl_sf_legendre_deriv2_alt_array (const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   These functions are similar to their :code:`arrayx` counterparts above,
   except they compute the recurrence factors by calling the function
   :func:`gsl_sf_legendre_precompute` on each call. These functions
   omit the Condon-Shortley phase factor, and store the output arrays
   in L-major order. These functions are provided for backward compatibility,
   but it is recommended to use instead
   the :code:`arrayx` functions which are more efficient when calculating
   ALFs for multiple input points :math:`x`.

.. function:: size_t gsl_sf_legendre_nlm(const size_t lmax)

   This function returns the total number of associated Legendre
   functions :math:`P_l^m(x)` for a given :data:`lmax`. The number is
   :code:`(lmax+1) * (lmax+2) / 2`.

.. function:: size_t gsl_sf_legendre_array_n (const size_t lmax)

   This function returns the minimum array size for maximum degree :data:`lmax`
   needed for the array versions of the associated Legendre functions.
   Size is calculated as the total number of :math:`P_l^m(x)` functions
   (see :func:`gsl_sf_legendre_nlm`),
   plus extra space for precomputing multiplicative factors used in the
   recurrence relations.

.. function:: size_t gsl_sf_legendre_array_index (const size_t l, const size_t m)

   This function returns the index into :data:`result_array`,
   :data:`result_deriv_array`, or :data:`result_deriv2_array` corresponding
   to :math:`P_l^m(x)`, :math:`P_l^{'m}(x)`, or :math:`P_l^{''m}(x)`. The
   index is given by :math:`l(l+1)/2 + m`.

   An inline version of this function is used if :macro:`HAVE_INLINE` is
   defined.

.. function:: double gsl_sf_legendre_Plm (int l, int m, double x)
              int gsl_sf_legendre_Plm_e (int l, int m, double x, gsl_sf_result * result)

   These routines compute the associated Legendre polynomial
   :math:`P_l^m(x)` for :math:`m \ge 0`,
   :math:`l \ge m`, and :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: double gsl_sf_legendre_sphPlm (int l, int m, double x)
              int gsl_sf_legendre_sphPlm_e (int l, int m, double x, gsl_sf_result * result)

   These routines compute the normalized associated Legendre polynomial
   (with Condon-Shortley phase)
   :math:`(-1)^m \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)` suitable
   for use in spherical harmonics.  The parameters must satisfy :math:`m \ge 0`,
   :math:`l \ge m`, and :math:`|x| \le 1`.
   These routines avoid the overflows
   that occur for the standard normalization of :math:`P_l^m(x)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[])
              int gsl_sf_legendre_deriv_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[])
              int gsl_sf_legendre_deriv_alt_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[])
              int gsl_sf_legendre_deriv2_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[], double result_deriv2_array[])
              int gsl_sf_legendre_deriv2_alt_array_e (const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   These functions are deprecated and will be removed in a future release.

Conical Functions
-----------------

The Conical Functions :math:`P^\mu_{-(1/2)+i\lambda}(x)`
and :math:`Q^\mu_{-(1/2)+i\lambda}`
are described in Abramowitz & Stegun, Section 8.12.

.. function:: double gsl_sf_conicalP_half (double lambda, double x)
              int gsl_sf_conicalP_half_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the irregular Spherical Conical Function
   :math:`P^{1/2}_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_mhalf (double lambda, double x)
              int gsl_sf_conicalP_mhalf_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the regular Spherical Conical Function
   :math:`P^{-1/2}_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_0 (double lambda, double x)
              int gsl_sf_conicalP_0_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the conical function
   :math:`P^0_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_1 (double lambda, double x)
              int gsl_sf_conicalP_1_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the conical function 
   :math:`P^1_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_sph_reg (int l, double lambda, double x)
              int gsl_sf_conicalP_sph_reg_e (int l, double lambda, double x, gsl_sf_result * result)

   These routines compute the Regular Spherical Conical Function
   :math:`P^{-1/2-l}_{-1/2 + i \lambda}(x)`
   for :math:`x > -1` and :math:`l \ge -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_cyl_reg (int m, double lambda, double x)
              int gsl_sf_conicalP_cyl_reg_e (int m, double lambda, double x, gsl_sf_result * result)

   These routines compute the Regular Cylindrical Conical Function
   :math:`P^{-m}_{-1/2 + i \lambda}(x)`
   for :math:`x > -1` and :math:`m \ge -1`.
.. Exceptional Return Values: GSL_EDOM

Radial Functions for Hyperbolic Space
-------------------------------------

The following spherical functions are specializations of Legendre
functions which give the regular eigenfunctions of the Laplacian on a
3-dimensional hyperbolic space :math:`H^3`.  Of particular interest is
the flat limit, :math:`\lambda \to \infty`, :math:`\eta \to 0`,
:math:`\lambda\eta` fixed.
  
.. function:: double gsl_sf_legendre_H3d_0 (double lambda, double eta)
              int gsl_sf_legendre_H3d_0_e (double lambda, double eta, gsl_sf_result * result)

   These routines compute the zeroth radial eigenfunction of the Laplacian on the
   3-dimensional hyperbolic space,

   .. math:: L^{H3d}_0(\lambda,\eta) := {\sin(\lambda\eta) \over \lambda\sinh(\eta)}

   for :math:`\eta \ge 0`.
   In the flat limit this takes the form
   :math:`L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_H3d_1 (double lambda, double eta)
              int gsl_sf_legendre_H3d_1_e (double lambda, double eta, gsl_sf_result * result)

   These routines compute the first radial eigenfunction of the Laplacian on
   the 3-dimensional hyperbolic space,

   .. math:: L^{H3d}_1(\lambda,\eta) := {1\over\sqrt{\lambda^2 + 1}} {\left(\sin(\lambda \eta)\over \lambda \sinh(\eta)\right)} \left(\coth(\eta) - \lambda \cot(\lambda\eta)\right)

   for :math:`\eta \ge 0`
   In the flat limit this takes the form 
   :math:`L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_H3d (int l, double lambda, double eta)
              int gsl_sf_legendre_H3d_e (int l, double lambda, double eta, gsl_sf_result * result)

   These routines compute the :data:`l`-th radial eigenfunction of the
   Laplacian on the 3-dimensional hyperbolic space :math:`\eta \ge 0` and
   :math:`l \ge 0`.
   In the flat limit this takes the form
   :math:`L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_H3d_array (int lmax, double lambda, double eta, double result_array[])

   This function computes an array of radial eigenfunctions
   :math:`L^{H3d}_l( \lambda, \eta)`
   for :math:`0 \le l \le lmax`.
.. Exceptional Return Values:
