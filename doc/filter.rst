*****************
Digital Filtering
*****************

Introduction
============

The filters discussed in this chapter are based on the following moving data
window which is centered on :math:`i`-th sample:

.. math:: W_i^H = \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+H} \right\}

Here, :math:`H` is a non-negative integer called the *window half-length*, which
represents the number of samples before and after sample :math:`i`.
The total window length is :math:`K = 2 H + 1`.

Nonlinear Digital Filters
=========================

The nonlinear digital filters described below are based on the window median, which is given
by

.. math:: m_i = \textrm{median} \left\{ W_i^H \right\} = \textrm{median} \left\{ x_{i-H}, \dots, x_i, \dots, x_{i+H} \right\}

The median is considered robust to local outliers, unlike the mean.
Median filters can preserve sharp edges while at the same removing signal noise, and are used
in a wide range of applications.

Standard Median Filter
----------------------

The *standard median filter* (SMF) simply replaces the sample :math:`x_i` by the median
:math:`m_i` of the window :math:`W_i^H`: This filter has one tuning parameter given
by :math:`H`. The standard median filter is considered highly resistant to
local outliers and local noise in the data sequence :math:`\{x_i\}`.

.. function:: gsl_filter_median_workspace * gsl_filter_median_alloc(const size_t K)

   This function initializes a workspace for standard median filtering using a symmetric centered moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(7K)`.

.. function:: void gsl_filter_median_free(gsl_filter_median_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_median(const gsl_filter_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_filter_median_workspace * w)

   This function applies a standard median filter to the input :data:`x`, storing the output in :data:`y`.
   The parameter :data:`endtype` specifies how the signal end points are handled. It
   is allowed to have :data:`x` = :data:`y` for an in-place filter.

Recursive Median Filter
-----------------------

The *recursive median filter* (RMF) is a modification of the SMF to include previous filter outputs
in the window before computing the median. The filter's response is

.. math:: y_i = \textrm{median} \left( y_{i-H}, \dots, y_{i-1}, x_i, x_{i+1}, \dots, x_{i+H} \right)

Sometimes, the SMF must be applied several times in a row to achieve adequate smoothing (i.e. a cascade filter).
The RMF, on the other hand, converges to a *root sequence* in one pass,
and can sometimes provide a smoother result than several passes of the SMF. A root sequence is an input which is
left unchanged by the filter.  So there is no need to apply a recursive median filter twice to an input vector.

.. function:: gsl_filter_rmedian_workspace * gsl_filter_rmedian_alloc(const size_t K)

   This function initializes a workspace for recursive median filtering using a symmetric centered moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(K)`.

.. function:: void gsl_filter_rmedian_free(gsl_filter_rmedian_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_rmedian(const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w)

   This function applies a recursive median filter to the input :data:`x`, storing the output in :data:`y`. It
   is allowed to have :data:`x` = :data:`y` for an in-place filter.

Impulse Detection Filter
------------------------

Impulsive noise is characterized by short sequences of data points distinct from those in the
surrounding neighborhood. This section describes a powerful class of filters, also known as
*impulse rejection filters* and *decision-based filters*, designed to detect and remove such outliers from data.
The filter's response is given by

.. math:: y_i = \left\{
                  \begin{array}{cc}
                    x_i, & |x_i - m_i| \le t S_i \\
                    m_i, & |x_i - m_i| > t S_i
                  \end{array}
                \right.

where :math:`m_i` is the median value of the window :math:`W_i^H`, :math:`S_i` is a robust estimate
of the scatter or dispersion for the window :math:`W_i^H`, and :math:`t` is a tuning parameter specifying
the number of scale factors needed to determine that a point is an outlier. The main idea is that the median
:math:`m_i` will be unaffected by a small number of outliers in the window, and so a given
sample :math:`x_i` is tested to determine how far away it is from the median in terms of the local
scale estimate :math:`S_i`. Samples which are more than :math:`t` scale estimates away from the median
are labeled as outliers and replaced by the window median :math:`m_i`. Samples which are less than
:math:`t` scale estimates from the median are left unchanged by the filter.

Note that when :math:`t = 0`, the impulse detection filter is equivalent to the standard median filter. When
:math:`t \rightarrow \infty`, it becomes the identity filter. This means the impulse detection filter can
be viewed as a "less aggressive" version of the standard median filter, becoming less aggressive as :math:`t` is
increased. Note that this filter modifies only samples identified as outliers, while the standard median
filter changes all samples to the local median, regardless of whether they are outliers. This fact, plus
the additional flexibility offered by the additional tuning parameter :math:`t` can make the impulse detection filter
a better choice for some applications.

It is important to have a robust and accurate scale estimate :math:`S_i` in order to
detect impulse outliers even in the presence of noise. The window standard deviation is not
typically a good choice, as it can be significantly perturbed by the presence of even one outlier.
GSL offers the following choices (specified by a parameter of type :type:`gsl_filter_scale_t`) for
computing the scale estimate :math:`S_i`, all of which are robust to the presence of impulse outliers.

.. type:: gsl_filter_scale_t

   This type specifies how the scale estimate :math:`S_i` of the window :math:`W_i^H` is calculated.

   .. macro:: GSL_FILTER_SCALE_MAD

      This option specifies the median absolute deviation (MAD) scale estimate, defined by

      .. math:: S_i = 1.4826 \times \textrm{median} \left\{ | W_i^H - m_i | \right\}

      This choice of scale estimate is also known as the *Hampel filter* in the statistical literature.
      See :ref:`here <sec_mad-statistic>` for more information.

   .. macro:: GSL_FILTER_SCALE_IQR

      This option specifies the interquartile range (IQR) scale estimate, defined as the difference between
      the 75th and 25th percentiles of the window :math:`W_i^H`,

      .. math:: S_i = 0.7413 \left( Q_{0.75} - Q_{0.25} \right)

      where :math:`Q_p` is the p-quantile of the window :math:`W_i^H`. The idea is to throw away the largest
      and smallest 25% of the window samples (where the outliers would be), and estimate a scale from the middle 50%.
      The factor :math:`0.7413` provides an unbiased estimate of the standard deviation for Gaussian data.

   .. macro:: GSL_FILTER_SCALE_SN

      This option specifies the so-called :math:`S_n` statistic proposed by Croux and Rousseeuw.
      See :ref:`here <sec_Sn-statistic>` for more information.

   .. macro:: GSL_FILTER_SCALE_QN

      This option specifies the so-called :math:`Q_n` statistic proposed by Croux and Rousseeuw.
      See :ref:`here <sec_Qn-statistic>` for more information.

.. warning::

   While the scale estimates defined above are much less sensitive to outliers than the standard deviation,
   they can suffer from an effect called *implosion*. The standard deviation of a window :math:`W_i^H` will be zero
   if and only if all samples in the window are equal. However, it is possible for the MAD of a window
   to be zero even if all the samples in the window are not equal. For example, if :math:`K/2 + 1` or more
   of the :math:`K` samples in the window are equal to some value :math:`x^{*}`, then the window median will
   be equal to :math:`x^{*}`. Consequently, at least :math:`K/2 + 1` of the absolute deviations
   :math:`|x_j - x^{*}|` will be zero, and so the MAD will be zero. In such a case, the Hampel
   filter will act like the standard median filter regardless of the value of :math:`t`. Caution should also
   be exercised if dividing by :math:`S_i`.

Because of the possibility of scale implosion, GSL offers a routine :func:`gsl_filter_impulse2` where
the user can input an additional parameter :data:`epsilon`. This parameter is used as a lower bound
on the :math:`S_i`. So for this function, the filter's response is

.. math:: y_i = \left\{
                  \begin{array}{cc}
                    x_i, & |x_i - m_i| \le t S_i \textrm{ or } S_i < \epsilon \\
                    m_i, & |x_i - m_i| > t S_i
                  \end{array}
                \right.

The function :func:`gsl_filter_impulse` sets :math:`\epsilon = 0`.

.. function:: gsl_filter_impulse_workspace * gsl_filter_impulse_alloc(const size_t K)

   This function initializes a workspace for impulse detection filtering using a symmetric moving window of
   size :data:`K`. Here, :math:`H = K / 2`. If :math:`K` is even, it is rounded up to the next
   odd integer to ensure a symmetric window. The size of the workspace is :math:`O(6K)`.

.. function:: void gsl_filter_impulse_free(gsl_filter_impulse_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_filter_impulse(const gsl_filter_end_t endtype, const gsl_filter_scale_t scale_type, const double t, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_impulse_workspace * w)
.. function:: int gsl_filter_impulse2(const gsl_filter_end_t endtype, const gsl_filter_scale_t scale_type, const double epsilon, const double t, const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier, gsl_filter_impulse_workspace * w)

   These functions apply an impulse detection filter to the input vector :data:`x`, storing the filtered output
   in :data:`y`. The tuning parameter :math:`t` is provided in :data:`t`. The lower
   bound :math:`\epsilon` for the scale estimates :math:`S_i` is provided in :data:`epsilon`.
   The window medians :math:`m_i` are stored in :data:`xmedian` and the :math:`S_i` are stored in :data:`xsigma` on output.
   The number of outliers detected is stored in :data:`noutlier` on output, while
   the locations of flagged outliers are stored in the boolean array :data:`ioutlier`. The input
   :data:`ioutlier` may be :code:`NULL` if not desired. It  is allowed to have :data:`x` = :data:`y` for an
   in-place filter.

Examples
========

Square Wave Signal Example
--------------------------

The following example program illustrates the median filters on a noisy
square wave signal. Median filters are well known for preserving sharp
edges in the input signal while reducing noise. The program constructs
a 5 Hz square wave signal with Gaussian noise added. Then the signal is
filtered with a standard median filter and recursive median filter using
a symmetric window of length :math:`K = 7`. The results are shown in
:numref:`fig_filt-edge`.

.. _fig_filt-edge:

.. figure:: /images/filt_edge.png
   :scale: 60%

   Original time series is in gray. The standard median filter output is in
   green and the recursive median filter output is in red.

Both filters preserve the sharp signal edges while reducing the noise. The
recursive median filter achieves a smoother result than the standard median
filter. The "blocky" nature of the output is characteristic of all median
filters. The program is given below.

.. include:: examples/filt_edge.c
   :code:

Impulse Detection Example
-------------------------

The following example program illustrates the impulse detection filter. First,
it constructs a sinusoid signal of length :math:`N = 1000` with Gaussian noise
added. Then, about 1% of the data are perturbed to represent large outliers. An
impulse detecting filter is applied with a window size :math:`K = 25` and
tuning parameter :math:`t = 4`, using the :math:`Q_n` statistic as the robust
measure of scale. The results are plotted in :numref:`fig_impulse`.

.. _fig_impulse:

.. figure:: /images/impulse.png
   :scale: 60%

   Original time series is in blue, filter output is in green, upper and
   lower intervals for detecting outliers are in red and yellow respectively.
   Detected outliers are marked with squares.

The program is given below.

.. include:: examples/impulse.c
   :code:
