/*********************************************************
*
*  Martin Eckstein, 2010
*  Contour Green funtion class
*
*********************************************************/
#ifndef GREEN_CNTR_FULL
#define GREEN_CNTR_FULL

/* #######################################################################################
#
#                single-particle Greenfunctions G(t,t')
#
##########################################################################################

     currently there exist three closely related Greenfunction types ...
	 
	 -  type T is double or float
     -  timesteps 0... ntau on imag contour,
	    0...nt on upper/lower real contour.
		for nt=-1, only the matsubara sector is stored
     -  each element has element_size_=size2_*size1_ complex<T> numbers
	    (only square matices allowed (size1_=size2_)
     -  strored elements:
        gtr(t,t'),les(t,t'),tv(t,t'),vt(t,t')
		mat(tau),mat(-tau)
		ret(t,t')=gtr-les  (stored for all t,t')
		
####################################################################################### */

/*!
    \mainpage Welcome to NESSi — Non-Equilibrium Systems Simulation

    \section S1 What is NESSi?

    \c NESSi is an open-source software package for the manipulation of nonequilibrium
    Green's functions defined on the Kadanoff–Baym contour.
    The Green's function method in its time-dependent formulation is a versatile framework
    for solving interacting many-body problems out of equilibrium.

    \c NESSi provides classes representing the various types of Green's functions,
    implements the basic operations on these functions, and allows one to solve the
    corresponding equations of motion. The library is aimed at the study of transient
    dynamics from an initial equilibrium state, induced by time-dependent model parameters.

    Key features:
    - Tools for constructing Feynman diagrams and solving equations of motion for
      nonequilibrium Green's functions on the Kadanoff–Baym contour.
    - High-order quadrature rules: for \f$N\f$ time slices, the error scales as
      \f$\mathcal{O}(N^{-p})\f$ with \f$p\f$ up to 7.
    - Efficient distributed-memory parallelization over reciprocal space for
      large-scale calculations on extended systems.
    - Memory-truncated time propagation and steady-state calculations (v2.0.0).

    \section S2 Full Documentation

    The full manual including background, usage examples, and API reference
    is available online at: [link to be added]
*/


#include "cntr_decl.hpp"
#include "cntr_impl.hpp"
#ifndef CNTR_NO_EXTERN_TEMPLATES
#include "cntr_extern_templates.hpp"
#endif

#ifdef CNTR_USE_NESS
#include "ness2.hpp"
#endif




#endif // GREEN_CNTR_FULL
