About CauseMap
=============

.. raw:: html

    <p>
    <table>
    <tr>

:raw-html:`<td>` |image1| :raw-html:`</td>`
:raw-html:`<td>` |image2| :raw-html:`</td>`
:raw-html:`<td>` |image3| :raw-html:`</td>`

.. raw:: html

    </tr>
    </table>
    </p>


Convergent cross mapping (CCM) provides a model-free approach to detecting dependencies
and establishing causality in complex non-linear systems, even in the presence of feedback loops and
unmeasured confounding [1]_. CCM derives this power from
explicitly capturing time-dependent dynamics through a technique known as
state-space reconstruction (SSR). In practice, this analysis
typically requires 25 or more time points, measured with sufficient density
to capture system dynamics.

|br|

How CauseMap works
------------------
To get a feel for the algorithm and the intuition behind it, we highly recommend
the following two videos from `George Sugihara <http://scrippsscholars.ucsd.edu/gsugihara/biocv>`_'s group at UCSD:


Time series as projections of full system dynamics
++++++++++++++++++++++++++++++++++++++++++++++++++
The following video will explain this concept using the Lorenz attractor (or butterfly manifold) as an example system.

.. raw:: html

   <iframe width="610" height="320" src="http://www.youtube.com/embed/6i57udsPKms" frameborder="1" allowfullscreen></iframe>

|

Reconstructing full system dynamics using mathemagics!
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. raw:: html

   <iframe width="610" height="320" src="http://www.youtube.com/embed/NrFdIz-D2yM" frameborder="1" allowfullscreen></iframe>


|

In summary, the first animation explains how time series can be viewed as projections of
higher-dimensional system dynamics. In this light, one can see that time series
of individual variables each contain information about their full causal system.
The second animation explains a corollary to this observation: shadows of full
causal systems can be constructed using individual time series. For further details,
please feel free to consult the original reference (below), or our description in
the `CauseMap manuscript <https://peerj.com/articles/824/>`_.

|

.. _params:

Model Parameters
----------------
Library size
++++++++++++
Libraries are the collection of data points used for cross mapping.
Convergent cross mapping gets its name because, if a relationship is truly causal,
predictive skill should converge with library size (read: amount of data used).
Library-relevant parameters are:

- ``libsizemin``: The minimum library size.
- ``libsizemax``: The maximum library size. For long time series, for example, you may not want to use the largest possible library.
- ``nboots``: This allows you to switch from the traditional "slide window" approach to bootstrap library selection. As an example, ``nboots=20`` will randomly select 20 libraries at each library size.

Note: the bootstrapping approach generally performs better than
sliding windows because it isn't as sensitive to secular trends.

**In future releases, this will be the default.**

Predicted data
+++++++++++++++
- ``predstart``: The first point to be predicted from the target time series.
- ``npred``: The number of points to be predicted by each library.
- ``tau_p``: The assumed time lag for the causal effect.

Manifold reconstruction
++++++++++++++++++++++++
- ``E_vals``: This sets the range of manifold dimensions to test. The higher dimensional reconstructions are necessary for causal systems with a larger number of components.
- ``tau_s``: The time lag to use for manifold reconstruction. A value of one is customary, but data may be thinned if autocorrelation is strong.


.. [1] Sugihara,G. *et al.* (2012) Detecting causality in complex ecosystems. *Science*, 338, 496–500.

.. |br| raw:: html

   <br />
