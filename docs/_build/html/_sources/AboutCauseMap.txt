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

.. _intro:
Introduction
-------------
Convergent cross mapping (CCM) provides a model-free approach to detecting dependencies
and establishing causality in complex non-linear systems, even in the presence of feedback loops and
unmeasured confounding [1]_. CCM derives this power from
explicitly capturing time-dependent dynamics through a technique known as
state-space reconstruction (SSR). In practice, this analysis
typically requires 25 or more time points, measured with sufficient density
to capture system dynamics.


How CauseMap works
------------------
To get a feel for the algorithm and the intuition behind it, we highly recommend
the following two videos from George Sugihara's group at UCSD.
|

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
the latest version of our `manuscript <https://github.com/cyrusmaher/CauseMap/CauseMap_latest.docx>`_.

|


.. [1] Sugihara,G. *et al.* (2012) Detecting causality in complex ecosystems. *Science*, 338, 496–500.
