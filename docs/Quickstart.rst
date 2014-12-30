Quickstart
==========

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

Installation
------------
CauseMap is an official Julia package. |br|
To install, the user can simply type ``Pkg.add("CauseMap")`` from the Julia REPL. |br|
CauseMap can then be imported with ``using CauseMap``.


Example script
--------------
CauseMap is designed for rapid analysis. After setting a few tuning parameters (see :ref:`params`),
the user is ready to analyze and plot the results.

.. literalinclude:: ../examples/CCM_example_para_didi.jl
   :language: julia


Result
++++++
This produces the plot below:

.. figure:: _images/ParaDidi_optim.jpeg
   :height: 250px

   In the left panel, we see convergence for both causal directions, suggesting a bi-directional causal effect.

|br|

Source data
~~~~~~~~~~~~

For a look at the data behind this analysis,
check out the interactive plots :download:`here <ParaDidiExample_plot.html>`.


- *The first two panels are the time series for the predator-prey system.*
- *The third and forth panels show the manifold reconstructions from each time series.*
- *The fifth panel gives a view of the full manifold.*


Note that corresponding points can be simultaneously highlighted across plots using the ``select`` tool.

|


.. |br| raw:: html

   <br />
