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
   :height: 350px

Interpretation
~~~~~~~~~~~~
In the left panel, we see convergence for both causal directions. The interpretation of this result is that the causal relationship between P. aurelia and D. nasutum is bi-directional. That is, the number of predators influences the number of prey, and vice-versa. 

The right panels show the dependence of the max ρ\ :sub:`ccm` on E (the dimensionality of the reconstructed system), and the supposed time lag of the causal effect (τ\ :sub:`p`). Overall, the patterning of these heatmaps demonstrates that max ρ\ :sub:`ccm` has a reasonable dependence on causal parameters of interest. While the max ρccm is relatively insensitive to the assumed dimensionality, the best-performing τp values correspond to either immediate causal effects, or those delayed by five days. 

Note that τ\ :sub:`p`\ =5 corresponds to the principal frequency of the Paramecium aurelia and Didinium nasutum time series. This suggests that the peak at τ\ :sub:`p`\ =5 is artifactual. Therefore, we are able to infer from these results that, as we would expect, predator and prey populations exert bidirectional effects in real-time. 


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
