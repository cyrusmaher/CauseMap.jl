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
To install, the user can simply type ``Pkg.add("CauseMap")`` from the Julia REPL.
CauseMap can then be imported with ``using CauseMap``.


Example script
--------------
CauseMap is designed for rapid analysis. After setting a few tuning parameters,
the user is ready to analyze and plot the results.

.. literalinclude:: ../examples/CCM_example_para_didi.jl
   :language: julia
   :lines: 1-10
The tuning parameters above are described in greater depth elsewhere. Now we were ready to get results!

.. literalinclude:: ../examples/CCM_example_para_didi.jl
   :language: julia
   :lines: 17-20

Result
++++++
This produces the plot below:

.. figure:: ../images/ParaDidi_optim.jpeg
   :height: 250px

   In the left panel, we see convergence for both causal directions, suggesting a bi-directional causal effect.



For a look at the data behind this analysis,
check out the interactive plots :doc:`here <ParaDidiExample>`.

|


.. |br| raw:: html

   <br />
