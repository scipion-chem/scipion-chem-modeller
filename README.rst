=======================
Modeller plugin
=======================

**Documentation under development, sorry for the inconvenience**

This is a **Scipion** plugin that offers different **Modeller tools**.
These tools will make it possible to carry out different functions for modelling protein structures
(e.g: introducing mutations to structures)

Therefore, this plugin allows to use programs from the Modeller software suite
within the Scipion framework.

==========================
Install this plugin
==========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.


1. **Install the plugin in Scipion**

**The installation is automatic, but you need to register into the
modeller website (https://salilab.org/modeller/registration.html) in order to obtain a license key,
which must be edited in the modeller/modlib/modeller/config.py file**.

Modeller is installed as a python module, however the module is only callable by using the bash file
modeller/bin/modpy.sh. For more information: https://salilab.org/modeller/download_installation.html

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or

.. code-block::

    scipion3 installp -p scipion-chem-modeller


- **Developer's version**

    1. **Download repository**:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-modeller.git

    2. **Switch to the desired branch** (master or devel):

    Scipion-chem-modeller is constantly under development and including new features.
    If you want a relatively older an more stable version, use master branch (default).
    If you want the latest changes and developments, user devel branch.

    .. code-block::

            cd scipion-chem-modeller
            git checkout devel

    3. **Install**:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-modeller --devel


