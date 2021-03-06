============
Contributing
============

Contributions to Rhodonite are very welcome, and always appreciated!

Rhodonite is a part of Nesta's work towards an `Open Science of Science 
<https://www.nesta.org.uk/blog/towards-an-open-science-science/>`_.

Every type of contribution makes a huge addition to the project. Tools like
Rhodonite can only improve and be sustained by the community, and we will always
give credit for any help recieved.

Bug reports
===========

When `reporting a bug <https://github.com/nestauk/rhodonite/issues>`_ please include:

    * Your operating system name and version.
    * Any details about your local setup that might be helpful in troubleshooting.
    * Detailed steps to reproduce the bug.

Documentation improvements
==========================

Rhodonite's documentation could always benefit from further additions, 
improvements and clarifications, whether as part of the official Rhodonite docs,
in docstrings, or even on the web in blog posts, articles, and such.

Documentation contributions are among the most important, so never feel like 
this is an unworthy task!

Feature requests and feedback
=============================

Rhodonite is a toolbox of network based algorithms for analysing and measuring activity and outputs
in science, innovation and research. We are always interested in incorporating new methods, so if you
have been involved in developing one (or several), please get in touch!

The best way to send feedback is to file an issue at https://github.com/nestauk/rhodonite/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

Development
===========

To set up `rhodonite` for local development:

1. Fork `rhodonite <https://github.com/nestauk/rhodonite>`_
   (look for the "Fork" button).
2. Clone your fork locally::

    git clone git@github.com:your_name_here/rhodonite.git

3. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

4. When you're done making changes, run all the checks, doc builder and spell checker with `tox <http://tox.readthedocs.io/en/latest/install.html>`_ one command::

    tox

5. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

6. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

If you need some code review or feedback while you're developing the code just make the pull request.

For merging, you should:

1. Include passing tests (run ``tox``) [1]_.
2. Update documentation when there's new API, functionality etc.
3. Add a note to ``CHANGELOG.rst`` about the changes.
4. Add yourself to ``AUTHORS.rst``.

.. [1] If you don't have all the necessary python versions available locally you can rely on Travis - it will
       `run the tests <https://travis-ci.org/georgerichardson/rhodonite/pull_requests>`_ for each change you add in the pull request.

       It will be slower though ...

.. Tips
.. ----
.. 
.. To run a subset of tests::
.. 
..     tox -e envname -- pytest -k test_myfeature
.. 
.. To run all the test environments in *parallel* (you need to ``pip install detox``)::
.. 
..     detox
