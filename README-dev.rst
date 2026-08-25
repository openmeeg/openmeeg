OpenMEEG developer notes
========================

.. highlight:: console

Notes for OpenMEEG maintainers. See `README.rst <README.rst>`_ for user-facing
build and installation instructions.

Versioning
----------

The version is **derived from git tags**; no file in the repository records it.
There is therefore nothing to bump before or after a release:

- The C++ side gets it in ``cmake/Version.cmake`` (included before ``project()``),
  which runs ``git describe`` and sets ``OPENMEEG_VERSION`` (the full string, used
  for ``OpenMEEGConfigure.h`` and CPack) and ``OPENMEEG_VERSION_NUMERIC`` (the
  ``X.Y.Z`` part, which is all ``project(VERSION ...)`` accepts).
- The Python side gets it from ``setuptools_scm`` (configured in
  ``wrapping/python/pyproject.toml``), which is what ``openmeeg.__version__``
  reports.

Both use the same rules, so the two versions always agree:

- On a tagged commit, the version is exactly the tag, e.g. ``2.5.16``.
- On any later commit, it is the *next* patch release plus the number of commits
  since the tag, e.g. ``2.5.17.dev68`` 68 commits after ``2.5.16``.

Two consequences worth knowing:

- Builds need the tags: clone normally (not ``--depth``), and CI checkouts must
  use ``fetch-depth: 0``. GitHub source archives work, too, because
  ``.git_archival.txt`` is filled in by ``export-subst`` on export.
- If neither is available (e.g. a tarball made by hand), configure with
  ``-DOPENMEEG_VERSION=X.Y.Z`` to set it explicitly.

Cutting a new release
---------------------

Releasing is meant to be done using the UI on GitHub, which can also handle the
tagging step simultaneously (there is an option for it in the UI). 
Because the version follows the tags, there is no version-bump PR before the
release and no restore-the-dev-version PR after it.

1. Check that ``main`` is in a releasable state: everything intended for the
   release is merged, and the ``build_and_test``, ``cibuildwheel`` and
   ``cibuildwheel-apps`` workflows are green on the commit you want to ship.

2. Visit https://github.com/openmeeg/openmeeg/releases/new to cut a new release,
   creating a new tag with the version number and using the "Generate release notes"
   option (tweaking the result as needed). These can be edited after the fact.
   Click the green button "Publish release" to push the tag and trigger the release
   workflows.

3. Wait for the release workflows to finish. Wheels will automatically be uploaded
   to PyPI. Wheels and app bundles will be attached to the GitHub release, which will
   be put back in draft mode. Once they're all attached, click "Publish release" to
   make it public. 

4. The `conda-forge feedstock <https://github.com/conda-forge/openmeeg-feedstock>`_
   is updated by the regro autotick bot, which opens a PR a few hours after the
   PyPI/GitHub release appears. Review and merge it with any necessary dependency
   or build changes.
