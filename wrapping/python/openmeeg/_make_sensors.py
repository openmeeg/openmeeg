# Build a Sensors object from plain arrays.
import numpy as np

from ._openmeeg_wrapper import Matrix, Sensors, Vector


def make_sensors(
    positions,
    orientations=None,
    *,
    labels=None,
    weights=None,
    radii=None,
    geometry=None,
):
    """Build a :class:`Sensors` from arrays.

    This wraps the (many) :class:`Sensors` constructors with a single, forgiving
    entry point: it accepts arrays in any memory layout, fills in sensible
    defaults, and validates shapes (raising instead of crashing).

    MEG (magnetometers/gradiometers), from positions and orientations::

        make_sensors(positions, orientations)

    EEG electrodes, projected onto a head ``geometry``::

        make_sensors(positions, geometry=geom)

    Parameters
    ----------
    positions : array-like, shape (n_sensors, 3)
        Sensor (integration point) positions.
    orientations : array-like, shape (n_sensors, 3) | None
        Sensor orientations. Required for MEG sensors. Leave as None for EEG
        electrodes, which instead require ``geometry``.
    labels : list of str | None
        One name per sensor. Defaults to ``["0", "1", ...]``.
    weights : array-like, shape (n_sensors,) | None
        Integration weights (MEG only). Defaults to ones.
    radii : array-like, shape (n_sensors,) | None
        Integration radii (MEG only). Defaults to zeros (point sensors).
    geometry : Geometry | None
        The head geometry. Required for EEG electrodes; optional for MEG.

    Returns
    -------
    sensors : Sensors
        The sensors object, usable in e.g. :func:`Head2MEGMat`,
        :func:`Head2EEGMat`, :func:`DipSource2MEGMat`.
    """
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(
            f"positions must have shape (n_sensors, 3), got {positions.shape}"
        )
    n_sensors = positions.shape[0]
    positions = Matrix(np.asfortranarray(positions))

    if orientations is None:
        if geometry is None:
            raise ValueError(
                "orientations=None builds EEG electrodes, which require a "
                "geometry; pass geometry=... (or pass orientations for MEG)."
            )
        return Sensors(positions, geometry)

    orientations = np.asarray(orientations, dtype=np.float64)
    if orientations.shape != (n_sensors, 3):
        raise ValueError(
            f"orientations must have the same shape as positions "
            f"({n_sensors}, 3), got {orientations.shape}"
        )
    orientations = Matrix(np.asfortranarray(orientations))

    if labels is None:
        labels = [str(idx) for idx in range(n_sensors)]
    else:
        labels = [str(label) for label in labels]
        if len(labels) != n_sensors:
            raise ValueError(f"labels must have length {n_sensors}, got {len(labels)}")

    def _to_vector(values, name, default):
        values = default if values is None else values
        values = np.asarray(values, dtype=np.float64).ravel()
        if values.shape != (n_sensors,):
            raise ValueError(
                f"{name} must have shape ({n_sensors},), got {values.shape}"
            )
        return Vector(np.ascontiguousarray(values))

    weights = _to_vector(weights, "weights", np.ones(n_sensors))
    radii = _to_vector(radii, "radii", np.zeros(n_sensors))

    args = [labels, positions, orientations, weights, radii]
    if geometry is not None:
        args.append(geometry)
    return Sensors(*args)
