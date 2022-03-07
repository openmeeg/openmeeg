import inspect
import openmeeg as om

# Make sure that doc is generated
assert inspect.getdoc(om.HeadMat) is not None

# Check docstring content
headmat_expected_docstring = \
    ("`HeadMat(const Geometry &geo, const Integrator &integrator=Integrator(3, 0,"
     "\n    0.005)) -> SymMatrix`  ")
assert inspect.getdoc(om.HeadMat) == headmat_expected_docstring
