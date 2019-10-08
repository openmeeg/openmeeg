import inspect
import openmeeg as om

# Make sure that doc is generated
assert inspect.getdoc(om.HeadMat) is not None

# Check docstring content
headmat_expected_docstring = "`HeadMat(const Geometry &geo, const unsigned gauss_order=3)`  \n\nConstructors\n------------\n* `HeadMat(const Geometry &geo, const unsigned gauss_order=3)`  \n\nC++ includes: assemble.h"
assert inspect.getdoc(om.HeadMat) == headmat_expected_docstring
