import inspect
import openmeeg as om


def test_doc():
    # Make sure that doc is generated
    doc = inspect.getdoc(om.HeadMat)
    assert doc is not None

    headmat_expected_docstring = """\
HeadMat(Geometry geo, Integrator const & integrator=Integrator(3,0,0.005)) \
-> SymMatrix"""
    assert (
        doc == headmat_expected_docstring
    ), f"got: {repr(doc)} != expected: {repr(headmat_expected_docstring)}"
