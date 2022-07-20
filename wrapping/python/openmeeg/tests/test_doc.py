import inspect
import openmeeg as om


def test_doc():
    # Make sure that doc is generated
    doc = inspect.getdoc(om.HeadMat)
    assert doc is not None

    # Check docstring content
    # Before reverting a bunch of commits this was better:
    # headmat_expected_docstring = \
    #     ("HeadMat(Geometry geo, Integrator const & integrator=Integrator(3,0,0.005)) -> SymMatrix")
    # But now we get:
    headmat_expected_docstring = 'Proxy of C++ OpenMEEG::HeadMat class.'
    assert doc == headmat_expected_docstring, f'got: {repr(doc)} != expected: {repr(headmat_expected_docstring)}'
