# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 2013

this provides two functions om_display_vtp, and om_display_vtk.
It can be called as:
from om_display import om_display_vtp # visualiation with VTK

@author: - E. Olivi

This work has been done for the
CNRS, Laboratoire de Neurosciences Cognitives, UMR 7291, 13331, Marseille, France
under the supervision of Boris Burle

"""

import vtk

# TODO: it only works for nested geometry: find a better alternative than
# connectivity filter

# Common functions ##########################

def UpdateColorBar(colorbar,mapper):
    if mapper.GetScalarMode() == 1:
        srange = mapper.GetInput().GetPointData().GetScalars().GetRange()
    else:
        srange = mapper.GetInput().GetCellData().GetScalars().GetRange()
    tc = vtk.vtkColorTransferFunction(); tc.SetColorSpaceToDiverging()
    tc.AddRGBPoint(srange[0],0,0,1); tc.AddRGBPoint(sum(srange)/len(srange),1,1,1)
    tc.AddRGBPoint(srange[1],1,0,0); mapper.SetLookupTable( tc ); colorbar.SetLookupTable( tc );
    colorbar.SetNumberOfLabels(3); colorbar.SetLabelFormat('%3.3e');

# This callback function does plot the selected points
def PickData(object, event, selactor, state, view, text_init):
    picker = object.GetPicker()
    hsel = vtk.vtkHardwareSelector()
    ren = picker.GetRenderer(); hsel.SetRenderer(ren)
    hsel.SetArea(ren.GetPickX1(),ren.GetPickY1(),ren.GetPickX2(),ren.GetPickY2())
    if state:
        hsel.SetFieldAssociation(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS)
    else:
        hsel.SetFieldAssociation(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
    sel = hsel.Select(); ex = vtk.vtkExtractSelection()
    if picker.GetMapper():
        if (not sel.GetNumberOfNodes()==0):
            ex.SetInputConnection(picker.GetMapper().GetInputConnection(0,0))
            ex.SetSelectionConnection(sel.GetProducerPort()); ex.Update()
            selmapper = vtk.vtkDataSetMapper(); data = ex.GetOutput()
            selmapper.SetInput(data); selmapper.ScalarVisibilityOff()
            selactor.SetMapper(selmapper); selactor.PickableOff(); 
            selactor.GetProperty().SetColor(0.0, 1.0, 0.0)
            selactor.GetProperty().SetOpacity(0.5); selactor.GetProperty().SetPointSize(5)
            ren.AddActor(selactor)
            PlotSelectedData(data, state, view, text_init)

def PlotSelectedData(data, state, view, text_init):
    nb_sources = 0
    for i in range(data.GetPointData().GetNumberOfArrays()):
        if data.GetPointData().GetGlobalIds('Potentials-'+str(i)):
            nb_sources += 1
    view.GetRenderer().RemoveActor2D(text_init)
    chart = vtk.vtkChartXY()
    view.GetScene().RemoveItem(0)
    view.GetScene().AddItem(chart)
    chart.SetShowLegend(True)
    table = vtk.vtkTable()
    X = vtk.vtkDoubleArray(); X.SetName("X")
    if state:
        numPoints = data.GetNumberOfCells()
        chart.GetAxis(0).SetTitle('Current'); 
    else:
        numPoints = data.GetNumberOfPoints()
        chart.GetAxis(0).SetTitle('Potential'); 
    chart.GetAxis(0).SetRange(0,nb_sources)
    for j in range(nb_sources):
        X.InsertNextValue(j)
    table.AddColumn(X)
    for i in range(numPoints):
        Y = vtk.vtkDoubleArray()
        if state:
            Y.SetName("id"+str(data.GetCellData().GetGlobalIds('Indices').GetValue(i)))
        else:
            Y.SetName("id"+str(data.GetPointData().GetGlobalIds('Indices').GetValue(i)))
        for j in range(nb_sources):
            if state:
                Y.InsertNextValue(data.GetCellData().GetGlobalIds('Currents-'+str(j)).GetValue(i))
            else:
                Y.InsertNextValue(data.GetPointData().GetGlobalIds('Potentials-'+str(j)).GetValue(i))
        table.AddColumn(Y)
        # Now add the line plots
        line = chart.AddPlot(0)
        line.SetInput(table,0,i+1)
        line.SetColor(0,1,0)
        line.SetWidth(1.0)
    view.GetRenderWindow().Render()

##################################################

def om_display_vtp(f, n = 0):
    """
    This function displays a VTK::vtp file generated with OpenMEEG.
    Such a file defines a polydata, containing points and triangles of several
    meshes which are labelled through a vtkAbstractArray (Strings) associated to
    the cells (mesh names).
    Results of the forward problem (or a cortical mapping) can be seen thanks to
    arrays associated to points and cells (respectively potentials and normals
    currents).
    """
    welcome = """Welcome\n\n
    Switch the button: To either see Potentials (on points) or Currents (on triangles)\n
    Move the slider to see all sources (columns of the input matrix)\n
    Press 'r': To select points/cells.\n"""

    # This callback function does updates the mappers for where n is the slider value
    def CleanPickData(object, event):
        for i in range(4):
            rens[i].RemoveActor(selactor)
        if buttonWidget.GetRepresentation().GetState():
            PickData(object, event, selactor, 1, view, text_init)
        else:
            PickData(object, event, selactor, 0, view, text_init)
    def SelectSource(object, event): # object will be the slider2D
        slidervalue = int(round(object.GetRepresentation().GetValue()))
        for i in range(4):
            mappers[i].GetInput().GetPointData().SetActiveScalars("Potentials-"+str(slidervalue))
            mappers[i].GetInput().GetCellData().SetActiveScalars("Currents-"+str(slidervalue))
            renWin.SetWindowName(renWin.GetWindowName()[0:(renWin.GetWindowName().find('-')+1)]+str(slidervalue))
            UpdateColorBar(colorBars[i], mappers[i])

    # This callback function does updates the Scalar Mode To Use
    def SelectMode(object, event):
        # object will be the buttonWidget
        for i in range(4):
            if (object.GetRepresentation().GetState()):
                mappers[i].SetScalarModeToUseCellData()
                renWin.SetWindowName(renWin.GetWindowName().replace('Potentials','Currents'))
            else:
                mappers[i].SetScalarModeToUsePointData()
                renWin.SetWindowName(renWin.GetWindowName().replace('Currents','Potentials'))
            UpdateColorBar(colorBars[i], mappers[i])

    # A window with an interactor
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(600, 600)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    iren.SetInteractorStyle(vtk.vtkInteractorStyleRubberBandPick())
    # A picker (to pick points/cells)
    picker = vtk.vtkRenderedAreaPicker()
    iren.SetPicker(picker)
    # Read the input file
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(f); reader.Update()
    poly = reader.GetOutput()
    renWin.SetWindowName(f+' Potentials-'+str(n))
    # determine the number of sources
    nb_sources = 0
    for i in range(poly.GetPointData().GetNumberOfArrays()):
        if poly.GetPointData().GetGlobalIds('Potentials-'+str(i)):
            nb_sources += 1
    if n < nb_sources:
        poly.GetPointData().SetActiveScalars('Potentials-'+str(n))
        poly.GetCellData().SetActiveScalars('Currents-'+str(n))
    # Get the mesh names
    cell_labels = poly.GetCellData().GetAbstractArray(0)
    assert(cell_labels.GetName()=='Names')
    s = set(); nb_meshes = 0; cell_ids = list()
    for i in range(cell_labels.GetNumberOfValues()):
        s.add(cell_labels.GetValue(i))
        if len(s)>nb_meshes:
            # if a label is added, store the ID for the connectivity filter
            cell_ids.append(i)
            nb_meshes += 1
    # Number of meshes
    assert(nb_meshes<=4)
    # Multiple viewports: 4
    xmins = [0,.5,0,.5]; xmaxs = [0.5,1,0.5,1]; ymins = [0,0,.5,.5]; ymaxs = [0.5,0.5,1,1]

    mappers   = [vtk.vtkPolyDataMapper() for i in range(4)]
    colorBars = [vtk.vtkScalarBarActor() for i in range(4)]
    actors    = [vtk.vtkActor() for i in range(4)]
    rens      = [vtk.vtkRenderer() for i in range(4)]

    for i in range(4):
        rens[i].SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i]);
        # Display the meshes
        if (i < nb_meshes):
            # Create a connectivity filter based on cell seeded region (to display
            # only one mesh per viewport)
            conn = vtk.vtkPolyDataConnectivityFilter()
            conn.SetInput(poly)
            conn.SetExtractionModeToCellSeededRegions()
            conn.AddSeed(cell_ids[i]); conn.Update()
            actor_meshname = vtk.vtkTextActor();
            actor_meshname.SetInput(cell_labels.GetValue(cell_ids[i]));
            actor_meshname.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport();
            actor_meshname.SetPosition(0.5, 0.85); tprop = actor_meshname.GetTextProperty(); tprop.SetFontSize(30)
            tprop.SetFontFamilyToArial(); tprop.SetColor(1, 1, 1); tprop.SetJustificationToCentered()
            mappers[i].SetInputConnection(conn.GetOutputPort())
            mappers[i].SetScalarModeToUsePointData(); mappers[i].Update()
            if nb_sources:
                rens[i].AddActor2D(colorBars[i])
            actors[i].SetMapper(mappers[i])
            rens[i].AddActor2D(actor_meshname)
            rens[i].AddActor(actors[i])
            if (i == 0):
                cam = rens[i].GetActiveCamera()
                rens[i].ResetCamera()
        else:
            # Create a plane to cut
            plane = vtk.vtkPlane(); plane.SetOrigin(0,0,0); plane.SetNormal(1,0,0);
            # Create cutter
            extract = vtk.vtkExtractPolyDataGeometry(); extract.SetInput(poly)
            extract.SetImplicitFunction(plane); extract.ExtractBoundaryCellsOff()
            mappers[i].SetInputConnection(extract.GetOutputPort())
            mappers[i].SetScalarModeToUsePointData(); mappers[i].Update()
            # Create plane actor
            actors[i].SetMapper(mappers[i])
            rens[i].AddActor(actors[i])
        rens[i].SetActiveCamera(cam)
        if nb_sources:
            UpdateColorBar(colorBars[i], mappers[i])
        renWin.AddRenderer(rens[i])
        renWin.Render();

    if nb_sources > 1:
        # Slider
        sliderWidget = vtk.vtkSliderWidget()
        slider = vtk.vtkSliderRepresentation2D(); slider.SetMaximumValue(nb_sources-1)
        slider.SetValue(n); slider.SetEndCapLength(0.01); slider.SetLabelFormat('%1.0f')
        slider.SetSliderWidth(0.05); slider.SetSliderLength(1./nb_sources)
        slider.GetPoint1Coordinate().SetCoordinateSystemToNormalizedViewport()
        slider.GetPoint1Coordinate().SetValue(.0 ,0.02)
        slider.GetPoint2Coordinate().SetCoordinateSystemToNormalizedViewport()
        slider.GetPoint2Coordinate().SetValue(1. ,0.02);
        sliderWidget.SetInteractor(iren); sliderWidget.SetRepresentation(slider);
        sliderWidget.SetAnimationModeToAnimate(); sliderWidget.EnabledOn();
        sliderWidget.AddObserver("InteractionEvent", SelectSource);
    if not nb_sources == 0:
        # The button for choosing Potentials/Currents
        buttonWidget = vtk.vtkButtonWidget()
        button = vtk.vtkTexturedButtonRepresentation2D(); button.SetNumberOfStates(2)
        tex1r = vtk.vtkImageData(); tex2r = vtk.vtkImageData();
        prop  = vtk.vtkTextProperty(); prop.SetFontSize(24);
        prop.SetColor(1,0,0); prop.SetBold(2); prop.SetShadow(2); 
        str2im = vtk.vtkFreeTypeStringToImage()
        str2im.RenderString(prop,'Potentials',tex1r)
        str2im.RenderString(prop,'Currents',tex2r)
        button.SetButtonTexture(0, tex1r)
        button.SetButtonTexture(1, tex2r)
        buttonWidget.SetInteractor(iren);
        buttonWidget.SetRepresentation(button);
        button.SetPlaceFactor(1);
        button.PlaceWidget([0., 100, 50, 500, 0, 0]);
        buttonWidget.On()
        buttonWidget.AddObserver(vtk.vtkCommand.StateChangedEvent,SelectMode);
        # Selection
        selactor = vtk.vtkActor()
        view = vtk.vtkContextView(); view.GetRenderWindow().SetWindowName('Plot')
        view.GetRenderWindow().SetPosition(600, 0); view.GetRenderWindow().SetSize(600, 600)
        # Welcome text
        text_init = vtk.vtkTextActor()
        text_init.SetPosition(10, 300)
        text_init.SetInput(welcome)
        text_init.GetTextProperty().SetColor(1.0, 0.0, 0.0)
        view.GetRenderer().AddActor2D(text_init)
        view.GetInteractor().Initialize()
        iren.AddObserver(vtk.vtkCommand.EndPickEvent,CleanPickData)
    iren.Initialize()
    iren.Start()

def om_display_vtk(f,d = 0,n = 0):
    """
    This function displays a VTK::vtk file generated with OpenMEEG.
    Such a file defines a polydata, containing points and triangles of a single
    mesh. Most often a EEG helmet mesh and associated leadfield.
    """
    welcome = """Welcome\n\n
    Move the slider to see all sources (columns of the input matrix)\n
    Press 'r': To select points/cells.\n"""

    # This callback function does updates the mappers for where n is the slider value
    def CleanPickData(object, event):
        ren.RemoveActor(selactor)
        PickData(object, event, selactor, 0, view, text_init)

    def SelectSource(object, event): # object will be the slider2D
        slidervalue = int(round(object.GetRepresentation().GetValue()))
        mapper.GetInput().GetPointData().SetActiveScalars("Potentials-"+str(slidervalue))
        renWin.SetWindowName(renWin.GetWindowName()[0:(renWin.GetWindowName().find('-')+1)]+str(slidervalue))
        UpdateColorBar(colorBar, mapper)

    # A window with an interactor
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(600, 600)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    iren.SetInteractorStyle(vtk.vtkInteractorStyleRubberBandPick())
    # A picker (to pick points/cells)
    picker = vtk.vtkRenderedAreaPicker()
    iren.SetPicker(picker)
    # Read the input file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(f); reader.Update()
    poly = reader.GetOutput()
    renWin.SetWindowName(f+' Potentials-'+str(n))
    # determine the number of sources
    nb_sources = 0
    for i in range(poly.GetPointData().GetNumberOfArrays()):
        if poly.GetPointData().GetGlobalIds('Potentials-'+str(i)):
            nb_sources += 1
    if nb_sources == 0: #the file doesn't provide potentials
        if not d.__class__ == int:
            assert(d.shape[0] == poly.GetNumberOfPoints())
            nb_sources = d.shape[1]
            pot = [vtk.vtkDoubleArray() for j in range(d.shape[1])]
            for j in range(d.shape[1]):
                pot[j].SetName('Potentials-'+str(j))
                for i in range(d.shape[0]):
                    pot[j].InsertNextValue(d[i,j])
                poly.GetPointData().AddArray(pot[j])
            poly.Update()
        if not poly.GetPointData().GetGlobalIds('Indices'):
            ind = vtk.vtkUnsignedIntArray()
            ind.SetName('Indices')
            for i in range(poly.GetNumberOfPoints()):
                ind.InsertNextValue(i)
            poly.GetPointData().AddArray(ind)

    poly.GetPointData().SetActiveScalars('Potentials-'+str(n))

    mapper   = vtk.vtkPolyDataMapper()
    colorBar = vtk.vtkScalarBarActor()
    actor    = vtk.vtkActor()
    ren      = vtk.vtkRenderer()
    mapper.SetInput(poly)
    mapper.SetScalarModeToUsePointData(); mapper.Update()
    actor.SetMapper(mapper)
    ren.AddActor(actor)
    if nb_sources:
        ren.AddActor2D(colorBar)
        UpdateColorBar(colorBar, mapper)
    renWin.AddRenderer(ren)
    renWin.Render()

    if nb_sources > 1:
        # Slider
        sliderWidget = vtk.vtkSliderWidget()
        slider = vtk.vtkSliderRepresentation2D(); slider.SetMaximumValue(nb_sources-1)
        slider.SetValue(n); slider.SetEndCapLength(0.01); slider.SetLabelFormat('%1.0f')
        slider.SetSliderWidth(0.05); slider.SetSliderLength(1./nb_sources)
        slider.GetPoint1Coordinate().SetCoordinateSystemToNormalizedViewport()
        slider.GetPoint1Coordinate().SetValue(.0, 0.02)
        slider.GetPoint2Coordinate().SetCoordinateSystemToNormalizedViewport()
        slider.GetPoint2Coordinate().SetValue(1., 0.02);
        sliderWidget.SetInteractor(iren); sliderWidget.SetRepresentation(slider);
        sliderWidget.SetAnimationModeToAnimate(); sliderWidget.EnabledOn();
        sliderWidget.AddObserver("InteractionEvent", SelectSource);
        # Selection
        selactor = vtk.vtkActor()
        view = vtk.vtkContextView(); view.GetRenderWindow().SetWindowName('Plot')
        view.GetRenderWindow().SetPosition(600, 0); view.GetRenderWindow().SetSize(600, 600)
        # Welcome text
        text_init = vtk.vtkTextActor()
        text_init.SetPosition(10, 300)
        text_init.SetInput(welcome)
        text_init.GetTextProperty().SetColor(1.0, 0.0, 0.0)
        view.GetRenderer().AddActor2D(text_init)
        view.GetInteractor().Initialize()
        iren.AddObserver(vtk.vtkCommand.EndPickEvent,CleanPickData)
    iren.Initialize()
    iren.Start()
