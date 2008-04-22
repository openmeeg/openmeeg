MODEL:=nerve1
CUTPLANE:=cutz0
# root of openmeeg build directory 
OPENMEEGBUILDDIR:=../../../../
# model directory (where nerve parameters, geometry, etc are stored)
MDIR:=$(OPENMEEGBUILDDIR)HeadModels/$(MODEL)/
# directory containing the evaluation positions
PDIR:=$(OPENMEEGBUILDDIR)Computations/$(MODEL)/positions/
# directory containing the results
RDIR:=$(OPENMEEGBUILDDIR)Computations/$(MODEL)/results/
ASSEMBLE=$(OPENMEEGBUILDDIR)src/om_assemble
MINVERSER=$(OPENMEEGBUILDDIR)src/om_minverser
GAIN=$(OPENMEEGBUILDDIR)src/om_gain
MAKE_NERVE=$(OPENMEEGBUILDDIR)tools/om_make_nerve 

$(MDIR)$(MODEL).LHS: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(ASSEMBLE) -LHS $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).LHS

$(MDIR)$(MODEL).LHSinv: $(MDIR)$(MODEL).LHS
	$(MINVERSER) $(MDIR)$(MODEL).LHS $(MDIR)$(MODEL).LHSinv

$(MDIR)$(MODEL).EITsource: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(ASSEMBLE) -EITsource $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).areas

$(MDIR)$(MODEL).EITstim: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).stimelec
	($(ASSEMBLE) -EITstim $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond\
	  $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).stimelec $(MDIR)$(MODEL).EITstim)

$(RDIR)$(MODEL).$(CUTPLANE).SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)
	$(ASSEMBLE) -SurfToVol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE) $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dx
	$(ASSEMBLE) -SurfToVol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dx $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dy
	$(ASSEMBLE) -SurfToVol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dy $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dz
	$(ASSEMBLE) -SurfToVol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dz $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE).gain: $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolPotEIT $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim $(RDIR)$(MODEL).$(CUTPLANE).gain

$(RDIR)$(MODEL).$(CUTPLANE)dx.gain: $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol $(MDIR)$(MODEL).LHSinv  $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolPotEIT $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol $(MDIR)$(MODEL).LHSinv  $(MDIR)$(MODEL).EITstim  $(RDIR)$(MODEL).$(CUTPLANE)dx.gain

$(RDIR)$(MODEL).$(CUTPLANE)dy.gain: $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolPotEIT $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol $(MDIR)$(MODEL).LHSinv  $(MDIR)$(MODEL).EITstim  $(RDIR)$(MODEL).$(CUTPLANE)dy.gain

$(RDIR)$(MODEL).$(CUTPLANE)dz.gain: $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolPotEIT $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim  $(RDIR)$(MODEL).$(CUTPLANE)dz.gain

$(MDIR)$(MODEL).parameters:
	$(MAKE_NERVE) -makeparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

nerve:
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

SurfToVol: $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol

gain: $(RDIR)$(MODEL).$(CUTPLANE).gain $(RDIR)$(MODEL).$(CUTPLANE)dx.gain $(RDIR)$(MODEL).$(CUTPLANE)dy.gain $(RDIR)$(MODEL).$(CUTPLANE)dz.gain

all:
	(make nerve;\
	make gain)

clean:
	(-@rm $(RDIR)$(MODEL).*SurfToVol* $(MDIR)$(MODEL).*LHS* $(MDIR)$(MODEL).stimelec\
	   $(MDIR)$(MODEL).EIT*)

purge:
	(clean\
	 -@rm $(RDIR)$(MODEL).*.gain)
