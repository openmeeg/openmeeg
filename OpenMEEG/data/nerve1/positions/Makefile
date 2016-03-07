MODEL:=nerve1
CUTPLANE:=cutz0
# root of openmeeg build directory 
OPENMEEGBUILDDIR:=../../../../
# model directory (where nerve parameters, geometry, etc are stored)
MDIR:=$(OPENMEEGBUILDDIR)Models/$(MODEL)/
# directory containing the evaluation positions
PDIR:=$(OPENMEEGBUILDDIR)Computations/$(MODEL)/positions/
# directory containing the results
RDIR:=$(OPENMEEGBUILDDIR)Computations/$(MODEL)/results/
ASSEMBLE=$(OPENMEEGBUILDDIR)src/om_assemble
MINVERSER=$(OPENMEEGBUILDDIR)src/om_minverser
GAIN=$(OPENMEEGBUILDDIR)src/om_gain
MAKE_NERVE=$(OPENMEEGBUILDDIR)tools/om_make_nerve 

$(MDIR)$(MODEL).HM: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(ASSEMBLE) -HeadMat $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).HM

$(MDIR)$(MODEL).HMinv: $(MDIR)$(MODEL).HM
	$(MINVERSER) $(MDIR)$(MODEL).HM $(MDIR)$(MODEL).HMinv

$(MDIR)$(MODEL).EITsource: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(ASSEMBLE) -EITsource $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).areas

$(MDIR)$(MODEL).EITstim: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).stimelec
	($(ASSEMBLE) -EITstim $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond\
	  $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).stimelec $(MDIR)$(MODEL).EITstim)

$(RDIR)$(MODEL).$(CUTPLANE).Surf2Vol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE) $(RDIR)$(MODEL).$(CUTPLANE).Surf2Vol

$(RDIR)$(MODEL).$(CUTPLANE)dx.Surf2Vol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dx
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dx $(RDIR)$(MODEL).$(CUTPLANE)dx.Surf2Vol

$(RDIR)$(MODEL).$(CUTPLANE)dy.Surf2Vol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dy
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dy $(RDIR)$(MODEL).$(CUTPLANE)dy.Surf2Vol

$(RDIR)$(MODEL).$(CUTPLANE)dz.Surf2Vol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dz
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dz $(RDIR)$(MODEL).$(CUTPLANE)dz.Surf2Vol

$(RDIR)$(MODEL).$(CUTPLANE).gain: $(RDIR)$(MODEL).$(CUTPLANE).Surf2Vol $(MDIR)$(MODEL).HMinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).HMinv $(MDIR)$(MODEL).EITstim  $(RDIR)$(MODEL).$(CUTPLANE).Surf2Vol $(RDIR)$(MODEL).$(CUTPLANE).gain

$(RDIR)$(MODEL).$(CUTPLANE)dx.gain: $(RDIR)$(MODEL).$(CUTPLANE)dx.Surf2Vol $(MDIR)$(MODEL).LHSinv  $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).HMinv  $(MDIR)$(MODEL).EITstim  $(RDIR)$(MODEL).$(CUTPLANE)dx.Surf2Vol $(RDIR)$(MODEL).$(CUTPLANE)dx.gain

$(RDIR)$(MODEL).$(CUTPLANE)dy.gain: $(RDIR)$(MODEL).$(CUTPLANE)dy.Surf2Vol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).HMinv  $(MDIR)$(MODEL).EITstim $(RDIR)$(MODEL).$(CUTPLANE)dy.Surf2Vol  $(RDIR)$(MODEL).$(CUTPLANE)dy.gain

$(RDIR)$(MODEL).$(CUTPLANE)dz.gain: $(RDIR)$(MODEL).$(CUTPLANE)dz.Surf2Vol $(MDIR)$(MODEL).LHSinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).HMinv $(MDIR)$(MODEL).EITstim  $(RDIR)$(MODEL).$(CUTPLANE)dz.Surf2Vol $(RDIR)$(MODEL).$(CUTPLANE)dz.gain

$(MDIR)$(MODEL).parameters:
	$(MAKE_NERVE) -makeparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

$(MDIR)$(MODEL).geom:
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

$(MDIR)$(MODEL).cond:
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

nerve:
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

Surf2Vol: $(RDIR)$(MODEL).$(CUTPLANE).Surf2Vol $(RDIR)$(MODEL).$(CUTPLANE)dx.Surf2Vol $(RDIR)$(MODEL).$(CUTPLANE)dy.Surf2Vol $(RDIR)$(MODEL).$(CUTPLANE)dz.Surf2Vol

gain: $(RDIR)$(MODEL).$(CUTPLANE).gain $(RDIR)$(MODEL).$(CUTPLANE)dx.gain $(RDIR)$(MODEL).$(CUTPLANE)dy.gain $(RDIR)$(MODEL).$(CUTPLANE)dz.gain

all:
	make gain

clean:
	(-@rm $(RDIR)$(MODEL).*Surf2Vol* $(MDIR)$(MODEL).*HM* $(MDIR)$(MODEL).stimelec\
	   $(MDIR)$(MODEL).EIT*)

purge:
	(clean\
	 -@rm $(RDIR)$(MODEL).*.gain)
