# set paths
OPENMEEGDATAPATH=../..
OPENMEEGBINPATH=../../../src
OPENMEEGTOOLSPATH=../../../tools

MODEL:=nerve1
CUTPLANE:=test0
# model directory (where nerve parameters, geometry, etc are stored)
MDIR:=$(OPENMEEGDATAPATH)/Models/$(MODEL)/
# directory containing the evaluation positions
PDIR:=$(OPENMEEGDATAPATH)/Computations/$(MODEL)/positions/
# directory containing the results
RDIR:=$(OPENMEEGDATAPATH)/Computations/$(MODEL)/results/
ASSEMBLE:=$(OPENMEEGBINPATH)/om_assemble
MINVERSER:=$(OPENMEEGBINPATH)/om_minverser
GAIN:=$(OPENMEEGBINPATH)/om_gain
MAKE_NERVE:=$(OPENMEEGTOOLSPATH)/om_make_nerve

$(MDIR)$(MODEL).lhs: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(ASSEMBLE) -HeadMat $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).lhs

$(MDIR)$(MODEL).lhsinv: $(MDIR)$(MODEL).lhs
	$(MINVERSER) $(MDIR)$(MODEL).lhs $(MDIR)$(MODEL).lhsinv

$(MDIR)$(MODEL).EITsource: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(ASSEMBLE) -EITsource $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).areas

$(MDIR)$(MODEL).EITstim: $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).stimelec
	$(ASSEMBLE) -EITstim $(MDIR)$(MODEL).EITsource $(MDIR)$(MODEL).stimelec $(MDIR)$(MODEL).EITstim

$(RDIR)$(MODEL).$(CUTPLANE).SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE) $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dx
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dx $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dy
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dy $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol: $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dz
	$(ASSEMBLE) -Surf2Vol $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(PDIR)$(CUTPLANE)dz $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol

$(RDIR)$(MODEL).$(CUTPLANE).gain: $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol $(RDIR)$(MODEL).$(CUTPLANE).gain

$(RDIR)$(MODEL).$(CUTPLANE)dx.gain: $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol $(MDIR)$(MODEL).lhsinv  $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim   $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dx.gain

$(RDIR)$(MODEL).$(CUTPLANE)dy.gain: $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol  $(RDIR)$(MODEL).$(CUTPLANE)dy.gain

$(RDIR)$(MODEL).$(CUTPLANE)dz.gain: $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim
	$(GAIN) -VolEEG $(MDIR)$(MODEL).lhsinv $(MDIR)$(MODEL).EITstim $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol  $(RDIR)$(MODEL).$(CUTPLANE)dz.gain

$(MDIR)$(MODEL).parameters:
	$(MAKE_NERVE) -makeparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

$(MDIR)$(MODEL).geom:
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

$(MDIR)$(MODEL).cond:
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

nerve: $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond
	$(MAKE_NERVE) -useparameters $(MDIR)$(MODEL).parameters $(MDIR)$(MODEL).geom $(MDIR)$(MODEL).cond $(MDIR)$(MODEL).int.tri $(MDIR)$(MODEL).ext.tri $(MDIR)$(MODEL).stimelec

SurfToVol: $(RDIR)$(MODEL).$(CUTPLANE).SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dx.SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dy.SurfToVol $(RDIR)$(MODEL).$(CUTPLANE)dz.SurfToVol

gain: $(RDIR)$(MODEL).$(CUTPLANE).gain $(RDIR)$(MODEL).$(CUTPLANE)dx.gain $(RDIR)$(MODEL).$(CUTPLANE)dy.gain $(RDIR)$(MODEL).$(CUTPLANE)dz.gain

all:
	make gain

clean:
	(rm $(RDIR)$(MODEL).*SurfToVol* $(MDIR)$(MODEL).*lhs* $(MDIR)$(MODEL).stimelec\
	   $(MDIR)$(MODEL).EIT*)

purge:
	(clean\
	 -@rm $(RDIR)$(MODEL).*.gain)
