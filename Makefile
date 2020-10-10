FC = gfortran

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized

.PHONY: all libs clean

all: libs lbfgsb_77_1 lbfgsb_77_2 lbfgsb_77_3 lbfgsb_90_1 lbfgsb_90_2 lbfgsb_90_3

lbfgsb_77_1: driver_77_1.o liblbfgsb.a
	$(FC) $(FFLAGS) driver_77_1.o -L. -llbfgsb -o x.lbfgsb_77_1

lbfgsb_77_2: driver_77_2.o liblbfgsb.a
	$(FC) $(FFLAGS) driver_77_2.o -L. -llbfgsb -o x.lbfgsb_77_2

lbfgsb_77_3: driver_77_3.o liblbfgsb.a
	$(FC) $(FFLAGS) driver_77_3.o -L. -llbfgsb -o x.lbfgsb_77_3

lbfgsb_90_1: driver_90_1.o liblbfgsb.a
	$(FC) $(FFLAGS) driver_90_1.o -L. -llbfgsb -o x.lbfgsb_90_1

lbfgsb_90_2: driver_90_2.o liblbfgsb.a
	$(FC) $(FFLAGS) driver_90_2.o -L. -llbfgsb -o x.lbfgsb_90_2

lbfgsb_90_3: driver_90_3.o liblbfgsb.a
	$(FC) $(FFLAGS) driver_90_3.o -L. -llbfgsb -o x.lbfgsb_90_3

driver_77_1.o: driver1.f
	$(FC) $(FFLAGS) driver1.f   -c -o driver_77_1.o

driver_77_2.o: driver2.f
	$(FC) $(FFLAGS) driver2.f   -c -o driver_77_2.o

driver_77_3.o: driver3.f
	$(FC) $(FFLAGS) driver3.f   -c -o driver_77_3.o

driver_90_1.o: driver1.f90
	$(FC) $(FFLAGS) driver1.f90 -c -o driver_90_1.o

driver_90_2.o: driver2.f90
	$(FC) $(FFLAGS) driver2.f90 -c -o driver_90_2.o

driver_90_3.o: driver3.f90
	$(FC) $(FFLAGS) driver3.f90 -c -o driver_90_3.o

## L-BFGS-B as static library

libs: liblbfgsb.a

liblbfgsb.a: setulb.o mainlb.o timer.o errclb.o prn3lb.o prn1lb.o active.o dcopy.o projgr.o cauchy.o freev.o formk.o cmprlb.o subsm.o lnsrlb.o prn2lb.o ddot.o dscal.o matupd.o formt.o bmv.o hpsolb.o daxpy.o dpofa.o dtrsl.o dcsrch.o dcstep.o
	ar rcs liblbfgsb.a setulb.o mainlb.o timer.o errclb.o prn3lb.o prn1lb.o active.o dcopy.o projgr.o cauchy.o freev.o formk.o cmprlb.o subsm.o lnsrlb.o prn2lb.o ddot.o dscal.o matupd.o formt.o bmv.o hpsolb.o daxpy.o dpofa.o dtrsl.o dcsrch.o dcstep.o

setulb.o: setulb.f
	$(FC) $(FFLAGS) setulb.f -c

mainlb.o: mainlb.f
	$(FC) $(FFLAGS) mainlb.f -c

timer.o: timer.f
	$(FC) $(FFLAGS) timer.f -c

errclb.o: errclb.f
	$(FC) $(FFLAGS) errclb.f -c

prn3lb.o: prn3lb.f
	$(FC) $(FFLAGS) prn3lb.f -c

prn1lb.o: prn1lb.f
	$(FC) $(FFLAGS) prn1lb.f -c

active.o: active.f
	$(FC) $(FFLAGS) active.f -c

dcopy.o: dcopy.f
	$(FC) $(FFLAGS) dcopy.f -c

projgr.o: projgr.f
	$(FC) $(FFLAGS) projgr.f -c

cauchy.o: cauchy.f
	$(FC) $(FFLAGS) cauchy.f -c

freev.o: freev.f
	$(FC) $(FFLAGS) freev.f -c

formk.o: formk.f
	$(FC) $(FFLAGS) formk.f -c

cmprlb.o: cmprlb.f
	$(FC) $(FFLAGS) cmprlb.f -c

subsm.o: subsm.f
	$(FC) $(FFLAGS) subsm.f -c

lnsrlb.o: lnsrlb.f
	$(FC) $(FFLAGS) lnsrlb.f -c

prn2lb.o: prn2lb.f
	$(FC) $(FFLAGS) prn2lb.f -c

ddot.o: ddot.f
	$(FC) $(FFLAGS) ddot.f -c

dscal.o: dscal.f
	$(FC) $(FFLAGS) dscal.f -c

matupd.o: matupd.f
	$(FC) $(FFLAGS) matupd.f -c

formt.o: formt.f
	$(FC) $(FFLAGS) formt.f -c

bmv.o: bmv.f
	$(FC) $(FFLAGS) bmv.f -c

hpsolb.o: hpsolb.f
	$(FC) $(FFLAGS) hpsolb.f -c

daxpy.o: daxpy.f
	$(FC) $(FFLAGS) daxpy.f -c

dpofa.o: dpofa.f
	$(FC) $(FFLAGS) dpofa.f -c

dtrsl.o: dtrsl.f
	$(FC) $(FFLAGS) dtrsl.f -c

dcsrch.o: dcsrch.f
	$(FC) $(FFLAGS) dcsrch.f -c

dcstep.o: dcstep.f
	$(FC) $(FFLAGS) dcstep.f -c



clean:
	rm -rf *.o *.a x.* docs/

