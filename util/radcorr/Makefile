FC := gfortran
FFLAGS ?= -O -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8

test_get_radcorr_3he: ../../src/get_radcorr_3he.f test_get_radcorr_3he.f
	$(FC) $(FFLAGS) ../../src/get_radcorr_3he.f test_get_radcorr_3he.f -o $@

check_born_consistency: ../../src/get_radcorr_3he.f \
	../../src/getsf_f1f2fit.f check_born_consistency.f
	$(FC) $(FFLAGS) ../../src/get_radcorr_3he.f \
	../../src/getsf_f1f2fit.f check_born_consistency.f -o $@

clean:
	rm -f test_get_radcorr_3he check_born_consistency
