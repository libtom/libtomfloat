#GCC makefile for LibTomFloat
#
#Tom St Denis

default: libtomfloat.a

CFLAGS += -O2 -g3 -Wall -W -I./

VERSION=0.02

#default files to install
LIBNAME=libtomfloat.a
HEADERS=tomfloat.h

#LIBPATH-The directory for libtomfloat to be installed to.
#INCPATH-The directory to install the header files for libtomfloat.
#DATAPATH-The directory to install the pdf docs.
DESTDIR=
LIBPATH=/usr/lib
INCPATH=/usr/include
DATAPATH=/usr/share/doc/libtomfloat/pdf


OBJECTS = \
mpf_init.o mpf_clear.o mpf_init_multi.o mpf_clear_multi.o mpf_init_copy.o \
\
mpf_copy.o mpf_exch.o mpf_abs.o mpf_neg.o \
\
mpf_cmp.o mpf_cmp_d.o \
\
mpf_normalize.o mpf_normalize_to.o mpf_iterations.o \
mpf_normalize_to_multi.o \
mpf_const_0.o    mpf_const_1r2.o  mpf_const_2rpi.o  mpf_const_e.o     \
mpf_const_l2e.o  mpf_const_pi.o   mpf_const_pi4.o   mpf_const_1pi.o   \
mpf_const_2pi.o  mpf_const_d.o    mpf_const_l10e.o  mpf_const_le2.o   \
mpf_const_le10.o\
mpf_const_pi2.o  mpf_const_r2.o   mpf_const_ln_d.o                    \
mpf_const_gamma.o\
\
mpf_mul_2.o mpf_div_2.o mpf_add.o mpf_sub.o mpf_mul.o mpf_sqr.o mpf_div.o \
mpf_add_d.o mpf_sub_d.o mpf_mul_d.o mpf_div_d.o \
\
mpf_invsqrt.o mpf_inv.o mpf_exp.o mpf_sqrt.o mpf_pow.o mpf_ln.o \
\
mpf_cos.o mpf_sin.o mpf_tan.o mpf_acos.o mpf_asin.o mpf_atan.o \
mpf_pow_d.o\
mpf_set_str.o mpf_get_str.o\
mpf_const_nan.o mpf_const_inf.o\
mpf_from_mp_int.o\
mpf_global_variables.o\
mpf_getdecimalexponent.o  mpf_getdecimalprecision.o\
mpf_set_double.o mpf_get_double.o\
mpf_ldexp.o mpf_frexp.o\
mpf_inv_old.o\
mpf_digits.o\
mpf_nthroot.o\
mpf_agm.o\
mpf_sincos.o\
mpf_floor.o mpf_ceil.o mpf_round.o\
mpf_trig_arg_reduct.o\
mpf_sincos.o\
mpf_const_eps.o\
mpf_dump.o\
mpf_sinh.o mpf_cosh.o mpf_tanh.o\
mpf_atanh.o mpf_kernel_atan.o\
mpf_atan2.o\
mpf_asinh.o mpf_acosh.o\
mp_acoth_binary_splitting.o mp_acot_binary_splitting.o\
mpf_lambertw.o



#mpf_trig_arg_reduct.o




libtomfloat.a: $(OBJECTS)
	$(AR) $(ARFLAGS) libtomfloat.a $(OBJECTS)
	ranlib libtomfloat.a

ex1: libtomfloat.a demos/ex1.o
	$(CC) demos/ex1.o libtomfloat.a -ltommath -lm -o ex1

#LTF user manual
mandvi: float.tex
	echo "hello" > float.ind
	latex float > /dev/null
	latex float > /dev/null
	makeindex float
	latex float > /dev/null

#LTF user manual [pdf]
manual:	mandvi
	pdflatex float >/dev/null
	rm -f float.aux float.dvi float.log float.idx float.lof float.out float.toc

install: libtomfloat.a
	install -d -g root -o root $(DESTDIR)$(LIBPATH)
	install -d -g root -o root $(DESTDIR)$(INCPATH)
	install -g root -o root $(LIBNAME) $(DESTDIR)$(LIBPATH)
	install -g root -o root $(HEADERS) $(DESTDIR)$(INCPATH)

clean:
	rm -f $(OBJECTS) libtomfloat.a *~ demos/*.o demos/*~ ex1
	rm -f float.aux float.dvi float.log float.idx float.lof float.out float.toc float.ilg float.ind float.pdf

zipup: clean manual
	cd .. ; rm -rf ltf* libtomfloat-$(VERSION) ; mkdir libtomfloat-$(VERSION) ; \
	cp -R ./libtomfloat/* ./libtomfloat-$(VERSION)/ ; \
	tar -c libtomfloat-$(VERSION)/* | bzip2 -9vvc > ltf-$(VERSION).tar.bz2 ; \
	zip -9 -r ltf-$(VERSION).zip libtomfloat-$(VERSION)/*
