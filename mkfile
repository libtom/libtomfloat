</$objtype/mkfile
LIB=/$objtype/lib/ape/libtomfloat.a
OFILES=mpf_abs.$O mpf_acos.$O mpf_add.$O mpf_add_d.$O mpf_asin.$O mpf_atan.$O mpf_clear.$O \
	mpf_clear_multi.$O mpf_cmp.$O mpf_cmp_d.$O mpf_const_0.$O mpf_const_1pi.$O \
	mpf_const_1r2.$O mpf_const_2pi.$O mpf_const_2rpi.$O mpf_const_d.$O mpf_const_e.$O \
	mpf_const_l10e.$O mpf_const_l2e.$O mpf_const_le2.$O mpf_const_ln_d.$O mpf_const_pi.$O \
	mpf_const_pi2.$O mpf_const_pi4.$O mpf_const_r2.$O mpf_const_sqrt_d.$O mpf_copy.$O \
	mpf_cos.$O mpf_div.$O mpf_div_2.$O mpf_div_d.$O mpf_exch.$O mpf_exp.$O \
	mpf_init.$O mpf_init_copy.$O mpf_init_multi.$O mpf_inv.$O mpf_invsqrt.$O \
	mpf_iterations.$O mpf_ln.$O mpf_mul.$O mpf_mul_2.$O mpf_mul_d.$O mpf_neg.$O \
	mpf_normalize.$O mpf_normalize_to.$O mpf_pow.$O mpf_sin.$O mpf_sqr.$O mpf_sqrt.$O \
	mpf_sub.$O mpf_sub_d.$O mpf_tan.$O 
HFILES=tomfloat.h /sys/include/ape/tommath_class.h /sys/include/ape/tommath_superclass.h \
	/sys/include/ape/tommath_class.h /sys/include/ape/tommath_superclass.h \
	/sys/include/ape/tommath_class.h /sys/include/ape/tommath_superclass.h \
	/sys/include/ape/tommath_class.h /sys/include/ape/limits.h /sys/include/ape/ctype.h \
	/sys/include/ape/stddef.h /sys/include/ape/stdlib.h /sys/include/ape/bsd.h \
	/sys/include/ape/stddef.h /sys/include/ape/string.h /sys/include/ape/sys/types.h \
	/sys/include/ape/stddef.h /$objtype/include/ape/stdarg.h /sys/include/ape/stdio.h \
	/sys/include/ape/tommath.h 

</sys/src/cmd/mksyslib

CC=pcc
LD=pcc
CFLAGS=-c -I. -I/sys/include/ape -I/$objtype/include/ape -B -D_POSIX_SOURCE -D_SUSV2_SOURCE \
	-D_BSD_EXTENSION 
LDFLAGS=-ltommath

install:V:
	cp tomfloat.h /sys/include/ape/

nuke:V:
	mk clean
	rm -f $LIB
	rm -f /sys/include/ape/tomfloat.h
