#include "stella_cedar.h"
#ifdef WITH_CEDAR


void cedar_copyto(const double *src, cedar_vec dst)
{
	cedar_real *dst_arr;
	cedar_vec_baseptr(dst, &dst_arr);
	int nd;
	cedar_vec_getdim(dst, &nd);
	if (nd == 2) {
		cedar_len ilen, jlen;
		cedar_vec_len2d(dst, &ilen, &jlen);

		cedar_len curr = 0;
		for (cedar_len j = 1; j < jlen - 1; j++) {
			for (cedar_len i = 1; i < ilen - 1; i++) {
				dst_arr[j*ilen + i] = src[curr];
				curr++;
			}
		}
	} else {
		cedar_len ilen, jlen, klen;
		cedar_vec_len3d(dst, &ilen, &jlen, &klen);

		cedar_len curr = 0;
		for (cedar_len k = 1; k < klen - 1; k++) {
			for (cedar_len j = 1; j < jlen - 1; j++) {
				for (cedar_len i = 1; i < ilen - 1; i++) {
					dst_arr[k*(ilen*jlen) + j*ilen + i] = src[curr];
					curr++;
				}
			}
		}
	}
}


void cedar_copyfrom(cedar_vec src, double *dst)
{
	cedar_real *src_arr;
	cedar_vec_baseptr(src, &src_arr);
	int nd;
	cedar_vec_getdim(src, &nd);
	if (nd == 2) {
		cedar_len ilen, jlen;
		cedar_vec_len2d(src, &ilen, &jlen);

		cedar_len curr = 0;
		for (cedar_len j = 1; j < jlen - 1; j++) {
			for (cedar_len i = 1; i < ilen - 1; i++) {
				dst[curr] = src_arr[j*ilen + i];
				curr++;
			}
		}
	} else {
		cedar_len ilen, jlen, klen;
		cedar_vec_len3d(src, &ilen, &jlen, &klen);

		cedar_len curr = 0;
		for (cedar_len k = 1; k < klen - 1; k++) {
			for (cedar_len j = 1; j < jlen - 1; j++) {
				for (cedar_len i = 1; i < ilen - 1; i++) {
					dst[curr] = src_arr[k*(ilen*jlen) + j*ilen + i];
					curr++;
				}
			}
		}
	}
}

#endif
