/* spm_datatypes.h */

/* constants etc defining analyze / spm image data types */

#define SPM_UNSIGNED_CHAR     2
#define SPM_SIGNED_SHORT      4
#define SPM_SIGNED_INT        8
#define SPM_FLOAT             16
#define SPM_DOUBLE            64
#define SPM_SIGNED_CHAR       (SPM_UNSIGNED_CHAR+128) 
#define SPM_UNSIGNED_SHORT    (SPM_SIGNED_SHORT+128) 
#define SPM_UNSIGNED_INT      (SPM_SIGNED_INT+128)

/* byte swapped types */
#define SPM_SIGNED_SHORT_S    (SPM_SIGNED_SHORT<<8)
#define SPM_SIGNED_INT_S      (SPM_SIGNED_INT<<8)
#define SPM_FLOAT_S           (SPM_FLOAT<<8)
#define SPM_DOUBLE_S          (SPM_DOUBLE<<8)
#define SPM_UNSIGNED_SHORT_S  (SPM_UNSIGNED_SHORT<<8)
#define SPM_UNSIGNED_INT_S    (SPM_UNSIGNED_INT<<8)

#define ISSWAPPED(A) 	(A > 256)
#define DESWAP(A)		(A >> 8)
