/*
 * John Ashburner
 * Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging
 */

/* SPM image data types */

#ifndef _SPM_DATATYPES_H_
#define _SPM_DATATYPES_H_

#define SPM_UNSIGNED_CHAR     2
#define SPM_SIGNED_SHORT      4
#define SPM_SIGNED_INT        8
#define SPM_FLOAT             16
#define SPM_DOUBLE            64
#define SPM_SIGNED_LONG_LONG  1024
#define SPM_UNSIGNED_LONG_LONG 1280
#define SPM_SIGNED_CHAR       (SPM_UNSIGNED_CHAR+128) 
#define SPM_UNSIGNED_SHORT    (SPM_SIGNED_SHORT+128) 
#define SPM_UNSIGNED_INT      (SPM_SIGNED_INT+128)

/* codes for byte swapped types */
#define SPM_SIGNED_SHORT_S       (32768 | SPM_SIGNED_SHORT)
#define SPM_SIGNED_INT_S         (32768 | SPM_SIGNED_INT)
#define SPM_FLOAT_S              (32768 | SPM_FLOAT)
#define SPM_DOUBLE_S             (32768 | SPM_DOUBLE)
#define SPM_SIGNED_LONG_LONG_S   (32768 | SPM_SIGNED_LONG_LONG) 
#define SPM_UNSIGNED_LONG_LONG_S (32768 | SPM_UNSIGNED_LONG_LONG)
#define SPM_UNSIGNED_SHORT_S     (32768 | SPM_UNSIGNED_SHORT)
#define SPM_UNSIGNED_INT_S       (32768 | SPM_UNSIGNED_INT)

#endif /* _SPM_DATATYPES_H_ */
