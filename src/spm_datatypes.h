/* 
 * $Id: spm_datatypes.h 938 2007-10-12 19:09:31Z john $
 * John Ashburner
 */

/* SPM image data types */

#ifndef _SPM_DATATYPES_H_
#define _SPM_DATATYPES_H_

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

#endif /* _SPM_DATATYPES_H_ */
