#ifndef lint
static char sccsid_connex[]="@(#)connex.h	1.2 JB Poline 96/03/04";
#endif
 


#define EPS2 	4.4409e-16

typedef struct cc_position_3d_str
	{
	float	t, sum, label,  *ptr_flt;
	int	x, y, z, dx, dy, dz, size;
	}
	cc_position_3d;


static	cc_position_3d	pos3;



/*------------------------------------ 3D ---------------------------*/

void comp_c_c6()		/* 6 connexity version */
{
	if(*(pos3.ptr_flt) > pos3.t)
	{
		/* pos3.sum += *(pos3.ptr_flt); */ /* for future use ? */
		(pos3.size)++;
		*(pos3.ptr_flt) = pos3.label;
	
		if(pos3.y > 0) {
			(pos3.y)--;
			pos3.ptr_flt -= pos3.dx;
			comp_c_c6();
			(pos3.y)++;
			pos3.ptr_flt += pos3.dx;
		}
		if(pos3.y < pos3.dy - 1) {
			(pos3.y)++;
			pos3.ptr_flt = pos3.ptr_flt + pos3.dx;
			comp_c_c6();
			(pos3.y)--;
			pos3.ptr_flt = pos3.ptr_flt - pos3.dx;
		}
		if(pos3.x < pos3.dx - 1) {
			(pos3.x)++;
			pos3.ptr_flt = pos3.ptr_flt + 1;
			comp_c_c6();
			(pos3.x)--;
			pos3.ptr_flt = pos3.ptr_flt - 1;
		}
		if(pos3.x > 0) {
			(pos3.x) --;
			pos3.ptr_flt = pos3.ptr_flt - 1;
			comp_c_c6();
			(pos3.x) ++;
			pos3.ptr_flt = pos3.ptr_flt + 1;
		}
		if(pos3.z < pos3.dz - 1) {
			(pos3.z)++;
			pos3.ptr_flt += (pos3.dx * pos3.dy);
			comp_c_c6();
			(pos3.z)--;
			pos3.ptr_flt -= (pos3.dx * pos3.dy);
		}
		if(pos3.z > 0) {
			(pos3.z) --;
			pos3.ptr_flt -= (pos3.dx * pos3.dy);
			comp_c_c6();
			(pos3.z) ++;
			pos3.ptr_flt += (pos3.dx * pos3.dy);
		}
	}
	return;
}




void comp_c_c18( )
{
	if(*(pos3.ptr_flt) > pos3.t)
	{
		/* pos3.sum += *(pos3.ptr_flt); */ /* for future use ? */
		(pos3.size)++;
		*(pos3.ptr_flt) = pos3.label;
	
		/* North */
		if(pos3.y > 0) {
			(pos3.y)--; pos3.ptr_flt -= pos3.dx;
			 comp_c_c18();
			(pos3.y)++; pos3.ptr_flt += pos3.dx;
		}
		/* South */
		if(pos3.y < pos3.dy - 1) {
			(pos3.y)++; pos3.ptr_flt += pos3.dx;
			 comp_c_c18();
			(pos3.y)--; pos3.ptr_flt -= pos3.dx;
		}
		/* East */
		if(pos3.x < pos3.dx - 1) {
			(pos3.x)++; pos3.ptr_flt += 1; 
			comp_c_c18();
			(pos3.x)--; pos3.ptr_flt -= 1; 
		}
		/* West */
		if(pos3.x > 0) {
			(pos3.x)--; pos3.ptr_flt -= 1; 
			comp_c_c18();
			(pos3.x)++; pos3.ptr_flt += 1; 
		}



		/* NEast */
		if(pos3.x < pos3.dx - 1 && pos3.y > 0) {
			(pos3.x)++; (pos3.y)--;
			pos3.ptr_flt = pos3.ptr_flt - pos3.dx + 1;
			comp_c_c18();
			(pos3.x)--; (pos3.y)++;
			pos3.ptr_flt = pos3.ptr_flt + pos3.dx - 1;
		}
		/* NWest */
		if(pos3.x > 0 && pos3.y > 0) {
			(pos3.x) --; (pos3.y)--;
			pos3.ptr_flt = pos3.ptr_flt - pos3.dx - 1;
			comp_c_c18();
			(pos3.x) ++; (pos3.y)++;
			pos3.ptr_flt = pos3.ptr_flt + pos3.dx + 1;
		}
		/* SEast */
		if(pos3.x < pos3.dx - 1 && pos3.y < pos3.dy - 1) {
			(pos3.x)++; (pos3.y)++;
			pos3.ptr_flt = pos3.ptr_flt + pos3.dx + 1;
			comp_c_c18();
			(pos3.x)--; (pos3.y)--;
			pos3.ptr_flt = pos3.ptr_flt - pos3.dx - 1;
		}
		/* SWest */
		if(pos3.x > 0 && pos3.y < pos3.dy - 1) {
			(pos3.x) --; (pos3.y)++;
			pos3.ptr_flt = pos3.ptr_flt + pos3.dx - 1;
			comp_c_c18();
			(pos3.x) ++; (pos3.y)--;
			pos3.ptr_flt = pos3.ptr_flt - pos3.dx + 1;
		}



		/* Lower */
		if(pos3.z < pos3.dz - 1) {
			(pos3.z)++; 
			pos3.ptr_flt += (pos3.dx * pos3.dy);
			comp_c_c18();
			(pos3.z)--; 
			pos3.ptr_flt -= (pos3.dx * pos3.dy);
		}
		/* Lower N */
		if(pos3.z < pos3.dz - 1 && pos3.y > 0) {
			(pos3.z)++; (pos3.y)--; 
			pos3.ptr_flt += ((pos3.dx * pos3.dy) - pos3.dx);
			comp_c_c18();
			(pos3.z)--; (pos3.y)++; 
			pos3.ptr_flt -= ((pos3.dx * pos3.dy) - pos3.dx);
		}
		/* Lower S */
		if(pos3.z < pos3.dz - 1 && pos3.y < pos3.dy - 1) {
			(pos3.z)++; (pos3.y)++; 
			pos3.ptr_flt += ((pos3.dx * pos3.dy) + pos3.dx);
			comp_c_c18();
			(pos3.z)--; (pos3.y)--; 
			pos3.ptr_flt -= ((pos3.dx * pos3.dy) + pos3.dx);
		}
		/* Lower E */
		if(pos3.z < pos3.dz - 1 && pos3.x < pos3.dx - 1) {
			(pos3.z)++; (pos3.x)++; 
			pos3.ptr_flt += ((pos3.dx * pos3.dy) + 1);
			comp_c_c18();
			(pos3.z)--; (pos3.x)--; 
			pos3.ptr_flt -= ((pos3.dx * pos3.dy) + 1);
		}
		/* Lower W */
		if(pos3.z < pos3.dz - 1 && pos3.x > 0) {
			(pos3.z)++; (pos3.x)--; 
			pos3.ptr_flt += ((pos3.dx * pos3.dy) - 1);
			comp_c_c18();
			(pos3.z)--; (pos3.x)++; 
			pos3.ptr_flt -= ((pos3.dx * pos3.dy) - 1);
		}


		/* Upper */
		if(pos3.z > 0) {
			(pos3.z)--; 
			pos3.ptr_flt += -(pos3.dx * pos3.dy);
			comp_c_c18();
			(pos3.z)++; 
			pos3.ptr_flt -= -(pos3.dx * pos3.dy);
		}
		/* Upper N */
		if(pos3.z > 0 && pos3.y > 0) {
			(pos3.z)--; (pos3.y)--; 
			pos3.ptr_flt += (-(pos3.dx * pos3.dy) - pos3.dx);
			comp_c_c18();
			(pos3.z)++; (pos3.y)++; 
			pos3.ptr_flt -= (-(pos3.dx * pos3.dy) - pos3.dx);
		}
		/* Upper S */
		if(pos3.z > 0 && pos3.y < pos3.dy - 1) {
			(pos3.z)--; (pos3.y)++; 
			pos3.ptr_flt += (-(pos3.dx * pos3.dy) + pos3.dx);
			comp_c_c18();
			(pos3.z)++; (pos3.y)--; 
			pos3.ptr_flt -= (-(pos3.dx * pos3.dy) + pos3.dx);
		}
		/* Upper E */
		if(pos3.z > 0 && pos3.x < pos3.dx - 1) {
			(pos3.z)--; (pos3.x)++; 
			pos3.ptr_flt += (-(pos3.dx * pos3.dy) + 1);
			comp_c_c18();
			(pos3.z)++; (pos3.x)--; 
			pos3.ptr_flt -= (-(pos3.dx * pos3.dy) + 1);
		}
		/* Upper W */
		if(pos3.z > 0 && pos3.x > 0) {
			(pos3.z)--; (pos3.x)--; 
			pos3.ptr_flt += (-(pos3.dx * pos3.dy) - 1);
			comp_c_c18();
			(pos3.z)++; (pos3.x)++; 
			pos3.ptr_flt -= (-(pos3.dx * pos3.dy) - 1);
		}


	}
	return;
}


/*---------------------------------------------------------------*/

