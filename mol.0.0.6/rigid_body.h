#ifndef _MOL_RIGID_H_
#define _MOL_RIGID_H_

struct rigidbody {
	double center[3];
	double *origin;
};

void ag2rigidbody( struct rigidbody *rigidbody, struct atomgrp *ag);

void rigidbody2ag(double *change, struct atomgrp *ag,
		  struct rigidbody *rigidbody);

void mol_rigidbody_grad(double *grad, struct atomgrp *ag, double *inp, double *origin);

#endif
