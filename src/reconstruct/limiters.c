
double slope_limiter(const double dqm, const double dqp) {
 
  double w = SLOPEW;
  
  double rq,dq,xi,xiR;
  
  dq = .5*(1 + w)*dqm + .5*(1-w)*dqp;
  
  xiR = 2/(1 - w + (1+w)*rq);
  
  if (dqp != 0 ) {
    rq = dqm/dqp;
  }
  else {
    rq = -1;
  }
  
  xi = user_limiter(rq,xiR);
  
  return xi*dq;
 
  
  
}

double minmod_limiter(const double rq, const double xiR) {
  
  if (rq <= 0) {
    return 0;
    
  }
  else {
    if (rq <= 1) {
      return rq;
    }
    else {
      return fmin(1,xiR);
    }
  }
 
}
double minmod_limiter(const double rq, const double xiR) {
  
  if (rq <= 0) {
    return 0;  
  }
  else if (rq <= 1) {
    return rq;
  else {
    return fmin(1,xiR);
  }

 
}
double superbee_limiter(const double rq, const double xiR) {
  
  if (rq <= 0) {
    return 0;
    
  }
  else if (rq <= .5) {
    return rq;
  }
  else if (rq <= 1) {
    return 1.0;
  }
  else {
    return fmin(2.0, fmin(xiR,rq));
  }
}

double vanleer_limiter(const double rq, const double xiR) {
  
  if (rq <= 0) {
    return 0;
    
  }
  else {
    return fmin(2*rq/(1+rq), xiR);
  }
}
double vanalbada_limiter(const double rq, const double xiR) {
  
  if (rq <= 0) {
    return 0;
    
  }
  else {
    return fmin(rq*(1+rq)/(1+rq*rq), xiR);
  }
}