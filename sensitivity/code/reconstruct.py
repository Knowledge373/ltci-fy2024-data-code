"""Independent arithmetic on fitted coefficients; no model is fitted here."""
import math
import numpy as np

def rk(A):
 s=np.linalg.svd(A,compute_uv=False);tol=max(A.shape)*np.finfo(float).eps*s[0] if len(s) else 0
 return {'rank':int(np.sum(s>tol)),'n':A.shape[0],'p':A.shape[1],'tolerance':float(tol),'singular_values':s.tolist()}

def expit(z):
 if z>=0:return 1/(1+math.exp(-z))
 e=math.exp(z);return e/(1+e)

def linpred(X,b):return np.array([math.fsum(float(v)*float(c) for v,c in zip(row,b)) for row in X])
