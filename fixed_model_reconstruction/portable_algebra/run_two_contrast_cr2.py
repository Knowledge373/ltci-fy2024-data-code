"""Unchanged existing Python translation of clubSandwich 0.7.0 CR2 algebra.
This is not execution of the R package. Source attribution: Pustejovsky and Tipton.
"""
import numpy as np
from scipy import linalg, stats

def matrix_power(A,p):
 ev,U=linalg.eigh(A,lower=True,driver='evr');w=np.zeros_like(ev);keep=ev>1e-12;w[keep]=ev[keep]**p
 return (U*w)@U.T,ev,U

def calc_model(X,y,pref,names,label):
 # Full rank QR OLS and bread, without regularization.
 Q,R=linalg.qr(X,mode='economic',pivoting=False)
 Ri=linalg.solve_triangular(R,np.eye(R.shape[0]),lower=False)
 M=Ri@Ri.T;beta=linalg.solve_triangular(R,Q.T@y,lower=False);e=y-X@beta
 groups=[np.flatnonzero(pref==g) for g in range(1,48)];XX=[X[ix] for ix in groups]
 MU,xeig,_=matrix_power(sum(x.T@x for x in XX),-1)
 L=linalg.cholesky(M,lower=True)
 G=[];H=[];scores=[];eigrows=[];A=[];evs=[];vecs=[]
 for g,(ix,x) in enumerate(zip(groups,XX),1):
  # Algebraically identical CR2 target I-H, stably evaluated via QR.
  a,ev,u=matrix_power(np.eye(len(ix))-Q[ix]@Q[ix].T,-.5)
  gg=M@x.T@a;hh=gg@x@L
  A.append(a);G.append(gg);H.append(hh);scores.append(gg@e[ix]);evs.append(ev);vecs.append(u)
  evraw=linalg.eigh(np.eye(len(ix))-x@MU@x.T,lower=True,driver='evr')[0]
  eigrows.append({'model':label,'prefecture':g,'n':len(ix),'minimum_adjustment_QR':float(ev.min()),'minimum_raw_official_MU_diagnostic':float(evraw.min()),'raw_MU_suppressed1e12_diagnostic':int((evraw<=1e-12).sum()),'material_negative_count':int((ev < -1e-10).sum()),'suppressed_1e12':int((ev<=1e-12).sum()),'retained_1e14':int((ev>1e-14).sum()),'retained_1e12':int((ev>1e-12).sum()),'retained_1e10':int((ev>1e-10).sum()),'min_retained':float(ev[ev>1e-12].min())})
 scores=np.stack(scores,axis=1);V=scores@scores.T
 dfs=[];Ps=[]
 for k in range(X.shape[1]):
  hh=np.stack([h[k] for h in H],axis=1)
  P=np.diag([np.dot(g[k],g[k]) for g in G])-hh.T@hh
  dfs.append(float(np.trace(P)**2/np.sum(P*P)));Ps.append(P)
 se=np.sqrt(np.diag(V));dfs=np.array(dfs);crit=stats.t.ppf(.975,dfs)
 return dict(X=X,y=y,names=names,Q=Q,R=R,M=M,L=L,MU=MU,beta=beta,e=e,groups=groups,XX=XX,G=G,H=H,A=A,V=V,df=dfs,Ps=Ps,se=se,lo=beta-crit*se,hi=beta+crit*se,eig=eigrows,evs=evs,vecs=vecs)

def contrast(model,c):
 b=float(c@model['beta']);v=float(c@model['V']@c);se=np.sqrt(v)
 gg=[c@g for g in model['G']];hh=np.stack([c@h for h in model['H']],axis=1)
 P=np.diag([g@g for g in gg])-hh.T@hh
 df=float(np.trace(P)**2/np.sum(P*P));crit=float(stats.t.ppf(.975,df))
 # Same contrast, direct independent residual-projector P check.
 block=np.zeros((model['X'].shape[0],47))
 for j,ix in enumerate(model['groups']):block[ix,j]=gg[j]
 q=block-model['Q']@(model['Q'].T@block);Pdirect=q.T@q
 return {'estimate_fraction':b,'variance_fraction_squared':v,'SE_fraction':float(se),'df':df,'t975':crit,'CI_lower_fraction':float(b-crit*se),'CI_upper_fraction':float(b+crit*se),'P':P,'Pdirect':Pdirect}
