import numpy as np
import ppr
import sklearn.preprocessing

def supsmu(x,y,span=0.0,alpha=0.0):
    '''wrapper to super smoother function
    '''
    n = len(x)
    w = np.ones_like(x)
    smo = np.ones_like(x)
    sc = np.empty((n,7),order="F")
    
    ppr.supsmu(x,y,w,1,span,alpha,smo,sc)
    
    #print sc
    
    return sc


class ProjectionPursuitRegression(object):
    def __init__(self,ml,mu,verbose=0):
        assert mu >= 1
        assert ml >= 1
        assert mu <= ml
        
        self.mu = mu
        self.ml = ml
        self.verbose = verbose
        self.scaler = sklearn.preprocessing.MinMaxScaler()
        
    def get_smoother(self,p):
        return self.t[p,:],self.f[p,:]*self.std + self.mean        

    def fit(self,_x,y):
        '''
        x is a matrix of dimension (n,p): n number of samples, p number of features
        y is a vector of size n: n number of samples or a matrix of (n,q), q number of responses
        '''
        assert len(y.shape) == 1

        #one dimensional
        ny = len(y)
        q = 1
        y = y.reshape((ny,1)).T

        nx,p = _x.shape
            
        assert  nx == ny
        
        n = nx

        self.features = p

        
        #apply scaler
        x = self.scaler.fit_transform(_x).T
        
        n = nx

        ml = self.ml
        mu = self.mu

        nsmod = ml * (p + q + 2 * n) + q + 7 + ml + 1
        nsp = n * (q + 15) + q + 3 * p
        ndp = p * (p + 1)/2 + 6 * p

        w = np.ones(n)
        ww = np.ones(q)
        smod = np.empty(nsmod,dtype="float32")
        sp = np.empty(nsp,dtype="float32")
        dp = np.empty(ndp,dtype="float64")

        x = np.asfortranarray(x,dtype=np.float32)
        y = np.asfortranarray(y,dtype=np.float32)
        ww = np.asfortranarray(ww,dtype=np.float32)
        w = np.asfortranarray(w,dtype=np.float32)
        smod = np.asfortranarray(smod,dtype=np.float32)

        ppr.parms.ifl = self.verbose
        ppr.smartr(ml,mu,w,x,y,ww,smod,sp,dp)

        #print "ml",smod[0],ml
        #print "p",smod[1],p
        #print "q",smod[2],q
        #print "n",smod[3],n
        #print "mu",smod[4],mu
        #print "nsmod",nsmod
        #print "nsp",nsp
        #print "ndp",ndp
        #print "mean",smod[5],np.mean(y)
        #print "std",smod[6],np.std(y)
        
        self.std = smod[6]
        self.mean = smod[5]

        alpha_slice_start = q + 6
        alpha_slice_size = (p * ml)
        alpha_slice = slice(alpha_slice_start,alpha_slice_start+alpha_slice_size)

        alpha = smod[alpha_slice]
        #print alpha
        alpha = alpha.reshape((ml,p))
        #print alpha,alpha_slice,alpha_slice_start

        beta_slice_start = alpha_slice_start + alpha_slice_size
        beta_slice_size = q*ml
        beta_slice = slice(beta_slice_start,beta_slice_start + beta_slice_size)
        beta = smod[beta_slice]
        beta = beta.reshape((q,ml))

        #print "beta",beta_slice,beta,beta_slice_start

        f_slice_start = beta_slice_start + beta_slice_size
        f_slice_size = (n * ml)
        f_slice = slice(f_slice_start,f_slice_start+f_slice_size)
        f = smod[f_slice]
        f = f.reshape((ml,n))

        #print "f",f,f_slice,f_slice_start

        t_slice_start = f_slice_start + f_slice_size
        t_slice_size = (n * ml)
        t_slice = slice(t_slice_start,t_slice_start+t_slice_size)
        t = smod[t_slice]
        t = t.reshape((ml,n))
        
        #for k in range(ml):
        #    l = t[k]
        #    print all(l[i] <= l[i+1] for i in xrange(len(l)-1))

        #print "t",t,t_slice,t_slice_start,t_slice_start+t_slice_size,len(smod)
        
        self.alpha = alpha
        self.beta = beta
        self.t = t
        self.f = f 
        
        return self
        
    def predict(self,_x):
        #project
        n,p = _x.shape

        assert p == self.features
        #p = np.zeros(n)
        
        #apply scaler
        x = self.scaler.transform(_x).T

        prediction = np.empty(n)

        for i in range(self.mu):
            p = np.dot(self.alpha[i],x)
            
            #print i,np.min(p),np.max(p)
            #print i,np.min(self.t[i]),np.max(self.t[i])
            #print i,np.min(self.f[i]),np.max(self.f[i])
            #interpolate in the smoother

            sp = np.interp(p,self.t[i],self.f[i])
            
            #print "sp",self.beta[0,i],sp
            prediction += sp*float(self.beta[0,i])
            
            #add beta
            #prediction *= float(self.beta[0,i])
        
        #rescale
        prediction = prediction * self.std + self.mean
        return prediction
        
    def score(self,x,y):
        prediction = self.predict(x)
        from sklearn.metrics import explained_variance_score
        
        return explained_variance_score(y,prediction)
