"""Line widths and conditional measurement covariance for W1/W2/W3.

The stored model fits forward-fold intrinsic profiles through a Gaussian LSF
and pixel integration. F_broad in units 1e-18 erg s^-1 cm^-2 is emitted flux
before absorption; width in km s^-1. Double-Gaussian width is of the sum.
Nonparametric covariance conditions on the stored subtraction and CV alpha.
"""
from common import *
from lrdwidth import linefit as core
from scipy.optimize import nnls
from astropy.cosmology import Planck18

def log_measurement(f):
    names=f['names'];p=np.array(f['params']);cov=np.array(f['covariance_conditional'])
    flux_idx=[names.index(n) for n in ['F_broad1','F_broad2'] if n in names]
    flux=p[flux_idx].sum();J=np.zeros((2,len(p)));J[0,flux_idx]=1/(np.log(10)*flux)
    if f['kind']!='double_gaussian':
        k=names.index('W_broad1');width=p[k];J[1,k]=1/(np.log(10)*width)
    else:
        # Exact same 10 km/s grid and outermost half-maximum crossings as T2.
        v=np.arange(-15000,15000.01,10.)
        def width_of(par):
            t=dict(zip(names,par));shape=sum(t['F_broad'+str(i)]*core.gauss(v,t['v_broad'+str(i)],t['W_broad'+str(i)]) for i in [1,2])
            return core.widths(v,shape)['FWHM_kms']
        width=width_of(p)
        for i,name in enumerate(names):
            if name not in ['F_broad1','F_broad2','v_broad1','v_broad2','W_broad1','W_broad2']:continue
            h=max(abs(p[i])*1e-5,1e-5);up=p.copy();dn=p.copy();up[i]+=h;dn[i]-=h
            J[1,i]=(np.log10(width_of(up))-np.log10(width_of(dn)))/(2*h)
    assert np.isclose(width,f['FWHM_intrinsic_kms'],rtol=1e-9)
    C=J@cov@J.T
    assert np.isfinite(C).all() and np.linalg.eigvalsh(C).min()>-1e-10
    return dict(width_kms=width,flux_cgs=flux*core.FUNIT,cov_log_flux_width=C,jacobian=J)

def np_bootstrap(d,npfit,nboot=200,seed=0):
    knots=np.array(npfit['knots_kms']);sol=np.array(npfit['intrinsic_knots']);fine=np.array(npfit['v_fine'])
    basis=np.maximum(1-np.abs(d['grid'][:,None]-knots)/125.,0);A=np.asarray(d['K']@basis)
    D2=np.diff(np.eye(len(knots)),n=2,axis=0);alpha=npfit['cv_alpha_grid'][npfit['selected_alpha_index']]
    mat=np.vstack([A/d['err'][:,None],np.sqrt(alpha)*D2]);pred=A@sol
    assert np.allclose(pred,npfit['nonparametric_model_on_pixels'],rtol=1e-9,atol=1e-9)
    rng=np.random.default_rng(seed);samples=[]
    for i in range(nboot):
        yy=pred+rng.normal(0,d['err']);p=nnls(mat,np.r_[yy/d['err'],np.zeros(len(D2))],maxiter=3000)[0]
        m=core.widths(fine,np.interp(fine,knots,p))
        samples.append([m['flux_window_erg_s_cm2'],m['W50_kms'],m['FWHM_kms'] or np.nan])
    a=np.array(samples);out={}
    for i,label in [(1,'nonparametric_W50'),(2,'nonparametric_FWHM')]:
        mask=np.isfinite(a[:,[0,i]]).all(1)&(a[:,[0,i]]>0).all(1)
        vals=a[mask][:,[0,i]];out[label]=dict(N_usable=int(mask.sum()),N_requested=nboot,cov_log_flux_width=np.cov(np.log10(vals),rowvar=False,ddof=1),flux_p16_p50_p84=np.quantile(vals[:,0],[.16,.5,.84]),width_p16_p50_p84=np.quantile(vals[:,1],[.16,.5,.84]))
        assert mask.sum()>=max(20,.8*nboot),('insufficient conditional bootstrap',label,mask.sum())
    return out,a

def make_row(meta,m):
    L=m['flux_cgs']*4*np.pi*Planck18.luminosity_distance(float(meta['z_fit'])).to_value('cm')**2
    return dict(source_id=meta['source_id'],file_id=meta['file_id'],cls=meta['lrd_class_stage2'],x=np.log10(L)-43,y=np.log10(m['width_kms']),cov=m['cov_log_flux_width'],width_kms=m['width_kms'],luminosity_erg_s=L,flux_erg_s_cm2=m['flux_cgs'])
