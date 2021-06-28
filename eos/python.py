from scipy.io import FortranFile
import numpy as np
class stagger_eos():
    """ Legacy stagger EOS tables. Not square, indirect indexing.
    """
    def __init__(self, data='table.dat'):
        fp=FortranFile(data)
        ((self.md,self.iupdate,self.nvar,self.mbox,self.eosxmin),
         (self.dbox,self.ul,self.ut,self.ur,self.eps,self.tff,self.grv,self.abnd)
        )=fp.read_record('(5)i4,(8)<f4')[0]
        ((self.tmean,self.tamp,self.rhm,self.xcorr,self.thmin,self.thmax,
          self.dth,self.eemin,self.eemax,self.deetab),(self.itab,self.mtab)
        )=fp.read_record('(10,{0:})f4,(2,{0:})i4'.format(self.md))[0]
        self.mtable=fp.read_ints()[0]
        self.tab=fp.read_record('i4')
        fp.close()
        self.njon = self.nvar-2.*self.mbox
        self.rhm1 = np.log(self.rhm[0])
        self.rhm2 = np.log(self.rhm[-1])
        self.drhm=(self.rhm2-self.rhm1)/(self.md-1)
        if (self.njon==0):
            self.njon=self.nvar-self.mbox

    def lookup(self, dim=[20,20,20],e=None,ee=None,d=None,lnd=None,pg=None,dpg1=None,dpg2=None,tt=None,ne=None,rk=None,scr=None,ss=None):
        mx=dim[0]
        my=dim[1]
        mz=dim[2]
        if type(lnd) is np.ndarray:
            lnd_loc=lnd
        elif type(d) is np.ndarray:
            lnd_loc=np.log(d)
        else:
            raise NameError ('no density given to eos!')
        if type(ee) is np.ndarray:
            ee_loc=ee
        elif e:
            if type(d) is np.ndarray:
                ee_loc=e/d
            else:
                ee_loc=e/np.exp(lnd_loc)
        else:
            raise NameError ('no energy given to eos!')
        for iz in range(0,mz):
            for iy in range(0,my):
                nx=np.zeros((mx), dtype=int)
                px=np.zeros((mx))
                for ix in range(0,mx):
                    algrk=(lnd_loc[ix,iy,iz]-self.rhm1)/self.drhm
                    nx[ix]=max(0,min(self.md-2,int(algrk)))
                    px[ix]=algrk-nx[ix]
                ntab  = np.zeros((mx), dtype=int)
                ntab1 = np.zeros((mx), dtype=int)
                py    = np.zeros((mx))
                ik    = np.zeros(mx, dtype=int)
                ik1   = np.zeros(mx, dtype=int)
                for ix in range(0,mx):
                    ntab[ix] = self.mtab[nx[ix]]
                    eek      = 1.+(ee_loc[ix,iy,iz]-self.eemin[nx[ix]])/self.deetab[nx[ix]]
                    eek      = max(eek,1.0)
                    kee      = min(ntab[ix]-1,max(1,int(eek)))
                    py[ix]   = eek - kee
                    ik[ix]   = self.itab[nx[ix]]+kee-1

                    ntab1[ix]= self.mtab[nx[ix]+1]
                    kee1     = kee - np.rint((self.eemin[nx[ix]+1]-self.eemin[nx[ix]])/self.deetab[nx[ix]])
                    kee1     = min(ntab1[ix]-1,max(1,kee1))
                    ik1[ix]  = self.itab[nx[ix]+1]+kee1-1
                f00 =np.zeros((mx))
                f01 =np.zeros((mx))
                f10 =np.zeros((mx))
                f11 =np.zeros((mx))
                fx00=np.zeros((mx))
                fx01=np.zeros((mx))
                fx10=np.zeros((mx))
                fx11=np.zeros((mx))
                fy00=np.zeros((mx))
                fy01=np.zeros((mx))
                fy10=np.zeros((mx))
                fy11=np.zeros((mx))
                for ix in range(0,mx):
                    qx   = 1. - px[ix]
                    pxqx = px[ix] * qx
                    pxpx = px[ix] * px[ix]
                    qxqx = qx * qx

                    qy   = 1. - py[ix]
                    pyqy = py[ix] * qy
                    pypy = py[ix] * py[ix]
                    qyqy = qy * qy

                    pxqy = px[ix] * qy
                    pxpy = px[ix] * py[ix]
                    qxqy = qx * qy
                    qxpy = qx * py[ix]

                    f00[ix] = qxqy * (1. + pxqx - pxpx + pyqy - pypy)
                    f01[ix] = qxpy * (1. + pxqx - pxpx + pyqy - qyqy)
                    f10[ix] = pxqy * (1. - qxqx + pxqx + pyqy - pypy)
                    f11[ix] = pxpy * (1. - qxqx + pxqx + pyqy - qyqy)

                    fx00[ix] =    qxqy * pxqx
                    fx01[ix] =    qxpy * pxqx
                    fx10[ix] =  - pxqy * pxqx
                    fx11[ix] =  - pxpy * pxqx

                    fy00[ix] =    qxqy * pyqy
                    fy01[ix] =  - qxpy * pyqy
                    fy10[ix] =    pxqy * pyqy
                    fy11[ix] =  - pxpy * pyqy

                if type(pg) is np.ndarray:
                    print ('working on pg')
                    for ix in range(0,mx):
                        pg[ix,iy,iz] = np.exp( \
                         f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                      +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                      + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                      + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                      + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                      + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                      + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                      + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                      + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                      + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                             )
                ik  = ik  + ntab
                ik1 = ik1 + ntab1

                if type(dpg1) is np.ndarray:
                    print ('working on dpg1')
                    for ix in range(0,mx):
                        dpg1[ix,iy,iz] = np.exp( \
                         f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                      +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                      + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                      + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                      + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                      + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                      + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                      + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                      + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                      + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                             )

                ik  = ik  + ntab
                ik1 = ik1 + ntab1

                if type(dpg2) is np.ndarray:
                    print ('working on dpg2')
                    for ix in range(0,mx):
                        dpg2[ix,iy,iz] = np.exp( \
                         f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                      +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                      + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                      + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                      + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                      + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                      + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                      + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                      + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                      + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                             )
                ik  = ik  + ntab
                ik1 = ik1 + ntab1

                if type(rk) is np.ndarray:
                    print ('working on rk')
                    for ix in range(0,mx):
                        rk[ix,iy,iz,0]   = np.exp( \
                             f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                          +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                          + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                          + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                          + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                          + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                          + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                          + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                          + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                          + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                                 )
                    if (rk.shape[3]>1):
                        for ij in range(1,rk.shape[3]):
                            rk[:,iy,iz,ij] = rk[:,iy,iz,ij-1]*10.0

                ik  = ik  + 3*ntab
                ik1 = ik1 + 3*ntab1

                if type(tt) is np.ndarray:
                    print ('working on tt')
                    for ix in range(0,mx):
                        tt[ix,iy,iz]   = ( \
                             f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                          +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                          + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                          + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                          + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                          + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                          + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                          + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                          + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                          + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                                 )

                ik  = ik  + 3*ntab
                ik1 = ik1 + 3*ntab1

                if type(ne) is np.ndarray:
                    print ('working on ne')
                    for ix in range(0,mx):
                        ne[ix,iy,iz]   = np.exp( \
                             f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                          +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                          + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                          + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                          + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                          + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                          + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                          + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                          + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                          + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                                 )

                if type(scr) is np.ndarray:
                    print ('working on src')
                    for ij in range(0,src.shape[3]):
                        for ix in range(0,mx):
                            ik[ix]  = ik[ix]  + 3*ntab[ix]
                            ik1[ix] = ik1[ix] + 3*ntab1[ix]
                            src[ix,iy,iz,ij] = np.exp( \
                             f00[ix] * self.tab[ik [ix]-1]+  f01[ix] * self.tab[ik [ix] + 1 -1] \
                          +  f10[ix] * self.tab[ik1[ix]-1]+  f11[ix] * self.tab[ik1[ix] + 1 -1] \
                          + fx00[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix]    ] \
                          + fx01[ix] * self.tab[ik [ix]-1 + 2 * ntab [ix] + 1] \
                          + fx10[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix]    ] \
                          + fx11[ix] * self.tab[ik1[ix]-1 + 2 * ntab1[ix] + 1] \
                          + fy00[ix] * self.tab[ik [ix]-1 +     ntab [ix]    ] \
                          + fy01[ix] * self.tab[ik [ix]-1 +     ntab [ix] + 1] \
                          + fy10[ix] * self.tab[ik1[ix]-1 +     ntab1[ix]    ] \
                          + fy11[ix] * self.tab[ik1[ix]-1 +     ntab1[ix] + 1] \
                                                 )
