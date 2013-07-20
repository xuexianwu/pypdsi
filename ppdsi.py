#!/usr/bin/env python
from sys import argv
import numpy as np

#Script of using a python wrapper to driver long-term pdsi calculation.
#2013-02-01 10:37:43

#$BDIR/pypdsi.py $climate $prefix $years $dims $awcfile $latfile
odir, prefix, years, dims, awc, lat = argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]

def pdsi_wrapper(vd='.', odir='.', dims=[1812,621,1405], yrs='1950 2100', bin='pdsicpp/pdsi', prdata=[], tadata=[], awc='', lat=''):
    """pdsi calculate configure includes:
    vd = ', ramdisk to speed up, otherwise eq odir
    odir = '', output dir 
    dims = [times, cols, rows], dimension of time series data
    yrs = '1950 2010', stat & end of year, string delimited with space
    prdata = [], list of sorted precipitation data files
    tadata = [], list of sorted temperature data files
    awc = '', filename of awc
    lat = '', filename of latitude
    bin = 'bin/pdsi' no change
    """
    import numpy as np
    rows, cols = range(dims[1]), range(dims[2])
    yrs = yrs.split()
    years = range(int(yrs[0]), int(yrs[1])+1)
    n_month = len(years) * 12
    pdsi_bip, pdsi_bip2 = odir+'/pdsi-bip', odir+'/pdsi-bip-cal'
    out = open(pdsi_bip,'w')
    out1 = open(pdsi_bip2,'w')
    fvalue = np.ones(n_month, 'f4') * (-99.99)
    pr_fn, ta_fn, ta_normal_fn = vd+'/monthly_P', vd+'/monthly_T', vd+'/mon_T_normal'
    command_line = '%s -i%s -o%s -m' % (bin, vd, vd) # using metric unit
    fn_parameter, out_pdsi_fn = vd+'/parameter', vd+'/monthly/original/PDSI.tbl'
    out_pdsi_fncal = vd+'/monthly/self_cal/PDSI.tbl'
    savefmt='%4i %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f'
    if not os.path.exists(odir): os.path.mkdir(odir)
    def rd(fn): return np.fromfile(fn, dtype='f4').reshape(dims[1], dims[2])
    def pdsi_core(pr, ta, awc, lat, years): # metrix unit is default. 
        n_yrs = len(pr) / 12
        pdsi = pr * 0
        pr_vec = np.zeros((n_yrs, 13), dtype='f4')
        ta_vec = np.zeros((n_yrs, 13), dtype='f4')
        ta_normal = np.zeros(12, dtype='f4')
        pr_vec[:, 0] = ta_vec[:, 0] = years 
        # create parameter file
        fpar = open(fn_parameter, 'w')
        fpar.write('%f %f' % (awc, lat))
        fpar.close()
        for i, v in enumerate(years):
            pr_vec[i, 1:] = pr[i * 12:i * 12 + 12]# / 25.4
            ta_vec[i, 1:] = ta[i * 12:i * 12 + 12]# * (9.0 / 5.0) + 32;
        np.savetxt(pr_fn, pr_vec, fmt=savefmt)
        np.savetxt(ta_fn, ta_vec, fmt=savefmt)
        fout = open(ta_normal_fn, 'w')
        for i in range(12):
            ta_normal[i] = np.mean(ta_vec[:, i + 1])
            fout.write('  %.2f' % ta_normal[i])
        fout.close()
        os.system(command_line)
        tmp = np.loadtxt(out_pdsi_fn, dtype='f4')
        tmp1 = np.loadtxt(out_pdsi_fncal, dtype='f4')
        return tmp[:, 1:].flatten(), tmp1[:, 1:].flatten()
    def pdsi_series():
        pr = np.load(exdir+'../exchange/climate.prism.prcp.timeseries')
        pr1 = np.load(exdir+'../exchange/climate.ea_rcp60.prcp.timeseries')
        pr = np.array(list(pr) + list(pr1)[60:])
        ta = np.load(exdir+'../exchange/climate.prism.tmax.timeseries')
        ta1 = np.load(exdir+'../exchange/climate.ea_rcp60.tmax.timeseries')
        ta = np.array(list(ta) + list(ta1)[60:])
        tm = np.load(exdir+'../exchange/climate.prism.tmin.timeseries')
        tm1 = np.load(exdir+'../exchange/climate.ea_rcp60.tmin.timeseries')
        tm = np.array(list(tm) + list(tm1)[60:])
        ta = (ta + tm) / 2.
        pdsi, cal = pdsi_core(pr, ta, 12.8, 40., years)
        pdsi.dump(exdir+'pdsi.origin.prismrcp60')
        cal.dump(exdir+'pdsi.cal.prismrcp60')
    def pdsi_area():
        Awc, Lat = rd(awc), rd(lat)
        #start this big loop
        for row in rows:
            pr = conv_bip(prdata, odir, dims, row)
            ta = conv_bip(tadata, odir, dims, row)
            for col in cols:
                if Awc[row, col] == 0 or Awc[row, col] < -9900.: 
                    out.write(fvalue.data)
                    out1.write(fvalue.data)
                    continue
                pdsi, cal = pdsi_core(pr[:,col], ta[:, col], Awc[row, col], Lat[row, col], years)
                out.write(pdsi.data)
                out1.write(cal.data)
            logger.info('row %d/%d' % (row, dims[1]))
        out.close()
        out1.close()
    pdsi_series()

def conv_bip(files, odir, dims, row):
    bipfn = odir + '.bmp-bip'
    import cStringIO
    mm = cStringIO.StringIO()
    for i, fn in enumerate(files):
        f = open(fn)
        offset = dims[2] * row * 4
        f.seek(offset)
        mm.write(f.read(dims[2] * 4))
    out = open(bipfn, 'w')
    out.write(mm.getvalue())
    out.close()
    mm.close()
    return np.fromfile(bipfn,dtype='f4').reshape(dims[0], dims[2])

def test_gen():
    for i in range(25):
        data = np.ones((100,75),'f4') * i
        oname = 'test/%02d' % i
        data.tofile(oname)

def test_bip():
    files = glob('prism.1900-2010/us_tdmean_2*.??')
    dims = [132, 621, 1405]
    rows, cols = dims[1], dims[2]
    for row in range(200, rows):
        data = conv_bip(files, dims, row)
        for v in data[:,200]:
            print '%.4f' % v
        break

def disp_bip():
    import matplotlib.pyplot as plt
    dims = [621,132,1405]
    dims = [132,621,1405]
    data = np.fromfile('bip', dtype='f4').reshape(dims[1], dims[0], dims[2])
    data = data[:,0,:]
    plt.figure()
    im = plt.imshow(histeq(data))
    plt.savefig('figure.png', dpi=100)

def ustest_disp():
    import matplotlib.pyplot as plt
    dims = [1332, 124, 281]
    data = np.fromfile('ustest/pdsi-bip','f4').reshape(dims[1],dims[2],dims[0])
    print data[50,50,:]
    #data = data[:,:,0]
    #plt.figure()
    #im = plt.imshow(histeq(data))
    #plt.savefig('figure.png', dpi=100)

def prism004():
    """US 0.04 degree data, 1950 - 2010"""
    src = '/nobackup/jxiong1/climate-data/us004/prism/'
    vd = 'tmp-pdsi-prism/'
    odir = '/nobackup/jxiong1/model-out/pdsi-prism/'
    dims = [360, 621, 1405]
    yrs = '1950 1979'
    bin = 'bin/pdsi'
    prdata = glob(src + 'prism.prcp.*')[:360]
    tadata = glob(src + 'prism.tdmean.*')[:360]
    prdata.sort()
    tadata.sort()
    awc = '/nobackup/jxiong1/climate-data/us004/us004.f4.awc'
    lat = '/nobackup/jxiong1/climate-data/us004/us004.f4.latitude'
    pdsi_wrapper(vd, odir, dims, yrs, bin, prdata, tadata, awc, lat)
def us80045():
    """US 0.04 degree data, 1950 - 2100 prism+rcp45"""
    src = '/nobackup/jxiong1/climate-data/us004/prism/'
    src1 = '/nobackup/jxiong1/climate-data/us004/ea_rcp45/'
    vd = '../tmp/pdsi45/'
    odir = vd
    dims = [1812, 621, 1405]
    yrs = '1950 2100'
    bin = 'pdsicpp/pdsi'
    prdata = fsearch(src + 'prism.prcp.*')
    tadata = fsearch(src + 'prism.tdmean.*')
    prdata += fsearch(src1 + 'ea_rcp45.prcp.*')[60:]
    tadata += fsearch(src1 + 'ea_rcp45.tdmean.*')[60:]
    #print prdata, tadata, len(prdata), len(tadata)
    awc = '/nobackup/jxiong1/climate-data/us004/us004.f4.awc'
    lat = '/nobackup/jxiong1/climate-data/us004/us004.f4.latitude'
    pdsi_wrapper(vd, odir, dims, yrs, bin, prdata, tadata, awc, lat)

def us80060():
    """US 0.04 degree data, 1950 - 2100 prism+rcp60 """
    src = '/nobackup/jxiong1/climate-data/us004/prism/'
    src1 = '/nobackup/jxiong1/climate-data/us004/ea_rcp60/'
    vd = '../tmp/pdsi60/'
    odir = vd
    dims = [1812, 621, 1405]
    yrs = '1950 2100'
    bin = 'pdsicpp60/pdsi'
    prdata = fsearch(src + 'prism.prcp.*')
    tadata = fsearch(src + 'prism.tdmean.*')
    prdata += fsearch(src1 + 'ea_rcp60.prcp.*')[60:]
    tadata += fsearch(src1 + 'ea_rcp60.tdmean.*')[60:]
    #print prdata, tadata, len(prdata), len(tadata)
    awc = '/nobackup/jxiong1/climate-data/us004/us004.f4.awc'
    lat = '/nobackup/jxiong1/climate-data/us004/us004.f4.latitude'
    pdsi_wrapper(vd, odir, dims, yrs, bin, prdata, tadata, awc, lat)

def resample():
    import scipy.ndimage as nd
    dims = [621, 1405]
    files = glob('prism.1900-2010/*.??')
    #files = 'us004.f4.awc us004.f4.latitude'.split()
    odir = 'small/'
    factor = [.2,.2]
    for i, f in enumerate(files):
        oname = odir + os.path.basename(f)
        data = np.fromfile(f,dtype='f4').reshape(dims[0], dims[1])
        data = nd.interpolation.zoom(data, factor)
        data /= 100.
        data.tofile(oname)

def scale():
    dims = [621, 1405]
    files = glob('prism.1900-2010/*.??')
    odir = 'scale/'
    for i, f in enumerate(files):
        oname = odir + os.path.basename(f)
        data = np.fromfile(f,dtype='f4')
        data /= 100.
        data.tofile(oname)

def asc2bin():
    files = glob('prism.1900-2010.asc/*tdmean_????.??')
    for f in files:
        oname = 'prism.1900-2010/' + f.split('/')[-1]
        if os.path.isfile(oname): continue
        data = np.loadtxt(f, dtype='f4', skiprows=6).reshape(621,1405)
        data = data[::-1] / 100.
        data.tofile(oname)
        print oname

def tm():
    src = '/nobackup/jxiong1/climate-data/us004/ea_rcp60'
    prefix = 'ea_rcp60'
    tdmean(src, prefix)

def tdmean(src, prefix):
    from glob import glob
    print(src + '/%s*.tmax.*')
    tmaxs = glob(src + '/%s*.tmax.*' % prefix)
    tmaxs.sort()
    tmins = glob(src + '/%s*.tmin.*' % prefix)
    tmins.sort()
    for i, f in enumerate(tmaxs):
        oname = f.replace('tmax', 'tdmean')
        try:
            data = np.fromfile(f, dtype='f4')
        except:
            print f
        try:
            data1 = np.fromfile(tmins[i], dtype='f4')
        except:
            print tmins[i]
        data = (data + data1) / 2
        data.tofile(oname)
        logger.info(os.path.basename(oname))

if __name__ == '__main__':
    if arg1 != '':
        for opt in argv[1:]:
            if opt == 'e': mai_edit()
            elif opt == 'd':
                mai_doc()
                exit()
        exec('%s()' % arg1)
        mai_done()
    else:
        mai_about()
