#!/usr/bin/env python
from sys import argv
import numpy as np, os

#Script of using a python wrapper to driver long-term pdsi calculation.
#2013-02-01 10:37:43
def pdsi_wrapper(vd='.', odir='.', dims=[1812,621,1405], years=range(1950,2100), bin='pdsicpp/pdsi', prdata=[], tadata=[], awc='', lat=''):
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
        Awc, Lat = rd(lat), rd(lat)
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
        out.close()
        out1.close()
    pdsi_area()

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

def asc2bin():
    files = glob('prism.1900-2010.asc/*tdmean_????.??')
    for f in files:
        oname = 'prism.1900-2010/' + f.split('/')[-1]
        if os.path.isfile(oname): continue
        data = np.loadtxt(f, dtype='f4', skiprows=6).reshape(621,1405)
        data = data[::-1] / 100.
        data.tofile(oname)
        print oname


#$BDIR/pypdsi.py $climate $prefix $years $dims $awcfile $latfile
climate, odir, prefix, years, dims, awc, lat = argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]
yrs = years.split('-')
yrs = range(int(yrs[0]), int(yrs[1])+1)
yms = ['%d%02d' % (y,m) for y in yrs for m in range(1,13)]
prdata = [ '%s/%s.prcp.%s.flt32' % (climate, prefix, ym) for ym in yms]
tadata = [ '%s/%s.tmax.%s.flt32' % (climate, prefix, ym) for ym in yms]
dims = dims.split('x')
dims = [int(dims[0]), int(dims[1]), int(dims[2])]

vd = 'vd/'
if not os.path.exists(vd): os.mkdir(vd)
bin = os.path.dirname(__file__) + '/pdsi'
print prdata, tadata, len(prdata), len(tadata)
pdsi_wrapper(vd, odir, dims, yrs, bin, prdata, tadata, awc, lat)
