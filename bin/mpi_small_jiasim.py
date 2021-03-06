import numpy as np
import sys, os, glob
from orphics.analysis.pipeline import mpi_distribute, MPIStats
import orphics.tools.stats as stats
import alhazen.io as aio
import orphics.tools.io as io
import orphics.analysis.flatMaps as fmaps
import warnings
import logging
logger = logging.getLogger()
with io.nostdout():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from enlib import enmap, lensing, resample
from alhazen.quadraticEstimator import Estimator
import alhazen.lensTools as lt
from ConfigParser import SafeConfigParser 
import enlib.fft as fftfast
import argparse
from mpi4py import MPI

# Add LiuConvergence
import astropy.io.fits as fits
class LiuConvergence(object):

    def __init__(self,root_dir="/global/cscratch1/sd/djbard/JiaSims/convergence_maps/Maps11000/"):
        self.root = root_dir
        size_deg = 3.5
        Npix = 2048.
        px = size_deg*60./Npix
        self.px = px
        self.shape, self.wcs = fmaps.rect_geometry(width_deg = size_deg,px_res_arcmin=px,proj="CAR",pol=False)
        self.modlmap = enmap.modlmap(self.shape,self.wcs)
        #print ('JIA print self.shape',self.shape)
        
    def get_kappa(self,index,z=1100):
        zstr = "{:.2f}".format(z)
        kappa_file = self.root+"Maps%02d"%(z*10)+"/WLconv_z"+zstr+"_"+str(index).zfill(4)+"r.fits"
        
        my_map = fits.open(kappa_file)[0]
        my_map = my_map.data
        #print ('JIA print my_map.shape',my_map.shape)

        assert my_map.shape == self.shape
        low_pass_ell = 40000
        retmap = enmap.ndmap(my_map,self.wcs)
        kmask = fmaps.mask_kspace(self.shape,self.wcs,lmax=low_pass_ell)
        retmap = enmap.ndmap(fmaps.filter_map_new(retmap,kmask),self.wcs)
        #print ('JIA print retmap.shape', retmap.shape)
        return retmap


lc = LiuConvergence(root_dir="/global/cscratch1/sd/djbard/JiaSims/convergence_maps/")
# forget the ini, let's get shape and wcs from Jia's sims themselves
shape,wcs = lc.shape,lc.wcs
shape_sim,wcs_sim = shape,wcs
shape_dat,wcs_dat = shape,wcs

# Runtime params that should be moved to command line
cosmology_section = "cc_nam" # this should ideally by Jia's sims cosmology


# Parse command line
parser = argparse.ArgumentParser(description='Verify lensing reconstruction.')
parser.add_argument("Exp", type=str,help='Experiment name.')
parser.add_argument("-N", "--nsim",     type=int,  default=None)
args = parser.parse_args()
Ntot = args.nsim


expf_name = args.Exp #"experiment_small"

# Get MPI comm
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numcores = comm.Get_size()    

gradCut = 2000 #None #2000

# i/o directories
#out_dir = os.environ['WWW']+"plots/"+expf_name+"_"+str(gradCut)+"_"  # for plots
out_dir = "output/jia_"+expf_name+"_"+str(gradCut)+"_"  # for plots



# Efficiently distribute sims over MPI cores
num_each,each_tasks = mpi_distribute(Ntot,numcores)
# Initialize a container for stats and stacks
mpibox = MPIStats(comm,num_each,tag_start=333)

if rank==0: print "At most ", max(num_each) , " tasks..."

# What am I doing?
my_tasks = each_tasks[rank]


# Read config
iniFile = "input/recon.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)

pol = False
#shape_sim, wcs_sim, shape_dat, wcs_dat = aio.enmaps_from_config(Config,sim_section,analysis_section,pol=pol)
#analysis_resolution =  np.min(enmap.extent(shape_dat,wcs_dat)/shape_dat[-2:])*60.*180./np.pi
#print 'analysis resolution',analysis_resolution
min_ell = fmaps.minimum_ell(shape_dat,wcs_dat)
lb = aio.ellbounds_from_config(Config,"reconstruction_small",min_ell)
tellmin = lb['tellminY']
tellmax = lb['tellmaxY']
pellmin = lb['pellminY']
pellmax = lb['pellmaxY']
kellmin = lb['kellmin']
kellmax = lb['kellmax']
parray_dat = aio.patch_array_from_config(Config,expf_name,shape_dat,wcs_dat,dimensionless=True)
parray_sim = aio.patch_array_from_config(Config,expf_name,shape_sim,wcs_sim,dimensionless=True)
lxmap_dat,lymap_dat,modlmap_dat,angmap_dat,lx_dat,ly_dat = fmaps.get_ft_attributes_enmap(shape_dat,wcs_dat)
lxmap_sim,lymap_sim,modlmap_sim,angmap_sim,lx_sim,ly_sim = fmaps.get_ft_attributes_enmap(shape_sim,wcs_sim)
lbin_edges = np.arange(kellmin,kellmax,300)
lbinner_dat = stats.bin2D(modlmap_dat,lbin_edges)
lbinner_sim = stats.bin2D(modlmap_sim,lbin_edges)

# === COSMOLOGY ===
theory, cc, lmax = aio.theory_from_config(Config,cosmology_section)

# switch to Jia's theory kk 
Clkk = np.loadtxt('output/inputkk.csv')
theory.loadGenericCls(Clkk[:,0],Clkk[:,1],'kk',lpad=41000)

parray_dat.add_theory(None,theory,lmax)
template_dat = fmaps.simple_flipper_template_from_enmap(shape_dat,wcs_dat)
nT = parray_dat.nT
nP = parray_dat.nP
if rank==0: io.quickPlot2d(nT,out_dir+"nt.png")
kbeam_dat = parray_dat.lbeam
kbeampass = kbeam_dat
if rank==0: io.quickPlot2d(kbeampass,out_dir+"kbeam.png")
fMaskCMB_T = fmaps.fourierMask(lx_dat,ly_dat,modlmap_dat,lmin=tellmin,lmax=tellmax)
fMaskCMB_P = fmaps.fourierMask(lx_dat,ly_dat,modlmap_dat,lmin=pellmin,lmax=pellmax)
fMask = fmaps.fourierMask(lx_dat,ly_dat,modlmap_dat,lmin=kellmin,lmax=kellmax)

with io.nostdout():
    qest = Estimator(template_dat,
                     theory,
                     theorySpectraForNorm=None,
                     noiseX2dTEB=[nT,nP,nP],
                     noiseY2dTEB=[nT,nP,nP],
                     fmaskX2dTEB=[fMaskCMB_T,fMaskCMB_P,fMaskCMB_P],
                     fmaskY2dTEB=[fMaskCMB_T,fMaskCMB_P,fMaskCMB_P],
                     fmaskKappa=fMask,
                     kBeamX = kbeampass,
                     kBeamY = kbeampass,
                     doCurl=False,
                     TOnly=not(pol),
                     halo=True,
                     uEqualsL=True,
                     gradCut=gradCut,verbose=False,
                     bigell=lmax)

    
#pixratio = analysis_resolution/Config.getfloat(sim_section,"pixel_arcmin")
#px_dat = analysis_resolution
lens_order = 5 #Config.getint(sim_section,"lens_order")
parray_sim = aio.patch_array_from_config(Config,expf_name,shape_sim,wcs_sim,dimensionless=True)
parray_sim.add_theory(None,theory,lmax)


k = -1
for index in my_tasks:
    kappa = lc.get_kappa(index+1,z=1100)
    #kappa = parray_sim.get_grf_kappa(seed=index+1000000)

    phi, fphi = lt.kappa_to_phi(kappa,parray_sim.modlmap,return_fphi=True)
    #alpha_pix = enmap.grad_pixf(fphi)
    grad_phi = enmap.grad(phi)
            
    if rank==0: print "Generating unlensed CMB..."
    unlensed = parray_sim.get_unlensed_cmb(seed=index)
    if rank==0: print "Lensing..."
    #lensed = lensing.lens_map_flat_pix(unlensed.copy(), alpha_pix.copy(),order=lens_order)
    lensed = lensing.lens_map(unlensed.copy(), grad_phi, order=lens_order, mode="spline", border="cyclic", trans=False, deriv=False, h=1e-7)

    
    # === ADD NOISE BEFORE DOWNSAMPLE
    # if rank==0: print "Beam convolving..."
    # olensed = enmap.ndmap(lensed.copy() if abs(pixratio-1.)<1.e-3 else resample.resample_fft(lensed.copy(),shape_dat),wcs_dat)
    # flensed = fftfast.fft(lensed,axes=[-2,-1])
    # flensed *= parray_sim.lbeam
    # lensed = fftfast.ifft(flensed,axes=[-2,-1],normalize=True).real
    # if rank==0: print "Adding noise..."    
    # noise = parray_sim.get_noise_sim(seed=index+20000)
    # lensed += noise    
    # if rank==0: print "Downsampling..."
    # cmb = lensed if abs(pixratio-1.)<1.e-3 else resample.resample_fft(lensed,shape_dat)

    
    # === ADD NOISE AFTER DOWNSAMPLE
    if rank==0: print "Beam convolving..."
    olensed = enmap.ndmap(lensed.copy(),wcs_sim) #if abs(pixratio-1.)<1.e-3 else resample.resample_fft(lensed.copy(),shape_dat),wcs_dat)
    flensed = fftfast.fft(olensed,axes=[-2,-1])
    flensed *= parray_dat.lbeam
    lensed = fftfast.ifft(flensed,axes=[-2,-1],normalize=True).real
    if rank==0: print "Adding noise..."    
    noise = parray_dat.get_noise_sim(seed=index+20000)

    lcents, noise1d = lbinner_dat.bin(fmaps.get_simple_power_enmap(noise))
    mpibox.add_to_stats('noisett',noise1d)        

    lensed += noise    
    if rank==0: print "Downsampling..."
    cmb = lensed

    
    
    cmb = enmap.ndmap(cmb,wcs_dat)
    if rank==0: print "Calculating powers for diagnostics..."
    utt2d = fmaps.get_simple_power_enmap(enmap.ndmap(unlensed,wcs_sim))
    ltt2d = fmaps.get_simple_power_enmap(olensed)
    ccents,utt = lbinner_dat.bin(utt2d)
    ccents,ltt = lbinner_dat.bin(ltt2d)
    mpibox.add_to_stats("ucl",utt)
    mpibox.add_to_stats("lcl",ltt)
            

    if rank==0: print "Reconstructing..."
    measured = cmb
    fkmaps = fftfast.fft(measured,axes=[-2,-1])
    qest.updateTEB_X(fkmaps,alreadyFTed=True)
    qest.updateTEB_Y()
    with io.nostdout():
        rawkappa = qest.getKappa("TT").real



        
    kappa_recon = enmap.ndmap(rawkappa,wcs_dat)
    apower = fmaps.get_simple_power_enmap(enmap1=kappa_recon)

    
    data_power_2d_TT = fmaps.get_simple_power_enmap(measured)
    sd = qest.N.super_dumb_N0_TTTT(data_power_2d_TT)
    lcents,sdp = lbinner_dat.bin(sd)
    
    mpibox.add_to_stats("superdumbs",sdp)
    n0subbed = apower - sd
    lcents,rclkk = lbinner_dat.bin(n0subbed)
    mpibox.add_to_stats("auto_n0subbed",rclkk)


    
    if rank==0: print "Downsampling input kappa..."
    downk = enmap.ndmap(kappa,wcs_dat)
    if rank==0: print "Calculating kappa powers and binning..."
    cpower = fmaps.get_simple_power_enmap(enmap1=kappa_recon,enmap2=downk)
    ipower = fmaps.get_simple_power_enmap(enmap1=downk)
    lcents, cclkk = lbinner_dat.bin(cpower)
    lcents, aclkk = lbinner_dat.bin(apower)
    lcents, iclkk = lbinner_dat.bin(ipower)

    mpibox.add_to_stats("cross",cclkk)
    mpibox.add_to_stats("ipower",iclkk)
    mpibox.add_to_stats("auto",aclkk)

    if rank==0 and index==0:
        io.quickPlot2d(cmb,out_dir+"cmb.png")
        io.quickPlot2d(measured,out_dir+"mcmb.png")
        io.quickPlot2d(kappa,out_dir+"inpkappa.png")
        io.quickPlot2d(kappa_recon,out_dir+"reconkappa.png")


mpibox.get_stacks()
mpibox.get_stats()



if rank==0:



    cstats = mpibox.stats['cross']
    istats = mpibox.stats['ipower']
    astats = mpibox.stats['auto']
    rstats = mpibox.stats['auto_n0subbed']
    nstats = mpibox.stats['superdumbs']

    area = unlensed.area()*(180./np.pi)**2.
    print "area: ", area, " sq.deg."
    fsky = area/41250.
    print "fsky: ",fsky
    diag = np.sqrt(np.diagonal(astats['cov'])*lcents*np.diff(lbin_edges)*fsky)
    diagr = np.sqrt(np.diagonal(rstats['cov'])*lcents*np.diff(lbin_edges)*fsky)

    pl = io.Plotter(scaleY='log',scaleX='log')
    pl.addErr(lcents,cstats['mean'],yerr=cstats['errmean'],marker="o",label="recon x cross")
    pl.add(lcents,istats['mean'],marker="x",ls="none",label="input")
    pl.add(lcents,diag,ls="-.",lw=2,label="diag no n0sub")
    pl.add(lcents,diagr,ls="-.",lw=2,label="diag n0sub")
    lcents,nlkk = lbinner_dat.bin(qest.N.Nlkk['TT'])
    ellrange = np.arange(2,kellmax,1)
    clkk = theory.gCl("kk",ellrange)
    clkk2d = theory.gCl("kk",modlmap_dat)
    ccents,clkk1d = lbinner_dat.bin(clkk2d)
    
    pl.addErr(lcents,rstats['mean']-clkk1d,yerr=rstats['errmean'],marker="o",alpha=0.5,label="auto n0subbed - clkk")
    pl.addErr(lcents,astats['mean'],yerr=astats['errmean'],marker="o",alpha=0.5,label="raw")
    pl.addErr(lcents,rstats['mean'],yerr=rstats['errmean'],marker="o",alpha=0.5,label="auto n0subbed")
    pl.add(lcents,nlkk,ls="--",label="theory n0")
    pl.add(lcents,nstats['mean'],ls="--",label="superdumb n0")
    pl.add(lcents,nstats['mean']+clkk1d,ls="--",label="superdumb n0 + clkk")
    pl.add(ellrange,clkk,color="k")
    pl.legendOn(loc="lower left",labsize=9)
    pl._ax.set_xlim(30,1e5)
    pl._ax.set_ylim(1e-12,1e-5)
    pl.done(out_dir+"cpower.png")

    np.savetxt(out_dir+'ell_recon_raw_autoN0subbed_jia.csv',np.vstack([lcents,cstats['mean'],astats['mean'],rstats['mean']]).T)


    io.quickPlot2d(stats.cov2corr(astats['covmean']),out_dir+"corr.png")
    io.quickPlot2d(stats.cov2corr(rstats['covmean']),out_dir+"rcorr.png")


    np.save(out_dir+str(area)+"sqdeg_covmat_dl300.npy",rstats['cov'])
    np.save(out_dir+str(area)+"sqdeg_autocovmat_dl300.npy",astats['cov'])
    np.save(out_dir+str(area)+"sqdeg_lbin_edges_dl300.npy",lbin_edges)
    import cPickle as pickle
    pickle.dump((lcents,mpibox.stats['noisett']['mean']),open(out_dir+"noise_mpismall.pkl",'wb'))

    pl = io.Plotter()
    ldiff = (cstats['mean']-istats['mean'])*100./istats['mean']
    lerr = cstats['errmean']*100./istats['mean']
    pl.addErr(lcents,ldiff,yerr=lerr,marker="o",ls="-")
    pl._ax.axhline(y=0.,ls="--",color="k")
    pl.done(out_dir+"powerdiff.png")
    
    iutt2d = theory.uCl("TT",parray_dat.modlmap)
    iltt2d = theory.lCl("TT",parray_dat.modlmap)
    ccents,iutt = lbinner_dat.bin(iutt2d)
    ccents,iltt = lbinner_dat.bin(iltt2d)
    uclstats = mpibox.stats["ucl"]
    lclstats = mpibox.stats["lcl"]

    utt = uclstats['mean']
    ltt = lclstats['mean']
    utterr = uclstats['errmean']
    ltterr = lclstats['errmean']


    pl = io.Plotter()




    pdiff = (utt-iutt)*100./iutt
    perr = 100.*utterr/iutt

    pl.addErr(ccents+25,pdiff,yerr=perr,marker="x",ls="none",label="unlensed")

    pdiff = (ltt-iltt)*100./iltt
    perr = 100.*ltterr/iltt

    pl.addErr(ccents+50,pdiff,yerr=perr,marker="o",ls="none",label="lensed")
    pl.legendOn(labsize=10,loc="lower left")
    pl._ax.axhline(y=0.,ls="--",color="k")
    pl._ax.set_ylim(-5.,5.)
    pl.done(out_dir+"clttpdiff.png")



    pl = io.Plotter(scaleY='log',scaleX='log')

    pl.add(ccents,iutt*ccents**2.)
    pl.addErr(ccents,utt*ccents**2.,yerr=utterr*ccents**2.,marker="x",ls="none",label="unlensed")

    pl.add(ccents,iltt*ccents**2.)
    pl.addErr(ccents,ltt*ccents**2.,yerr=ltterr*ccents**2.,marker="o",ls="none",label="lensed")

    pl.legendOn(labsize=10)
    pl.done(out_dir+"clttp.png")
