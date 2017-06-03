import numpy as np

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

files = ['data/powtable_DMONLY_WMAP7','data/powtable_AGN_WMAP7']
for f in files:
    data = np.loadtxt(f+'.dat')

    z = data[:,0]
    Pk = data[:,2]
    k = data[:,1]
    for i in range(len(z)):
        if z[i]>z[0]:
            nrow = i
            break
    z = np.unique(z)
    k = np.unique(k)
    Pk.shape = (Pk.size//nrow,nrow)
    np.savetxt(f+'_ready.dat',np.hstack([k.reshape(len(k),1),Pk.T]))
    
    line = '####\t'+'\t'.join([str(x) for x in z])
    line_prepender(f+'_ready.dat',line)
