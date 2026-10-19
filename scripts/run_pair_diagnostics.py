"""Reproducible synthetic benchmark; no SDSS/DESI data are loaded."""
import csv
import json
from pathlib import Path
import numpy as np
from research.pair_diagnostics import landy_szalay, spatial_jackknife, random_catalogue_sensitivity


def main(output: Path):
    output.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261019)
    data=rng.uniform(0,100,(160,3))
    random=rng.uniform(0,100,(500,3))
    edges=np.array([0,20,35,50,65,85,110,180])
    result=landy_szalay(data,random,edges)
    replicates,cov=spatial_jackknife(data,random,edges,box_size=100)
    sensitivity=random_catalogue_sensitivity(data,edges,box_size=100,random_sizes=(250,500),repeats=6)
    with (output/'pair_diagnostics.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['r_lo_Mpc_h','r_hi_Mpc_h','DD','DR','RR','xi','jackknife_sigma','random250_sigma','random500_sigma'])
        for i in range(len(result.xi)):
            writer.writerow([edges[i],edges[i+1],int(result.dd[i]),int(result.dr[i]),int(result.rr[i]),
                             float(result.xi[i]),float(np.sqrt(cov[i,i])),float(sensitivity[250]['std'][i]),float(sensitivity[500]['std'][i])])
    manifest={'synthetic':True,'seed':20261019,'data_count':len(data),'random_count':len(random),
              'jackknife_regions':len(replicates),'covariance':cov.tolist(),
              'warnings':['Illustrative uniform cubic window, not a measured survey correlation function.',
                          'Random catalogue dispersion excludes sample variance and survey-systematic uncertainty.',
                          'Small RR and limited jackknife regions can produce noisy or singular covariance.']}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


if __name__=='__main__':
    main(Path('results/pair_diagnostics'))
