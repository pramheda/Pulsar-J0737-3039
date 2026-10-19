"""Digitize published vector paths and test geometric consistency, not author code."""
from pathlib import Path
import sys,os,re,json,hashlib,xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];CAMPAIGN=ROOT.parent
sys.path.insert(0,str(CAMPAIGN/'trajectory/code'))
from galactic_trace import sky,integrate,coordinates,TIME_MYR,BIN,timing_binary,CLUSTERS
import numpy as np
from scipy.optimize import minimize_scalar,least_squares

def extract(path):
    root=ET.parse(path).getroot();lines=[];grid=[]
    number=r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?'
    for e in root.iter():
        if not e.tag.endswith('path') or e.get('fill')!='none':continue
        d=e.get('d','')
        if any(c in d for c in 'CQZ'):continue
        points=np.array([float(x) for x in re.findall(number,d)]).reshape(-1,2)
        if e.get('stroke')=='rgb(69.018555%, 69.018555%, 69.018555%)' and len(points)==2:
            grid.append(points)
        if e.get('stroke') in ['rgb(100%, 0%, 0%)','rgb(0%, 0%, 0%)','rgb(0%, 0%, 100%)'] and len(points)>10:
            lines.append((e.get('stroke'),points))
    out={};cal=[]
    for panel in range(2):
        gr=[g for g in grid if (g[:,0].mean()>450)==bool(panel)]
        xs=sorted(g[0,0] for g in gr if g[0,0]==g[1,0])
        ys=sorted(g[0,1] for g in gr if g[0,1]==g[1,1])
        assert len(xs)==5 and len(ys)==5
        cx=np.polyfit(xs,np.arange(5)*.04+7.92,1)
        cy=np.polyfit(ys,np.arange(5)*(.025 if panel else .15)+(-.025 if panel else -2.55),1)
        cal.append(dict(panel=panel,x_ticks=xs,y_ticks=ys,x_transform=cx.tolist(),y_transform=cy.tolist()))
        for color,points in lines:
            if (points[:,0].mean()>450)!=bool(panel):continue
            name={'rgb(100%, 0%, 0%)':'binary','rgb(0%, 0%, 0%)':'Theia_58','rgb(0%, 0%, 100%)':'OC_0450'}[color]
            out.setdefault(name,{})['xz' if panel else 'xy']=np.stack([np.polyval(cx,points[:,0]),np.polyval(cy,points[:,1])],axis=1).tolist()
    return dict(source=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                calibration=cal,paths=out,method='Original vector centerline vertices calibrated from five visually verified tick labels per axis; no raster tracing.',
                limitations=['The points describe a published plotted curve, not independent observations or recovered inputs.','No numerical uncertainty model can turn this comparison into a formal likelihood ratio.','Path coordinates carry PDF/SVG quantization; graph labels and shared point order are checked explicitly.'])

def check(data):
    lookback=np.linspace(0,15,1501);tt=-lookback
    paths=data['paths'];output=[]
    for name,path in paths.items():
        xy=np.array(path['xy']);xz=np.array(path['xz'])
        assert xy.shape==xz.shape
        xgap=np.max(abs(xy[:,0]-xz[:,0]));assert xgap<1e-5
        target=np.stack([(xy[:,0]+xz[:,0])/2,xy[:,1],xz[:,1]],axis=1)
        if name!='binary':
            c=next(c for c in CLUSTERS if c['name']==name)
            coord=sky(c['ra_deg'],c['dec_deg'],c['distance_pc']/1000,c['pm_ra_cosdec_masyr'],c['pm_dec_masyr'],c['heliocentric_radial_velocity_kms'])
            xyz=coordinates(integrate(coord,tt),tt/TIME_MYR)
            interp=np.stack([np.interp(target[:,0],xyz[::-1,0],xyz[::-1,j]) for j in range(3)],axis=1)
            residual=(interp-target)*1000
            output.append(dict(name=name,max_panel_x_discrepancy_pc=xgap*1000,rms_distance_residual_pc=float(np.sqrt(np.mean(np.sum(residual**2,axis=1)))),coordinate_rms_pc=np.sqrt(np.mean(residual**2,axis=0)).tolist()))
            continue
        for branch,b in [('printed',BIN),('timing',timing_binary())]:
            def residual(params):
                d,rv=params
                xyz=coordinates(integrate(sky(b['ra_deg'],b['dec_deg'],d,b['pmra'],b['pmdec'],rv),tt),tt/TIME_MYR)
                if np.any(np.diff(xyz[:,0])>=0):raise ValueError('Non-monotonic x')
                interp=np.stack([np.interp(target[:,0],xyz[::-1,0],xyz[::-1,j]) for j in range(3)],axis=1)
                return ((interp-target)*1000).ravel()
            nominal=residual([.74,66.25]).reshape(-1,3)
            opt=minimize_scalar(lambda d:np.sum(residual([d,66.25])**2),bounds=(.38,1.1),method='bounded',options={'xatol':1e-9})
            fit2=least_squares(residual,[opt.x,66.25],bounds=([.38,-200],[1.1,200]),xtol=1e-10,ftol=1e-10,gtol=1e-10)
            output.append(dict(name=name,branch=branch,max_panel_x_discrepancy_pc=xgap*1000,
                nominal_rms_distance_pc=float(np.sqrt(np.mean(np.sum(nominal**2,axis=1)))),
                fixed_published_rv_kms=66.25,best_distance_kpc=float(opt.x),distance_fit_rms_pc=float(np.sqrt(opt.fun/len(target))),
                free_distance_radial_fit=fit2.x.tolist(),free_fit_rms_pc=float(np.sqrt(np.sum(fit2.fun**2)/len(target))),
                nominal_residual_pc=nominal.tolist()))
    return output

if __name__=='__main__':
    data=extract(ROOT/'sources/j0737_oc_orbits.svg')
    data['comparison']=check(data)
    (ROOT/'results/vector_figure.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data['comparison'],indent=2))
