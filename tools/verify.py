import math

# ---- synthetic subject (units as in the data dictionary) --------------------
i000, i015, i030, i060, i090, i120 = 12., 70., 90., 85., 70., 55.   # uU/mL
g000, g015, g030, g060, g090, g120 = 100., 155., 160., 150., 140., 130.  # mg/dL
cp000, cp015, cp030, cp060, cp090, cp120 = .6, 1.0, 1.5, 1.8, 1.6, 1.4  # nmol/L
glcg000, glcg030 = 15., 9.       # pmol/L
trig, chdl = 150., 45.           # mg/dL
bmi, weight, waist, hip, height = 30., 85., 100., 105., 170.  # kg/m2,kg,cm,cm,cm
whr = waist/hip; whtr = waist/height
age, htn, hba1c = 55., 1., 7.0
adipon, leptin, pin = 8., 20., 15.
male = 1

# conversions
def g_mmol(x): return x/18.0
def tg_mmol(x): return x/88.57
def hdl_mmol(x): return x/38.67
def i_pmol(x): return x*6.945
g000m,g015m,g030m,g060m,g090m,g120m = map(g_mmol,(g000,g015,g030,g060,g090,g120))
trig_m, chdl_m = tg_mmol(trig), hdl_mmol(chdl)
i000p,i030p,i090p,i120p = map(i_pmol,(i000,i030,i090,i120))

def auc(times, vals):
    return sum((times[k+1]-times[k])*(vals[k]+vals[k+1])/2 for k in range(len(times)-1))
def iauc(times, vals):
    base=vals[0]; return auc(times,[v-base for v in vals])

R={}
# fasting sensitivity
R['raynaud']=40/i000
R['homa_ir']=i000*g000m/22.5
R['isi']=1/R['homa_ir']
R['firi']=i000*g000m/25
R['quicki']=1/(math.log10(i000)+math.log10(g000))
R['igr']=i000/g000m
R['gir']=g000m/i000
R['isi_bennett']=1/(math.log(i000)*math.log(g000m))
R['mcauley']=math.exp(2.63-0.28*math.log(i000)-0.31*math.log(trig_m))
R['ohkura']=20/(cp000*g000m)
R['homa_ad']=i000*g000m/(22.5*adipon)
R['belfiore_basal']=2/(((i000/9.46)*(g000m/5.08))+1)
R['egdr_bmi']=19.02-0.22*bmi-3.26*htn-0.61*hba1c
R['egdr_wc']=21.158-0.09*waist-3.407*htn-0.551*hba1c
R['egdr_whr']=24.31-12.22*whr-3.29*htn-0.57*hba1c
# lipid/anthro
R['tg_hdl']=trig/chdl
R['tyg']=math.log(trig*g000/2)
R['tyg_bmi']=R['tyg']*bmi
R['mets_ir']=(math.log(2*g000+trig)*bmi)/math.log(chdl)
R['mets_vf']=4.466+0.011*math.log(R['mets_ir'])**3+3.239*math.log(whtr)**3+0.319*male+0.594*math.log(age)
R['spise']=(600*chdl**0.185)/(trig**0.2*bmi**1.338)
R['vai']=(waist/(39.68+1.88*bmi))*(trig_m/1.03)*(1.31/chdl_m)
R['lap']=(waist-65)*trig_m
R['lar']=leptin/adipon
# OGTT sensitivity
gmean=auc([0,30,60,90,120],[g000,g030,g060,g090,g120])/120
imean=auc([0,30,60,90,120],[i000,i030,i060,i090,i120])/120
gmean_mmol=gmean/18.0
R['gmean']=gmean; R['imean']=imean
R['matsuda']=10000/math.sqrt(g000*i000*gmean*imean)
R['matsuda_30']=10000/math.sqrt(g000*i000*(auc([0,30],[g000,g030])/30)*(auc([0,30],[i000,i030])/30))
R['matsuda_60']=10000/math.sqrt(g000*i000*(auc([0,30,60],[g000,g030,g060])/60)*(auc([0,30,60],[i000,i030,i060])/60))
R['gutt_isi0120']=(75000+(g000-g120)*0.19*weight)/(120*gmean*math.log10(imean))
R['cederholm']=(75000+(g000m-g120m)*1.15*180*0.19*weight)/(120*gmean_mmol*math.log10(imean))
R['avignon_sib']=1e8/(i000*g000m*(150*weight))
R['avignon_si2h']=1e8/(i120*g120m*(150*weight))
R['avignon_sim']=(0.137*R['avignon_sib']+R['avignon_si2h'])/2
R['stumvoll_isi']=0.226-0.0032*bmi-0.0000645*i120p-0.00375*g090m
R['stumvoll_dem']=0.222-0.00333*bmi-0.0000779*i120p-0.000422*age
R['stumvoll_modi']=0.156-0.0000459*i120p-0.000321*i000p-0.00541*g120m
R['stumvoll_mcr']=18.8-0.271*bmi-0.0052*i120p-0.27*g090m
# OGIS
P1,P2,P3,P4,P5,P6=6.5,1951,4514,792,0.0118,173; V=10000; Gcl=5; T=30; dose=75
bsa=0.1640443958298*weight**0.515*(0.01*height)**0.422
d0=5.551*dose/bsa
cl=P4*((P1*d0-V*(g120m-g090m)/T)/g090m+P3/g000m)/(i090p-i000p+P2)
b=(P5*(g090m-Gcl)+1)*cl
R['ogis']=(b+math.sqrt(b*b+4*P5*P6*(g090m-Gcl)*cl))/2
# Belfiore ogtt
iarea=0.5*i000+i060+0.5*i120; garea=0.5*g000m+g060m+0.5*g120m
R['belfiore_ogtt']=2/(((iarea/91.87)*(garea/11.36))+1)
# secretion
R['homa_b']=20*i000/(g000m-3.5)
R['igi_30']=(i030-i000)/(g030m-g000m)
R['igi_120']=(i120-i000)/(g120m-g000m)
R['cir30']=i030/(g030*(g030-70)) if g030>70 else float('nan')
R['kadowaki']=(i030-i000)/g030
R['stumvoll_1st']=1283+1.829*i030p-138.7*g030m+3.772*i000p
R['stumvoll_2nd']=287+0.416*i030p-26.07*g030m+0.926*i000p
R['bigtt_30_120']=math.exp(8.20+0.00178*i000+0.00168*i030-0.000383*i120-0.314*g000m-0.109*g030m+0.0781*g120m+0.180*male+0.032*bmi)
R['bigtt_60_120']=math.exp(8.19+0.00339*i000+0.00152*i060-0.000959*i120-0.389*g000m-0.142*g060m+0.164*g120m+0.256*male+0.038*bmi)
# c-peptide
R['cpep_gluc_0']=cp000/g000m
R['cpep_gluc_30']=cp030/g030m
R['cpi_30']=(cp030-cp000)/(g030m-g000m)
R['icpr']=i000p/(cp000*1000)   # corrected molar ratio (script *1000 is a 1e6 slip)
aCP120=auc([0,30,60,90,120],[cp000,cp030,cp060,cp090,cp120])
aI120=auc([0,30,60,90,120],[i000,i030,i060,i090,i120])
R['cpep_ins_molar_auc']=(aCP120*1000)/(aI120*6.945)
# AUC + ratios
aG30=auc([0,30],[g000,g030]); aI30=auc([0,30],[i000,i030])
R['auc_gluc_0_30']=aG30; R['auc_ins_0_30']=aI30
R['auc_ins_gluc_0_30']=aI30/aG30
R['auc_gluc_0_120']=auc([0,30,60,90,120],[g000,g030,g060,g090,g120])
R['auc_ins_gluc_0_120']=aI120/R['auc_gluc_0_120']
R['iauc_gluc_0_30']=iauc([0,30],[g000,g030])
R['iauc_ins_0_30']=iauc([0,30],[i000,i030])
R['iauc_ins_gluc_0_30']=R['iauc_ins_0_30']/R['iauc_gluc_0_30'] if R['iauc_gluc_0_30']>0 else float('nan')
# disposition
R['odi']=R['igi_30']*(1/i000)
R['odi_matsuda']=R['igi_30']*R['matsuda']
R['odi_matsuda_30']=R['igi_30']*R['matsuda_30']
R['odi_quicki']=R['igi_30']*R['quicki']
R['odi_isi']=R['igi_30']*R['isi']
R['odi_cp']=R['cpi_30']*R['matsuda']
R['issi2']=R['auc_ins_gluc_0_120']*R['matsuda']
# glucagon / proinsulin
R['glucagon_ins_ratio']=glcg000/i000p
R['glucagon_suppr_abs']=glcg000-glcg030
R['glucagon_suppr_pct']=(glcg000-glcg030)/glcg000*100
R['proinsulin_insulin_ratio']=pin/i000p

for k in sorted(R): print(f"{k} = {R[k]:.6g}")
