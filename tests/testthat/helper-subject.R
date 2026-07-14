# A single synthetic subject used across tests. Expected index values were
# computed independently in Python from the GRADE/DPP data dictionary
# (see tools/verify.py in the source repo) and are hard-coded here so the R
# implementation is checked against an external reference.

subj <- list(
  # raw inputs
  i000 = 12, i015 = 70, i030 = 90, i060 = 85, i090 = 70, i120 = 55,   # uU/mL
  g000 = 100, g015 = 155, g030 = 160, g060 = 150, g090 = 140, g120 = 130, # mg/dL
  cp000 = 0.6, cp015 = 1.0, cp030 = 1.5, cp060 = 1.8, cp090 = 1.6, cp120 = 1.4, # nmol/L
  glcg000 = 15, glcg030 = 9,           # pmol/L
  trig = 150, chdl = 45,               # mg/dL
  bmi = 30, weight = 85, waist = 100, hip = 105, height = 170,
  age = 55, htn = 1, hba1c = 7.0,
  adipon = 8, leptin = 20, pin = 15, male = 1
)
subj$whr  <- subj$waist / subj$hip
subj$whtr <- subj$waist / subj$height

# derived unit conversions (independent arithmetic)
g <- function(x) x / 18                 # mg/dL -> mmol/L glucose
subj$g000m <- g(subj$g000); subj$g030m <- g(subj$g030); subj$g060m <- g(subj$g060)
subj$g090m <- g(subj$g090); subj$g120m <- g(subj$g120); subj$g015m <- g(subj$g015)
subj$trig_m <- subj$trig / 88.57
subj$chdl_m <- subj$chdl / 38.67
subj$i000p <- subj$i000 * 6.945; subj$i030p <- subj$i030 * 6.945
subj$i090p <- subj$i090 * 6.945; subj$i120p <- subj$i120 * 6.945

# externally computed expected values
EXP <- list(
  raynaud = 3.33333, homa_ir = 2.96296, isi = 0.3375, firi = 2.66667,
  quicki = 0.324762, igr = 2.16, gir = 0.462963, isi_bennett = 0.23468,
  mcauley = 5.87617, ohkura = 6, homa_ad = 0.37037, belfiore_basal = 0.837785,
  egdr_bmi = 4.89, egdr_wc = 4.894, egdr_whr = 5.3919,
  tg_hdl = 3.33333, tyg = 8.92266, tyg_bmi = 267.68, mets_ir = 46.1659,
  mets_vf = 7.30051, spise = 4.70315, vai = 1.92649, lap = 59.2751, lar = 2.5,
  matsuda = 2.91094, matsuda_30 = 3.5453, matsuda_60 = 2.90598,
  gutt_isi0120 = 2.38565, cederholm = 40.01,
  avignon_sib = 117.647, avignon_si2h = 19.745, avignon_sim = 17.9313,
  stumvoll_isi = 0.0761959, stumvoll_dem = 0.0691341, stumvoll_modi = 0.072643,
  stumvoll_mcr = 6.58373, ogis = 355.82, belfiore_ogtt = 0.748607,
  homa_b = 116.757, igi_30 = 23.4, igi_120 = 25.8, cir30 = 0.00625,
  kadowaki = 0.4875, stumvoll_1st = 1507.69, stumvoll_2nd = 392.46,
  bigtt_30_120 = 1544.13, bigtt_60_120 = 1887.74,
  cpep_gluc_0 = 0.108, cpep_gluc_30 = 0.16875, cpi_30 = 0.27,
  icpr = 0.1389, cpep_ins_molar_auc = 3.05038,
  auc_gluc_0_30 = 3900, auc_ins_0_30 = 1530, auc_ins_gluc_0_30 = 0.392308,
  auc_gluc_0_120 = 16950, auc_ins_gluc_0_120 = 0.49292,
  iauc_gluc_0_30 = 900, iauc_ins_0_30 = 1170, iauc_ins_gluc_0_30 = 1.3,
  odi = 1.95, odi_matsuda = 68.1159, odi_matsuda_30 = 82.9599,
  odi_quicki = 7.59942, odi_isi = 7.8975, odi_cp = 0.785953, issi2 = 1.43486,
  glucagon_ins_ratio = 0.179986, glucagon_suppr_abs = 6, glucagon_suppr_pct = 40,
  proinsulin_insulin_ratio = 0.179986
)

TOL <- 1e-4
