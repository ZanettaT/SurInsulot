import re, numpy as np

src = open("/sessions/brave-dazzling-ramanujan/mnt/outputs/homa2calc_src/homa2calc-main/R/tables.R").read()

def get_vec(name):
    m = re.search(name + r"\s*<-\s*c\(([^)]*)\)", src)
    return [float(x) for x in m.group(1).replace("\n"," ").split(",") if x.strip()]

def get_mat(name):
    m = re.search(name + r"\s*<-\s*matrix\(c\((.*?)\),\s*nrow=21", src, re.S)
    nums = [float(x) for x in m.group(1).replace("\n"," ").split(",") if x.strip()]
    assert len(nums) == 441, (name, len(nums))
    return np.array(nums).reshape((21, 21), order="F")   # R fills column-major

TAB = {}
for mode in ["ins", "spec", "cpep"]:
    TAB[mode] = dict(
        g=get_vec(f"homa2_{mode}_gluc_axis"),
        h=get_vec(f"homa2_{mode}_horm_axis"),
        b=get_mat(f"homa2_{mode}_b_mat"),
        s=get_mat(f"homa2_{mode}_s_mat"),
    )

def interp(glucose, hormone, mode):
    ga, ha, bm, sm = TAB[mode]["g"], TAB[mode]["h"], TAB[mode]["b"], TAB[mode]["s"]
    if not (ga[0] <= glucose <= ga[-1] and ha[0] <= hormone <= ha[-1]):
        return (np.nan, np.nan, np.nan)
    g1 = max(i for i,v in enumerate(ga) if v <= glucose)
    g2 = min(i for i,v in enumerate(ga) if v >= glucose)
    h1 = max(i for i,v in enumerate(ha) if v <= hormone)
    h2 = min(i for i,v in enumerate(ha) if v >= hormone)
    tx = 0 if g2==g1 else (glucose-ga[g1])/(ga[g2]-ga[g1])
    ty = 0 if h2==h1 else (hormone-ha[h1])/(ha[h2]-ha[h1])
    def bl(m):
        return ((1-tx)*(1-ty)*m[g1,h1] + tx*(1-ty)*m[g2,h1]
                + (1-tx)*ty*m[g1,h2] + tx*ty*m[g2,h2])
    b, s = bl(bm), bl(sm)
    return (round(b,1), round(s,1), 100/s)

checks = [
    ("ins", 3.0, 20, 139.8, 308.9),
    ("ins", 5.2, 58, 93.0, 91.4),
    ("ins", 14.0, 96, 21.7, 42.7),
    ("ins", 5.0, 60, 131.9, 92.8),     # interpolated (README anchor)
    ("cpep", 3.0, 0.2, 151.8, 272.6),
    ("cpep", 14.0, 3.5, 97.2, 8.9),
    ("spec", 3.0, 20, 150.3, 276.9),
]
ok = True
for mode, g, h, eb, es in checks:
    b, s, ir = interp(g, h, mode)
    good = (abs(b-eb) < 0.05 and abs(s-es) < 0.05)
    ok = ok and good
    print(f"{mode:5s}({g},{h}) -> b={b} s={s} ir={ir:.4f}  expect b={eb} s={es}  {'OK' if good else 'MISMATCH'}")

# out of range
print("oor ins(2.0,50) ->", interp(2.0, 50, "ins"), "(expect nan)")
print("oor ins(5.0,500) ->", interp(5.0, 500, "ins"), "(expect nan)")
# also record the wrapper-relevant value: DPP uses insulin uU->pmol x6.945
print("ALL ANCHORS", "PASS" if ok else "FAIL")
