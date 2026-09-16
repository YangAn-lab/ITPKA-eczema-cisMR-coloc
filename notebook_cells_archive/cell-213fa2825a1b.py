import numpy as np
from scipy import stats

print("="*72)
print("P2-02 严谨性审计 A：移交文档数值自洽性（b / OR / 95%CI / p 互相重算）")
print("="*72)

def check(name, b, or_, lo, hi, p_stated):
    or_calc = np.exp(b)
    se_from_ci = (np.log(hi) - np.log(lo)) / (2*1.959964)
    z_from_ci = b / se_from_ci
    p_from_ci = 2*stats.norm.sf(abs(z_from_ci))
    z_from_p = stats.norm.isf(p_stated/2)
    se_from_p = b / z_from_p
    ci_from_p = (np.exp(b-1.959964*se_from_p), np.exp(b+1.959964*se_from_p))
    print(f"\n[{name}]")
    print(f"  声称: b={b}, OR={or_} ({lo}-{hi}), p={p_stated}")
    print(f"  exp(b)={or_calc:.3f}  vs 声称OR={or_}  -> {'OK' if abs(or_calc-or_)<0.011 else 'MISMATCH'}")
    print(f"  由CI反推: SE={se_from_ci:.4f}, z={z_from_ci:.2f}, p={p_from_ci:.2e}  (声称p={p_stated})")
    print(f"  由p反推:  SE={se_from_p:.4f}, 对应95%CI=({ci_from_p[0]:.3f},{ci_from_p[1]:.3f})  vs 声称({lo},{hi})")
    print(f"  判定: {'自洽（差异在舍入范围内）' if (abs(np.log10(p_from_ci)-np.log10(p_stated))<0.35 and abs(ci_from_p[0]-lo)<0.02 and abs(ci_from_p[1]-hi)<0.02) else '需复核'}")

check("湿疹端点 L12_DERMATITISECZEMA", 0.237, 1.27, 1.17, 1.38, 2.9e-08)
check("特应性皮炎端点 L12_ATOPIC", 0.323, 1.38, 1.22, 1.56, 1.4e-07)

print("\n" + "="*72)
print("多重检验与工具变量强度")
print("="*72)
n_tests = 498  # scan_sasp_skin.csv 行数-表头 = 498 对
bonf = 0.05/n_tests
print(f"Bonferroni阈值: 0.05/{n_tests} = {bonf:.2e}")
print(f"湿疹 p=2.9e-08  通过: {2.9e-08 < bonf}   AD p=1.4e-07 通过: {1.4e-07 < bonf}")
z_exp = stats.norm.isf(7.4e-29/2)
print(f"暴露侧(eQTLGen) p=7.4e-29 -> z={z_exp:.2f}, F≈z²={z_exp**2:.0f}  (>>10, 强工具变量)")
print(f"coloc阈值: PP.H4=0.962 > 0.75 通过; AD PP.H4=0.561 未过（移交文档如实标注为支持性）")
print(f"FinnGen R12 样本量核对: 67474+432874={67474+432874} (R12总样本≈500,348, 一致)")