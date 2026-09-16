rows = [
 # 类别, 数据集, 组织, 位点, 基因, NES/beta, p, 备注
 ["位点身份","GTEx v8 (hg38)","-","rs11635906 = chr15:41487062 A/G","-","-","-",
  "b37=15:41779260 A/G；位于ITPKA转录起始位点上游6,331 bp（两种build一致）；移交文档“ITPKA基因内”表述需更正为“启动子近端上游”"],
 ["皮肤eQTL验证","GTEx v8","Skin_Sun_Exposed_Lower_leg (n=605)","rs11635906","ITPKA","-","不显著",
  "ITPKA在该组织无任何显著cis-eQTL（显著对列表为空）"],
 ["皮肤eQTL验证","GTEx v8","Skin_Not_Sun_Exposed_Suprapubic (n=517)","rs11635906","ITPKA","-","不显著","同上"],
 ["皮肤eQTL验证","GTEx v10","两个皮肤组织","rs11635906","ITPKA","-","不显著","v10中同样无ITPKA显著皮肤eQTL"],
 ["皮肤中该变异的其他靶基因","GTEx v8","Skin_Sun_Exposed_Lower_leg","rs11635906","OIP5-AS1","+0.384","2.6e-35","lncRNA；皮肤中该变异最强信号"],
 ["皮肤中该变异的其他靶基因","GTEx v8","Skin_Sun_Exposed_Lower_leg","rs11635906","CHP1","-0.118","3.0e-10","邻近基因"],
 ["皮肤中该变异的其他靶基因","GTEx v8","Skin_Sun_Exposed_Lower_leg","rs11635906","RPAP1","-0.143","1.3e-06","邻近基因"],
 ["皮肤中该变异的其他靶基因","GTEx v8","Skin_Not_Sun_Exposed_Suprapubic","rs11635906","OIP5-AS1","+0.380","8.5e-32",""],
 ["皮肤中该变异的其他靶基因","GTEx v8","Skin_Not_Sun_Exposed_Suprapubic","rs11635906","RPAP1","-0.150","1.0e-05",""],
 ["血液交叉验证（独立队列复现eQTLGen）","GTEx v8","Whole_Blood (n=670)","rs11635906","ITPKA","+0.216","6.8e-06","G等位上调ITPKA；与eQTLGen血液发现方向待与原始扫描结果核对（写作阶段）"],
 ["血液交叉验证","GTEx v10","Whole_Blood","rs11635906","ITPKA","+0.246","3.3e-08","v10中更强"],
 ["表达量背景","GTEx v8","Skin_Sun_Exposed_Lower_leg","-","ITPKA","-","-","中位表达0.828 TPM（高于全血0.155）——皮肤无eQTL非因不表达"],
 ["表达量背景","GTEx v8","Skin_Not_Sun_Exposed_Suprapubic","-","ITPKA","-","-","中位表达0.737 TPM"],
]
ev = pd.DataFrame(rows, columns=["类别","数据集","组织","位点","基因","NES","p值","备注"])
ev.to_csv("/mnt/results/P2-02_step1_GTEx皮肤eQTL验证_证据表.csv", index=False, encoding="utf-8-sig")
print(ev.to_string(index=False, max_colwidth=60))
print("\n已保存: P2-02_step1_GTEx皮肤eQTL验证_证据表.csv")