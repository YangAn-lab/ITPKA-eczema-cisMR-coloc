# 分析E: 重画图4 v2（修正布局后）
%run /workspace/fig4_redraw.py
import shutil
shutil.copy('/workspace/P2-02_Figure4_tissue_localization.png', '/mnt/results/P2-02_Figure4_tissue_localization.png')
shutil.copy('/workspace/P2-02_Figure4_tissue_localization.svg', '/mnt/results/P2-02_Figure4_tissue_localization.svg')
print('copied to /mnt/results/')