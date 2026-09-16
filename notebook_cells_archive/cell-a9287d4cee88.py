exec(open('/workspace/fig5_v5.py').read())
import shutil
for ext in ['png', 'svg']:
    shutil.copy(f'/workspace/P2-02_Figure5_mechanism_hypothesis_v5.{ext}',
                f'/mnt/results/P2-02_Figure5_mechanism_hypothesis_v5.{ext}')
print('copied to /mnt/results')