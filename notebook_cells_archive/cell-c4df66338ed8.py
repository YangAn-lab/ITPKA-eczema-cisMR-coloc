import importlib, subprocess, sys
try:
    import cellxgene_census
    print('cellxgene_census', cellxgene_census.__version__)
except ImportError:
    print('installing...')
    subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'cellxgene-census'], check=False)
    import cellxgene_census
    print('installed', cellxgene_census.__version__)