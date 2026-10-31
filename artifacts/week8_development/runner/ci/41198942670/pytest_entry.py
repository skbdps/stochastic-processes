import os,pathlib,sys,json
sys.path.insert(0,str(pathlib.Path.cwd()))
def main():
    import ou_irregular.runner
    import pytest
    os.environ.pop('MKL_NUM_THREADS',None);os.environ.pop('VECLIB_MAXIMUM_THREADS',None)
    from ou_irregular.experiments.provenance import environment
    out=pathlib.Path(os.environ['RUNNER_TEMP'])/'e4-checks'
    assert environment()==json.loads((out/'runtime.json').read_text())
    args=['-q','-ra','tests','--ignore=tests/week4/test_week3_isolation.py','--ignore=tests/week7/test_plots_integrity.py','--ignore=tests/week7/test_plots_reference.py','--ignore=tests/week7/test_selection.py','--junitxml='+str(out/'current.xml')]
    (out/'current-command.json').write_text(json.dumps(args)+'\n')
    return pytest.main(args)
if __name__=='__main__':raise SystemExit(main())
