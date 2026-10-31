"""E4 figure contracts from archived R=4 evidence; sampling/fitting forbidden."""
from copy import deepcopy
import ast
import inspect
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt

from ou_irregular.experiments.provenance import ROOT, CORE_FILES, canonical
from ou_irregular.week8 import audit, execution, plots, reporting, source
from ou_irregular.week8.design import validation_design
from .test_audit import mutate_table, mutate_partition


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path, value): Path(path).write_text(canonical(value)+'\n')
def snapshot(path):
    return {str(p.relative_to(path)): sha(p) for p in path.rglob('*') if p.is_file() and p.name != '.execution.lock'}
def forbidden(*a, **k): raise AssertionError('figure stage must not generate sources, refit, or replace the audit')


@pytest.fixture(autouse=True)
def no_science(monkeypatch):
    for name in ('sample','gaps'): monkeypatch.setattr(source.E4Source, name, forbidden)
    for name in ('generate_gaps','simulate_ou'): monkeypatch.setattr(source, name, forbidden)
    for name in ('audit_run','independent_source','reconstruct_schedule','source_records'): monkeypatch.setattr(audit, name, forbidden)
    monkeypatch.setattr(np.random, 'default_rng', forbidden)


@pytest.fixture(scope='module')
def saved(tmp_path_factory):
    # This is a copy of retained evidence, not a new invocation of a scientific run.
    parent = Path(os.environ['OU_E4_PLOTS_FIXTURE_ROOT']) if 'OU_E4_PLOTS_FIXTURE_ROOT' in os.environ else tmp_path_factory.mktemp('e4_plot')
    parent.mkdir(parents=True, exist_ok=True)
    out = parent/'saved_archive'
    assert not out.exists()
    original = ROOT/'artifacts/week8_development/audit/ci/native_fixture'
    shutil.copytree(original/'native_archive', out)
    shutil.copyfile(original/'native_audit.json', out/'audit.json')
    return out


@pytest.fixture
def archive(saved, tmp_path):
    out = tmp_path/'e4'
    shutil.copytree(saved, out)
    return out


@pytest.fixture(scope='module')
def tables(saved):
    return {name: pd.read_csv(saved/name, float_precision='round_trip') for name in reporting.TABLE_KEYS}


@pytest.fixture(scope='module')
def built(tables):
    figures = plots.build_figures(tables, validation_design())
    for fig in figures.values(): plt.close(fig)
    yield figures


@pytest.mark.parametrize('stem', plots.STEMS)
def test_each_planned_figure_is_single_axes_with_validation_label(built, stem):
    fig = built[stem]
    assert len(fig.axes) == 1
    assert any('VALIDATION ONLY' in t.get_text() for t in fig.texts)
    assert fig.axes[0].get_ylabel()
    assert len(plots.STEMS) == 18 and len(set(plots.STEMS)) == 18



@pytest.mark.parametrize('parameter', ['theta','mu','sigma'])
def test_coverage_contains_exact_saved_values_for_all_seven_cells(tables, built, parameter):
    fig = built['W8_coverage_'+parameter]
    for method, _, _ in plots.METHODS:
        data = tables['parameter_summary.csv'].query('estimator == @method and parameter == @parameter').set_index('cell_id')
        line = next(l for l in fig.axes[0].lines if l.get_gid() == method+'/operational_coverage')
        expected = [100*data.loc[c.identity,'coverage_operational'] for c in validation_design().cells]
        np.testing.assert_array_equal(line.get_ydata(), expected)
        assert len(line.get_xdata()) == 7


@pytest.mark.parametrize('metric,field,factor', [('bias','excess_bias',1),('coverage','operational_coverage_difference',100),('rmse_ratio','rmse_ratio',1)])
def test_matched_mu_is_absolute_not_relative_and_uses_saved_independent_values(tables, built, metric, field, factor):
    ax = built['W8_matched_'+metric+'_mu'].axes[0]
    cells = [c for c in validation_design().cells if c.mechanism == 'mmpp']
    for method, _, _ in plots.METHODS:
        data = tables['matched_iid_comparisons.csv'].query('estimator == @method and parameter == "mu"').set_index('cell_id')
        line = next(l for l in ax.lines if l.get_gid() == method+'/'+field)
        np.testing.assert_array_equal(line.get_ydata(), [factor*data.loc[c.identity,field] for c in cells])
    assert (ax.get_yscale() == 'log') == (metric == 'rmse_ratio')


def test_source_span_sd_uses_saved_across_path_mcse_and_not_gap_count(tables, built):
    data = tables['source_summary.csv'].set_index('cell_id')
    line = next(l for l in built['W8_source_span_sd'].axes[0].lines if l.get_gid() == 'realized/span_sd')
    np.testing.assert_array_equal(line.get_ydata(), [2*data.loc[c.identity,'realized_span_mcse'] for c in validation_design().cells])
    assert 'not the CV' in (ROOT/'docs/week8_plots.md').read_text() or 'not equated' in (ROOT/'docs/week8_plots.md').read_text()


def test_all_point_and_interval_availability_series_are_separate(built):
    gids = {l.get_gid() for l in built['W8_availability'].axes[0].lines if l.get_gid()}
    assert gids == {method+'/'+kind+'_availability' for method,_,_ in plots.METHODS for kind in ('point','interval')}


def test_missing_and_unavailable_uncertainty_are_not_deleted_or_filled():
    fig, ax = plt.subplots()
    try:
        plots._series(ax, [0,1,2], [2,np.nan,3], [1,np.nan,np.nan], [3,np.nan,np.nan], label='literal', marker='o', gid='test')
        assert list(next(l for l in ax.lines if l.get_gid() == 'test').get_ydata()) == [2,3]
        assert [t.get_text() for t in ax.texts] == ['NA','†']
    finally: plt.close(fig)


def test_zero_ratio_has_explicit_off_log_scale_label_not_false_missing():
    fig, ax = plt.subplots()
    try:
        plots._series(ax, [0,1], [0,1], [np.nan,.5], [np.nan,2], label='literal', marker='o', gid='test', ratio=True)
        assert any(t.get_text() == '0' for t in ax.texts)
        assert list(next(l for l in ax.lines if l.get_gid() == 'test').get_ydata()) == [1]
    finally: plt.close(fig)


@pytest.mark.parametrize('values,lower,upper', [([1],[2],[3]),([1],[0],[np.inf]),([-1],[-2],[0])])
def test_invalid_errorbars_or_negative_ratios_are_rejected(values, lower, upper):
    fig, ax = plt.subplots()
    try:
        with pytest.raises(ValueError): plots._series(ax,[0],values,lower,upper,label='bad',marker='o',gid='bad',ratio=True)
    finally: plt.close(fig)


def test_saved_readback_is_pure_and_binds_all_expected_inputs(archive):
    before = snapshot(archive)
    result = plots.load_inputs(archive, validation_design())
    assert snapshot(archive) == before
    assert result.binding['audit_sha256'] == sha(archive/'audit.json')
    assert set(result.binding['input_table_sha256']) == set(execution.TABLE_FILES)
    assert result.binding['raw_input_sha256'] == json.loads((archive/'table_binding.json').read_text())['raw_input_sha256']


@pytest.mark.parametrize('field,value', [('full_completion_audit',False),('replayed_fits',41),('source_paths_verified',27),('raw_input_sha256','0'*64),('schedule_trace_hashes_verified',0)])
def test_invalid_audit_declarations_block_before_any_output(archive, field, value):
    path=archive/'audit.json'; report=json.loads(path.read_bytes()); report[field]=value;write_json(path,report)
    before=snapshot(archive)
    with pytest.raises(ValueError): plots.write_figures(archive,validation_design())
    assert snapshot(archive)==before and not plots.renderer_state_path(archive).exists()


def test_rehashed_endpoint_descriptor_corruption_is_rejected(archive):
    path=archive/'audit.json'; report=json.loads(path.read_bytes());report['replay_records'][0]['times_sha256']='0'*64;write_json(path,report)
    with pytest.raises(ValueError): plots.load_inputs(archive,validation_design())


def test_rehashed_table_mutation_does_not_evade_saved_audit(archive):
    mutate_table(archive,'parameter_summary.csv',lambda cols,rows: rows[0].update(bias='123'))
    with pytest.raises(ValueError): plots.load_inputs(archive,validation_design())


def test_rehashed_raw_mutation_does_not_evade_saved_audit(archive):
    cell=validation_design().cells[0]
    def change(rows):
        for row in rows:
            row['raw_mean_gap'] += .01
    mutate_partition(archive,cell,change)
    with pytest.raises(ValueError): plots.load_inputs(archive,validation_design())


@pytest.mark.parametrize('kind',['missing_table','extra_file','symlink','unfrozen_full_alias'])
def test_invalid_input_layout_or_design_is_read_only(archive,kind):
    if kind=='missing_table': (archive/'source_summary.csv').unlink()
    elif kind=='extra_file': (archive/'unknown').write_text('x')
    elif kind=='symlink':
        p=archive/'source_summary.csv'; payload=p.read_bytes();p.unlink();q=archive.parent/'target';q.write_bytes(payload);p.symlink_to(q)
    else:
        p=archive/'run_metadata.json';m=json.loads(p.read_bytes());m['require_frozen']=True;write_json(p,m)
    before=snapshot(archive)
    with pytest.raises((ValueError,FileNotFoundError)): plots.write_figures(archive,validation_design())
    assert snapshot(archive)==before


@pytest.fixture(scope='module')
def rendered(saved):
    # Retains actual PNG/PDF bytes for separate visual inspection after tests.
    out=saved.parent/'rendered_archive';shutil.copytree(saved,out)
    result=plots.write_figures(out,validation_design())
    return out,result


def test_native_render_has_all_pairs_and_repeat_is_pure(rendered):
    out,result=rendered; before=snapshot(out)
    assert set(result['figures'])==plots.FILENAMES and len(result['figures'])==36
    assert result['new_sources_or_fits'] is False and result['visual_review_performed'] is False
    for stem in plots.STEMS:
        assert (out/'figures'/(stem+'.png')).read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
        assert (out/'figures'/(stem+'.pdf')).read_bytes().startswith(b'%PDF')
    assert plots.write_figures(out,validation_design())==result
    assert plots.verify_figures(out,validation_design())==result
    assert snapshot(out)==before


@pytest.fixture
def rendered_copy(rendered,tmp_path):
    original,_=rendered;out=tmp_path/'rendered';shutil.copytree(original,out)
    src=plots.renderer_state_path(original);dst=plots.renderer_state_path(out);dst.parent.mkdir();shutil.copyfile(src,dst)
    return out


def test_identical_partial_adoption_recreates_only_missing_pair(rendered_copy):
    out=rendered_copy;before=snapshot(out)
    (out/'figures/figure_manifest.json').unlink()
    (out/'figures'/next(iter(sorted(plots.FILENAMES)))).unlink()
    plots.write_figures(out,validation_design())
    assert snapshot(out)==before


def test_conflicting_late_partial_blocks_every_missing_output(rendered_copy):
    out=rendered_copy;names=sorted(plots.FILENAMES)
    (out/'figures/figure_manifest.json').unlink();(out/'figures'/names[0]).unlink();(out/'figures'/names[-1]).write_bytes(b'changed')
    before=snapshot(out)
    with pytest.raises(ValueError,match='partial figure'): plots.write_figures(out,validation_design())
    assert snapshot(out)==before and not (out/'figures'/names[0]).exists()


def test_incomplete_pair_fails_before_state_or_figure_publication(archive,monkeypatch):
    monkeypatch.setattr(plots,'figure_bytes',lambda fig,stem:{stem+'.png':b'png'})
    before=snapshot(archive)
    with pytest.raises(ValueError,match='pair incomplete'): plots.write_figures(archive,validation_design())
    assert snapshot(archive)==before and not plots.renderer_state_path(archive).exists()


def test_missing_original_renderer_state_blocks_partial_adoption(rendered_copy):
    out=rendered_copy;(out/'figures/figure_manifest.json').unlink();plots.renderer_state_path(out).unlink();plots.renderer_state_path(out).parent.rmdir()
    before=snapshot(out)
    with pytest.raises(ValueError,match='original renderer'): plots.write_figures(out,validation_design())
    assert snapshot(out)==before


def test_complete_embedded_renderer_allows_readback_without_sibling(rendered_copy):
    out=rendered_copy;plots.renderer_state_path(out).unlink();plots.renderer_state_path(out).parent.rmdir()
    before=snapshot(out)
    plots.verify_figures(out,validation_design())
    assert snapshot(out)==before


@pytest.mark.parametrize('kind',['figure_bytes','manifest_binding','renderer_settings','renderer_dependency','missing_figure','unknown_figure'])
def test_complete_figure_tampering_is_not_repaired(rendered_copy,kind):
    out=rendered_copy;p=out/'figures/figure_manifest.json';m=json.loads(p.read_bytes());name=sorted(plots.FILENAMES)[0]
    if kind=='figure_bytes': (out/'figures'/name).write_bytes(b'bad')
    elif kind=='missing_figure': (out/'figures'/name).unlink()
    elif kind=='unknown_figure': (out/'figures/unknown.png').write_bytes(b'bad')
    elif kind=='manifest_binding': m['input_binding']['audit_sha256']='0'*64;write_json(p,m)
    elif kind=='renderer_settings': m['renderer']['descriptor']['renderer_environment']['png_dpi']=1;write_json(p,m)
    else: m['renderer']['provenance']['source_files']['ou_irregular/week8/plots.py']='0'*64;write_json(p,m)
    before=snapshot(out)
    with pytest.raises(ValueError): plots.write_figures(out,validation_design())
    assert snapshot(out)==before


def test_root_drift_during_rendering_preserves_outside_write(archive,monkeypatch):
    original=plots.figure_bytes;called=[]
    def change(fig,stem):
        if not called:
            p=archive/'run_metadata.json';m=json.loads(p.read_bytes());m['outside_change']='retained';write_json(p,m);called.append(1)
        return original(fig,stem)
    monkeypatch.setattr(plots,'figure_bytes',change)
    with pytest.raises(execution.ManifestDriftError): plots.write_figures(archive,validation_design())
    assert json.loads((archive/'run_metadata.json').read_bytes())['outside_change']=='retained'
    assert not (archive/'figures').exists()
    assert not plots.renderer_state_path(archive).exists()


def test_completed_seal_is_read_only_and_wrong_artifact_map_is_rejected(rendered_copy):
    out=rendered_copy
    review={'reviewer':'literal unit fixture; not actual visual inspection','reviewed_pngs':sorted(n for n in plots.FILENAMES if n.endswith('.png')),
            'all_panels_legible':True,'figure_manifest_sha256':sha(out/'figures/figure_manifest.json')}
    write_json(out/'visual_review.json',review)
    m=json.loads((out/'run_metadata.json').read_bytes());m.update(status='completed',audit_sha256=sha(out/'audit.json'),visual_review_sha256=sha(out/'visual_review.json'))
    m['artifact_sha256']={k:v for k,v in snapshot(out).items() if k!='run_metadata.json'};write_json(out/'run_metadata.json',m)
    before=snapshot(out);plots.verify_figures(out,validation_design());assert snapshot(out)==before
    m['artifact_sha256']['source_summary.csv']='0'*64;write_json(out/'run_metadata.json',m)
    with pytest.raises(ValueError): plots.verify_figures(out,validation_design())


def test_dependency_closure_includes_all_imported_local_modules():
    probe = """import json, sys
from pathlib import Path
from ou_irregular.experiments.provenance import ROOT, CORE_FILES
from ou_irregular.week8 import plots
modules = {str(Path(m.__file__).resolve().relative_to(ROOT))
           for n,m in list(sys.modules.items()) if n.startswith('ou_irregular.')
           and getattr(m,'__file__',None) and Path(m.__file__).resolve().is_relative_to(ROOT)}
print(json.dumps(sorted(modules - (set(CORE_FILES)|set(plots.FILES)))))
"""
    result=subprocess.run([sys.executable,'-c',probe],cwd=ROOT,text=True,capture_output=True,check=True)
    assert json.loads(result.stdout)==[]


def test_no_explicit_palette_or_multiaxes_subplot_in_builders():
    tree=ast.parse(inspect.getsource(plots))
    for node in ast.walk(tree):
        if not isinstance(node,ast.Call): continue
        assert not any(k.arg in {'color','c','cmap','facecolor','edgecolor','ecolor'} for k in node.keywords)
        if isinstance(node.func,ast.Attribute) and node.func.attr=='subplots':
            assert not node.args
            assert not any(k.arg in {'nrows','ncols'} for k in node.keywords)


@pytest.mark.parametrize('bad',['abc',b'abc'])
def test_dependency_strings_rejected_without_reading_or_creating_outputs(tmp_path,bad):
    with pytest.raises(ValueError): plots.write_figures(tmp_path/'missing',validation_design(),extra_files=bad)
    assert not list(tmp_path.iterdir())


def test_zero_clipping_display_keeps_a_nonnegative_axis(built):
    ax=built['W8_source_clipping'].axes[0]
    assert ax.get_ylim()[0]==0 and ax.get_ylim()[1]>=1
    line=next(l for l in ax.lines if l.get_gid()=='source/clipping')
    assert np.all(line.get_ydata()==0)  # retained fixture has zero clips, not imputed zeros


def test_external_replot_preserves_entire_source_and_reuses_saved_values(saved,tmp_path):
    before=snapshot(saved);out=tmp_path/'external'
    result=plots.write_figures(saved,validation_design(),destination=out)
    assert snapshot(saved)==before
    assert {p.name for p in out.iterdir()}=={'figures','.execution.lock'}
    assert set(result['figures'])==plots.FILENAMES
    dest_before=snapshot(out)
    assert plots.verify_figures(saved,validation_design(),destination=out)==result
    assert plots.write_figures(saved,validation_design(),destination=out)==result
    assert snapshot(out)==dest_before and snapshot(saved)==before
    assert result['visual_review_performed'] is False


@pytest.mark.parametrize('relation',['same','child','parent','unrelated_files'])
def test_external_replot_refuses_unsafe_destination_without_changes(archive,tmp_path,relation):
    target={'same':archive,'child':archive/'derived','parent':archive.parent,'unrelated_files':tmp_path/'occupied'}[relation]
    if relation=='unrelated_files': target.mkdir();(target/'notes.txt').write_text('retain')
    before=snapshot(tmp_path)
    with pytest.raises(ValueError): plots.write_figures(archive,validation_design(),destination=target)
    assert snapshot(tmp_path)==before


def test_external_replot_accepts_a_valid_sealed_source_without_mutation(rendered_copy,tmp_path):
    out=rendered_copy
    review={'reviewer':'unit-test assertion, not actual inspection','reviewed_pngs':sorted(n for n in plots.FILENAMES if n.endswith('.png')),
            'all_panels_legible':True,'figure_manifest_sha256':sha(out/'figures/figure_manifest.json')}
    write_json(out/'visual_review.json',review)
    m=json.loads((out/'run_metadata.json').read_bytes());m.update(status='completed',audit_sha256=sha(out/'audit.json'),visual_review_sha256=sha(out/'visual_review.json'))
    m['artifact_sha256']={k:v for k,v in snapshot(out).items() if k!='run_metadata.json'};write_json(out/'run_metadata.json',m)
    before=snapshot(out);target=tmp_path/'separate'
    result=plots.write_figures(out,validation_design(),destination=target)
    assert snapshot(out)==before
    assert result['figures']==json.loads((out/'figures/figure_manifest.json').read_bytes())['figures']
