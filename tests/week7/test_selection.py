"""Independent recovery checks for outcome-informed selection and input seals.

Constructed archives below are disposable metadata fixtures, not numerical
experiments. The one published-input test reads existing evidence only.
"""
from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path

import pandas as pd
import pytest

from ou_irregular.experiments.audit import artifact_fingerprint, semantic_metadata_digest
from ou_irregular.experiments.design import Cell
from ou_irregular.experiments.provenance import ROOT, canonical, capture, file_digest
from ou_irregular.week7 import selection


def row(cid='a'*64, rep=0, method='exact', *, point=True, interval=True,
        white=None, score=None, error='', point_status='valid', interval_status='valid'):
    return {'cell_id':cid, 'rep':rep, 'estimator':method, 'point_valid':point,
            'interval_valid':interval, 'pipeline_error':error,
            'point_status':point_status, 'interval_status':interval_status,
            'hessian_whitened_difference':white, 'standardized_score':score}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(canonical(value)+'\n')


def bytes_map(root):
    return {str(p.relative_to(root)):p.read_bytes() for p in root.rglob('*') if p.is_file()}


def sealed_fixture(path):
    """Minimal bound metadata fixture; deliberately contains no simulated data."""
    path.mkdir()
    design=selection.followup_design().validation(1)
    write_json(path/'resolved_design.json', design.to_dict())
    (path/'data.csv').write_text('fixture\n1\n')
    table={'data.csv':{'sha256':file_digest(path/'data.csv'),'rows':1}}
    write_json(path/'table_manifest.json',table)
    (path/'figures').mkdir()
    (path/'figures/panel.png').write_bytes(b'fixture image bytes, no scientific claim')
    figures={'design_sha256':design.digest,
             'input_table_sha256':{'data.csv':file_digest(path/'data.csv')},
             'figures':{'panel.png':file_digest(path/'figures/panel.png')}}
    write_json(path/'figures/figure_manifest.json',figures)
    metadata={'status':'completed','provenance':capture(design),
              'completed_cells':{},'pipeline_errors':0}
    audit={'status':'passed','full_completion_audit':True,
           'audited_run_metadata_sha256':semantic_metadata_digest(metadata),
           'audited_artifact_sha256':artifact_fingerprint(path)}
    write_json(path/'audit.json',audit)
    metadata['audit_sha256']=file_digest(path/'audit.json')
    review={'reviewer':'unit-test fixture, not a scientific review',
            'all_panels_legible':True,'reviewed_pngs':['panel.png'],
            'figure_manifest_sha256':file_digest(path/'figures/figure_manifest.json')}
    write_json(path/'visual_review.json',review)
    rebind_completion(path,metadata)
    return design, metadata, figures


def rebind_completion(path,metadata):
    metadata['visual_review_sha256']=file_digest(path/'visual_review.json')
    metadata['artifact_sha256']={str(p.relative_to(path)):file_digest(p)
        for p in sorted(path.rglob('*')) if p.is_file() and p.name!='.execution.lock'
        and p!=path/'run_metadata.json'}
    write_json(path/'run_metadata.json',metadata)


def test_all_unavailable_and_pipeline_sources_are_retained_without_duplicate_sources():
    records=[row(rep=3,point=False,interval=False,point_status='boundary',interval_status='point_unavailable'),
             row(rep=3,method='pfml',point=False,interval=False,point_status='boundary',interval_status='point_unavailable'),
             row(rep=4,point=True,interval=False,interval_status='unstable_hessian'),
             row(rep=5,method='euler',point=False,interval=False,error='source exception',
                 point_status='source_failure',interval_status='source_failure'),
             row(rep=6)]
    sources,failures=selection.select_records(records)
    assert [source['rep'] for source in sources]==[3,4,5]
    assert len(failures)==4
    assert len(sources[0]['reasons'])==2
    assert failures[2]['point_status']=='valid' and failures[2]['interval_status']=='unstable_hessian'
    assert all(reason.startswith('unavailable:') for source in sources for reason in source['reasons'])


def test_accepted_extrema_are_per_method_metric_and_ties_are_input_order_independent():
    a,b,c='a'*64,'b'*64,'c'*64
    records=[row(b,0,white=.5,score=.1),row(a,9,white=.5,score=.1),
             row(a,2,white=.5,score=.1),row(c,1,white=.2,score=.9),
             row(b,3,'pfml',white=.7,score=.8),row(c,4,'euler',white=.2,score=.3)]
    sources,failures=selection.select_records(records)
    assert failures==[]
    expected={(a,2):{'accepted_maximum:exact:hessian_whitened_difference'},
              (c,1):{'accepted_maximum:exact:standardized_score'},
              (b,3):{'accepted_maximum:pfml:hessian_whitened_difference','accepted_maximum:pfml:standardized_score'},
              (c,4):{'accepted_maximum:euler:hessian_whitened_difference','accepted_maximum:euler:standardized_score'}}
    assert {(s['cell_id'],s['rep']):set(s['reasons']) for s in sources}==expected
    assert selection.select_records(reversed(records))==(sources,failures)


def test_nonfinite_missing_and_interval_invalid_values_cannot_win_accepted_extrema():
    records=[row(rep=0,white=.1,score=.2),row(rep=1,white=float('inf'),score=float('nan')),
             row(rep=2,white=None,score=None),
             row(rep=3,interval=False,white=1e6,score=1e6,interval_status='unstable_hessian')]
    sources,failures=selection.select_records(records)
    assert [s['rep'] for s in sources]==[0,3]
    assert len(failures)==1
    assert all(reason.startswith('unavailable:') for reason in sources[1]['reasons'])


def test_duplicate_fit_keys_are_rejected_before_selection_can_double_count():
    repeated=row(rep=7)
    with pytest.raises(ValueError,match='duplicate'):
        selection.select_records([repeated,dict(repeated)])


def test_dated_exploratory_designs_preserve_counts_truth_and_physical_floor_ids():
    follow=selection.followup_design();floor=selection.floor_design()
    assert (follow.name,follow.scope,follow.seed,follow.replications)==('week7_e3_followup','exploratory',20261007,10000)
    assert [(c.n,c.q,c.cv) for c in follow.cells]==[(250,2.,0.),(250,2.,.5)]
    assert (floor.name,floor.scope,floor.seed,floor.replications)==('week7_floor_sensitivity','exploratory',20261009,2000)
    assert {(c.n,c.q,c.cv,c.floor_fraction) for c in floor.cells}=={
        (1000,.5,cv,epsilon) for cv in (0.,1.5,2.,3.) for epsilon in (1e-8,1e-6,1e-4)}
    assert len({c.identity for c in floor.cells})==12
    assert len({replace(c,floor_fraction=1e-6).identity for c in floor.cells})==4
    assert all((c.theta,c.mu,c.sigma,c.floor_cv_threshold)==(1.,0.,.5,1.5) for c in (*follow.cells,*floor.cells))
    assert len(follow.cells)*follow.replications==20000
    assert len(floor.cells)*floor.replications==24000


@pytest.mark.parametrize('factory',[selection.followup_design,selection.floor_design])
def test_validation_namespace_changes_only_scope_name_seed_and_replications(factory):
    full=factory();development=full.validation(3)
    assert development.cells==full.cells
    assert development.scope=='validation' and development.name==full.name+'_validation'
    assert development.seed==full.seed+1_000_000 and development.replications==3
    assert development.digest!=full.digest


def test_valid_synthetic_seal_is_read_only(tmp_path):
    path=tmp_path/'fixture';design,metadata,figures=sealed_fixture(path)
    before=bytes_map(path)
    selected_design,identity=selection.sealed_input(path)
    assert selected_design==design
    assert identity['audit_sha256']==metadata['audit_sha256']
    assert identity['metadata_sha256']==file_digest(path/'run_metadata.json')
    selection.verify_completion_seal(path,metadata,figures)
    assert bytes_map(path)==before


@pytest.mark.parametrize('state',['running_incomplete','raw_complete_pending_audit','audited_pending_visual_review'])
def test_incomplete_input_states_are_rejected_without_writes(tmp_path,state):
    path=tmp_path/'fixture';_,metadata,_=sealed_fixture(path)
    metadata['status']=state;write_json(path/'run_metadata.json',metadata)
    before=bytes_map(path)
    with pytest.raises(ValueError,match='completed, fully audited'):
        selection.sealed_input(path)
    assert bytes_map(path)==before


@pytest.mark.parametrize('field,value',[('status','review_required'),('full_completion_audit',False)])
def test_unpassed_or_partial_audit_cannot_supply_selection(tmp_path,field,value):
    path=tmp_path/'fixture';sealed_fixture(path)
    audit=json.loads((path/'audit.json').read_text());audit[field]=value;write_json(path/'audit.json',audit)
    with pytest.raises(ValueError,match='completed, fully audited'):
        selection.sealed_input(path)


@pytest.mark.parametrize('damage',['scientific_metadata','raw_bytes','audit_bytes','provenance'])
def test_scientific_and_audit_binding_changes_are_detected(tmp_path,damage):
    path=tmp_path/'fixture';_,metadata,_=sealed_fixture(path)
    if damage=='scientific_metadata':
        metadata['pipeline_errors']=1;write_json(path/'run_metadata.json',metadata)
    elif damage=='raw_bytes': (path/'data.csv').write_text('changed\n2\n')
    elif damage=='audit_bytes': (path/'audit.json').write_text((path/'audit.json').read_text()+' ')
    else:
        metadata['provenance']['design_sha256']='0'*64;write_json(path/'run_metadata.json',metadata)
    before=bytes_map(path)
    with pytest.raises(ValueError):selection.sealed_input(path)
    assert bytes_map(path)==before


@pytest.mark.parametrize('damage',['extra_file','removed_file','changed_figure'])
def test_completion_seal_requires_the_entire_archive_map(tmp_path,damage):
    path=tmp_path/'fixture';_,metadata,figures=sealed_fixture(path)
    if damage=='extra_file': (path/'unexpected.txt').write_text('added')
    elif damage=='removed_file': (path/'data.csv').unlink()
    else: (path/'figures/panel.png').write_bytes(b'changed')
    with pytest.raises(ValueError,match='artifact map changed'):
        selection.verify_completion_seal(path,metadata,figures)


@pytest.mark.parametrize('field,value',[('all_panels_legible',False),('reviewed_pngs',[]),
    ('figure_manifest_sha256','0'*64),('reviewer','   ')])
def test_visual_review_content_is_checked_even_with_rebound_file_hashes(tmp_path,field,value):
    path=tmp_path/'fixture';_,metadata,figures=sealed_fixture(path)
    review=json.loads((path/'visual_review.json').read_text());review[field]=value
    write_json(path/'visual_review.json',review);rebind_completion(path,metadata)
    with pytest.raises(ValueError,match='visual review content'):
        selection.verify_completion_seal(path,metadata,figures)


def test_build_selection_rejects_nonregistered_headline_design_before_reading_rows(monkeypatch,tmp_path):
    monkeypatch.setattr(selection,'sealed_input',lambda p:(selection.followup_design(),{}))
    monkeypatch.setattr(selection,'saved_rows',lambda p:pytest.fail('read rows before registered-design check'))
    with pytest.raises(ValueError,match='exact registered E1/E2'):
        selection.build_selection(tmp_path/'e1',tmp_path/'e2')


def test_build_selection_preserves_experiment_labels_and_fixed_failure_cells(monkeypatch,tmp_path):
    designs={'e1':selection.load_e1(ROOT/'configs/week5_e1.yaml'),
             'e2':selection.load_e2(ROOT/'configs/week6_e2.yaml')}
    common=Cell(n=1000,q=.5,cv=.5)
    failed=[Cell(n=250,q=2,cv=cv) for cv in (0,.5)]
    records={'e1':[row(common.identity,4,white=.1,score=.2)],
             'e2':[row(common.identity,4,white=.3,score=.4),
                   *[row(c.identity,10,point=False,interval=False,point_status='boundary',
                         interval_status='point_unavailable') for c in failed]]}
    monkeypatch.setattr(selection,'sealed_input',lambda p:(designs[Path(p).name],{'fixture':Path(p).name}))
    monkeypatch.setattr(selection,'saved_rows',lambda p:iter(records[Path(p).name]))
    unresolved={'cell_id':common.identity,'n':1000,'theta_mean_gap':.5,'cv':.5,
                'estimator':'pfml','parameter':'sigma','practical_classification':'unresolved'}
    monkeypatch.setattr(selection.pd,'read_csv',lambda p:pd.DataFrame([unresolved]))
    output=selection.build_selection(tmp_path/'e1',tmp_path/'e2')
    assert output['scope']=='exploratory_selection_after_complete_headlines'
    assert {(s['experiment'],s['cell_id'],s['rep']) for s in output['sources']}=={
        ('E1',common.identity,4),('E2',common.identity,4),
        *[('E2',c.identity,10) for c in failed]}
    assert output['followup_cells']==[c.to_dict() for c in failed]
    assert len(output['failed_fit_records'])==2
    assert [r['experiment'] for r in output['unresolved_parameter_rows']]==['E1','E2']
    records['e2'].append(row(Cell(n=1000,q=1,cv=1).identity,0,interval=False,interval_status='unexpected'))
    with pytest.raises(ValueError,match='dated two-cell'):
        selection.build_selection(tmp_path/'e1',tmp_path/'e2')


@pytest.mark.parametrize('name',['week5_e1','week6_e2'])
def test_published_completed_inputs_verify_without_regeneration_or_mutation(name,monkeypatch):
    path=ROOT/'artifacts'/name
    from ou_irregular.week4.simulation import OuSampleSource
    monkeypatch.setattr(OuSampleSource,'sample',lambda *a:pytest.fail('input verification regenerated a source'))
    before=file_digest(path/'run_metadata.json')
    design,identity=selection.sealed_input(path)
    assert design.scope=='headline'
    assert identity['metadata_sha256']==before==file_digest(path/'run_metadata.json')
