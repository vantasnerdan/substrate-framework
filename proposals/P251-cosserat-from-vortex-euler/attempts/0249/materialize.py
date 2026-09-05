from pathlib import Path
from datetime import datetime, timezone
import hashlib
import copy
import yaml

root=Path(__file__).resolve().parents[4]
attempt=root/'proposals/P251-cosserat-from-vortex-euler/attempts/0249'
registry_path=root/'governance/claims.yaml'
payload_path=attempt/'claim-payload.yaml'
payload=yaml.safe_load(payload_path.read_text())
old_registry=yaml.safe_load(registry_path.read_text())
assert not any(c['id']==payload['id'] for c in old_registry['claims'])
assert len(old_registry['claims'])==269
for p in payload['evidence']:
    assert (root/p).is_file(), p
release_path=root/'governance/releases/current.yaml'
release=yaml.safe_load(release_path.read_text())
assert release['release']=='v0.181.0'
stamp=datetime.now(timezone.utc).isoformat(timespec='seconds')
base='1edefee'
import subprocess
base=subprocess.check_output(['git','rev-parse',base],cwd=root,text=True).strip()
review=attempt/'materialization-review.md'
assert review.is_file()
sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
campaign=root/'campaigns/P251-periodic-joint-continuum'
campaign.mkdir()
adjudication={
 'schema_version':1,'campaign':'P251','canonical_issue':'https://github.com/vantasnerdan/substrate-framework/issues/200',
 'base_release':'v0.181.0','source_baseline':base,'adjudicated_at':stamp,
 'adjudication_scope':'One reviewed prepared periodic actual Euler history/action/current continuum and its same-cell finite invariant tube; compact Euclidean parent remains active.',
 'claims_promoted':['C-CST-017'],'release':'v0.182.0',
 'frozen_statement':{'artifact':str(payload_path.relative_to(root)),'sha256':sha(payload_path)},
 'review':{'reviewer':'herdr optical-review pane w3:p3','non_author_of_claim_delta':True,
   'transaction':'proposals/P251-cosserat-from-vortex-euler/attempts/0247/review.md',
   'sha256':sha(root/'proposals/P251-cosserat-from-vortex-euler/attempts/0247/review.md'),
   'materialization':str(review.relative_to(root)), 'materialization_sha256':sha(review),
   'geometry_reviewer':'herdr geometry-review pane w3:p4',
   'geometry':'proposals/P251-cosserat-from-vortex-euler/attempts/0245/README.md',
   'geometry_sha256':sha(root/'proposals/P251-cosserat-from-vortex-euler/attempts/0245/README.md'),
   'independence_scope':'Separate non-author processes reviewed the new suppliers and joint implication; root authored the source delta and additive extraction.'},
 'positive_result':'Actual common prepared periodic fields satisfy the coupled physical equations, inherited action and complete bulk constitutive virtual work through second spatial order, with exact same-cell geometry and density.',
 'evidence_scope':'Analytic source construction and reviewed error ordering; exact residual and material-current identities; eight direct importable API tests. No all-time or compact Euclidean completion is inferred.',
 'archive_policy':'All269 prior claim objects and prior immutable campaigns retain their statements; attempts and review hashes preserve the source history.',
 'validation':{'implementation':'proposals/P251-cosserat-from-vortex-euler/attempts/0249/receipt.md','promotion':'proposals/P251-cosserat-from-vortex-euler/attempts/0249/receipt.md'},
 'in_boundary_debt':[],'full_coupled_continuum':'active','goal_verdict':'active','terminal_pr_eligible':False}
(campaign/'adjudication.yaml').write_text(yaml.safe_dump(adjudication,sort_keys=False,width=100))
lines=payload_path.read_text().splitlines()
addition='  - '+lines[0]+'\n'+''.join('    '+line+'\n' if line else '\n' for line in lines[1:])
newtext=registry_path.read_text().rstrip()+'\n'+addition
new_registry=yaml.safe_load(newtext)
assert new_registry['claims'][:-1]==old_registry['claims']
assert new_registry['claims'][-1]==payload
registry_path.write_text(newtext)
release=copy.deepcopy(release)
release.update(release='v0.182.0',source_baseline='substrate-framework@'+base,released_at=stamp)
release['accepted_claims'].append(payload['id'])
out=yaml.safe_dump(release,sort_keys=False)
(root/'governance/releases/v0.182.0.yaml').write_text(out)
release_path.write_text(out)
print('Materialized C-CST-017 / v0.182.0; all269 prior claim objects identical; parent active.')
