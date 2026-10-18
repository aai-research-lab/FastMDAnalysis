"""Review observed preparation changes without inventing a correspondence."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest


DATA=Path(__file__).parent/"data"/"preparation_changes"


def _api():
    from fastmdxplora.gui import preparation_review
    return preparation_review


@pytest.fixture
def study(tmp_path):
    root=tmp_path/"study"
    shutil.copytree(DATA/"setup",root/"setup")
    return root


def test_the_producers_observations_and_notice_remain_distinct_from_choices(study):
    payload=_api().preparation_payload(study)
    assert payload["available"] and payload["status"]=="complete"
    assert "do not independently" in payload["notice"]
    assert [e["operation"] for e in payload["events"]]==[
        "pdbfixer_input","add_missing_hydrogens","prepared_system","recorded_setup_choices"]
    assert payload["events"][1]["details"]["ph"]==7.4
    assert payload["events"][3]["details"]["requested_settings"]["ph"]==7.4


def test_added_hydrogen_is_compared_by_author_identity_in_the_actual_snapshots(study):
    api=_api();payload=api.preparation_payload(study)
    change=api.preparation_change(study,"event-0002",payload["version"])
    assert change["available"] and change["compared"]
    assert len(change["before"]["atoms"])==4 and len(change["after"]["atoms"])==5
    assert change["added"]==[["A","1","","GLY","H"]]
    assert change["removed"]==[] and change["unchanged"]==4
    hydrogen=next(a for a in change["after"]["atoms"] if a["identity"][-1]=="H")
    assert hydrogen["xyz"]==pytest.approx([.7,2,3])


def test_atom_serials_and_line_order_are_not_used_as_correspondence():
    api=_api();pdb=(DATA/"setup"/"audit"/"source-0001.pdb").read_text()
    atoms,reason=api.pdb_atoms(pdb)
    reordered="\n".join(line[:6]+"99999"+line[11:] for line in reversed(pdb.splitlines()) if line.startswith("ATOM"))
    other,error=api.pdb_atoms(reordered)
    assert not reason and not error
    assert {tuple(a["identity"]) for a in atoms}=={tuple(a["identity"]) for a in other}


@pytest.mark.parametrize("kind",["alternate","duplicate","nonfinite","models"])
def test_ambiguous_or_invalid_atom_correspondence_is_withheld(kind):
    pdb=(DATA/"setup"/"audit"/"source-0001.pdb").read_text()
    atom=next(line for line in pdb.splitlines() if line.startswith("ATOM"))
    if kind=="alternate":pdb=pdb.replace(atom,atom[:16]+"A"+atom[17:])
    if kind=="duplicate":pdb+="\n"+atom
    if kind=="nonfinite":pdb=pdb.replace(atom,atom[:30]+"     nan"+atom[38:])
    if kind=="models":pdb="MODEL        1\n"+pdb+"MODEL        2\n"+pdb
    atoms,reason=_api().pdb_atoms(pdb)
    assert not atoms and reason


def test_a_changed_audit_cannot_reuse_a_previous_review_version(study):
    api=_api();payload=api.preparation_payload(study)
    target=study/"setup"/"preparation_audit.json";target.write_text(target.read_text()+"\n")
    changed=api.preparation_change(study,"event-0002",payload["version"])
    assert not changed["available"] and "changed" in changed["reason"]


def test_changed_snapshot_coordinates_cannot_be_presented_as_the_saved_change(study):
    api=_api();payload=api.preparation_payload(study)
    target=study/"setup"/"audit"/"source-0002.pdb";target.write_text(target.read_text()+"\n")
    change=api.preparation_change(study,"event-0002",payload["version"])
    assert not change["compared"] and change["after"] is None
    assert "changed" in change["after_reason"] and "added" not in change


def test_snapshot_checksums_keep_windows_line_endings(study):
    api=_api();target=study/"setup"/"audit"/"source-0002.pdb"
    raw=target.read_bytes().replace(b"\n",b"\r\n");target.write_bytes(raw)
    audit=study/"setup"/"preparation_audit.json";value=json.loads(audit.read_text())
    value["sources"]["source-0002"]["sha256"]=hashlib.sha256(raw).hexdigest();audit.write_text(json.dumps(value))
    change=api.preparation_change(study,"event-0002",api.preparation_payload(study)["version"])
    assert change["compared"] and len(change["added"])==1


def test_an_audit_cannot_name_an_arbitrary_file_as_its_snapshot(study):
    api=_api();audit=study/"setup"/"preparation_audit.json";value=json.loads(audit.read_text())
    value["sources"]["source-0002"]["snapshot"]="../../private.txt";audit.write_text(json.dumps(value))
    change=api.preparation_change(study,"event-0002",api.preparation_payload(study)["version"])
    assert not change["compared"] and change["after"] is None
    assert "not retained" in change["after_reason"]


def test_a_link_outside_the_study_is_not_followed(study,tmp_path):
    api=_api();target=study/"setup"/"audit"/"source-0002.pdb";target.unlink()
    private=tmp_path/"private.pdb";private.write_text("private marker")
    try:target.symlink_to(private)
    except OSError:pytest.skip("symlinks unavailable")
    change=api.preparation_change(study,"event-0002",api.preparation_payload(study)["version"])
    assert not change["available"] and "private marker" not in str(change)
