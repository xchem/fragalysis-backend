"""Delete the latest ExperimentUpload of a target, restoring the previous state.

This is the reusable core of the upload-deletion feature, modelled on
:mod:`viewer.target_delete`: it is deliberately free of any HTTP or management-command
concerns, so the *caller* owns authorization and the deployment-mode guard.

Only the *latest* upload can be deleted, but the operation is repeatable: removing upload
N makes N-1 the latest, which this can then remove in turn, down to the target's first
upload - for which :func:`viewer.target_delete.delete_target` is the right tool.

What the upload added is read from its own ``meta_aligner.yaml`` rather than inferred from
the database. The bundle is the input that created the rows, so it cannot drift from them:
the ``aligned_files`` nesting ``crystal -> chain -> ligand -> altloc -> version`` composes
to exactly the ``longcode`` the loader stores (``target_loader.process_site_observation``),
which is a 1:1 handle on the row. ``SiteObservation.version`` records the same thing and is
used as a cross-check - the two agreed on every target tested, so a disagreement means one
of the assumptions is wrong on that data and the run refuses rather than guesses.

Three properties make repeated deletion safe, and this module is organised around them:

* ``Versionable.version`` is immutable, so row identification does not drift as uploads go;
* experiments are re-pointed to ``min(surviving observation version)``, a fixed point -
  "the previous upload" is correct once but collapses under iteration;
* ``superseded`` is recomputed from the surviving rows only, so it converges rather than
  accumulating.

Not everything a load did is reversible, and the plan reports rather than papers over it:
file paths on carried-over experiments that this upload overwrote have no recorded previous
value, and compound reconciliation merges are not undone (the merged rows are duplicates,
so the surviving compound carries the structure - see ``_apply_compound_supersessions``).
Compounds the upload created are removed only once nothing anywhere points at them.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import yaml
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.db import transaction
from django.db.models import Count, Exists, Min, Model, OuterRef, Q

from viewer.cache import clear_view_cache
from viewer.compound_dedup import _through_compound_fk_attname
from viewer.media_watcher import deletions_expected
from viewer.models import (
    CanonSite,
    CanonSiteConf,
    Compound,
    Experiment,
    ExperimentUpload,
    Pose,
    QuatAssembly,
    SiteObservation,
    Xtalform,
    XtalformQuatAssembly,
    XtalformSite,
)
from viewer.utils import longcode_from_tag, strip_version

logger = logging.getLogger(__name__)

#: The metadata file that records what a bundle contained.
METADATA_FILE = "meta_aligner.yaml"

#: ``superseded`` is recomputed for these. CanonSite is excluded on purpose: its supersede
#: identity includes ``version`` (see :attr:`viewer.models.CanonSite.SUPERSEDE_FIELDS`), so
#: the invariant does not describe the rows the loader actually flags. Repairing those is a
#: separate job; a deletion must not quietly change rows it did not touch.
RESTORED_MODELS = (SiteObservation, CanonSiteConf, XtalformSite)

#: Models the metadata declares, by the name used in the metadata-derived sets.
DECLARED_MODELS = (SiteObservation, CanonSite, CanonSiteConf, XtalformSite)

#: Experiment file fields, and the ``crystallographic_files`` key that supplies each.
#: Mirrors ``target_loader.process_experiment``, which stores
#: ``<target dir>/<the yaml's 'file' value>`` in the field.
EXPERIMENT_FILE_KEYS = {
    "pdb_info": "xtal_pdb",
    "mtz_info": "xtal_mtz",
    "cif_info": "ligand_cif",
}

#: The list-valued one: ``map_info`` is an ArrayField of the event-map files.
EXPERIMENT_MAP_FIELD = "map_info"
EXPERIMENT_MAP_KEY = "ligand_binding_events"


class NotDeletable(Exception):
    """The upload cannot be deleted: not the latest, the only one, or unreadable.

    ``whole_target_instead`` marks the one refusal that deleting the entire target would
    actually resolve, so a caller can offer that without having to match on the message.
    """

    def __init__(self, *args, whole_target_instead: bool = False):
        super().__init__(*args)
        self.whole_target_instead = whole_target_instead


@dataclass
class UploadDeletionPlan:
    """What a deletion did, or (from :func:`plan_upload_deletion`) would do."""

    upload_id: int
    upload_version: int
    target_title: str
    doomed: dict[str, list[int]] = field(default_factory=dict)
    orphaned: dict[str, list[int]] = field(default_factory=dict)
    poses_repointed: list[tuple[int, int, int]] = field(default_factory=list)
    poses_removed: list[int] = field(default_factory=list)
    experiments_repointed: list[tuple[int, int]] = field(default_factory=list)
    experiment_files_restored: list[tuple[int, list]] = field(default_factory=list)
    experiment_files_cleared: list[tuple[int, list]] = field(default_factory=list)
    children_repointed: list[tuple[str, int, int, int]] = field(default_factory=list)
    refs_cleared: list[tuple[str, int]] = field(default_factory=list)
    restored: dict[str, int] = field(default_factory=dict)
    files_removed: list[Path] = field(default_factory=list)
    files_kept: list[Path] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """A human-readable digest, used by the management command."""
        lines = [
            f"Target {self.target_title}: upload {self.upload_version} "
            f"(ExperimentUpload pk={self.upload_id})",
        ]
        for name, pks in sorted(self.doomed.items()):
            lines.append(f"  delete {len(pks):>6} {name}")
        for name, pks in sorted(self.orphaned.items()):
            if pks:
                lines.append(f"  orphan {len(pks):>6} {name}")
        lines.extend(
            [
                f"  poses  {len(self.poses_repointed):>6} re-pointed to an earlier main",
                f"  poses  {len(self.poses_removed):>6} removed (no members left)",
                f"  exps   {len(self.experiments_repointed):>6} re-pointed to an earlier upload",
                f"  exps   {len(self.experiment_files_restored):>6} had file paths restored",
                f"  exps   {len(self.experiment_files_cleared):>6} had file paths cleared (this upload first supplied them)",
                f"  refs   {len(self.refs_cleared):>6} reference FKs cleared",
                f"  files  {len(self.files_removed):>6} removed, "
                f"{len(self.files_kept)} kept (still referenced)",
            ]
        )
        if self.restored:
            flags = ", ".join(f"{k}={v}" for k, v in sorted(self.restored.items()))
            lines.append(f"  superseded flags changed: {flags}")
        for warning in self.warnings:
            lines.append(f"  WARNING: {warning}")
        return "\n".join(lines)


# --------------------------------------------------------------------------------------
# eligibility
# --------------------------------------------------------------------------------------


def _check_deletable(upload: ExperimentUpload) -> None:
    """Refuse anything but the latest upload of a multi-upload target.

    The first upload is not deletable because there is no previous state to restore to -
    ``viewer.target_delete.delete_target`` removes the whole graph and its media instead.
    """
    siblings = ExperimentUpload.objects.filter(target=upload.target)
    if siblings.count() < 2:
        raise NotDeletable(
            f"{upload.target.title} has only this upload, so deleting it would mean "
            "deleting the target itself.",
            whole_target_instead=True,
        )
    latest = siblings.order_by("upload_version").last()
    if latest is None or latest.pk != upload.pk:
        raise NotDeletable(
            f"Only the latest upload of a target can be deleted; "
            f"{upload.target.title} is at upload "
            f"{latest.upload_version if latest else '?'}."
        )


# --------------------------------------------------------------------------------------
# what the bundle says it added
# --------------------------------------------------------------------------------------


def read_upload_metadata(upload: ExperimentUpload) -> dict[str, set[str]]:
    """The natural keys of everything ``meta_aligner.yaml`` says this upload created.

    The file is cumulative: it lists what the *target* contains after the load, each entry
    carrying the version it was created at. The entries whose version equals this upload's
    own version number are the ones it added; the rest were re-found, not created.

    Site observations come back as longcodes, composed from the ``aligned_files`` path in
    exactly the way :mod:`viewer.target_loader` composes them. The other three sections key
    on ``<name><separator><version>``, split the same way the loader splits them.
    """
    path = Path(upload.get_upload_path()) / METADATA_FILE
    if not path.is_file():
        raise NotDeletable(
            f"{path} is missing; cannot establish what this upload added. "
            "The bundle's own metadata is the authority for the deletion."
        )

    meta = yaml.safe_load(path.read_text(encoding="utf-8")) or {}

    declared_format = str(meta.get("data_format_version", ""))
    supported_major = settings.XCA_DATA_FORMAT_VERSION.split(".")[0]
    if declared_format.split(".")[0] != supported_major:
        raise NotDeletable(
            f"{path} declares data format version '{declared_format}'; only major "
            f"version {supported_major} is supported. Older bundles nest aligned_files "
            "without an altloc level, so their entries cannot be mapped onto rows."
        )

    try:
        version_number = int(meta["version_number"])
    except (KeyError, TypeError, ValueError) as exc:
        raise NotDeletable(f"{path} has no usable 'version_number'.") from exc

    if version_number != upload.upload_version:
        raise NotDeletable(
            f"{path} declares version_number {version_number} but the upload is recorded "
            f"as version {upload.upload_version}. Refusing to guess which is right."
        )

    longcodes = set()
    for crystal, entry in (meta.get("crystals") or {}).items():
        for chain, by_ligand in ((entry or {}).get("aligned_files") or {}).items():
            for ligand, by_altloc in (by_ligand or {}).items():
                for altloc, by_version in (by_altloc or {}).items():
                    for version in by_version or {}:
                        if int(version) != version_number:
                            continue
                        longcodes.add(
                            longcode_from_tag(
                                f"{crystal}/{chain}/{ligand}/{altloc}/{version}"
                            )
                        )

    def created(section: str, separator: str) -> set[str]:
        names = set()
        for key in meta.get(section) or {}:
            try:
                name, version = strip_version(str(key), separator=separator)
            except ValueError:
                # Not versioned - it cannot have been created by this upload alone.
                continue
            if version == version_number:
                names.add(name)
        return names

    return {
        SiteObservation.__name__: longcodes,
        CanonSite.__name__: created("canon_sites", "+"),
        CanonSiteConf.__name__: created("conformer_sites", "+"),
        XtalformSite.__name__: created("xtalform_sites", "/"),
    }


def _target_media_prefix(target) -> str:
    """The media-relative directory holding everything the loader wrote for a target."""
    return f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/{target.zip_archive}"


def _crystal_files(upload) -> dict[str, dict]:
    """``{crystal code: {yaml key: path or [paths]}}`` from an upload's metadata.

    Paths are exactly as the yaml spells them - relative to the target directory and
    carrying their own ``upload_N/`` prefix - which is what
    ``process_experiment`` joins onto the target directory to build the stored field.
    """
    path = Path(upload.get_upload_path()) / METADATA_FILE
    if not path.is_file():
        return {}
    meta = yaml.safe_load(path.read_text(encoding="utf-8")) or {}

    def one(value):
        """``(path, source_file)`` from either a file dict or a bare path string."""
        if isinstance(value, dict):
            return value.get("file"), value.get("source_file")
        return (value, None) if isinstance(value, str) else (None, None)

    crystals = {}
    for code, entry in (meta.get("crystals") or {}).items():
        files = (entry or {}).get("crystallographic_files") or {}
        found: dict = {}
        for key in EXPERIMENT_FILE_KEYS.values():
            name, source = one(files.get(key))
            if name:
                found[key] = (name, source)
        events = files.get(EXPERIMENT_MAP_KEY)
        if isinstance(events, list):
            pairs = [one(item) for item in events]
            usable = [(n, src) for n, src in pairs if n]
            if usable:
                found[EXPERIMENT_MAP_KEY] = usable
        if found:
            crystals[code] = found
    return crystals


def _restore_experiment_files(upload: ExperimentUpload, plan) -> None:
    """Point experiments back at the files they used before this upload overwrote them.

    A re-upload re-supplies a carried-over crystal's files and rewrites
    ``pdb_info``/``mtz_info``/``cif_info``/``map_info`` to its own directory, recording
    nothing about the previous value. Left alone, deleting the upload would either strand
    those rows on files inside a directory that should be gone, or force the sweep to
    keep the directory - which is only half a deletion.

    The previous value is recoverable because ``meta_aligner.yaml`` is cumulative: the
    newest surviving upload's metadata still lists every carried-over crystal with the
    path that was current before this upload. Uploads are searched newest first, and a
    candidate is accepted only when it points outside the doomed directory *and* the file
    is actually on disk.
    """
    target = upload.target
    prefix = f"{_target_media_prefix(target)}/{upload.upload_data_dir}/"
    media_root = Path(settings.MEDIA_ROOT)
    target_dir = _target_media_prefix(target)

    earlier = list(
        ExperimentUpload.objects.filter(target=target)
        .exclude(pk=upload.pk)
        .order_by("-upload_version")
    )
    if not earlier:
        return
    metadata: dict[int, dict] = {}

    def candidate(code, key):
        """The newest surviving ``(paths, sources)`` declared for this crystal's file.

        ``None`` means no upload that remains *declares* it - which is not a failure but
        a fact about the previous state: the row had no such file before this upload.

        Deliberately does NOT require the file to be on disk. What is being restored is
        the value the field held before this upload, and the metadata is the record of
        that value. An earlier version also demanded the file exist, which turned a
        restore into a clear whenever a file had gone missing for unrelated reasons -
        silently dropping a path the database should still carry. A missing file is
        reported instead, so the real problem stays visible.
        """
        for other in earlier:
            if other.pk not in metadata:
                metadata[other.pk] = _crystal_files(other)
            entry = (metadata[other.pk].get(code) or {}).get(key)
            if not entry:
                continue
            pairs = entry if isinstance(entry, list) else [entry]
            paths = [f"{target_dir}/{name}" for name, _ in pairs]
            if any(path.startswith(prefix) for path in paths):
                # Still inside the doomed directory: keep looking further back.
                continue
            sources = [source for _, source in pairs]
            if isinstance(entry, list):
                return paths, [x for x in sources if x]
            return paths[0], sources[0]
        return None

    def report_missing(exp, file_field, paths):
        absent = [p for p in paths if not (media_root / p).is_file()]
        if absent:
            plan.warnings.append(
                f"Experiment {exp.pk} ({exp.code}): {file_field} restored to "
                f"{absent[0]}, which is not on disk. The path is what the previous "
                "upload's metadata records; the missing file is a separate problem."
            )

    # Only experiments that survive: the rest go with the upload.
    survivors = Experiment.objects.filter(experiment_upload__target=target).exclude(
        experiment_upload=upload
    )

    for exp in survivors.iterator():
        restored, cleared = [], []

        for file_field, key in EXPERIMENT_FILE_KEYS.items():
            if not str(getattr(exp, file_field) or "").startswith(prefix):
                continue
            source_field = f"{file_field}_source_file"
            found = candidate(exp.code, key)
            if found is None:
                # No upload that remains declares this crystal's file, so before this
                # upload the row simply had none - clearing it restores that exactly.
                # Leaving it would strand the row on a directory that is being removed.
                setattr(exp, file_field, None)
                setattr(exp, source_field, None)
                cleared.extend([file_field, source_field])
                continue
            path, source = found
            setattr(exp, file_field, path)
            setattr(exp, source_field, source)
            restored.extend([file_field, source_field])
            report_missing(exp, file_field, [path])

        maps = [str(m) for m in (getattr(exp, EXPERIMENT_MAP_FIELD) or [])]
        if any(m.startswith(prefix) for m in maps):
            source_field = f"{EXPERIMENT_MAP_FIELD}_source_files"
            found = candidate(exp.code, EXPERIMENT_MAP_KEY)
            if found is None:
                setattr(exp, EXPERIMENT_MAP_FIELD, [])
                setattr(exp, source_field, [])
                cleared.extend([EXPERIMENT_MAP_FIELD, source_field])
            else:
                paths, sources = found
                setattr(exp, EXPERIMENT_MAP_FIELD, paths)
                setattr(exp, source_field, sorted(set(sources)))
                restored.extend([EXPERIMENT_MAP_FIELD, source_field])
                report_missing(exp, EXPERIMENT_MAP_FIELD, paths)

        if restored or cleared:
            exp.save(update_fields=sorted(set(restored + cleared)))
            if restored:
                plan.experiment_files_restored.append((exp.pk, sorted(set(restored))))
            if cleared:
                plan.experiment_files_cleared.append((exp.pk, sorted(set(cleared))))


def _scoped_ids(model, target) -> list[int]:
    """Primary keys of ``model`` belonging to ``target``.

    ``filter_manager.by_target()`` returns an *annotated* queryset and Django refuses
    ``.delete()`` on one, so everything downstream works from a plain list of pks.
    """
    return list(model.filter_manager.by_target(target).values_list("pk", flat=True))


def _resolve_declared(
    target, declared: dict[str, set[str]], version_number: int, plan
) -> dict[str, list[int]]:
    """Turn the metadata's natural keys into primary keys, scoped to the target.

    A declared key that matches no row is a warning rather than an error: the loader skips
    entries it cannot process (an observation with no parseable ligand is logged and
    dropped), so the bundle can legitimately declare more than the database holds.
    """
    lookups = {
        SiteObservation.__name__: (SiteObservation, "longcode__in"),
        CanonSite.__name__: (CanonSite, "name__in"),
        CanonSiteConf.__name__: (CanonSiteConf, "name__in"),
        XtalformSite.__name__: (XtalformSite, "xtalform_site_id__in"),
    }
    resolved: dict[str, list[int]] = {}
    for name, keys in declared.items():
        model, lookup = lookups[name]
        scoped = model.objects.filter(pk__in=_scoped_ids(model, target))
        if model is SiteObservation:
            # The longcode already carries the version, so it identifies the row alone.
            qs = scoped.filter(**{lookup: sorted(keys)})
        else:
            qs = scoped.filter(**{lookup: sorted(keys)}, version=version_number)
        resolved[name] = list(qs.values_list("pk", flat=True))
        missing = len(keys) - len(resolved[name])
        if missing > 0:
            plan.warnings.append(
                f"{missing} of {len(keys)} {name} entries declared by the bundle matched "
                "no row; the loader may have skipped them."
            )
    return resolved


def _cross_check_against_version(target, version_number: int, resolved) -> None:
    """Compare the metadata-derived sets against ``version == N``.

    Two independent derivations of the same thing: the bundle that created the rows, and
    the version stamp the loader wrote on them. They agreed on every target tested, so a
    disagreement means one of the two is wrong on this data - refuse rather than guess.
    """
    for model in DECLARED_MODELS:
        name = model.__name__
        by_version = set(
            model.objects.filter(
                pk__in=_scoped_ids(model, target), version=version_number
            ).values_list("pk", flat=True)
        )
        by_metadata = set(resolved[name])
        if by_version == by_metadata:
            continue
        only_meta = sorted(by_metadata - by_version)
        only_version = sorted(by_version - by_metadata)
        raise NotDeletable(
            f"{name}: the bundle's metadata and the stored version numbers disagree about "
            f"what upload {version_number} created. "
            f"Declared by metadata only: {only_meta or 'none'}. "
            f"Stamped version {version_number} only: {only_version or 'none'}. "
            "Refusing to delete on an ambiguous set."
        )


# --------------------------------------------------------------------------------------
# defusing the cascades
# --------------------------------------------------------------------------------------


def _chain_filter(obj) -> dict:
    """The supersede identity of ``obj``, as filter kwargs.

    Relations are compared by id so no extra rows are fetched. A NULL component (``cmpd``
    is nullable) yields ``IS NULL``, which pairs nothing - deliberately, so orphaned rows
    never chain to each other.
    """
    model = type(obj)
    kwargs = {}
    for name in model.SUPERSEDE_FIELDS:
        django_field = model._meta.get_field(name)  # pylint: disable=protected-access
        attname = django_field.attname if django_field.is_relation else name
        kwargs[attname] = getattr(obj, attname)
    return kwargs


def _predecessor(obj, exclude_ids, queryset=None):
    """The highest-versioned surviving row sharing ``obj``'s supersede identity."""
    model = type(obj)
    qs = model.objects.all() if queryset is None else queryset
    return (
        qs.filter(version__lt=obj.version, **_chain_filter(obj))
        .exclude(pk__in=exclude_ids)
        .order_by("-version")
        .first()
    )


def _repoint_pose_mains(doomed_obs: list[int], plan, commit: bool) -> None:
    """Give every pose whose main is being deleted its predecessor back.

    ``Pose.main_site_observation`` is ``OneToOne(CASCADE)``, so deleting the main would
    take the pose with it and ``SET_NULL`` the pose off its surviving members. Measured on
    a production dump, 421 poses have a doomed main; 240 of those empty out entirely and
    are removed in the orphan sweep, and every one of the 181 that survive has the doomed
    main's chain predecessor still among its own members.
    """
    poses = Pose.objects.filter(main_site_observation_id__in=doomed_obs).select_related(
        "main_site_observation"
    )
    for pose in poses:
        old = pose.main_site_observation
        survivors = SiteObservation.objects.filter(pose=pose).exclude(pk__in=doomed_obs)
        new_main = _predecessor(old, doomed_obs, queryset=survivors)
        if new_main is None:
            # No predecessor among the members. If nothing survives, the pose is empty and
            # the orphan sweep removes it - the cascade would have been harmless. If
            # something does survive, letting the cascade run would delete the pose and
            # orphan its members, so fall back to any surviving member.
            new_main = survivors.order_by("-version").first()
            if new_main is None:
                plan.poses_removed.append(pose.pk)
                continue
            plan.warnings.append(
                f"Pose {pose.pk} ({pose.display_name}): main {old.pk} has no surviving "
                f"predecessor in the pose; fell back to member {new_main.pk}."
            )
        plan.poses_repointed.append((pose.pk, old.pk, new_main.pk))
        if commit:
            pose.main_site_observation = new_main
            pose.save(update_fields=["main_site_observation"])


def _repoint_relinked_children(
    doomed: dict[str, list[int]], plan, commit: bool
) -> None:
    """Move surviving observations off doomed site rows, back onto the predecessors.

    The loader recomputes ``canon_site_conf`` and ``xtalform_site`` on every load, so an
    observation created by an *earlier* upload can end up pointing at a site row this
    upload created. Both FKs are ``CASCADE``, so deleting that row would destroy a
    surviving observation. Re-point it at the chain predecessor of the doomed row instead;
    where there is none the doomed row is kept (see :func:`_delete_declared`) rather than
    taking the observation with it.
    """
    for model, fk in (
        (CanonSiteConf, "canon_site_conf"),
        (XtalformSite, "xtalform_site"),
    ):
        doomed_ids = doomed[model.__name__]
        if not doomed_ids:
            continue
        stranded = SiteObservation.objects.filter(
            **{f"{fk}_id__in": doomed_ids}
        ).exclude(pk__in=doomed[SiteObservation.__name__])
        for obs in stranded.select_related(fk):
            old = getattr(obs, fk)
            new = _predecessor(old, doomed_ids)
            if new is None:
                plan.warnings.append(
                    f"SiteObservation {obs.pk} ({obs.longcode}) was re-linked by this "
                    f"upload to {model.__name__} {old.pk}, which has no earlier version; "
                    f"that {model.__name__} is kept so the observation survives."
                )
                continue
            plan.children_repointed.append((model.__name__, obs.pk, old.pk, new.pk))
            if commit:
                setattr(obs, fk, new)
                obs.save(update_fields=[fk])


def _clear_reference_fks(doomed: dict[str, list[int]], plan, commit: bool) -> None:
    """NULL the reference FKs of parents that are not themselves being deleted.

    ``CanonSiteConf.ref_site_observation`` and ``CanonSite.ref_conf_site`` are
    ``OneToOne(CASCADE)`` pointing *down* the graph, so deleting the referent destroys the
    parent - and for a CanonSite that cascades on into XtalformSite and into observations
    from earlier uploads.
    """
    for model, fk, referent in (
        (CanonSiteConf, "ref_site_observation", SiteObservation),
        (CanonSite, "ref_conf_site", CanonSiteConf),
    ):
        qs = model.objects.filter(
            **{f"{fk}_id__in": doomed[referent.__name__]}
        ).exclude(pk__in=doomed[model.__name__])
        for pk in qs.values_list("pk", flat=True):
            plan.refs_cleared.append((model.__name__, pk))
        if commit:
            qs.update(**{f"{fk}_id": None})


def _repoint_experiments(
    upload: ExperimentUpload, doomed_obs: list[int], plan, commit: bool
) -> None:
    """Re-point experiments that keep surviving observations, rather than deleting them.

    The new owner is the upload matching ``min(version)`` of the surviving observations,
    *not* simply the previous upload: "previous" is right for a single deletion but
    collapses under iteration, since the experiment would then be deleted one peel later.

    On bundles in scope this is a no-op - ``experiment_upload`` genuinely means "created
    by", so an experiment created by the latest upload can only hold that upload's rows. It
    is kept as a guard because the FK is not trustworthy on every row already in the
    database: before the fix in 6314c86b (2025-06-09) it was re-stamped wholesale on every
    load, and on such a target deleting the latest upload would otherwise take hundreds of
    experiments holding much older observations with it.
    """
    by_version = {
        other.upload_version: other
        for other in ExperimentUpload.objects.filter(target=upload.target).exclude(
            pk=upload.pk
        )
    }
    for exp in Experiment.objects.filter(experiment_upload=upload):
        earliest = (
            SiteObservation.objects.filter(experiment=exp)
            .exclude(pk__in=doomed_obs)
            .aggregate(v=Min("version"))["v"]
        )
        if earliest is None:
            # Nothing survives on this experiment: the cascade correctly takes it.
            continue
        owner = by_version.get(earliest)
        if owner is None:
            plan.warnings.append(
                f"Experiment {exp.pk} ({exp.code}) keeps observations at version "
                f"{earliest} but no upload has that version; it is left on the deleted "
                "upload and will be removed with it."
            )
            continue
        plan.experiments_repointed.append((exp.pk, owner.pk))
        if commit:
            exp.experiment_upload = owner
            exp.save(update_fields=["experiment_upload"])


# --------------------------------------------------------------------------------------
# deletion and the sweep
# --------------------------------------------------------------------------------------


def _delete_declared(doomed: dict[str, list[int]], plan) -> None:
    """Delete the declared rows, observations first.

    The site models are deleted only when nothing surviving still points at them.
    :func:`_repoint_relinked_children` moves observations off them first, so anything left
    holding on is a row with no earlier version to fall back to - keeping it is the lesser
    evil against cascading away a survivor.
    """
    SiteObservation.objects.filter(pk__in=doomed[SiteObservation.__name__]).delete()

    for model, probe in (
        (CanonSiteConf, "siteobservation"),
        (XtalformSite, "siteobservation"),
        (CanonSite, "canonsiteconf"),
    ):
        held = set(
            model.objects.filter(pk__in=doomed[model.__name__])
            .annotate(n=Count(probe))
            .filter(n__gt=0)
            .values_list("pk", flat=True)
        )
        if held:
            plan.warnings.append(
                f"{model.__name__} {sorted(held)} were created by this upload but still "
                "have surviving dependants; kept rather than cascading them away."
            )
        removable = [pk for pk in doomed[model.__name__] if pk not in held]
        model.objects.filter(pk__in=removable).delete()


def _collect_sweep_candidates(doomed: dict[str, list[int]]) -> dict[str, set[int]]:
    """The rows a deletion could plausibly empty: the parents of everything being deleted.

    Gathered *before* any row goes, because afterwards there is no way to tell a row the
    deletion emptied from one that was already empty. That distinction matters: a target
    can legitimately carry rows nothing points at - A71EV2A has 10 CanonSiteConfs with no
    site observation from its very first upload - and a deletion must leave those alone.
    """
    doomed_obs = doomed[SiteObservation.__name__]
    parents = SiteObservation.objects.filter(pk__in=doomed_obs).values_list(
        "canon_site_conf_id", "xtalform_site_id", "pose_id"
    )
    candidates: dict[str, set[int]] = {
        CanonSiteConf.__name__: set(),
        XtalformSite.__name__: set(),
        Pose.__name__: set(),
        CanonSite.__name__: set(),
        Xtalform.__name__: set(),
        QuatAssembly.__name__: set(),
        Compound.__name__: set(),
    }
    candidates[Compound.__name__].update(
        pk
        for pk in SiteObservation.objects.filter(pk__in=doomed_obs).values_list(
            "cmpd_id", flat=True
        )
        if pk
    )
    for conf_id, xtalform_site_id, pose_id in parents:
        if conf_id:
            candidates[CanonSiteConf.__name__].add(conf_id)
        if xtalform_site_id:
            candidates[XtalformSite.__name__].add(xtalform_site_id)
        if pose_id:
            candidates[Pose.__name__].add(pose_id)

    # A doomed site row's own parents can be emptied by its removal.
    conf_sites = (
        set(doomed[CanonSiteConf.__name__]) | candidates[CanonSiteConf.__name__]
    )
    candidates[CanonSite.__name__].update(
        CanonSiteConf.objects.filter(pk__in=conf_sites).values_list(
            "canon_site_id", flat=True
        )
    )
    xtalform_sites = (
        set(doomed[XtalformSite.__name__]) | candidates[XtalformSite.__name__]
    )
    for canon_site_id, xtalform_id in XtalformSite.objects.filter(
        pk__in=xtalform_sites
    ).values_list("canon_site_id", "xtalform_id"):
        if canon_site_id:
            candidates[CanonSite.__name__].add(canon_site_id)
        if xtalform_id:
            candidates[Xtalform.__name__].add(xtalform_id)
    return candidates


def _unreferenced_compounds(candidate_pks: set[int], project) -> list[int]:
    """Of ``candidate_pks``, the compounds nothing in the database points at any more.

    Every incoming relation is walked through Django's own metadata rather than a hand
    written list, so a model added later cannot quietly make this wrong; the M2M through
    tables are resolved the way :mod:`viewer.compound_dedup` resolves them, so the sweep
    and the dedup report can never disagree about what counts as a referrer.

    Restricted to ``project`` for the same reason the loader's own compound retirement is
    (``_apply_compound_supersessions``): a compound belonging to someone else is not ours
    to remove.
    """
    remaining = set(
        Compound.objects.filter(pk__in=candidate_pks, project=project).values_list(
            "pk", flat=True
        )
    )
    for rel in Compound._meta.related_objects:  # pylint: disable=protected-access
        if not remaining:
            break
        if rel.many_to_many:
            through = rel.through
            attname = _through_compound_fk_attname(through)
            used = through.objects.filter(**{f"{attname}__in": remaining}).values_list(
                attname, flat=True
            )
        else:
            related = rel.related_model
            if related._meta.auto_created:  # pylint: disable=protected-access
                # An auto-created through model: reached via its own m2m entry above.
                continue
            used = related.objects.filter(
                **{f"{rel.field.name}__in": remaining}
            ).values_list(rel.field.attname, flat=True)
        remaining -= set(used)
    return sorted(remaining)


def _sweep_orphans(target, candidates: dict[str, set[int]], plan) -> None:
    """Remove what the deletion emptied out - and nothing else.

    Nothing below SiteObservation carries upload provenance, and none is needed: a row goes
    when the deletion took away the last thing pointing at it. The sweep is restricted to
    ``candidates`` (see :func:`_collect_sweep_candidates`) so that rows which were already
    unreferenced before the deletion are left exactly as they were found. The order
    matters - each sweep can orphan the next one up, and compounds come last because the
    sweeps above are themselves referrers.
    """
    empty_poses = list(
        Pose.objects.filter(
            pk__in=set(_scoped_ids(Pose, target)) & candidates[Pose.__name__]
        )
        .annotate(n=Count("site_observations"))
        .filter(n=0)
        .values_list("pk", flat=True)
    )
    plan.poses_removed.extend(empty_poses)
    Pose.objects.filter(pk__in=empty_poses).delete()

    sweeps: tuple[tuple[type[Model], tuple[str, ...]], ...] = (
        (CanonSiteConf, ("siteobservation",)),
        (XtalformSite, ("siteobservation",)),
        (CanonSite, ("canonsiteconf", "xtalformsite", "pose")),
        (Xtalform, ("experiment", "xtalformsite")),
        (QuatAssembly, ("xtalformquatassembly",)),
    )
    for model, probes in sweeps:
        scoped = set(_scoped_ids(model, target)) & candidates[model.__name__]
        if not scoped:
            continue
        qs = model.objects.filter(pk__in=scoped)
        for i, probe in enumerate(probes):
            qs = qs.annotate(**{f"n{i}": Count(probe)})
        condition = Q()
        for i in range(len(probes)):
            condition &= Q(**{f"n{i}": 0})
        pks = list(qs.filter(condition).values_list("pk", flat=True))
        if not pks:
            continue
        plan.orphaned.setdefault(model.__name__, []).extend(pks)
        # Removing these can in turn empty their own parents.
        if model is CanonSiteConf:
            candidates[CanonSite.__name__].update(
                model.objects.filter(pk__in=pks).values_list("canon_site_id", flat=True)
            )
        elif model is XtalformSite:
            for canon_site_id, xtalform_id in model.objects.filter(
                pk__in=pks
            ).values_list("canon_site_id", "xtalform_id"):
                candidates[CanonSite.__name__].add(canon_site_id)
                candidates[Xtalform.__name__].add(xtalform_id)
        elif model is Xtalform:
            candidates[QuatAssembly.__name__].update(
                XtalformQuatAssembly.objects.filter(xtalform_id__in=pks).values_list(
                    "quat_assembly_id", flat=True
                )
            )
        model.objects.filter(pk__in=pks).delete()

    # Compounds last: the sweeps above (poses especially) are themselves referrers, so a
    # compound only becomes unreferenced once they have gone. Unlike the site models a
    # compound is project-scoped rather than target-scoped, so it gets its own check
    # across every relation instead of a per-target probe.
    compound_pks = _unreferenced_compounds(
        candidates[Compound.__name__], target.project
    )
    if compound_pks:
        plan.orphaned.setdefault(Compound.__name__, []).extend(compound_pks)
        Compound.objects.filter(pk__in=compound_pks).delete()


def restore_superseded(target) -> dict[str, int]:
    """Re-derive ``superseded`` from the rows that remain.

    A row is superseded exactly when another surviving row shares its supersede identity
    and carries a higher version. This is why no undo journal is needed: the loader's
    ``.update(superseded=True)`` records only a count (and, being a queryset update,
    bypasses simple_history), but the flag is a pure function of the rows that exist.
    Checked against a production dump: 0 mismatches over 11,033 SiteObservations, 757
    CanonSiteConfs and 1,023 XtalformSites.

    Returns the number of rows whose flag actually changed, per model.
    """
    changed = {}
    for model in RESTORED_MODELS:
        ids = _scoped_ids(model, target)
        scoped = model.objects.filter(pk__in=ids)
        later = model.objects.filter(
            pk__in=ids,
            version__gt=OuterRef("version"),
            **{f: OuterRef(f) for f in model.SUPERSEDE_FIELDS},
        )
        # Resolve to pks first: an UPDATE whose predicate is a correlated subquery on the
        # same table is needlessly fragile, and the row counts here are small.
        superseded_ids = list(scoped.filter(Exists(later)).values_list("pk", flat=True))
        n = model.objects.filter(pk__in=superseded_ids, superseded=False).update(
            superseded=True
        )
        n += (
            model.objects.filter(pk__in=ids, superseded=True)
            .exclude(pk__in=superseded_ids)
            .update(superseded=False)
        )
        changed[model.__name__] = n
    return changed


# --------------------------------------------------------------------------------------
# media
# --------------------------------------------------------------------------------------


def _file_values(instance) -> list[str]:
    """Every stored file path on ``instance``, including ArrayField-of-FileField."""
    values: list[str] = []
    for f in type(instance)._meta.get_fields():  # pylint: disable=protected-access
        if not getattr(f, "concrete", False):
            continue
        if isinstance(f, ArrayField) and hasattr(f.base_field, "upload_to"):
            values.extend(str(v) for v in (getattr(instance, f.name) or []) if v)
        elif hasattr(f, "upload_to"):
            value = getattr(instance, f.name)
            if value:
                values.append(str(value))
    return values


def _referenced_paths(target) -> set[Path]:
    """Absolute paths of every target-loader file still referenced by a surviving row.

    Computed after the database work, inside the same call but before any file is touched,
    so the answer describes the state the deletion actually leaves behind.
    """
    media_root = Path(settings.MEDIA_ROOT)
    referenced: set[Path] = set()
    querysets = (
        SiteObservation.objects.filter(pk__in=_scoped_ids(SiteObservation, target)),
        Experiment.objects.filter(experiment_upload__target=target),
    )
    for qs in querysets:
        for instance in qs.iterator():
            for value in _file_values(instance):
                path = Path(value)
                referenced.add(path if path.is_absolute() else media_root / path)
    return referenced


def _bundle_archive_paths(upload: ExperimentUpload) -> list[Path]:
    """Every place this upload's bundle archive could be sitting.

    ``ExperimentUpload.file`` records only the archive's basename
    (``target_loader.py``: ``file=self.data_bundle``, and ``data_bundle`` is
    ``Path(...).name``), while ``_move_and_save_target_experiment`` moves the archive
    into the *target's* directory. ``get_download_path()`` therefore resolves one
    directory too high, and a delete driven by it silently removes nothing - which is
    how six bundles totalling 8.7 GB survived a full peel.

    The same move renames the archive to ``<stem>_<version_dir><suffix>`` when a file of
    that name is already there, again without recording the new name, so that spelling
    has to be tried too.
    """
    candidates = [Path(upload.get_download_path())]

    name = Path(str(upload.file)).name
    if name and upload.target.zip_archive:
        target_dir = (
            Path(settings.MEDIA_ROOT)
            / settings.TARGET_LOADER_MEDIA_DIRECTORY
            / str(upload.target.zip_archive)
        )
        stem, suffix = Path(name).stem, Path(name).suffix
        candidates.append(target_dir / name)
        if upload.upload_data_dir:
            candidates.append(target_dir / f"{stem}_{upload.upload_data_dir}{suffix}")

    seen: dict[Path, None] = {}
    for path in candidates:
        seen.setdefault(path, None)
    return list(seen)


def _remove_upload_media(upload: ExperimentUpload, plan) -> None:
    """Remove the deleted upload's directory, and its bundle archive.

    Only the deleted upload's own directory is swept. An earlier version also swept every
    surviving upload's directory, on the theory that a file kept because some row
    referenced it would otherwise leak once that row went. That was over-reach, and it
    broke repeated deletion: a surviving upload's crystallographic file is unreferenced
    precisely while a later upload's copy is the one in use, and it is exactly the file
    :func:`_restore_experiment_files` must point back at when that later upload is
    deleted. Sweeping it left the next peel with nothing to restore to.

    A surviving upload's directory is that upload's payload; it is not this deletion's to
    tidy. Anything genuinely orphaned there is the ``cleanup_media`` command's business.

    Runs after the database transaction has committed, so a database failure can never
    leave deleted files behind live rows.
    """
    target = upload.target
    referenced = _referenced_paths(target)
    deleted_dir = (
        Path(settings.MEDIA_ROOT)
        / settings.TARGET_LOADER_MEDIA_DIRECTORY
        / str(target.zip_archive)
        / (upload.upload_data_dir or "")
        if target.zip_archive and upload.upload_data_dir
        else None
    )

    if deleted_dir is not None and deleted_dir.is_dir():
        for path in sorted(deleted_dir.rglob("*")):
            if not path.is_file():
                continue
            if path in referenced:
                # A row still needs it, so removing it would break that row. This is an
                # anomaly now that file paths are restored, and it is reported below.
                plan.files_kept.append(path)
                continue
            try:
                path.unlink()
                plan.files_removed.append(path)
            except OSError:
                plan.warnings.append(f"Could not remove {path}.")

    if deleted_dir is not None and deleted_dir.is_dir():
        # Prune what is now empty; anything left is a file a surviving row still needs.
        for path in sorted(deleted_dir.rglob("*"), reverse=True):
            if path.is_dir() and not any(path.iterdir()):
                path.rmdir()
        if not any(deleted_dir.iterdir()):
            deleted_dir.rmdir()
        else:
            left = sorted(p for p in deleted_dir.rglob("*") if p.is_file())
            plan.warnings.append(
                f"{deleted_dir} could not be removed: {len(left)} file(s) in it are "
                "still referenced by surviving rows. Experiment file paths are restored "
                "from the previous upload's metadata, so this means that metadata did "
                f"not supply a replacement. First: {left[:3]}"
            )

    for archive in _bundle_archive_paths(upload):
        if not archive.is_file():
            continue
        try:
            archive.unlink()
            plan.files_removed.append(archive)
        except OSError:
            plan.warnings.append(f"Could not remove bundle archive {archive}.")


# --------------------------------------------------------------------------------------
# orchestration
# --------------------------------------------------------------------------------------


def _resolve(upload: ExperimentUpload, plan) -> dict[str, list[int]]:
    """The doomed set, from the bundle's metadata, cross-checked against the database."""
    declared = read_upload_metadata(upload)
    resolved = _resolve_declared(upload.target, declared, upload.upload_version, plan)
    _cross_check_against_version(upload.target, upload.upload_version, resolved)
    return resolved


def plan_upload_deletion(upload: ExperimentUpload) -> UploadDeletionPlan:
    """Work out what deleting ``upload`` would do, without changing anything.

    Drives ``--dry-run``. Everything that inspects the database is shared with
    :func:`delete_latest_upload`; only the writes are skipped.
    """
    _check_deletable(upload)
    plan = UploadDeletionPlan(
        upload_id=upload.pk,
        upload_version=upload.upload_version,
        target_title=upload.target.title,
    )
    plan.doomed = _resolve(upload, plan)
    doomed_obs = plan.doomed[SiteObservation.__name__]
    _repoint_experiments(upload, doomed_obs, plan, commit=False)
    _repoint_pose_mains(doomed_obs, plan, commit=False)
    # Read-only: this only inspects metadata and the filesystem.
    with transaction.atomic():
        _restore_experiment_files(upload, plan)
        transaction.set_rollback(True)
    _repoint_relinked_children(plan.doomed, plan, commit=False)
    _clear_reference_fks(plan.doomed, plan, commit=False)
    return plan


def delete_latest_upload(upload: ExperimentUpload) -> UploadDeletionPlan:
    """Delete the target's latest upload and restore the state before it.

    The caller is responsible for authorization and for the deployment-mode check.

    :param upload: the ExperimentUpload to delete; must be the latest of a target that has
        more than one.
    :raises NotDeletable: if the upload is not the latest, is the only one, or its bundle
        metadata is missing, in an unsupported format, or disagrees with the database.
    """
    _check_deletable(upload)
    target = upload.target
    plan = UploadDeletionPlan(
        upload_id=upload.pk,
        upload_version=upload.upload_version,
        target_title=target.title,
    )
    logger.info(
        "Deleting upload %s (version %s) of target %s",
        upload.pk,
        upload.upload_version,
        target.title,
    )

    plan.doomed = _resolve(upload, plan)
    doomed_obs = plan.doomed[SiteObservation.__name__]
    # Gathered before anything is deleted: afterwards a row the deletion emptied is
    # indistinguishable from one that was already empty.
    candidates = _collect_sweep_candidates(plan.doomed)

    # Everything removed from here on was asked for, so keep it out of the media
    # deletion watcher's log - that watcher exists to catch files vanishing when
    # nothing should have touched them, and a peel removes thousands at a time.
    with deletions_expected(f"delete upload {upload.upload_version} of {target.title}"):
        with transaction.atomic():
            _repoint_experiments(upload, doomed_obs, plan, commit=True)
            _restore_experiment_files(upload, plan)
            _repoint_pose_mains(doomed_obs, plan, commit=True)
            _repoint_relinked_children(plan.doomed, plan, commit=True)
            _clear_reference_fks(plan.doomed, plan, commit=True)
            _delete_declared(plan.doomed, plan)
            _sweep_orphans(target, candidates, plan)
            plan.restored = restore_superseded(target)
            upload.delete()

        # post_delete signals are bypassed by cascade and queryset deletes alike, so
        # the view cache has to be cleared explicitly - the loader does the same after
        # a load.
        clear_view_cache("tag", "pose", "site-observation")

        # Files last: the database is already consistent, so a filesystem problem
        # degrades to debris rather than to rows pointing at things that are gone.
        # Inside the suppression, because this is where the bulk of the removals are.
        _remove_upload_media(upload, plan)

    logger.info(
        "Deleted upload %s of target %s: %s",
        plan.upload_id,
        plan.target_title,
        {k: len(v) for k, v in plan.doomed.items()},
    )
    return plan
