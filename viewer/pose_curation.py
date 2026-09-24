"""Find and fix poses damaged by superseding.

This is a *curation toolkit*, not a migration and not a management command. It
is meant to be driven from ``python manage.py shell`` alongside ad-hoc code
while the picture is still forming::

    from viewer.pose_curation import *

    print(survey())                      # whole database
    print(survey("A71EV2A"))             # one target

    for d in diagnose_superseded_mains("A71EV2A"):
        print(d.explain())

    plan = plan_fix_superseded_mains("A71EV2A")
    print(plan)                          # read-only: nothing written yet
    plan.apply()                         # now it writes, in one transaction

Every ``find_*``/``diagnose_*`` function is read-only. Every ``plan_*`` function
is read-only too: it returns a :class:`Plan`, which prints what it *would* do
and only touches the database when you call :meth:`Plan.apply`. That two-step
split is deliberate - there is no ``dry_run`` flag to forget.

Finders return querysets wherever they can, so you can keep chaining in the
shell (``find_superseded_members("A71EV2A").filter(code__startswith="x0")``).

Why this exists
---------------
``TargetLoader.supersedes`` records new->old pairings during an upload and
``_refresh_poses`` consumes them, but that map is per-upload and in memory: it
can do nothing for rows already in the database. Repairing the backlog means
*reconstructing* the pairing, and the only honest way to do that is the model's
own :attr:`SiteObservation.SUPERSEDE_FIELDS` - deliberately **not** the
conformer site or xtalform site, both of which are recomputed between uploads
and were what made the original lookup unrecoverable.

Reconstruction is a guess where the loader had a fact, so the planners here
refuse to act on any case they cannot call unambiguously. Those land in
``plan.skipped`` with a reason rather than being quietly repaired; reading that
list is the point of the exercise.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Iterable

from django.db import transaction
from django.db.models import Count, F, Max, Q, QuerySet

from viewer.models import Pose, SiteObservation, Target

logger = logging.getLogger(__name__)

__all__ = [
    "resolve_target",
    "identity_of",
    "find_successors",
    "find_poses_with_superseded_main",
    "find_superseded_members",
    "find_poses_with_superseded_members",
    "find_main_not_in_pose",
    "find_poses_without_main",
    "find_dead_poses",
    "find_orphan_observations",
    "find_split_groups",
    "find_anchor_violations",
    "find_display_name_mismatches",
    "find_duplicate_live_identities",
    "find_version_inversions",
    "find_wrongly_superseded",
    "survey",
    "diagnose_superseded_main",
    "diagnose_superseded_mains",
    "plan_fix_superseded_mains",
    "plan_sever_superseded_members",
    "plan_unsupersede",
    "plan_repair_display_names",
    "plan_delete_dead_poses",
    "Plan",
]


# --------------------------------------------------------------------------
# scoping helpers
# --------------------------------------------------------------------------


def resolve_target(target: Target | str | int | None) -> Target | None:
    """Accept a Target, a title, a pk or None (meaning "the whole database")."""
    if target is None or isinstance(target, Target):
        return target
    if isinstance(target, int):
        return Target.objects.get(pk=target)
    return Target.objects.get(title=target)


def _poses(target=None) -> QuerySet:
    target = resolve_target(target)
    if target is None:
        return Pose.objects.all()
    # the manager owns the definition of "belongs to this target"; don't
    # reinvent the join path here
    return Pose.filter_manager.by_target(target)


def _observations(target=None) -> QuerySet:
    target = resolve_target(target)
    if target is None:
        return SiteObservation.objects.all()
    return SiteObservation.filter_manager.by_target(target)


# --------------------------------------------------------------------------
# supersede identity
# --------------------------------------------------------------------------

#: Column names behind :attr:`SiteObservation.SUPERSEDE_FIELDS`, so a FK is
#: compared as ``experiment_id`` and never triggers a fetch.
IDENTITY_ATTNAMES: tuple[str, ...] = tuple(
    SiteObservation._meta.get_field(name).attname  # pylint: disable=protected-access
    for name in SiteObservation.SUPERSEDE_FIELDS
)

#: Identity fields whose NULL makes the identity meaningless rather than merely
#: unset. ``cmpd`` is the live case: the eight orphan observations have a NULL
#: compound, and ``filter(cmpd_id=None)`` compiles to ``cmpd_id IS NULL``, so
#: without this guard they would pair with each other. ``seq_id``/``chain_id``
#: NULLs are left alone - they match consistently and mean what they say.
IDENTITY_REQUIRED: tuple[str, ...] = ("experiment_id", "cmpd_id")


def identity_of(obs: SiteObservation) -> dict[str, Any]:
    """The supersede identity of ``obs`` as a filter kwargs dict."""
    return {attname: getattr(obs, attname) for attname in IDENTITY_ATTNAMES}


def _identity_gap(ident: dict[str, Any]) -> list[str]:
    return [name for name in IDENTITY_REQUIRED if ident.get(name) is None]


def find_successors(
    obs: SiteObservation,
    *,
    live_only: bool = True,
    allow_null_identity: bool = False,
) -> QuerySet:
    """Rows sharing ``obs``'s supersede identity, newest first.

    This is the reconstruction the loader does not have to do. It matches on
    ``SUPERSEDE_FIELDS`` only; anything derived from the site model is excluded
    on purpose.

    Ordered ``-version, -pk`` so ``.first()`` is the immediate successor - but
    see :func:`find_version_inversions`, which tells you whether that ordering
    is actually safe in your data.
    """
    ident = identity_of(obs)
    if _identity_gap(ident) and not allow_null_identity:
        return SiteObservation.objects.none()
    qs = SiteObservation.objects.filter(**ident).exclude(pk=obs.pk)
    if live_only:
        qs = qs.filter(superseded=False)
    return qs.order_by("-version", "-pk")


# --------------------------------------------------------------------------
# finders - all read-only
# --------------------------------------------------------------------------


def find_poses_with_superseded_main(target=None) -> QuerySet:
    """Poses whose main observation is flagged superseded (the headline bug)."""
    return _poses(target).filter(main_site_observation__superseded=True)


def find_superseded_members(target=None) -> QuerySet:
    """Superseded observations still sitting in a pose.

    These are what the frontend serves as pose members. ``_refresh_poses`` now
    severs them on upload; this finds the ones already in the database.
    """
    return _observations(target).filter(superseded=True, pose__isnull=False)


def find_poses_with_superseded_members(target=None) -> QuerySet:
    return _poses(target).filter(site_observations__superseded=True).distinct()


def find_main_not_in_pose(target=None) -> QuerySet:
    """Poses whose main observation is not a member of that same pose."""
    return (
        _poses(target)
        .exclude(main_site_observation__pose_id=F("pk"))
        .filter(main_site_observation__isnull=False)
    )


def find_poses_without_main(target=None) -> QuerySet:
    return _poses(target).filter(main_site_observation__isnull=True)


def find_dead_poses(target=None) -> QuerySet:
    """Poses with no live observation in them at all - invisible in the FE."""
    return (
        _poses(target)
        .annotate(
            live_members=Count(
                "site_observations", filter=Q(site_observations__superseded=False)
            ),
        )
        .filter(live_members=0)
    )


def find_orphan_observations(target=None, *, live_only: bool = True) -> QuerySet:
    """Observations with no pose.

    Restricted to rows that came from an experiment: computed-set and virtual
    observations legitimately have neither experiment nor pose.
    """
    qs = _observations(target).filter(pose__isnull=True, experiment__isnull=False)
    if live_only:
        qs = qs.filter(superseded=False)
    return qs


def find_split_groups(target=None) -> QuerySet:
    """``(canon_site, compound)`` groups holding more than one pose.

    ``_generate_poses`` handles a group one pose at a time, so existing poses in
    a split group are never revisited.
    """
    return (
        _poses(target)
        .values("canon_site_id", "compound_id")
        .annotate(
            poses=Count("id"),
        )
        .filter(poses__gt=1)
        .order_by("-poses")
    )


def find_anchor_violations(target=None) -> QuerySet:
    """Members whose canon site is not their pose's canon site.

    A pose is anchored on its canon site, so this should never happen. Members
    with no conformer site at all are included - they cannot satisfy the anchor
    either.
    """
    return (
        _observations(target)
        .filter(pose__isnull=False)
        .exclude(
            canon_site_conf__canon_site_id=F("pose__canon_site_id"),
        )
    )


def find_display_name_mismatches(target=None) -> list[dict[str, Any]]:
    """Poses whose ``display_name`` names an observation that is not in them.

    Walks the poses in Python - there is no join that expresses "the name is one
    of the member codes". Fine for a curation tool, slow if you aim it at the
    whole database.
    """
    out = []
    qs = _poses(target).prefetch_related("site_observations")
    for pose in qs:
        if not pose.display_name:
            continue
        codes = {so.code for so in pose.site_observations.all()}
        if pose.display_name not in codes:
            out.append(
                {
                    "pose_id": pose.pk,
                    "display_name": pose.display_name,
                    "member_codes": sorted(c for c in codes if c),
                }
            )
    return out


def find_duplicate_live_identities(target=None) -> QuerySet:
    """Two *live* rows sharing a full supersede identity.

    This should be impossible: it means superseding never fired for them.
    """
    return (
        _observations(target)
        .filter(superseded=False)
        .values(*IDENTITY_ATTNAMES)
        .annotate(live_rows=Count("id"))
        .filter(live_rows__gt=1)
        .order_by("-live_rows")
    )


def find_version_inversions(target=None) -> QuerySet:
    """Identity groups where a superseded row outranks the live one by version.

    If this returns anything, picking a survivor with ``order_by("-version")``
    is unsafe - ``supersede_fields`` carries no version constraint, so a
    lower-version upload can supersede a higher-version row.
    """
    return (
        _observations(target)
        .values(*IDENTITY_ATTNAMES)
        .annotate(
            max_superseded=Max("version", filter=Q(superseded=True)),
            max_live=Max("version", filter=Q(superseded=False)),
        )
        .filter(max_superseded__gt=F("max_live"))
    )


def find_wrongly_superseded(target=None) -> QuerySet:
    """Superseded rows that nothing live ever replaced.

    Identity groups where every row is flagged superseded. Either the flag was
    set in error, or the successor was later removed, or the crystal really was
    withdrawn - this cannot tell you which, which is why
    :func:`plan_unsupersede` takes explicit pks rather than a target.
    """
    dead_groups = (
        _observations(target)
        .values(*IDENTITY_ATTNAMES)
        .annotate(
            live_rows=Count("id", filter=Q(superseded=False)),
            dead_rows=Count("id", filter=Q(superseded=True)),
        )
        .filter(live_rows=0, dead_rows__gt=0)
    )

    q = Q()
    for group in dead_groups:
        q |= Q(**{name: group[name] for name in IDENTITY_ATTNAMES})
    if not q:
        return SiteObservation.objects.none()
    return _observations(target).filter(q, superseded=True).order_by("pk")


def survey(target=None) -> "Survey":
    """Count every defect class at once. Read-only."""
    return Survey(
        target=str(resolve_target(target) or "<all targets>"),
        counts={
            "poses": _poses(target).count(),
            "observations": _observations(target).count(),
            "poses with superseded main": find_poses_with_superseded_main(
                target
            ).count(),
            "superseded members (rows)": find_superseded_members(target).count(),
            "poses holding superseded members": find_poses_with_superseded_members(
                target
            ).count(),
            "main not in its own pose": find_main_not_in_pose(target).count(),
            "poses with no main": find_poses_without_main(target).count(),
            "dead poses (no live member)": find_dead_poses(target).count(),
            "live orphan observations": find_orphan_observations(target).count(),
            "split (canon_site, compound) groups": find_split_groups(target).count(),
            "anchor violations": find_anchor_violations(target).count(),
            "display_name mismatches": len(find_display_name_mismatches(target)),
            "duplicate live identities": find_duplicate_live_identities(target).count(),
            "version inversions": find_version_inversions(target).count(),
            "wrongly superseded (no live successor)": find_wrongly_superseded(
                target
            ).count(),
        },
    )


@dataclass
class Survey:
    target: str
    counts: dict[str, int]

    def __str__(self) -> str:
        width = max(len(k) for k in self.counts)
        lines = [f"pose curation survey - {self.target}", "=" * (width + 10)]
        lines.extend(f"{k.ljust(width)}  {v:>6}" for k, v in self.counts.items())
        return "\n".join(lines)

    __repr__ = __str__


# --------------------------------------------------------------------------
# diagnosis of a superseded main
# --------------------------------------------------------------------------

SUCCESSOR_SAME_POSE = "successor_same_pose"
SUCCESSOR_OTHER_POSE = "successor_other_pose"
SUCCESSOR_UNPOSED = "successor_unposed"
NO_SUCCESSOR_LIVE_MEMBERS = "no_successor_live_members"
NO_SUCCESSOR_DEAD = "no_successor_dead"
NULL_IDENTITY = "null_identity"

#: Kinds :func:`plan_fix_superseded_mains` will act on unaided.
UNAMBIGUOUS = frozenset(
    {SUCCESSOR_SAME_POSE, SUCCESSOR_UNPOSED, NO_SUCCESSOR_LIVE_MEMBERS}
)


@dataclass
class MainDiagnosis:
    """Why one pose's main is superseded, and what could replace it."""

    pose: Pose
    main: SiteObservation
    kind: str
    successor: SiteObservation | None
    live_members: list[SiteObservation]
    note: str = ""

    @property
    def actionable(self) -> bool:
        return self.kind in UNAMBIGUOUS

    def explain(self) -> str:
        lines = [
            f"pose {self.pose.pk} ({self.pose.display_name!r})  ->  {self.kind}",
            f"  main      {self.main.pk} {self.main.code!r} "
            f"v{self.main.version} superseded={self.main.superseded}",
        ]
        if self.successor is not None:
            lines.append(
                f"  successor {self.successor.pk} {self.successor.code!r} "
                f"v{self.successor.version} pose={self.successor.pose_id}"
            )
        else:
            lines.append("  successor none")
        lines.append(
            "  live members "
            + (
                ", ".join(f"{m.pk}:{m.code}" for m in self.live_members)
                if self.live_members
                else "none"
            )
        )
        if self.note:
            lines.append(f"  note      {self.note}")
        return "\n".join(lines)

    __str__ = explain
    __repr__ = explain


def diagnose_superseded_main(pose: Pose) -> MainDiagnosis:
    """Classify a single pose whose main is superseded."""
    main = pose.main_site_observation
    live_members = list(
        pose.site_observations.filter(superseded=False).order_by("-version", "pk")
    )

    gap = _identity_gap(identity_of(main))
    if gap:
        return MainDiagnosis(
            pose=pose,
            main=main,
            kind=NULL_IDENTITY,
            successor=None,
            live_members=live_members,
            note=f"identity unusable, NULL: {', '.join(gap)}",
        )

    successor = find_successors(main).first()

    if successor is None:
        if live_members:
            return MainDiagnosis(
                pose=pose,
                main=main,
                kind=NO_SUCCESSOR_LIVE_MEMBERS,
                successor=None,
                live_members=live_members,
                note="nothing replaced the main, but the pose still has live "
                "members to promote from",
            )
        return MainDiagnosis(
            pose=pose,
            main=main,
            kind=NO_SUCCESSOR_DEAD,
            successor=None,
            live_members=[],
            note="nothing live anywhere in this pose; either the main was "
            "wrongly superseded or the pose should go",
        )

    if successor.pose_id == pose.pk:
        kind, note = SUCCESSOR_SAME_POSE, ""
    elif successor.pose_id is None:
        kind, note = (
            SUCCESSOR_UNPOSED,
            "successor holds no pose; it can be adopted into this one",
        )
    else:
        kind, note = (
            SUCCESSOR_OTHER_POSE,
            f"successor already belongs to pose {successor.pose_id} - a canon-site "
            "hop or a split group. Needs a human decision.",
        )

    return MainDiagnosis(
        pose=pose,
        main=main,
        kind=kind,
        successor=successor,
        live_members=live_members,
        note=note,
    )


def diagnose_superseded_mains(target=None) -> list[MainDiagnosis]:
    """Diagnose every pose with a superseded main, worst cases last."""
    qs = find_poses_with_superseded_main(target).select_related("main_site_observation")
    return sorted(
        (diagnose_superseded_main(p) for p in qs),
        key=lambda d: (d.actionable, d.kind, d.pose.pk),
    )


# --------------------------------------------------------------------------
# actions and plans
# --------------------------------------------------------------------------


@dataclass
class SetPoseMain:
    pose_id: int
    new_main_id: int
    old_main_id: int | None

    destructive = False

    def describe(self) -> str:
        return f"pose {self.pose_id}: main {self.old_main_id} -> {self.new_main_id}"

    def apply(self) -> None:
        Pose.objects.filter(pk=self.pose_id).update(
            main_site_observation_id=self.new_main_id
        )


@dataclass
class SetPoseMember:
    obs_id: int
    pose_id: int | None
    was_pose_id: int | None

    destructive = False

    def describe(self) -> str:
        return f"observation {self.obs_id}: pose {self.was_pose_id} -> {self.pose_id}"

    def apply(self) -> None:
        SiteObservation.objects.filter(pk=self.obs_id).update(pose_id=self.pose_id)


@dataclass
class SetSuperseded:
    obs_id: int
    value: bool

    destructive = False

    def describe(self) -> str:
        return f"observation {self.obs_id}: superseded -> {self.value}"

    def apply(self) -> None:
        SiteObservation.objects.filter(pk=self.obs_id).update(superseded=self.value)


@dataclass
class SetDisplayName:
    pose_id: int
    value: str | None
    was: str | None

    destructive = False

    def describe(self) -> str:
        return f"pose {self.pose_id}: display_name {self.was!r} -> {self.value!r}"

    def apply(self) -> None:
        Pose.objects.filter(pk=self.pose_id).update(display_name=self.value)


@dataclass
class DeletePose:
    pose_id: int

    destructive = True

    def describe(self) -> str:
        return f"pose {self.pose_id}: DELETE"

    def apply(self) -> None:
        Pose.objects.filter(pk=self.pose_id).delete()


@dataclass
class Plan:
    """A set of intended writes. Printing it writes nothing."""

    name: str
    actions: list[Any] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)

    @property
    def destructive(self) -> bool:
        return any(a.destructive for a in self.actions)

    def __str__(self) -> str:
        lines = [
            f"Plan: {self.name}",
            f"  {len(self.actions)} action(s), {len(self.skipped)} skipped",
        ]
        if self.destructive:
            lines.append("  *** DESTRUCTIVE - apply(allow_destructive=True) ***")
        lines.append("")
        for action in self.actions:
            lines.append(f"  + {action.describe()}")
        if self.skipped:
            lines.append("")
            lines.append("  skipped:")
            for what, why in self.skipped:
                lines.append(f"  - {what}: {why}")
        return "\n".join(lines)

    __repr__ = __str__

    def apply(self, *, allow_destructive: bool = False) -> dict[str, int]:
        """Run every action in one transaction. This writes."""
        if self.destructive and not allow_destructive:
            raise RuntimeError(
                f"{self.name} contains destructive actions; "
                "call apply(allow_destructive=True) if that is what you want"
            )
        counts: dict[str, int] = {}
        with transaction.atomic():
            for action in self.actions:
                logger.info("pose curation [%s]: %s", self.name, action.describe())
                action.apply()
                key = type(action).__name__
                counts[key] = counts.get(key, 0) + 1
        return counts


def _mains_taken() -> set[int]:
    """Observation pks already serving as some pose's main.

    ``Pose.main_site_observation`` is a OneToOneField, so handing one
    observation to two poses is an IntegrityError rather than a bad row.
    """
    return set(
        Pose.objects.filter(main_site_observation__isnull=False).values_list(
            "main_site_observation_id", flat=True
        )
    )


# --------------------------------------------------------------------------
# planners
# --------------------------------------------------------------------------


def plan_fix_superseded_mains(
    target=None, *, kinds: Iterable[str] | None = None
) -> Plan:
    """Re-point poses whose main is superseded at the right live observation.

    Only the unambiguous diagnoses are acted on. ``SUCCESSOR_OTHER_POSE`` (the
    canon-site hops and split groups) and ``NO_SUCCESSOR_DEAD`` are reported in
    ``plan.skipped`` instead - both need a decision this code cannot make.

    Pass ``kinds`` to override that, once you have made the decision.
    """
    wanted = frozenset(kinds) if kinds is not None else UNAMBIGUOUS
    plan = Plan(name="fix superseded pose mains")
    taken = _mains_taken()

    for diag in diagnose_superseded_mains(target):
        label = f"pose {diag.pose.pk}"
        if diag.kind not in wanted:
            plan.skipped.append((label, f"{diag.kind}: {diag.note}"))
            continue

        # live_members is ordered newest-version-first by diagnose_*
        new_main = (
            diag.live_members[0]
            if diag.kind == NO_SUCCESSOR_LIVE_MEMBERS
            else diag.successor
        )
        if new_main is None:
            # reachable only via an explicit kinds= override that selects a
            # diagnosis with nothing to promote (NO_SUCCESSOR_DEAD, NULL_IDENTITY)
            plan.skipped.append(
                (label, f"{diag.kind}: no observation available to become main")
            )
            continue

        if new_main.pk in taken and new_main.pk != diag.pose.main_site_observation_id:
            plan.skipped.append(
                (
                    label,
                    f"observation {new_main.pk} is already the main of another "
                    "pose; main_site_observation is one-to-one",
                )
            )
            continue

        if new_main.pose_id != diag.pose.pk:
            plan.actions.append(
                SetPoseMember(
                    obs_id=new_main.pk,
                    pose_id=diag.pose.pk,
                    was_pose_id=new_main.pose_id,
                )
            )
        plan.actions.append(
            SetPoseMain(
                pose_id=diag.pose.pk,
                new_main_id=new_main.pk,
                old_main_id=diag.pose.main_site_observation_id,
            )
        )
        taken.discard(diag.pose.main_site_observation_id)
        taken.add(new_main.pk)

    return plan


def plan_sever_superseded_members(target=None) -> Plan:
    """Clear ``pose`` on superseded observations, so they stop being members.

    Skips any row that is still its pose's main - severing that would strand the
    main outside its own pose. Run :func:`plan_fix_superseded_mains` first and
    those rows become severable.
    """
    plan = Plan(name="sever superseded pose members")
    taken = _mains_taken()

    for obs in find_superseded_members(target).order_by("pk"):
        if obs.pk in taken:
            plan.skipped.append(
                (
                    f"observation {obs.pk}",
                    f"still the main of pose {obs.pose_id}; fix the main first",
                )
            )
            continue
        plan.actions.append(
            SetPoseMember(obs_id=obs.pk, pose_id=None, was_pose_id=obs.pose_id)
        )

    return plan


def plan_unsupersede(obs_ids: Iterable[int]) -> Plan:
    """Clear the superseded flag on specific observations.

    Deliberately takes explicit pks and never a target: "nothing live shares
    this identity" (:func:`find_wrongly_superseded`) is evidence, not a verdict.
    Look at the rows, decide, then pass the ones you mean.
    """
    plan = Plan(name="un-supersede observations")
    for obs in SiteObservation.objects.filter(pk__in=list(obs_ids)).order_by("pk"):
        if not obs.superseded:
            plan.skipped.append((f"observation {obs.pk}", "already live"))
            continue
        live = find_successors(obs).count()
        if live:
            plan.skipped.append(
                (
                    f"observation {obs.pk}",
                    f"{live} live row(s) already share its identity; "
                    "un-superseding would create a duplicate",
                )
            )
            continue
        plan.actions.append(SetSuperseded(obs_id=obs.pk, value=False))
    return plan


def plan_repair_display_names(target=None) -> Plan:
    """Rename poses whose ``display_name`` names no member, after the main."""
    plan = Plan(name="repair pose display names")
    mismatches = {m["pose_id"]: m for m in find_display_name_mismatches(target)}
    if not mismatches:
        return plan

    qs = Pose.objects.filter(pk__in=mismatches).select_related("main_site_observation")
    for pose in qs:
        main = pose.main_site_observation
        if main is None or not main.code:
            plan.skipped.append(
                (f"pose {pose.pk}", "no main observation to take a name from")
            )
            continue
        if main.superseded:
            plan.skipped.append(
                (f"pose {pose.pk}", "main is superseded; fix the main first")
            )
            continue
        plan.actions.append(
            SetDisplayName(pose_id=pose.pk, value=main.code, was=pose.display_name)
        )
    return plan


def plan_delete_dead_poses(target=None) -> Plan:
    """Delete poses with no live member at all.

    Destructive, and it destroys curation: a pose carries decisions a person
    made. Only poses that are *also* empty of any successor elsewhere are
    proposed; anything whose observations moved to another pose is reported so
    you can see the move first.
    """
    plan = Plan(name="delete dead poses")
    for pose in find_dead_poses(target).select_related("main_site_observation"):
        main = pose.main_site_observation
        if main is None:
            plan.actions.append(DeletePose(pose_id=pose.pk))
            continue
        successor = find_successors(main).first()
        if successor is not None:
            plan.skipped.append(
                (
                    f"pose {pose.pk}",
                    f"its main has a live successor ({successor.pk}) in pose "
                    f"{successor.pose_id}; this is a leftover shell, confirm the "
                    "move before deleting",
                )
            )
            continue
        plan.actions.append(DeletePose(pose_id=pose.pk))
    return plan
