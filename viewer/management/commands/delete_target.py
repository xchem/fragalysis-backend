"""Delete a target's latest upload, or the whole target.

Nothing is deleted unless an action is asked for explicitly. With no arguments at all
it prints this help; name a target without asking for an action and it shows that
target's uploads and the commands that would delete them::

    python manage.py delete_target --pk 1

    Uploads:

      --pk 1   A71EV2A   (proposal lb18145-1)
          upload 1  upload_1   2026-07-22    964 observations
          upload 2  upload_2   2026-07-23      6 observations

    To delete the latest upload:
        python manage.py delete_target --pk 1 --latest-upload --dry-run
        python manage.py delete_target --pk 1 --latest-upload

The target is named by ``--pk`` or by ``--title``. A title is unique only within a
proposal (Target has a unique constraint on (title, project)), so ``--title`` may need
``--proposal`` alongside it; ``--pk`` never does.

Deleting the latest upload restores the database to the state it was in after the
previous one: the rows the upload created go, ``superseded`` flags are recomputed, pose
mains are given their predecessors back, compounds the upload introduced are removed if
nothing else points at them, and the upload's media files are swept. It is repeatable -
run it again and the next upload down comes off - until only the first upload remains,
at which point ``--entire-target`` is the only thing left to do.

Both actions are disabled in production, matching ``TargetView.destroy``. The guard is
repeated here because a management command has no request to carry it.
"""

from django.core.management.base import BaseCommand, CommandError

from api.utils import deployment_mode_is_production
from viewer.models import ExperimentUpload, SiteObservation, Target
from viewer.target_delete import delete_target
from viewer.upload_delete import (
    NotDeletable,
    delete_latest_upload,
    plan_upload_deletion,
)


class Command(BaseCommand):
    help = (
        "Delete a target's latest experiment upload (--latest-upload) or the whole "
        "target (--entire-target). Name the target with --pk, or --title plus "
        "--proposal. Naming a target without an action lists its uploads. "
        "Disabled in production."
    )

    def add_arguments(self, parser):
        which = parser.add_argument_group("which target")
        which.add_argument(
            "--pk",
            type=int,
            help="Primary key of the target. Always unambiguous.",
        )
        which.add_argument(
            "--title",
            help="Title of the target, e.g. A71EV2A. May need --proposal.",
        )
        which.add_argument(
            "--proposal",
            help=(
                "Proposal / target-access string (the Project title). Only needed "
                "when --title names more than one target."
            ),
        )

        action = parser.add_argument_group("what to delete")
        action.add_argument(
            "--latest-upload",
            action="store_true",
            help="Delete the target's latest experiment upload.",
        )
        action.add_argument(
            "--entire-target",
            action="store_true",
            help=(
                "Delete the whole target - every upload, the full database graph and "
                "all its media."
            ),
        )

        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Report what would be deleted without changing anything.",
        )
        parser.add_argument(
            "--detail",
            action="store_true",
            help="List the primary keys of everything affected, not just the counts.",
        )

    def handle(self, *args, **options):
        del args
        # Never in production - the same rule target deletion already lives under.
        # There is no request here, so the check cannot be delegated to a permission
        # class.
        if deployment_mode_is_production():
            raise CommandError("Target and upload deletion are disabled in production.")

        if options["latest_upload"] and options["entire_target"]:
            raise CommandError(
                "--latest-upload and --entire-target are different operations; ask "
                "for one of them."
            )

        if not (options["latest_upload"] or options["entire_target"]):
            # Nothing named and nothing asked for: the caller wants to know how to use
            # this, not to be shown the whole database.
            if options["pk"] is None and not options["title"]:
                self.print_help("manage.py", "delete_target")
                return
            # A target was named but no action. Listing is not destructive, so an
            # ambiguous title is answered by showing every match with its pk and
            # proposal - which is exactly how a caller finds the one they want -
            # rather than by refusing as the delete paths do.
            self._list(self._matching_targets(options))
            return

        target = self._resolve_target(options)
        if target is None:
            raise CommandError(
                "Name the target to act on with --pk, or --title (plus --proposal if "
                "the title is ambiguous). Run with no arguments for help."
            )

        if options["entire_target"]:
            self._delete_entire_target(target, options)
        else:
            self._delete_latest_upload(target, options)

    # ------------------------------------------------------------------ #
    # Finding the target
    # ------------------------------------------------------------------ #

    def _target_queryset(self, options):
        """Targets matching --pk / --title / --proposal, newest constraint last."""
        queryset = Target.objects.select_related("project")
        if options["proposal"]:
            queryset = queryset.filter(project__title=options["proposal"])
        if options["pk"] is not None:
            queryset = queryset.filter(pk=options["pk"])
        if options["title"]:
            queryset = queryset.filter(title=options["title"])
        return queryset.order_by("pk")

    def _matching_targets(self, options):
        """Every target the options match, for the non-destructive listing path."""
        if options["pk"] is not None and options["title"]:
            raise CommandError("Give either --pk or --title, not both.")
        targets = list(self._target_queryset(options))
        if not targets:
            raise CommandError(self._nothing_matched(options))
        return targets

    def _nothing_matched(self, options) -> str:
        if options["pk"] is not None:
            return f"No target with --pk {options['pk']}."
        where = f" in proposal '{options['proposal']}'" if options["proposal"] else ""
        return f"No target titled '{options['title']}'{where}."

    def _resolve_target(self, options):
        """The single target the options name, or None when none were given.

        Returning None rather than erroring is what lets a bare invocation list every
        target: with nothing to act on, showing what exists is more useful than a
        usage message.
        """
        if options["pk"] is None and not options["title"]:
            return None

        if options["pk"] is not None and options["title"]:
            raise CommandError("Give either --pk or --title, not both.")

        queryset = Target.objects.select_related("project")
        if options["proposal"]:
            queryset = queryset.filter(project__title=options["proposal"])

        if options["pk"] is not None:
            target = queryset.filter(pk=options["pk"]).first()
            if target is None:
                raise CommandError(f"No target with --pk {options['pk']}.")
            return target

        candidates = list(queryset.filter(title=options["title"]))
        if not candidates:
            where = (
                f" in proposal '{options['proposal']}'" if options["proposal"] else ""
            )
            raise CommandError(f"No target titled '{options['title']}'{where}.")
        if len(candidates) > 1:
            listing = "\n".join(
                f"    --pk {t.pk}   {t.title}   (proposal {t.project.title})"
                for t in sorted(candidates, key=lambda t: t.pk)
            )
            raise CommandError(
                f"--title '{options['title']}' names {len(candidates)} targets - a "
                f"title is unique only within a proposal. Add --proposal, or use "
                f"--pk:\n{listing}"
            )
        return candidates[0]

    # ------------------------------------------------------------------ #
    # Listing
    # ------------------------------------------------------------------ #

    def _list(self, targets):
        """Show each matching target's uploads, and the commands that act on them."""
        if len(targets) > 1:
            self.stdout.write(
                f"{len(targets)} targets match - a title is unique only within a "
                f"proposal. Use --pk to pick one.\n"
            )
        self.stdout.write("Uploads:\n")
        for item in targets:
            self.stdout.write(
                f"  --pk {item.pk}   {item.title}   (proposal {item.project.title})"
            )
            uploads = list(
                ExperimentUpload.objects.filter(target=item).order_by("upload_version")
            )
            if not uploads:
                self.stdout.write("      (no uploads)")
                continue
            for upload in uploads:
                # Observations are counted by the version that created them, not
                # through the experiment FK: a carried-over experiment holds
                # observations from several uploads, so the FK would credit them all
                # to whichever upload created the experiment.
                created = SiteObservation.objects.filter(
                    experiment__experiment_upload__target=item,
                    version=upload.upload_version,
                ).count()
                when = (
                    upload.commit_datetime.date().isoformat()
                    if upload.commit_datetime
                    else "-"
                )
                self.stdout.write(
                    f"      upload {upload.upload_version}  "
                    f"{str(upload.upload_data_dir or '-'):10} "
                    f"{when}  {created:5} observations"
                )

        example = targets[0]
        self.stdout.write("\nTo delete the latest upload:")
        self.stdout.write(
            f"    python manage.py delete_target --pk {example.pk} "
            "--latest-upload --dry-run"
        )
        self.stdout.write(
            f"    python manage.py delete_target --pk {example.pk} --latest-upload"
        )
        self.stdout.write("\nTo delete the target outright:")
        self.stdout.write(
            f"    python manage.py delete_target --pk {example.pk} --entire-target"
        )

    # ------------------------------------------------------------------ #
    # Deleting
    # ------------------------------------------------------------------ #

    def _delete_entire_target(self, target, options):
        uploads = ExperimentUpload.objects.filter(target=target).count()
        if options["dry_run"]:
            self.stdout.write(
                f"Would delete target {target.title} (pk={target.pk}) entirely: "
                f"{uploads} upload(s), its whole database graph and all its media."
            )
            return
        delete_target(target)
        self.stdout.write(
            self.style.SUCCESS(
                f"Deleted target {target.title}: {uploads} upload(s) and all media."
            )
        )

    def _delete_latest_upload(self, target, options):
        upload = (
            ExperimentUpload.objects.filter(target=target)
            .order_by("upload_version")
            .last()
        )
        if upload is None:
            raise CommandError(f"Target '{target.title}' has no uploads.")

        try:
            if options["dry_run"]:
                plan = plan_upload_deletion(upload)
            else:
                plan = delete_latest_upload(upload)
        except NotDeletable as exc:
            # The service is surface-agnostic, so it cannot name this command's flag -
            # but only one refusal is actually answered by deleting the whole target.
            hint = (
                " Use --entire-target to do that." if exc.whole_target_instead else ""
            )
            raise CommandError(f"{exc}{hint}") from exc

        self.stdout.write(plan.summary())
        if options["detail"]:
            self._write_detail(plan)

        if options["dry_run"]:
            self.stdout.write(self.style.WARNING("Dry run - nothing was changed."))
        else:
            self.stdout.write(
                self.style.SUCCESS(
                    f"Deleted upload {plan.upload_version} of {plan.target_title}."
                )
            )

    def _write_detail(self, plan):
        for name, pks in sorted(plan.doomed.items()):
            if pks:
                self.stdout.write(f"  {name}: {sorted(pks)}")
        for name, pks in sorted(plan.orphaned.items()):
            if pks:
                self.stdout.write(f"  {name} (orphaned): {sorted(pks)}")
        for pose_pk, old, new in plan.poses_repointed:
            self.stdout.write(f"  Pose {pose_pk}: main {old} -> {new}")
        for model, obs_pk, old, new in plan.children_repointed:
            self.stdout.write(f"  SiteObservation {obs_pk}: {model} {old} -> {new}")
        for exp_pk, upload_pk in plan.experiments_repointed:
            self.stdout.write(f"  Experiment {exp_pk} -> ExperimentUpload {upload_pk}")
        for exp_pk, fields in plan.experiment_files_restored:
            self.stdout.write(f"  Experiment {exp_pk}: restored {', '.join(fields)}")
        for path in plan.files_kept:
            self.stdout.write(f"  kept (still referenced): {path}")
