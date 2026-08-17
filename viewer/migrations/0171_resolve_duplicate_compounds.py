# Generated migration to resolve duplicate compounds using auto-merge logic
# Created 2026-08-03

import logging
from django.db import migrations
from django.db.models import Count
from viewer.compound_curation import CONTENT_FIELDS

logger = logging.getLogger(__name__)


def _present(value) -> bool:
    """True if the value is real content (not NULL / not blank)."""
    if value is None:
        return False
    if isinstance(value, str) and value.strip() == "":
        return False
    return True


def resolve_duplicate_compounds(apps, schema_editor):
    """Auto-merge duplicate compounds with no genuine conflicts.

    For each inchi_key with multiple compounds:
    - NULL/empty values use existing value (not a conflict)
    - Description/comments concatenate with " | "
    - Genuine conflicts (both sides different) are logged but not resolved
    - Merge winner = first compound (lowest ID), others deleted
    - Relink all foreign key and ManyToMany references to point to winner

    Models updated:
    - ExperimentCompound, Pose, SiteObservation, CompoundIdentifier,
      ActivityPoint, Result (ForeignKey)
    - DesignSet_compounds (ManyToMany through table)
    """
    logger.info("="*70)
    logger.info("STARTING COMPOUND DUPLICATION RESOLUTION")
    logger.info("="*70)

    Compound = apps.get_model('viewer', 'Compound')
    Project = apps.get_model('viewer', 'Project')
    ExperimentCompound = apps.get_model('viewer', 'ExperimentCompound')
    Pose = apps.get_model('viewer', 'Pose')
    SiteObservation = apps.get_model('viewer', 'SiteObservation')
    CompoundIdentifier = apps.get_model('viewer', 'CompoundIdentifier')
    ActivityPoint = apps.get_model('viewer', 'ActivityPoint')
    Result = apps.get_model('viewer', 'Result')
    DesignSet = apps.get_model('viewer', 'DesignSet')
    DesignSetCompounds = apps.get_model('viewer', 'DesignSet_compounds')

    CONCATENABLE_FIELDS = {"description", "comments"}

    stats = {
        'projects_processed': 0,
        'duplicates_found': 0,
        'auto_merged': 0,
        'genuine_conflicts': 0,
        'references_relinked': 0,
        'merged_details': [],
    }

    projects = list(Project.objects.all())
    logger.info(f"Loading {len(projects)} project(s) to scan")

    for project in projects:
        stats['projects_processed'] += 1
        logger.info(f"[{stats['projects_processed']}/{len(projects)}] Scanning project: {project.title}")

        # Group compounds by (inchi_key, compound_code) to handle different codes separately
        duplicates = list(
            Compound.objects.filter(project=project, inchi_key__gt='')
            .values('inchi_key', 'compound_code')
            .annotate(count=Count('id'))
            .filter(count__gt=1)
        )

        if duplicates:
            logger.info(f"Found {len(duplicates)} duplicate group(s) in project {project.title}")

        for dup in duplicates:
            inchi_key = dup['inchi_key']
            compound_code = dup['compound_code']
            compounds = list(
                Compound.objects.filter(project=project, inchi_key=inchi_key, compound_code=compound_code)
                .order_by('id')
            )

            if len(compounds) < 2:
                continue

            stats['duplicates_found'] += 1
            winner = compounds[0]  # Keep the lowest ID
            losers = compounds[1:]


            # Check for genuine conflicts
            genuine_conflicts = {}
            merged_values = {}

            for field in CONTENT_FIELDS:
                winner_val = getattr(winner, field)
                loser_vals = [getattr(c, field) for c in losers]

                # Collect unique non-empty values
                unique_vals = set()
                for val in [winner_val] + loser_vals:
                    if _present(val):
                        unique_vals.add(str(val))

                if len(unique_vals) == 0:
                    # All empty -> no conflict, no change needed
                    continue
                elif len(unique_vals) == 1:
                    # All same value -> no conflict
                    merged_values[field] = list(unique_vals)[0]
                    continue

                # Multiple different values
                if field in CONCATENABLE_FIELDS:
                    # Concatenate description/comments
                    parts = []
                    if _present(winner_val):
                        parts.append(str(winner_val).strip())
                    for loser_val in loser_vals:
                        if _present(loser_val):
                            parts.append(str(loser_val).strip())
                    merged_values[field] = " | ".join(parts)
                else:
                    # Genuine conflict -> log but don't resolve
                    genuine_conflicts[field] = {
                        'winner': str(winner_val),
                        'losers': [str(v) for v in loser_vals],
                    }

            # If no genuine conflicts, proceed with merge
            if not genuine_conflicts:
                # Update winner with merged values
                for field, value in merged_values.items():
                    setattr(winner, field, value)
                winner.save()

                # Relink all foreign key references from losers to winner
                loser_ids = [c.id for c in losers]
                relinked = 0

                # ExperimentCompound.compound -> winner
                relinked += ExperimentCompound.objects.filter(compound__in=losers).update(compound=winner)

                # Pose.compound -> winner
                relinked += Pose.objects.filter(compound__in=losers).update(compound=winner)

                # SiteObservation.cmpd -> winner
                relinked += SiteObservation.objects.filter(cmpd__in=losers).update(cmpd=winner)

                # CompoundIdentifier.compound -> winner
                relinked += CompoundIdentifier.objects.filter(compound__in=losers).update(compound=winner)

                # ActivityPoint.cmpd_id -> winner
                relinked += ActivityPoint.objects.filter(cmpd_id__in=losers).update(cmpd_id=winner)

                # Result.compound -> winner (nullable)
                relinked += Result.objects.filter(compound__in=losers).update(compound=winner)

                # DesignSet_compounds.compound -> winner (ManyToMany through table)
                relinked += DesignSetCompounds.objects.filter(compound__in=losers).update(compound=winner)

                stats['references_relinked'] += relinked

                # Delete losers (cascade should handle any remaining references)
                for loser in losers:
                    loser.delete()

                stats['auto_merged'] += 1
                merged_field_names = list(merged_values.keys())
                logger.debug(f"AUTO-MERGED InChIKey {inchi_key[:16]}... (IDs: {[winner.id] + loser_ids}): Relinked {relinked} refs")
                if merged_field_names:
                    logger.debug(f"  Fields consolidated: {', '.join(merged_field_names)}")
                stats['merged_details'].append({
                    'inchi_key': inchi_key,
                    'winner_id': winner.id,
                    'loser_ids': loser_ids,
                    'merged_fields': merged_field_names,
                    'references_relinked': relinked,
                })
            else:
                # Genuine conflicts remain -> needs manual curation
                stats['genuine_conflicts'] += 1
                conflict_fields = list(genuine_conflicts.keys())
                logger.info(f"SKIPPED InChIKey {inchi_key[:16]}... (IDs: {loser_ids}): Has genuine conflicts in {conflict_fields}")
                # Log detailed conflict info for debugging
                for field, conflict_info in genuine_conflicts.items():
                    logger.debug(f"  {field}: winner={conflict_info['winner']}, losers={conflict_info['losers']}")

    # Log results
    logger.info("="*70)
    logger.info("COMPOUND DUPLICATION RESOLUTION SUMMARY")
    logger.info("="*70)
    logger.info(f"Projects processed: {stats['projects_processed']}")
    logger.info(f"Duplicate groups found: {stats['duplicates_found']}")
    logger.info(f"Auto-merged: {stats['auto_merged']}")
    logger.info(f"Genuine conflicts (not resolved): {stats['genuine_conflicts']}")
    logger.info(f"Foreign key references relinked: {stats['references_relinked']}")

    if stats['merged_details']:
        logger.info("Merged compounds:")
        for detail in stats['merged_details']:
            logger.info(f"  InChIKey {detail['inchi_key'][:16]}... :")
            logger.info(f"    Winner: {detail['winner_id']}, Merged losers: {detail['loser_ids']}")
            if detail['merged_fields']:
                logger.info(f"    Fields consolidated: {', '.join(detail['merged_fields'])}")
            logger.info(f"    References relinked: {detail['references_relinked']}")
    logger.info("="*70)


def reverse_merge(apps, schema_editor):
    """Reverse is not possible - merged compounds are deleted."""
    pass


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0170_merge_20260722_1450'),
    ]

    operations = [
        migrations.RunPython(resolve_duplicate_compounds, reverse_merge),
    ]
