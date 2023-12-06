import os

from django.conf import settings
from django.core.management.base import BaseCommand
from django.utils import timezone
from viewer.models import Target, SiteObservation, SiteObservationTag, TagCategory

# from scoring.models import MolGroup
from scoring.models import SiteObservationGroup


class Command(BaseCommand):
    help = 'Add moleculeTag record for existing mol_groups for a given target. This effectively adds molecule tags for all the sites for the Target'

    def add_arguments(self, parser):
        parser.add_argument('target', type=str, help='Target to be corrected')
        parser.add_argument(
            'update',
            type=str,
            help='Whether to update the target (yes) or display what will be updated (no)',
        )

    def handle(self, *args, **kwargs):
        # Unused args
        del args

        tags_existing = 0
        tags_to_be_created = 0
        tags_created = 0

        time_start = timezone.now().strftime('%X')
        self.stdout.write("Start %s" % time_start)
        target_name = kwargs['target']
        self.stdout.write("target_name: %s" % target_name)
        if kwargs['update'] == 'yes':
            update = True
        else:
            update = False
        self.stdout.write("update: %s" % update)

        target = Target.objects.filter(title=target_name)

        if not target:
            self.stdout.write("Target %s not found" % target_name)
            exit(1)
        else:
            self.stdout.write("Updating tags for Target %s" % target[0].title)

        # First, try sites file - e.g. /code/media/targets/mArh/sites.csv
        # If this is there, then the new sites functionality was used.
        sites_filepath = os.path.join(
            settings.MEDIA_ROOT, 'targets', target_name, 'sites.csv'
        )
        if os.path.isfile(sites_filepath):
            expected_sites = (
                sum(1 for line in open(sites_filepath, encoding='utf-8')) - 1
            )
            self.stdout.write("Expected number of sites: %s" % expected_sites)
            # These should correspond to the sites for the target held in sites.csv
            mol_groups = SiteObservationGroup.objects.filter(
                target_id__title=target_name, group_type="MC"
            )
            tag_type = 'site'
        else:
            # The sites should correspond to the centres of mass. The sites will be generated from them
            mol_groups = SiteObservationGroup.objects.filter(
                target_id__title=target_name, group_type="MC", description="c_of_m"
            )
            expected_sites = len(mol_groups)
            self.stdout.write("Expected number of sites: %s" % expected_sites)
            tag_type = 'c_of_e'

        if not mol_groups:
            self.stdout.write("No sites found for target")
            exit(1)

        for idx, mol_group in enumerate(mol_groups):
            self.stdout.write(
                "mol_group description: {}, index: {}".format(
                    mol_group.description, idx
                )
            )
            # A molecule tag record should not exist, but if it does go no further
            try:
                mol_tag = SiteObservationTag.objects.get(mol_group=mol_group)
            except:
                mol_tag = None

            if tag_type == 'site':
                tag_name = mol_group.description
            else:
                tag_name = 'c_of_m_{}'.format(idx)

            if mol_tag:
                self.stdout.write(
                    "Tag already exists for {}, index: {}".format(
                        mol_group.description, idx
                    )
                )
                tags_existing += 1
                continue
            else:
                self.stdout.write("Tag to be created for %s" % tag_name)
                self.stdout.write("    mol_tag.tag = %s" % tag_name)
                self.stdout.write(
                    "    mol_tag.category = %s"
                    % TagCategory.objects.get(category='Sites')
                )
                self.stdout.write("    mol_tag.target = %s" % target[0])
                self.stdout.write("    mol_tag.mol_group = %s" % mol_group)
                self.stdout.write(
                    "    mol_tag.molecules = %s"
                    % [mol['id'] for mol in mol_group.mol_id.values()]
                )
                tags_to_be_created += 1

            # If update flag is set then actually create molecule Tags.
            if update:
                mol_tag = SiteObservationTag()
                mol_tag.tag = tag_name
                mol_tag.category = TagCategory.objects.get(category='Sites')
                mol_tag.target = target[0]
                mol_tag.mol_group = mol_group
                mol_tag.save()
                for mol in mol_group.mol_id.values():
                    this_mol = SiteObservation.objects.get(id=mol['id'])
                    mol_tag.molecules.add(this_mol)
                tags_created += 1

        self.stdout.write("Expected number of sites: %s" % expected_sites)
        self.stdout.write("tags_existing %s" % tags_existing)
        self.stdout.write("tags_to_be_created %s" % tags_to_be_created)
        self.stdout.write("tags_created: %s" % tags_created)

        if tags_to_be_created == expected_sites:
            self.stdout.write('Looking good - tags_to_be_created = expected sites')

        if tags_created == expected_sites:
            self.stdout.write('Looking good - tags_created = expected sites')

        time_end = timezone.now().strftime('%X')
        self.stdout.write("End %s" % time_end)
