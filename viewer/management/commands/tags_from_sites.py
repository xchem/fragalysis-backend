import os

from django.conf import settings
from django.core.management.base import BaseCommand
from django.utils import timezone
from viewer.models import (
    Target,
    Molecule,
    MoleculeTag,
    TagCategory
)
from scoring.models import MolGroup

class Command(BaseCommand):
    help = 'Add moleculeTag record for existing mol_groups for a given target. This effectively adds molecule tags for all the sites for the Target'

    def add_arguments(self, parser):
        parser.add_argument('target', type=str, help='Target to be corrected')
        parser.add_argument('update', type=str, help='Whether to update the target (yes) or display what will be updated (no)')

    def handle(self, *args, **kwargs):
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

        # e.g. /code/media/targets/mArh
        sites_filepath = os.path.join(settings.MEDIA_ROOT, 'targets', target_name, 'sites.csv')
        expected_sites = sum(1 for line in open(sites_filepath)) - 1
        self.stdout.write("Expected number of sites: %s" % expected_sites)

        target = Target.objects.filter(title=target_name)

        if not target:
            self.stdout.write("Target %s not found" % target_name)
            exit(1)
        else:
            self.stdout.write("Updating tags for Target %s" % target[0].title)

        # These should correspond to the sites for the target held in sites.csv
        mol_groups = MolGroup.objects.filter(target_id__title=target_name, group_type = "MC", )

        if not mol_groups:
            self.stdout.write("No sites found for target")
            exit(1)

        for mol_group in mol_groups:
            self.stdout.write("site name %s" % mol_group.description)
            # A molecule tag record should not exist, but if it does go no further
            try:
                mol_tag = MoleculeTag.objects.get(mol_group=mol_group)
            except:
                mol_tag = None

            if mol_tag:
                self.stdout.write("Tag already exists for %s" % mol_group.description)
                tags_existing += 1
                continue
            else:
                self.stdout.write("Tag to be created for %s" % mol_group.description)
                self.stdout.write("    mol_tag.tag = %s" % mol_group.description)
                self.stdout.write("    mol_tag.category = %s" % TagCategory.objects.get(category='Sites'))
                self.stdout.write("    mol_tag.target = %s" % target[0])
                self.stdout.write("    mol_tag.mol_group = %s" % mol_group)
                self.stdout.write("    mol_tag.molecules = %s" % [mol['id'] for mol in mol_group.mol_id.values()])
                tags_to_be_created += 1


            # If update flag is set then actually create molecule Tags.
            if update:
                mol_tag = MoleculeTag()
                mol_tag.tag = mol_group.description
                mol_tag.category = TagCategory.objects.get(category='Sites')
                mol_tag.target = target[0]
                mol_tag.mol_group = mol_group
                mol_tag.save()
                for mol in mol_group.mol_id.values():
                    this_mol = Molecule.objects.get(id=mol['id'])
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
