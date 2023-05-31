import json
import os

from django.core.management.base import BaseCommand, CommandError
from viewer.models import JobOverride


class Command(BaseCommand):
    help = "Initialises the JobOverride table"

    def handle(self, *args, **options):
        
        day_1_override_file = os.path.join('/code', 'viewer', 'squonk', 'day-1-job-override.json')
        if not os.path.isfile(day_1_override_file):
            raise CommandError('Day-1 override file "%s" does not exist' % day_1_override_file)

        jo_count = JobOverride.objects.all().count()
        if jo_count == 0:           
            with open(day_1_override_file, 'r') as f:
                day_1_override = json.load(f)
            jo = JobOverride()
            jo.override = day_1_override
            jo.save()
            self.stdout.write(
                self.style.SUCCESS('Successfully added initial JobOverride record')
            )
        else:
            self.stdout.write(
                self.style.SUCCESS('Nothing to initialise in JobOverride - record count is %s' % jo_count)
            )
