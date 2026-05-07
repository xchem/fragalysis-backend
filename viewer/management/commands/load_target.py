import shutil
import sys
from pathlib import Path

from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand
from django.db import IntegrityError

from viewer.cache import clear_all_view_caches
from viewer.target_loader import load_target


class Command(BaseCommand):
    help = 'Load target directly from media/ directory.'

    def add_arguments(self, parser):
        parser.add_argument('tarball', type=str, help='Data archive to load')
        parser.add_argument('-u', '--username', type=str, help='User FedID')
        parser.add_argument(
            '-p',
            '--proposal_ref',
            type=str,
            required=True,
            help='Proposal number (target access string)',
        )

    def handle(self, *args, **kwargs):
        # Unused args
        del args

        # the loader, internally expects user ID, not username,
        # because that's how the view resolves it. translate to id
        username = kwargs.get('username', None)
        if not username:
            if settings.DEPLOYMENT_MODE == 'DEPLOYMENT':
                msg = 'Username is required'
                self.stdout.write(self.style.ERROR(msg))
                sys.exit()
            else:
                user_id = settings.ANONYMOUS_USER
        else:
            try:
                user_id = get_user_model().objects.get(username=username).pk
            except get_user_model().DoesNotExist:
                msg = (
                    f'User {username} does not exist in the database.'
                    + ' Log into fragalysis first to ensure your user is'
                    ' populated.'
                )
                self.stdout.write(self.style.ERROR(msg))
                sys.exit()

        # upload view uploads the tarball to /code/media/tmp/.
        # the user may not
        tarball_path = Path(kwargs['tarball'])
        temp_path = Path(settings.MEDIA_ROOT).joinpath('tmp')
        temp_path.mkdir(exist_ok=True)

        # shutil instead of path.rename because doesn't work in local
        # deployment with mounted volume
        shutil.move(str(tarball_path), str(temp_path.joinpath(tarball_path.name)))

        try:
            load_target(
                str(temp_path.joinpath(tarball_path.name)),
                proposal_ref=kwargs['proposal_ref'],
                user_id=user_id,
            )
            # SiteObservation, Pose and SiteObservationTag rows have just
            # been written; flush cached responses that read them.
            clear_all_view_caches()
            # self.stdout.write(self.style.SUCCESS('Data imported'))
        except KeyError as err:
            self.stdout.write(self.style.ERROR(err.args[0]))
        except IntegrityError as err:
            self.stdout.write(self.style.ERROR(err.args[0]))
        except FileNotFoundError as err:
            self.stdout.write(self.style.ERROR(err.args[0]))
