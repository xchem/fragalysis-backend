"""Create the media S3 bucket if it does not exist (issue 1038).

For compose stacks (development and the integration tests), which run an
S3-compatible store (SeaweedFS) that starts empty. A one-shot service runs this
before the backend starts, because the start-up check (viewer.E003) needs the
bucket. It uses MEDIA_S3_STORAGE_OPTIONS, and works with any S3 endpoint::

    python manage.py create_media_bucket --wait-s 60

A production bucket is created outside the stack (i.e. on STFC Echo).
"""
import time

from botocore.exceptions import BotoCoreError, ClientError
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured
from django.core.management.base import BaseCommand, CommandError

from frag_s3_storage import S3Storage

# The bucket does not exist yet, so the start-up checks (viewer.E003) would fail.
_NO_SYSTEM_CHECKS: list = []


class Command(BaseCommand):
    help = "Create the media S3 bucket (MEDIA_S3_BUCKET_NAME) if it does not exist."
    requires_system_checks = _NO_SYSTEM_CHECKS

    def add_arguments(self, parser):
        parser.add_argument(
            "--wait-s",
            type=int,
            default=60,
            help="Try to reach the S3 endpoint for this long, once a second"
            " (default: 60).",
        )

    def handle(self, *args, **kwargs):
        del args
        options = settings.MEDIA_S3_STORAGE_OPTIONS
        try:
            storage = S3Storage(**options)
        except ImproperlyConfigured as exc:
            raise CommandError(f"{exc} (set MEDIA_S3_BUCKET_NAME)") from exc

        if self._exists(storage, attempts=max(1, kwargs["wait_s"])):
            self.stdout.write(f"The bucket '{storage.bucket_name}' already exists.")
            return

        region = options.get("region_name")
        if region and region != "us-east-1":
            storage.client.create_bucket(
                Bucket=storage.bucket_name,
                CreateBucketConfiguration={"LocationConstraint": region},
            )
        else:
            storage.client.create_bucket(Bucket=storage.bucket_name)
        self.stdout.write(
            self.style.SUCCESS(f"Created the bucket '{storage.bucket_name}'.")
        )

    @staticmethod
    def _exists(storage: S3Storage, attempts: int) -> bool:
        """True if the bucket exists. Tries again while the endpoint is starting."""
        for attempt in range(attempts):
            try:
                storage.client.head_bucket(Bucket=storage.bucket_name)
                return True
            except ClientError as exc:
                if exc.response.get("Error", {}).get("Code") in ("404", "NoSuchBucket"):
                    return False
                raise CommandError(
                    f"The S3 endpoint refused the request: {exc}"
                ) from exc
            except BotoCoreError as exc:
                if attempt + 1 == attempts:
                    raise CommandError(
                        f"The S3 endpoint could not be reached: {exc}"
                    ) from exc
                time.sleep(1)
        return False
